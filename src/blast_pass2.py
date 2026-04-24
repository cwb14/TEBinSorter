"""
blast_pass2.py — MMseqs2-based pass-2 classification for HMM-unclassified
sequences.

Ported from TEBinSorter's original blastn-based pass-2 to mmseqs2, following
the edits on github.com/cwb14/TEsorter branch `my-new-idea2` (head b398509).
The filename, public symbols (`blast_pass2`, `store_blast_hits`,
`classify_from_blast`), and the SQLite `blast_hits` table are retained to keep
the rest of the pipeline (pipeline.py, tesorter_compat.py, classifier.py,
results.py) unchanged. mmseqs bits/fident/qcov are remapped to the BLAST-
shaped columns at insert time:

    pident   = fident * 100      (0..100)
    qcovs    = qcov   * 100      (0..100)
    length   = alnlen
    bitscore = bits
    evalue   = 0.0               (sentinel; unread downstream)
    slen     = 0                 (sentinel; unread downstream)

Sequences unclassified by HMM are searched against the classified pool with
a single multithreaded mmseqs2 easy-search (no chunking; mmseqs is internally
parallel). Split alignments are merged per (query, target) pair within a
500 bp gap window before best-hit selection.
"""

import logging
import os
import tempfile
import time
from collections import defaultdict

import pyfastx

import mmseqs
import pass2_external

log = logging.getLogger(__name__)


# Developer toggle, deliberately not exposed on the CLI.
#   True  -> require BOTH qcov AND tcov >= coverage threshold (strict).
#   False -> require AT LEAST ONE of qcov, tcov >= coverage threshold
#            (Wicker et al. 80-80-80 style: "candidate must cover ≥80% of
#            at least one of the elements being compared").
# mmseqs2 has no native "at-least-one" cov-mode, so either way this is
# enforced post-hoc in classify_from_blast's SQL WHERE clause.
REQUIRE_BOTH_COVERAGE = False


def _get_classified_ids(conn):
    """Get classified sequence IDs per database from classifier results.

    Returns:
        classified: {base_seq: set(databases)} — which databases classified each seq
    """
    classified = defaultdict(set)

    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}

    if "legacy_hits" in tables:
        for row in conn.execute(
            "SELECT DISTINCT base_seq, database FROM legacy_hits"
        ):
            classified[row[0]].add(row[1])

    return dict(classified)


def split_classified_unclassified(input_fasta, classified_ids, outdir,
                                  seq_type="nucl"):
    """Split input into classified-pool FASTA and unclassified-query FASTA.

    One pass through the input. When seq_type=="nucl", sequences written to
    the DB FASTA are uppercased and stripped of non-ATCG characters (mmseqs
    doesn't tolerate ambiguous bases as cleanly as blastn).
    """
    os.makedirs(outdir, exist_ok=True)

    db_fasta = os.path.join(outdir, "blast_db.fa")
    qry_fasta = os.path.join(outdir, "blast_query.fa")

    nucl = (seq_type == "nucl")
    fa = pyfastx.Fasta(input_fasta, build_index=True)

    n_classified = 0
    n_unclassified = 0
    db_seq_to_dbs = {}

    with open(db_fasta, "w") as dbh, open(qry_fasta, "w") as qh:
        for rec in fa:
            name = rec.name
            seq = str(rec.seq)
            if nucl:
                seq = "".join(c for c in seq.upper() if c in "ATCG")

            if name in classified_ids:
                dbh.write(f">{name}\n{seq}\n")
                db_seq_to_dbs[name] = classified_ids[name]
                n_classified += 1
            else:
                qh.write(f">{name}\n{seq}\n")
                n_unclassified += 1

    log.info(f"  Split: {n_classified} classified (DB), "
             f"{n_unclassified} unclassified (query)")

    return db_fasta, qry_fasta, db_seq_to_dbs


def run_mmseqs(query_fa, db_fa, m8_out, tmpdir,
               seqtype="nucl", ncpu=4,
               min_identity=80, min_coverage=80, min_length=80,
               cov_mode=0, sensitivity=None):
    """Run a single mmseqs2 easy-search.

    The thresholds (identity/coverage/length) are NOT passed as search-side
    filters. mmseqs2's prefilter + --cov-mode 0 + --min-seq-id interacts
    destructively on nucleotide data: hits that clearly pass the final
    thresholds post-alignment get dropped before they are scored. This mirrors
    TEBinSorter's blastn pass-2, which also runs unfiltered at search time and
    applies the 80-80-80 rule post-hoc. The final filter is enforced by
    classify_from_blast's SQL WHERE on the percent-scaled columns.

    cov_mode is still passed through — it shapes the geometric coverage
    definition mmseqs uses for `qcov`, which we then compare against the
    --pass2-rule coverage threshold in SQL.
    """
    mmseqs.mmseqs_easy_search(
        db_seq=db_fa,
        qry_seq=query_fa,
        out_m8=m8_out,
        tmpdir=tmpdir,
        seqtype=seqtype,
        ncpu=ncpu,
        min_seq_id=0.0,
        min_cov=0.0,
        cov_mode=cov_mode,
        min_aln_len=0,
        sensitivity=sensitivity,
    )
    return m8_out


def parse_mmseqs_output(m8_path):
    """Parse mmseqs2 M8, merge split alignments, return hit-dicts shaped for
    store_blast_hits. Percent columns (pident, qcovs, tcovs) are on 0-100;
    evalue is sentinel 0.0 because mmseqs's M8 format doesn't carry it and no
    downstream reader consumes it."""
    best = mmseqs.parse_mmseqs_m8_besthit(m8_path)

    hits = []
    for qid, rc in best.items():
        hits.append({
            "qseqid": rc.qseqid,
            "sseqid": rc.sseqid,
            "pident": rc.fident * 100.0,
            "length": rc.alnlen,
            "evalue": 0.0,
            "bitscore": rc.bits,
            "qlen": rc.qlen,
            "slen": rc.tlen,
            "qcovs": rc.qcov * 100.0,
            "tcovs": rc.tcov * 100.0,
        })
    return hits


def store_blast_hits(conn, hits, db_seq_to_dbs):
    """Store hits in SQLite blast_hits table.

    Adds a tcovs column so the pass-2 filter can enforce target-side coverage.
    classified_by records which databases flagged the target. For sequences
    that came from --pass2-classified-fasta, db_seq_to_dbs contains {"external"}.
    """
    conn.execute("""
        CREATE TABLE IF NOT EXISTS blast_hits (
            qseqid      TEXT NOT NULL,
            sseqid      TEXT NOT NULL,
            pident      REAL NOT NULL,
            length      INTEGER NOT NULL,
            evalue      REAL NOT NULL,
            bitscore    REAL NOT NULL,
            qlen        INTEGER NOT NULL,
            slen        INTEGER NOT NULL,
            qcovs       REAL NOT NULL,
            tcovs       REAL NOT NULL,
            classified_by TEXT NOT NULL
        )
    """)

    rows = []
    for h in hits:
        dbs = db_seq_to_dbs.get(h["sseqid"], set())
        classified_by = ",".join(sorted(dbs)) if dbs else "unknown"
        rows.append((
            h["qseqid"], h["sseqid"], h["pident"], h["length"],
            h["evalue"], h["bitscore"], h["qlen"], h["slen"],
            h["qcovs"], h["tcovs"], classified_by,
        ))

    conn.executemany(
        "INSERT INTO blast_hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    conn.commit()
    log.info(f"  Stored {len(rows)} mmseqs hits")


def classify_from_blast(conn, classifications, database=None,
                        min_identity=80, min_coverage=80, min_length=80):
    """Classify unclassified sequences from stored pass-2 hits.

    Filter requires pident ≥ min_identity, qcovs ≥ min_coverage,
    tcovs ≥ min_coverage, length ≥ min_length. (min_coverage is applied to
    BOTH query and target coverage.) Best hit per query is selected by
    bitscore. Targets inherit Order + Superfamily; Clade is set to 'unknown'.
    """
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
    if "blast_hits" not in tables:
        return []

    if REQUIRE_BOTH_COVERAGE:
        where = ("WHERE pident >= ? AND qcovs >= ? AND tcovs >= ? "
                 "AND length >= ?")
        params = [min_identity, min_coverage, min_coverage, min_length]
    else:
        # At-least-one-side coverage (Wicker et al. 80-80-80 interpretation)
        where = ("WHERE pident >= ? AND (qcovs >= ? OR tcovs >= ?) "
                 "AND length >= ?")
        params = [min_identity, min_coverage, min_coverage, min_length]
    log.info(f"  coverage mode: {'BOTH' if REQUIRE_BOTH_COVERAGE else 'AT-LEAST-ONE'} "
             f"qcov/tcov >= {min_coverage}")

    if database:
        where += " AND classified_by LIKE ?"
        params.append(f"%{database}%")

    rows = conn.execute(f"""
        SELECT qseqid, sseqid, pident, qcovs, tcovs, length, bitscore
        FROM blast_hits
        {where}
        ORDER BY bitscore DESC
    """, params).fetchall()

    classified_set = set(classifications.keys())

    best = {}
    for qid, sid, pident, qcovs, tcovs, length, bitscore in rows:
        if qid in classified_set:
            continue
        if qid not in best:
            best[qid] = (sid, pident, qcovs, tcovs, length, bitscore)

    new_classifications = []
    no_source = 0
    for qid, (sid, pident, qcovs, tcovs, length, bitscore) in best.items():
        if sid in classifications:
            source = classifications[sid]
            new_classifications.append({
                "id": qid,
                "order": source["order"],
                "superfamily": source["superfamily"],
                "clade": "unknown",
                "complete": "none",
                "strand": "?",
                "domains": "none",
                "blast_source": sid,
                "blast_pident": pident,
                "blast_qcovs": qcovs,
                "blast_tcovs": tcovs,
                "blast_bitscore": bitscore,
            })
        else:
            no_source += 1

    if no_source:
        log.info(f"    {no_source} pass-2 hits to unclassified targets (skipped)")

    log.info(f"  pass-2: {len(new_classifications)} sequences classified "
             f"(from {len(best)} hits passing filters)")
    return new_classifications


def blast_pass2(input_fasta, conn, hmm_classifications=None,
                seq_type="nucl", n_processors=4,
                min_identity=80, min_coverage=80, min_length=80,
                outdir=None,
                pass2_classified_fasta=None,
                sensitivity=None,
                cov_mode=0):
    """mmseqs2-based pass-2 pipeline.

    Args:
        input_fasta: path to input FASTA
        conn: sqlite3 connection with HMM results
        hmm_classifications: dict of {seq_id: {order, superfamily, clade, ...}}
                             from classifier.classify_sequences(). Mutated in
                             place if pass2_classified_fasta is supplied.
        seq_type: "nucl" (mmseqs --search-type 3) or "prot" (--search-type 1)
        n_processors: mmseqs --threads
        min_identity/min_coverage/min_length: percent thresholds (0-100, 0-100, bp)
        outdir: output directory; a subdirectory `mmseqs_pass2/` is created inside
        pass2_classified_fasta: optional FASTA of prior classifications to
                                augment the pass-2 target pool
        sensitivity: optional mmseqs -s value
        cov_mode: mmseqs --cov-mode (0=query, 1=target, 2=shorter). Default 0.

    Returns:
        list of new classification dicts for sequences rescued by pass-2.
    """
    t0 = time.time()
    mmseqs.check_mmseqs()

    classified_ids = _get_classified_ids(conn)
    if not classified_ids:
        log.info("  No classified sequences for mmseqs pass-2")
        return []

    log.info(f"  mmseqs pass-2: {len(classified_ids)} classified sequences as targets")

    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="tebinsorter_mmseqs_")
    work = os.path.join(outdir, "mmseqs_pass2")

    db_fasta, qry_fasta, db_seq_to_dbs = split_classified_unclassified(
        input_fasta, classified_ids, work, seq_type=seq_type)

    if os.path.getsize(qry_fasta) == 0:
        log.info("  No unclassified sequences to search")
        return []

    if hmm_classifications is None:
        hmm_classifications = {}

    if pass2_classified_fasta:
        updated = pass2_external.update_classified_fasta_headers(
            pass2_classified_fasta, hmm_classifications, work
        )
        src = updated or pass2_classified_fasta
        pass2_external.extend_hmm_classifications_from_fasta(
            hmm_classifications, src, db_seq_to_dbs
        )
        merged_db = os.path.join(work, "pass2_db_merged.fa")
        pass2_external.merge_classified_fastas(
            merged_db, db_fasta, src, clean_nucl=(seq_type == "nucl")
        )
        db_fasta = merged_db

    if os.path.getsize(db_fasta) == 0:
        log.info("  pass-2 target FASTA is empty; skipping mmseqs")
        return []

    m8_out = os.path.join(work, "pass2.m8")
    tmpdir = os.path.join(work, "mmseqs_tmp")

    log.info(f"  Running mmseqs easy-search with {n_processors} threads")
    t1 = time.time()
    run_mmseqs(
        qry_fasta, db_fasta, m8_out, tmpdir,
        seqtype=seq_type, ncpu=n_processors,
        min_identity=min_identity, min_coverage=min_coverage,
        min_length=min_length,
        cov_mode=cov_mode, sensitivity=sensitivity,
    )
    t2 = time.time()
    log.info(f"  mmseqs search: {t2 - t1:.1f}s")

    hits = parse_mmseqs_output(m8_out)
    log.info(f"  {len(hits)} mmseqs best-hit records after split-alignment merge")

    if hits:
        store_blast_hits(conn, hits, db_seq_to_dbs)

    log.info(f"  {len(hmm_classifications)} HMM classifications available for inheritance")

    new_cls = classify_from_blast(
        conn, hmm_classifications,
        min_identity=min_identity,
        min_coverage=min_coverage,
        min_length=min_length,
    )

    t3 = time.time()
    log.info(f"  mmseqs pass-2 total: {t3 - t0:.1f}s")

    return new_cls
