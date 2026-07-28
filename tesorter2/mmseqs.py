"""
mmseqs.py — MMseqs2 wrapper for pass-2 similarity search.

Ported from github.com/cwb14/TEsorter branch `my-new-idea2`
(TEsorter/modules/Mmseqs.py, head commit b398509). Public API and
split-alignment merging are kept byte-for-byte compatible with upstream to keep
future cherry-picks easy; the only substitutions are:

  * `subprocess.run` + `logging.getLogger` replace TEsorter's `RunCmdsMP.run_cmd`
  * default `cov_mode=0` (query coverage) instead of 2 — matches b398509's
    "Update cov_mode to 0 for improved performance" which overrode the upstream
    default at every caller
"""

import logging
import os
import shutil
import subprocess
from collections import OrderedDict, defaultdict

log = logging.getLogger(__name__)


def check_mmseqs(mmseqs_bin="mmseqs"):
    """Fail fast with a clear message if mmseqs is not on PATH."""
    if shutil.which(mmseqs_bin) is None:
        raise RuntimeError(
            f"{mmseqs_bin!r} not found on PATH. Install mmseqs2 "
            f"(conda: `mamba install -c bioconda mmseqs2`) and retry."
        )


def mmseqs_version(mmseqs_bin="mmseqs"):
    """Log the mmseqs version for reproducibility. Non-fatal if missing."""
    if shutil.which(mmseqs_bin) is None:
        log.warning(f"{mmseqs_bin!r} not found on PATH")
        return None
    result = subprocess.run(
        [mmseqs_bin, "version"], capture_output=True, text=True, check=False
    )
    v = (result.stdout or result.stderr).strip().splitlines()[0] if (
        result.stdout or result.stderr
    ) else "unknown"
    log.info(f"mmseqs version: {v}")
    return v


def mmseqs_easy_search(db_seq, qry_seq, out_m8, tmpdir,
                       seqtype="nucl", ncpu=4,
                       min_seq_id=0.0, min_cov=0.0, cov_mode=0, min_aln_len=0,
                       sensitivity=None,
                       mmseqs_bin="mmseqs"):
    """
    Run `mmseqs easy-search` with a custom output format that preserves the
    alignment coordinates needed for split-alignment merging.

    The query FASTA is pre-cleaned to ATCG-only (written to
    {tmpdir}/query.cleaned.fa) before invocation, mirroring upstream.

    seqtype:
        "nucl" -> --search-type 3
        "prot" -> --search-type 1
    """
    os.makedirs(tmpdir, exist_ok=True)

    cleaned_qry = os.path.join(tmpdir, "query.cleaned.fa")
    clean_fasta_atcg_only(qry_seq, cleaned_qry)
    qry_seq = cleaned_qry

    if seqtype == "nucl":
        search_type = 3
    elif seqtype == "prot":
        search_type = 1
    else:
        raise ValueError(f"Unknown seqtype {seqtype}")

    # qstart/qend/tstart/tend are needed for the split-alignment merge step.
    # tlen + tcov added so post-hoc filter can enforce target-side coverage
    # (cov-mode 0 + -c would enforce both sides in-search, but we keep -c at
    # 0 because mmseqs's k-mer prefilter drops valid hits at higher thresholds).
    fmt = ("query,target,fident,alnlen,qlen,qcov,bits,"
           "qstart,qend,tstart,tend,tlen,tcov")

    cmd = (
        f"{mmseqs_bin} easy-search "
        f"{qry_seq} {db_seq} {out_m8} {tmpdir} "
        f"--threads {ncpu} "
        f"--search-type {search_type} "
        f"--format-output \"{fmt}\" "
        f"--min-seq-id {min_seq_id} "
        f"-c {min_cov} --cov-mode {cov_mode} "
        f"--min-aln-len {min_aln_len} "
    )
    if sensitivity is not None:
        cmd += f"-s {float(sensitivity)} "

    log.info(f"mmseqs cmd: {cmd}")
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"mmseqs easy-search failed (exit {result.returncode}): "
            f"{(result.stderr or result.stdout)[:2000]}"
        )
    return out_m8


def clean_fasta_atcg_only(in_fa, out_fa):
    """
    Strip non-ATCG characters from every sequence in a FASTA; uppercase.
    Line-by-line so it scales to whole-genome inputs. Headers are preserved.
    """
    with open(in_fa) as fin, open(out_fa, "w") as fout:
        seq_buf = []
        header = [None]

        def flush():
            if header[0] is None:
                return
            seq = "".join(seq_buf).upper()
            seq = "".join(b for b in seq if b in ("A", "T", "C", "G"))
            fout.write(header[0] + "\n")
            fout.write(seq + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header[0] = line
                seq_buf.clear()
            else:
                seq_buf.append(line)
        flush()


class MmseqsM8Record:
    __slots__ = (
        "qseqid", "sseqid", "fident", "alnlen", "qlen", "qcov", "bits",
        "qstart", "qend", "tstart", "tend", "tlen", "tcov",
    )

    def __init__(self, line):
        # format: query,target,fident,alnlen,qlen,qcov,bits,
        #         qstart,qend,tstart,tend,tlen,tcov
        vals = line.rstrip("\n").split("\t")
        self.qseqid = vals[0]
        self.sseqid = vals[1]
        self.fident = float(vals[2])   # 0..1
        self.alnlen = int(vals[3])
        self.qlen = int(vals[4])
        self.qcov = float(vals[5])     # 0..1
        self.bits = float(vals[6])
        self.qstart = int(vals[7])     # 1-based
        self.qend = int(vals[8])       # 1-based
        self.tstart = int(vals[9])     # 1-based
        self.tend = int(vals[10])      # 1-based
        self.tlen = int(vals[11])
        self.tcov = float(vals[12])    # 0..1


# Maximum gap (bp) on both query and target sides within which two blocks are
# considered part of the same MMseqs2-split alignment rather than two
# genuinely independent local alignments. Hard-coded upstream; kept here too.
_MAX_SPLIT_GAP = 500


def _qlohi(r):
    """Normalize (qstart, qend) as (lo, hi). MMseqs reports reverse-strand
    nucleotide hits with qstart>qend; the upstream chain/merge code treats
    coords as 1D intervals, so we need lo<=hi."""
    return (r.qstart, r.qend) if r.qstart <= r.qend else (r.qend, r.qstart)


def _tlohi(r):
    return (r.tstart, r.tend) if r.tstart <= r.tend else (r.tend, r.tstart)


def _chain_blocks(records):
    """
    Partition records (all sharing the same query-target pair) into contiguous
    chains. Two consecutive blocks (sorted by query-interval lo) belong to the
    same chain when the gap between them is <= _MAX_SPLIT_GAP on *both* the
    query and the target. A negative gap (overlap) always satisfies the
    criterion. Strand is normalized via _qlohi / _tlohi so reverse-strand hits
    don't produce negative gaps.
    """
    blocks = sorted(records, key=lambda r: _qlohi(r))
    chains = [[blocks[0]]]
    for r in blocks[1:]:
        prev = chains[-1][-1]
        _, prev_qhi = _qlohi(prev)
        r_qlo, _ = _qlohi(r)
        _, prev_thi = _tlohi(prev)
        r_tlo, _ = _tlohi(r)
        qgap = r_qlo - prev_qhi - 1
        tgap = r_tlo - prev_thi - 1
        if qgap <= _MAX_SPLIT_GAP and tgap <= _MAX_SPLIT_GAP:
            chains[-1].append(r)
        else:
            chains.append([r])
    return chains


def _merge_chain(records):
    """
    Merge records forming a single contiguous chain:
      - union the aligned query intervals  -> qcov
      - union the aligned target intervals -> tcov
      - sum bit-scores
      - length-weighted average fident
    """
    if len(records) == 1:
        return records[0]

    blocks = sorted(records, key=lambda r: _qlohi(r))

    merged_q_intervals = []
    for r in blocks:
        lo, hi = _qlohi(r)
        if merged_q_intervals and lo <= merged_q_intervals[-1][1]:
            prev_lo, prev_hi = merged_q_intervals[-1]
            merged_q_intervals[-1] = (prev_lo, max(prev_hi, hi))
        else:
            merged_q_intervals.append((lo, hi))
    merged_q_alnlen = sum(hi - lo + 1 for lo, hi in merged_q_intervals)

    # Target intervals need independent sorting since query-sorted order
    # can leave target intervals out of order (especially reverse-strand).
    t_sorted = sorted([_tlohi(r) for r in blocks])
    merged_t_intervals = []
    for lo, hi in t_sorted:
        if merged_t_intervals and lo <= merged_t_intervals[-1][1]:
            prev_lo, prev_hi = merged_t_intervals[-1]
            merged_t_intervals[-1] = (prev_lo, max(prev_hi, hi))
        else:
            merged_t_intervals.append((lo, hi))
    merged_t_alnlen = sum(hi - lo + 1 for lo, hi in merged_t_intervals)

    total_weight = sum(r.alnlen for r in blocks)
    weighted_fident = (
        sum(r.fident * r.alnlen for r in blocks) / total_weight
        if total_weight else 0.0
    )
    merged_bits = sum(r.bits for r in blocks)

    ref = blocks[0]
    synthetic = MmseqsM8Record.__new__(MmseqsM8Record)
    synthetic.qseqid = ref.qseqid
    synthetic.sseqid = ref.sseqid
    synthetic.fident = weighted_fident
    synthetic.alnlen = merged_q_alnlen
    synthetic.qlen = ref.qlen
    synthetic.qcov = merged_q_alnlen / ref.qlen if ref.qlen else 0.0
    synthetic.bits = merged_bits
    synthetic.qstart = merged_q_intervals[0][0]
    synthetic.qend = merged_q_intervals[-1][1]
    synthetic.tstart = blocks[0].tstart
    synthetic.tend = blocks[-1].tend
    synthetic.tlen = ref.tlen
    synthetic.tcov = merged_t_alnlen / ref.tlen if ref.tlen else 0.0
    return synthetic


def _merge_split_alignments(records):
    """Group by (query, target) -> chain -> merge each chain -> flat list."""
    groups = defaultdict(list)
    for r in records:
        groups[(r.qseqid, r.sseqid)].append(r)

    merged = []
    for grp in groups.values():
        for chain in _chain_blocks(grp):
            merged.append(_merge_chain(chain))
    return merged


def parse_mmseqs_m8_besthit(m8_path):
    """
    Parse an MMseqs2 M8 (in the custom format written by mmseqs_easy_search),
    merge split alignments per query-target pair, then keep the best hit per
    query by bit-score (post-merge).
    """
    all_records = []
    if not os.path.exists(m8_path):
        return OrderedDict()
    with open(m8_path) as f:
        for line in f:
            if not line.strip():
                continue
            all_records.append(MmseqsM8Record(line))

    if not all_records:
        return OrderedDict()

    merged = _merge_split_alignments(all_records)

    best = OrderedDict()
    for rc in merged:
        if rc.qseqid not in best or rc.bits > best[rc.qseqid].bits:
            best[rc.qseqid] = rc
    return best
