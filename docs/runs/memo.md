# memo — TEBinSorter_mmseqs port

> **Naming note (added 2026-07-28):** this codebase was renamed **TEsorter2** upstream
> (`KGerhardt/TESorter2`) after this run. Names, paths and commands below are left exactly
> as they were at the time so the run stays reproducible — read "TEBinSorter" as "TEsorter2".

## Date
2026-04-23

## Purpose

Port the user's `blastn → mmseqs2` pass-2 edits (plus sequence cleaning and
external-FASTA pass-2 augmentation) from their TEsorter fork onto the
much-faster TEBinSorter codebase. TEBinSorter is a ~500× faster rewrite
of the upstream TEsorter HMM pass-1; its pass-2 still used `blastn`. This
port swaps the pass-2 engine for `mmseqs easy-search` while leaving the
HMM speedup untouched. End result: full mmseqs pass-2 semantics on top
of the fast pass-1.

## Source-of-truth / provenance

- Upstream TEsorter: `github.com/zhangrengang/TEsorter` at tag `v1.5.1`
  (commit `6c14898`). Clone at
  `/data/chris/wheat/ltrharvest/v2/v3/testing/TEsorter/`.
- User's TEsorter fork: `github.com/cwb14/TEsorter`. Clone at
  `/data/chris/wheat/ltrharvest/v2/v3/testing/TEsorter_mmseq/`.
  Authoritative mmseqs edits live on **branch `origin/my-new-idea2`**
  (head commit `b398509`, 10 commits on top of v1.5.1):
  ```
  b398509 Update cov_mode to 0 for improved performance
  5410719 Enhance MmseqsM8Record and alignment processing
  d857bbb Update app.py
  625b425 Enhance comments on mmseq2 tie-breaking and repeat identification
  8eaa64c Implement header update for pass2 classified fasta
  faeaeb8 Add sequence cleaning functions and update merging
  4c25c96 Implement clean_fasta_atcg_only function
  60b5b0f Refactor input ID extraction from fasta file
  d16f7ca Enhance pass-2 classification with external FASTA support
  15673ba Add Mmseqs module and update app integration
  ```
  The weaker 1-commit branch `my-new-idea` is fully subsumed by
  `my-new-idea2`; the unrelated branch `backup/dumb-edits` is a separate
  minimap2 experiment and was ignored.
- Port target: TEBinSorter at `/data/chris/wheat/ltrharvest/v2/v3/testing/TEBinSorter/`
  (most recent commit `7144950`).
- Port destination: `/data/chris/wheat/ltrharvest/v2/v3/testing/TEBinSorter_mmseqs/`.

## Environment / versions

- Active conda env during development: `synLTR`
  (`/home/chris/bin/mambaforge/envs/synLTR`).
- `mmseqs` binary: version 17.b804f (in `synLTR/bin/mmseqs`). Must be on
  `$PATH` at runtime.
- Python deps (unchanged from TEBinSorter): pyhmmer ≥ 0.10, pyfastx ≥ 2.0,
  numpy. Install as the TEBinSorter README describes (`pip install
  pyhmmer pyfastx numpy`).

## Files touched

- NEW `src/mmseqs.py` — near-literal port of
  `TEsorter_mmseq@my-new-idea2:TEsorter/modules/Mmseqs.py`. One runtime
  difference: uses `subprocess.run` + `logging` instead of TEsorter's
  `RunCmdsMP.run_cmd` (to match TEBinSorter's idiom). `cov_mode` default
  is 0 in the signature (b398509's settled tuning), not 2.
- NEW `src/pass2_external.py` — port of the five `--pass2-classified-fasta`
  helpers from `TEsorter_mmseq@my-new-idea2:TEsorter/app.py`. Emits dicts
  shaped like TEBinSorter's reconciled classifications dict (rather than
  the upstream `CommonClassification` namedtuple). Uses `pyfastx` instead
  of Biopython `SeqIO`.
- EDITED `src/blast_pass2.py` — filename intentionally kept (so
  `pipeline.py`, `tesorter_compat.py`, `classifier.py`, and the SQLite
  `blast_hits` table name stay unchanged). `blastn`/`makeblastdb` internals
  gone; `multiprocessing.Pool` chunking gone. Single `mmseqs easy-search`
  call. mmseqs bits/fident/qcov are remapped into the BLAST-shaped
  columns at insert time; `evalue` and `slen` hold sentinel `0.0` / `0`
  (no downstream reader consumes them).
- EDITED `src/pipeline.py` — added 5 CLI args (`-dp2/--disable-pass2`,
  `-rule/--pass2-rule`, `--pass2-classified-fasta`, `--mmseqs-sensitivity`,
  `--mmseqs-cov-mode`). Pass-2 block splits `--pass2-rule` into floats,
  calls `blast_pass2()` with the expanded signature.
- EDITED `src/tesorter_compat.py` — same 5 CLI args added; its
  `blast_pass2()` call forwards them.
- EDITED `README.md` — added a mmseqs section at the top and documented
  the new CLI options. Original TEBinSorter README preserved below.

## Commands run

Clone + scaffold:
```
cp -a /data/chris/wheat/ltrharvest/v2/v3/testing/TEBinSorter \
      /data/chris/wheat/ltrharvest/v2/v3/testing/TEBinSorter_mmseqs
rm -rf /data/chris/wheat/ltrharvest/v2/v3/testing/TEBinSorter_mmseqs/.git
```

Source extraction:
```
cd /data/chris/wheat/ltrharvest/v2/v3/testing/TEsorter_mmseq
git show origin/my-new-idea2:TEsorter/modules/Mmseqs.py   # source of src/mmseqs.py
git diff 6c14898..origin/my-new-idea2 -- TEsorter/app.py   # source of src/pass2_external.py helpers
```

Smoke tests (with stubs for pyhmmer/pyfastx since this sandbox lacks both):
```
PYTHONPATH=/tmp/stubs:TEBinSorter_mmseqs/src python TEBinSorter_mmseqs/src/pipeline.py --help
PYTHONPATH=/tmp/stubs:TEBinSorter_mmseqs/src python TEBinSorter_mmseqs/src/tesorter_compat.py --help
# Plus a Python unit-check that exercises parse_cls_from_fasta_header,
# _COORD_HEADER_RE, _merge_split_alignments, split-gap gating, and
# best-hit selection — all passed.
```

## Expected outputs

- `src/blast_pass2.py` still produces a SQLite `blast_hits` table with
  columns `qseqid, sseqid, pident, length, evalue, bitscore, qlen, slen,
  qcovs, classified_by`. `pident` and `qcovs` are on 0–100 (remapped from
  mmseqs 0–1). `evalue` and `slen` are sentinel `0.0` / `0`.
- `classifier.py:568` still reads `blast_source` from the per-hit dict;
  unchanged.
- `.cls.tsv` output for pass-1 rows will be byte-identical to stock
  TEBinSorter's. Pass-2 rows will differ on borderline similarity hits
  (expected — mmseqs and blastn disagree on marginal alignments) but the
  vast majority should agree.

## Verification — RUN (on `synLTR` env after `pip install pyhmmer pyfastx`)

End-to-end run completed successfully on `TEsorter/example_data/rice6.9.5.liban`
(2,431 TE sequences, REXdb v4+metazoa, 8 processors).

| Run | Total wall | HMM pass-1 | Pass-2 wall | Pass-2 hits | Pass-2 classifications |
|---|---|---|---|---|---|
| stock TEBinSorter (blastn) | 5.2 s | 2.7 s | 0.8 s | 249 | 14 |
| TEBinSorter_mmseqs | 12.8 s | 2.7 s | 8.7 s | 385 (merged best-hit) | 11 |

**HMM pass-1 output is byte-identical** between the two builds (verified
via `diff rice6.9.5.liban.rexdb.cls.tsv`). The only differences are in
pass-2 rescue classifications, which disagree on 9 borderline sequences
out of ~1800 unclassified queries:

  - 8 shared
  - 3 classified only by mmseqs (`Os0185#DNAnona/MULEtir`, `Os1639#DNAnona/CACTA`, `Os2155#DNAnona/PILE`)
  - 6 classified only by blastn (`Os0066#MITE/Tourist`, `Os0625_mPing#MITE/Tourist`, `Os0657#DNAnona/MULE`, `Os1435#DNAnona/hAT`, `Os1997#DNAnona/CACTA`, `Os2754#DNAnona/MLE`)

This is the expected "mmseqs and blastn disagree on marginal alignments"
outcome. Spot check: `Os2754#DNAnona/MLE` — mmseqs found 70.4% query
coverage (below the 80% rule); blastn's cumulative qcovs via multiple
HSPs pushed it above 80%. This is a real tool difference, not a port bug.

### Two upstream bugs fixed during verification

**Bug 1 — search-time filter kills valid hits.** Upstream's
`classify_by_mmseqs` passes `--min-seq-id 0.8 -c 0.8 --min-aln-len 80`
to `mmseqs easy-search`. mmseqs's k-mer prefilter + `--cov-mode 0`
interacts destructively: hits that clearly pass the thresholds
post-alignment (e.g. `Os3552→Os0203` at fident=0.977, qcov=1.000) get
dropped before scoring. The port now runs mmseqs unfiltered at search
time and applies the 80-80-80 rule post-merge via the SQL `WHERE` in
`classify_from_blast`. This matches TEBinSorter's existing blastn path,
which also runs unfiltered and filters post-hoc. On rice data this is
the difference between "0 pass-2 classifications" and "11".

**Bug 2 — split-alignment merger breaks on reverse-strand queries.**
Upstream's `_chain_blocks` / `_merge_chain` assume `qstart <= qend`.
mmseqs nucleotide search reports reverse-strand hits with `qstart > qend`.
The raw arithmetic `hi - lo + 1` then goes negative, poisoning the merged
`qcov`. The port adds a `_qlohi` / `_tlohi` normalization so lo/hi are
always (min, max). On rice data this recovered 2 classifications (e.g.
`Os3552#DNAnona/CACTG`).

Both fixes are documented inline in the ported modules with comments
explaining why the port diverges from `my-new-idea2`.

### Reproducing the verification

```bash
cd /data/chris/wheat/ltrharvest/v2/v3/testing

# 1. baseline (blastn pass-2)
python TEBinSorter/src/pipeline.py \
    TEsorter/example_data/rice6.9.5.liban \
    -d rexdb -p 8 -o /tmp/rice.tebinsorter

# 2. mmseqs port
python TEBinSorter_mmseqs/src/pipeline.py \
    TEsorter/example_data/rice6.9.5.liban \
    -d rexdb -p 8 --pass2-rule 80-80-80 \
    -o /tmp/rice.tebinsorter_mmseqs

# 3. diff the classification TSVs
diff <(sort /tmp/rice.tebinsorter/*.cls.tsv) \
     <(sort /tmp/rice.tebinsorter_mmseqs/*.cls.tsv) | head -60
# Expect: pass-1 rows byte-identical; pass-2 rows differ only on
# borderline similarity hits.

# 4. pass-1-only sanity (mmseqs path bypassed):
python TEBinSorter_mmseqs/src/pipeline.py \
    TEsorter/example_data/rice6.9.5.liban -d rexdb --disable-pass2 \
    -o /tmp/rice.mmseqs_p1only
diff <(sort /tmp/rice.tebinsorter/*.cls.tsv | grep -v "^#") \
     <(sort /tmp/rice.mmseqs_p1only/*.cls.tsv | grep -v "^#")
# Expect: differences only where pass-2 added rows in run 1.

# 5. missing-binary sanity:
PATH=/usr/bin python TEBinSorter_mmseqs/src/pipeline.py \
    TEsorter/example_data/rice6.9.5.liban -d rexdb
# Expect: clean RuntimeError "'mmseqs' not found on PATH."

# 6. SQLite schema spot-check:
sqlite3 /tmp/rice.tebinsorter_mmseqs/*.db \
    "SELECT qseqid, sseqid, pident, qcovs, bitscore FROM blast_hits LIMIT 3;"
# Expect: pident/qcovs in 0–100, bitscore populated.
```

## Notes / caveats

- The mmseq-branch upstream has typos `min_identtity` / `min_coverge`
  in `app.py:classify_by_mmseqs`. These are fixed to `min_identity` /
  `min_coverage` in the port (matching TEBinSorter's existing signature).
- `_MAX_SPLIT_GAP = 500` is hard-coded per upstream; not CLI-exposed.
- Sentinel `evalue=0.0` / `slen=0` in `blast_hits` is a deliberate choice
  to avoid a schema migration. No downstream code reads these columns
  (verified via grep across `src/`).
- Upstream `my-new-idea2` comment `b398509` confirms cov_mode=0 is the
  tuned choice; the port bakes this in as the signature default.
- `--pass2-classified-fasta` headers must be `>id#Order/Superfamily/Clade`;
  coordinate-shaped IDs `chr:start-end` additionally get their `unknown`
  fields upgraded from pass-1 results via `_COORD_HEADER_RE`.

## Reproducibility

```
claude --resume "port TEsorter-mmseqs edits onto TEBinSorter"
```

---

# 2026-04-24 update — dual-coverage filter + rape-library benchmarks

## What changed in the code

### 1. Target-side coverage is now captured and filterable

Previously the port only wrote `qcov` into mmseqs's `--format-output` string,
so the SQL pass-2 filter was query-only (`WHERE qcovs >= C`). Users pointed
out that the Wicker et al. 80-80-80 rule is meant to check coverage of the
candidate against an element, and that requires knowing target-side coverage
too.

Code changes:

- `src/mmseqs.py`
  - `--format-output` now includes `tlen,tcov` at the end:
    `"query,target,fident,alnlen,qlen,qcov,bits,qstart,qend,tstart,tend,tlen,tcov"`
  - `MmseqsM8Record.__slots__` adds `tlen`, `tcov`.
  - `_merge_chain()` unions target-side intervals (independently sorted,
    since query-sorted order can leave target intervals out of sequence for
    reverse-strand hits) and computes `tcov = merged_t_alnlen / tlen`.

- `src/blast_pass2.py`
  - `parse_mmseqs_output()` sets `tcovs = rec.tcov * 100` on the hit dict and
    `slen = rec.tlen` (no longer a sentinel — useful for downstream sanity).
  - `store_blast_hits()` schema now has `tcovs REAL NOT NULL`.
  - `classify_from_blast()` SQL filter includes tcov.

### 2. Developer toggle: `REQUIRE_BOTH_COVERAGE`

A module-level constant at the top of `src/blast_pass2.py`:

```python
# True  -> require BOTH qcov AND tcov >= coverage threshold (strict).
# False -> require AT LEAST ONE of qcov, tcov >= coverage threshold
#          (Wicker et al. 80-80-80: "candidate must cover ≥80% of at least
#          one of the elements being compared").
REQUIRE_BOTH_COVERAGE = False
```

- Not exposed on the CLI — edit the source to flip it.
- Defaults to `False` (OR) after reviewing Wicker 2007's original language.
- mmseqs has `--cov-mode 0,1,2` natively (both / target / query) but no
  "at-least-one" mode, so the OR semantics are implemented post-hoc in SQL:
  `WHERE pident >= ? AND (qcovs >= ? OR tcovs >= ?) AND length >= ?`
- `classify_from_blast` logs which mode is active at the start of pass-2:
  `coverage mode: AT-LEAST-ONE qcov/tcov >= 80`

### 3. Architectural limitation noted in the code

Spot-check during benchmarking: mmseqs nucleotide (`--search-type 3`) emits
**one best diagonal per (query, target) pair**, regardless of sensitivity,
zdrop, alt-ali, or prefilter mode. When two real homologous regions exist
between a query and target, mmseqs reports only the highest-scoring one.

This means:
- qcov and tcov for a given (q,t) are both derived from the same single
  alignment. They're coupled through the sequences' lengths: for a 10 kb
  query and 5 kb target with a 5 kb single alignment, qcov=50%, tcov=100%.
- Under OR mode, this collapses toward `max(qcov, tcov) >= C`, which on
  most pass-2 workloads behaves indistinguishably from qcov-only.
- A sibling port (`TEBinSorter_minimap2/`) was built because minimap2
  *does* emit multiple chains per (q,t) and can union across them for a
  more informative coverage calculation.

Confirmed empirically: `--diag-score 0`, `--alt-ali 3`, `-s 7.5`, `-k 11`,
`--exhaustive-search 1`, `--alignment-mode 2`, target-FASTA swap — none
recover both regions of the test pair Q=NW_026014909.1:83178-91907 vs
T=NW_026014943.1:54536-70563 on our rape data.

## Rape-library benchmarks (Brassica napus, LTR_retriever intact LTR-RTs)

Input: 180,556-sequence full library, subsampled to 5,000 with
`seqkit sample -s 42 -n 5000` → 4,892 unique sequences after dedup.
Environment: synLTR conda env; 16 threads.

### 5k subset, `--pass2-rule 80-80-80`, AT-LEAST-ONE (OR) mode — default

| Tool | Wall | Max RSS | Pass-2 rescues | Combined |
|---|---:|---:|---:|---:|
| blastn (stock TEBinSorter reference) | 17:37 | 76.6 GB | 1,005 | 2,522 |
| **mmseqs port** | **1:19** | **8.6 GB** | **397** | **1,914** |
| minimap2 port (reference)            | 0:45  | 6.4 GB  | 598   | 2,115 |

### 5k subset, `--pass2-rule 70-50-80`, AT-LEAST-ONE (OR) mode

| Tool | Wall | Pass-2 rescues | Combined |
|---|---:|---:|---:|
| blastn                    | 18:38 | 1,145 (re-filtered from 80-80-80 run) | — |
| **mmseqs port**           | **1:07** | **773** | **2,290** |
| minimap2 port (reference) | 0:46  | 939   | 2,456 |

### Archived, for reference: BOTH-coverage (AND) mode on same data

| Rule | mmseqs rescues | minimap2 rescues |
|---|---:|---:|
| 80-80-80 AND | 0   | 16  |
| 70-50-80 AND | 33  | 152 |

AND mode is pathologically strict on LTR-RT data because most candidate/element
pairs are length-asymmetric; OR is both biologically correct per Wicker and
recovers ~25× more rescues.

## Observed behavior in the mmseqs port specifically

- **mmseqs OR-80-80-80 (397 rescues) ≡ mmseqs qcov-only @ 80-80-80 (pre-update, 397).**
  Not coincidence: OR with one-diagonal search makes the effective rule
  `max(qcov, tcov) >= C`. For pairs where target is longer than query,
  `max` is usually qcov; for shorter-target pairs, it's tcov. Either way,
  the OR rule adds little to qcov-only on single-diagonal evidence.
- Therefore **the mmseqs port's practical advantage over blastn is wall time
  + memory, not biological precision**. Recall is capped by mmseqs's one-
  diagonal architecture.

## Notes on the new schema

The `blast_hits` table in this port now has an extra `tcovs` column (11
total columns, up from 10). Any external code that writes directly to this
DB will need updating. Downstream TEBinSorter code only reads via
`classify_from_blast`, which was updated here — no external change needed.

## What is NOT changed

- The four bug fixes from the 2026-04-23 port (reverse-strand `_qlohi`, etc)
  remain in place.
- Default `--mmseqs-cov-mode 0` on the CLI is unchanged; with our `-c 0.0`
  search-side threshold, cov-mode is effectively a no-op anyway (its filter
  never triggers). Kept for compatibility with the sibling minimap2 port's
  API shape.
- The `--pass2-classified-fasta` external-augmentation pipeline is unchanged.
- HMM pass-1 is untouched — verified byte-identical to stock TEBinSorter on
  both rice and rape.

## Bottom line

This port works, is ~13× faster and 9× lower memory than stock TEBinSorter's
blastn pass-2 at the same thresholds, and the new dual-coverage schema makes
its semantics explicit. On LTR-RT data specifically, the sibling minimap2
port consistently recovers more pass-2 rescues because minimap2 can produce
multiple chains per (query, target) pair and genuinely differentiate qcov
from tcov. If you need maximum recall while enforcing Wicker's 80-80-80
rule, prefer minimap2. If your use case is protein-DB-scale classification
or if throughput matters more than coverage nuance, mmseqs remains a sound
choice.

