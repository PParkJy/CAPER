# CAPER
"CAPER: Code-Aware Paired-End Read Merging in DNA Data Storage" was submitted to *IEEE Transactions on Computational Biology and Bioinformatics* (October 2026).

## Dataset
You can get the R1/R2 FASTQ files and checksum file from [CAPER-release-beta](https://github.com/PParkJy/CAPER/releases/tag/v1.0.0.-beta).  
Because of their size, the R1/R2 FASTQ files are compressed and split into two parts.  
Run the following in the `data/` folder:

```bash
cd data
cat paired_reads.tar.gz.part_00 paired_reads.tar.gz.part_01 > paired_reads.tar.gz
sha256sum -c paired_reads.tar.gz.sha256
tar -xzf paired_reads.tar.gz
```

1. Join the two parts into one `paired_reads.tar.gz`.
2. Verify the checksum. It should print `paired_reads.tar.gz: OK`.  
If it prints `FAILED`, the archive is corrupted; download the parts again.
3. Extract the archive. This produces the two FASTQ files:
   - `OP2-2_S1_L001_R1_RMother.fastq`: R1 reads
   - `OP2-2_S1_L001_R2_RMother.fastq`: R2 reads

## How to build
Requirements: g++ 7 or later (C++17). 

```bash
cd src
g++ -O2 -std=c++17 -pthread -o caper caper.cpp -lm
chmod +x sequence_analysis/*      
./caper --help                    # description of every parameter
```

Always run from `src/`. Input is read from `../data/` and results are written to `../result/`.

```
./caper <start_index> <end_index> <RSN> [alt_mode] [options]
```

| Argument | Meaning |
|---|---|
| `start_index`, `end_index` | Range of idx (trials) to run. The idx is the RNG seed of the read draw, so rerunning an idx draws exactly the same reads. |
| `RSN` | Number of reads (read_limit) drawn and decoded per idx. Each idx is decoded once at this read_limit; there is no retry with more reads. |
| `alt_mode` | `0` Baseline (no recovery), `1` EQ (flip equal-quality positions only), `2` DQ (flip any mismatch position) |

Main options:

| Option | Meaning |
|---|---|
| `--starcode 0/1` | Cluster-level recovery off (Naive) / on (NTL). |
| `--max-alt K --max-flips M --max-qsum Q` | Read-level recovery: K positions searched, at most M positions flipped per candidate, and an upper limit Q on the sum of the flipped positions' q_gap values. |
| `--cluster-tool starcode --star-dist D` | Starcode clustering with edit distance D. |
| `--cluster-tool cdhit --cdhit-id F` | CD-HIT clustering with identity F. |
| `--multilen L` | Reads of length 152 ± L also enter clustering. |
| `--ml-consensus weighted` | Use the quality-weighted consensus (for Proposed-NTL). |

### Commands per method
**Baseline-Naive**
```bash
./caper 1 100 80200 0 --merger pear --starcode 0
```

**Proposed-Naive (EQ)**
```bash
./caper 1 100 80200 1 --merger pear --starcode 0 \
    --max-alt 5 --max-flips 5 --max-qsum 20
```

**Proposed-Naive (DQ)**
```bash
./caper 1 100 80200 2 --merger pear --starcode 0 \
    --max-alt 5 --max-flips 5 --max-qsum 20
```

**Baseline-NTL**
```bash
./caper 1 100 80200 0 --merger pear --starcode 1 \
    --cluster-tool starcode --star-dist 7 --multilen 5
```

**Proposed-NTL (EQ)**
```bash
./caper 1 100 80200 1 --merger pear --starcode 1 \
    --cluster-tool starcode --star-dist 7 --multilen 5 --ml-consensus weighted \
    --max-alt 5 --max-flips 5 --max-qsum 20
```

**Proposed-NTL (DQ)**
```bash
./caper 1 100 80200 2 --merger pear --starcode 1 \
    --cluster-tool starcode --star-dist 7 --multilen 5 --ml-consensus weighted \
    --max-alt 5 --max-flips 5 --max-qsum 20
```

Notes:
- To use CD-HIT for NTL, replace `--cluster-tool starcode --star-dist 7` with `--cluster-tool cdhit --cdhit-id 0.94`.
- With `--multilen 0` (length-152 reads only), leave out `--ml-consensus`; otherwise the program stops with an error.
- For speed, add `--merge-jobs 5 --msa-jobs 3 --tool-threads 3 --tmp-dir /dev/shm`. The results do not change.
- To keep the console output, append `> log.txt 2>&1`.

### Mergers
Choose the merger with `--merger` and keep all other options as they are.

| Merger | Options | Settings |
|---|---|---|
| PEAR (default) | `--merger pear` | p-value 0.01 (`--pear-pvalue 0.01`) |
| NGmerge | `--merger ngmerge --ngmerge-mismatch 0.4` | mismatch fraction 0.4 |
| FLASH | `--merger flash --flash-max-overlap 150 --flash-mismatch 0.4` | max overlap 150, mismatch fraction 0.4 |

For example, Proposed-Naive (DQ) with FLASH:
```bash
./caper 1 100 80200 2 --merger flash \
    --flash-max-overlap 150 --flash-mismatch 0.4 --starcode 0 \
    --max-alt 5 --max-flips 5 --max-qsum 20
```

## Output
Each run writes two files to `../result/`, creating the folder if it does not exist.

`result_<start>_<end>_<suffix>.txt`

One line per idx.

| Column | Meaning |
|---|---|
| `idx` | trial number |
| `result` | `1` success (identical to the original), `0` decoding failure (not every chunk was recovered), `2` decoding finished but bit errors remain |
| `clusters` | strands handed to the LT decoder |
| `trial` | strands the decoder actually read |
| `LT_real_N` | strands stored with degree ≥ 2 |
| `observance` | strands that still held an unrecovered chunk when read |
| `inferred` | chunks recovered (16050 on success) |
| `base_error` | bit error rate against the original |
| `read_limit` | number of reads used (= `RSN`) |

`summary_<start>_<end>_<suffix>.txt`

A single line for the whole run.

| Column | Meaning |
|---|---|
| `read_limit` | number of reads used (= `RSN`) |
| `success_count` | number of idx decoded successfully (`result` = 1) |
| `total` | total number of idx |

### Filename suffix

The suffix is derived from the run settings, so runs with different settings never overwrite each other's results.

| Method | Suffix (PEAR) |
|---|---|
| Baseline-Naive | `org_p10` |
| Proposed-Naive (EQ) | `eqq_p10` |
| Proposed-Naive (DQ) | `allmm_p10` |
| Baseline-NTL | `org_rescue_d7_ml5_p10` |
| Proposed-NTL (EQ) | `eqq_rescue_d7_ml5w_p10` |
| Proposed-NTL (DQ) | `allmm_rescue_d7_ml5w_p10` |

What each part means:
- `org` / `eqq` / `allmm`: alt_mode 0 / 1 / 2
- `_rescue`: NTL (`--starcode 1`)
- `_d7`: Starcode distance
- `_cdhit940`: CD-HIT identity 0.94
- `_ml5` / `_ml5w`: `--multilen 5`; `w` means weighted consensus
- `_p10`: PEAR p-value 0.01
- `_ngmerge_p400`: NGmerge mismatch 0.4
- `_flash_M150_x400`: FLASH max overlap 150, mismatch 0.4
- `_k`, `_f`, `_qs`: added only when K, M or Q differ from the defaults (5, 5, 20)

## Contact
E-main: wldus8677@gmail.com  
Homepage: [Computing, Intelligence, and Communication Lab](https://cctl.jnu.ac.kr/)
