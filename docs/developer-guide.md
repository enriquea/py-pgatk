# Developer & Maintainer Guide

This guide is intended for developers who are new to the **pgatk** codebase and
need a quick yet thorough orientation.  It answers four questions:

1. What does the tool do?
2. Which improvements are feasible short-term?
3. Which extensions are worth focusing on in the short/mid-term?
4. Should we invest in refining the use cases with real data and results?

---

## 1. What the Tool Does

### One-sentence summary

**pgatk** (ProteoGenomics Analysis ToolKit) builds custom protein-sequence
databases from genomic and transcriptomic data so that mass-spectrometry
search engines can identify novel, variant, and non-canonical peptides.

### Why it matters

Standard proteomics searches use *reference* protein databases (e.g. UniProt).
Any peptide that differs from the reference — because of a somatic mutation, a
rare germline variant, a non-coding RNA translation, or an alternative reading
frame — will be missed.  **pgatk** bridges this gap by:

| Step | What pgatk does |
|------|-----------------|
| **Download** | Fetches transcript sequences, gene annotations, and variant calls from ENSEMBL, ClinVar, COSMIC, cBioPortal, or NCBI. |
| **Translate** | Applies genomic variants (SNPs, indels, frameshifts) to transcript sequences and translates the altered DNA into proteins. |
| **Discover ORFs** | Performs 3-frame or 6-frame translation of non-coding genes (lncRNAs, pseudogenes, antisense transcripts) to find novel open reading frames and micropeptides. |
| **Build databases** | Produces FASTA files compatible with all major search engines (MaxQuant, MSFragger, SearchGUI, Comet, DIA-NN, Proteome Discoverer). |
| **Generate decoys** | Creates target-decoy databases (4 methods) for false-discovery-rate estimation. |
| **Post-process** | Maps identified peptides back to genomic coordinates (GFF3), performs *in silico* digestion, and provides spectrum-level validation tools. |

### Architecture at a glance

```
pgatk/
├── cli.py                   # Click entry point  →  `pgatk <command>`
├── commands/                # 16 CLI sub-commands (thin wrappers)
├── ensembl/                 # ENSEMBL data download & VCF → protein translation
├── clinvar/                 # ClinVar variant processing (no VEP required)
├── cgenomes/                # COSMIC & cBioPortal cancer-genomics integration
├── proteomics/db/           # Decoy-database generation
├── proteogenomics/          # SpectrumAI (MS2 validation), BLAST search
├── db/                      # Peptide-to-genome mapping, in-silico digestion
├── toolbox/                 # Shared utilities: VCF parsing, REST calls, config
├── config/                  # YAML defaults + registry pattern
└── tests/                   # unittest + pytest suite with bundled test data
```

The data-flow pattern is the same for every command:

```
CLI option parsing (Click)
  → Configuration merging (YAML defaults + user overrides)
    → Service class instantiation
      → Data processing (download / translate / generate)
        → FASTA or GFF3 output
```

### CLI commands (grouped)

| Group | Commands |
|-------|----------|
| **Data downloaders** | `ensembl-downloader`, `ncbi-downloader`, `cosmic-downloader`, `cbioportal-downloader` |
| **Variant → protein** | `vcf-to-proteindb`, `clinvar-to-proteindb`, `cosmic-to-proteindb`, `cbioportal-to-proteindb` |
| **Sequence translation** | `dnaseq-to-proteindb`, `threeframe-translation` |
| **Database processing** | `generate-decoy`, `ensembl-check` |
| **Post-processing** | `digest-mutant-protein`, `map-peptide2genome`, `spectrumai`, `blast_get_position` |

### Key dependencies

| Library | Role |
|---------|------|
| Biopython | Biological sequence handling & translation |
| Click | CLI framework |
| gffutils | GTF/GFF annotation parsing |
| pandas | Tabular data manipulation (VCF processing) |
| pyteomics | Proteomics utilities (digestion, decoy generation) |
| pybedtools | BED/interval operations (requires system `bedtools`) |
| pyopenms | OpenMS bindings for mass-spectrometry data |
| pyahocorasick | Fast multi-pattern string matching |

External tools required: `bedtools`, `htslib` (tabix/bgzip).

---

## 2. Short-Term Feasible Improvements

These can be tackled independently in 1–2 week sprints with no architectural
changes.

### 2.1 Fix the variable-shadowing bug in decoy generation

In `pgatk/proteomics/db/protein_database_decoy.py`, a local variable named
`decoy_sequence` shadows the imported `pyteomics.fasta.decoy_sequence` function.
This is a latent correctness risk and should be renamed (e.g. to
`decoy_seq_record`).

### 2.2 Replace `print()` with structured logging

Several service classes mix `print()` statements and `logging.getLogger()`
calls.  Standardising on `logger.info()` / `logger.warning()` improves
debuggability and allows users to control verbosity via `--verbose` / `--quiet`
flags.

### 2.3 Modernise Python idioms

- **`pathlib`** instead of `os.path.join` — safer path handling.
- **`dataclasses`** or `NamedTuple` for configuration and model classes —
  less boilerplate, better IDE support.
- **Type hints** on public APIs — enables static analysis (mypy) and improves
  contributor onboarding.

### 2.4 Improve test coverage and CI

- The CI workflow (`pythonapp.yml`) only discovers files matching `*_tests.py`,
  but newer test files use the `test_*.py` convention (e.g. `test_ensembl_core.py`,
  `test_vcf_utils.py`, the entire `test_clinvar/` directory).  Updating the
  discovery pattern to `test*.py` will immediately increase coverage without
  writing new tests.
- Add `pytest` markers or `unittest.skip` for tests that require network
  access, so the offline suite runs faster.

### 2.5 Tighten exception handling

Replace broad `except Exception` blocks with specific exception types.
Unexpected errors should propagate rather than being silently swallowed, making
debugging far easier for users and maintainers.

### 2.6 Eliminate the `sys.argv` import-time side-effect

The scripts in `pgatk/db/` (e.g. `map_peptide2genome.py`,
`digest_mutant_protein.py`) parse `sys.argv` at module import time, which makes
them non-importable as library modules and breaks test isolation.  Convert them
to proper Click sub-commands consistent with the rest of the CLI.

---

## 3. Short/Mid-Term Extensions

These require more design work but deliver high user impact.  Detailed plans
already exist in `docs/plans/`.

### 3.1 Unified FASTA header format *(short-term, 2–4 weeks)*

**Design doc:** `docs/plans/2026-03-03-protein-accession-design.md`

Current FASTA headers differ across variant sources and cause parsing failures
in search engines like SearchGUI.  The approved `pgvar|` prefix scheme unifies
all variant proteins under a consistent, search-engine-compatible header.  This
is a prerequisite for any downstream tool integration and should be prioritised.

### 3.2 Graph-based transcript modelling *(mid-term, phased over months)*

**Design doc:** `docs/plans/2026-03-01-pgatk-graph-engine-design.md`

The current engine processes one variant at a time.  For a transcript with *N*
variants, this produces *N* single-variant sequences instead of the
combinatorially correct set.  A directed acyclic graph (DAG) per transcript
solves this and is how the competing tool *moPepGen* (Nature Biotechnology 2025)
works.  The design document defines a phased rollout:

| Phase | Scope | Adds |
|-------|-------|------|
| 0 | Infrastructure cleanup | Bug fixes, pathlib, type hints, dataclasses |
| 1 | Graph core | TranscriptGraph, SNP/indel co-occurrence |
| 2 | Cancer & clinical | ClinVar parser, graph-aware COSMIC/cBioPortal |
| 3 | Advanced events | Gene fusions, RNA editing, circular RNA |
| 4 | Performance | Profiling, optional Cython hot-paths |

Phase 0 overlaps with section 2 above and can begin immediately.

### 3.3 Performance optimisation of VCF processing *(short/mid-term)*

The core function `vcf_to_proteindb()` in `pgatk/ensembl/ensembl.py` iterates
over VCF rows with `pandas.iterrows()` and performs per-row gffutils database
lookups.  Two focused changes can yield large speedups:

1. **Vectorise** the pandas operations (group-by transcript, then process).
2. **Cache transcript features** from gffutils to avoid redundant queries.

These are independent of the graph engine and can be done now.

### 3.4 Support for additional variant sources

- **gnomAD v4** structural variants (SVs).
- **PharmGKB** pharmacogenomics variants.
- **dbNSFP** pre-computed functional annotations.

Each new source follows the established pattern: a downloader module + a
`*_to_proteindb` command + YAML config + tests.

---

## 4. Should We Invest in Refining Use Cases with Real Data?

**Yes — this is one of the highest-impact investments we can make.**

### Why

1. **Credibility & adoption.**  The existing `docs/use-cases.md` has 13
   end-to-end workflows, but they are *recipes* (command sequences) rather than
   *validated results*.  Users cannot tell whether a workflow produces
   biologically meaningful output until they run it themselves.  Adding concrete
   numbers (database size, peptide counts, FDR curves) turns recipes into
   evidence.

2. **Regression testing.**  Real-data runs serve as integration tests.  If a
   code change silently alters the output FASTA (e.g. different number of
   proteins or changed headers), a recorded baseline catches it immediately.

3. **Benchmarking against competitors.**  The field now has *moPepGen*,
   *vcf2prot*, *PG2*, and others.  Publishing side-by-side comparisons (database
   size, unique peptide yield, runtime) on the same datasets positions pgatk in
   the landscape and highlights where the graph engine (section 3.2) will close
   the gap.

4. **Reproducibility & publication.**  The original pgatk paper (Umer et al.
   2022) demonstrated 43,501 non-canonical peptides across 64 cell lines.
   Reproducing and extending this benchmark with current code validates that the
   tool still delivers and provides material for a methods-update publication.

### Recommended approach

| Action | Effort | Outcome |
|--------|--------|---------|
| Pick 2–3 public datasets (e.g. PRIDE PXD identifiers) and document exact commands + expected outputs | 1–2 weeks | Reproducible "gold standard" workflows |
| Record summary statistics (# proteins, # unique variant peptides, database size, runtime) per workflow | Minimal | Quantitative baselines for regression testing |
| Add a CI job that runs a small subset of these workflows on every PR | 1 week | Automated integration testing |
| Write a comparison table vs. moPepGen / vcf2prot on the same dataset | 2–3 weeks | Competitive positioning for the README and documentation |

### Suggested datasets

- **PRIDE PXD004452** — used in the original Umer et al. 2022 benchmark (64
  human cell lines).
- **CPTAC** — well-curated cancer proteogenomics cohorts with matched WGS/WES.
- **ENSEMBL GRCh38 + gnomAD v3** — population-scale variant database for
  germline-variant benchmarking.

---

## Quick-Start for New Maintainers

### Set up a development environment

```bash
# Clone and install in editable mode
git clone https://github.com/bigbio/pgatk.git
cd pgatk
pip install -e ".[dev]"

# System dependencies (on Ubuntu/Debian)
sudo apt-get install bedtools tabix

# Or via Conda
conda install -c conda-forge -c bioconda bedtools htslib pyopenms
```

### Run the tests

```bash
cd pgatk
# Current CI pattern (misses some test files — see section 2.4)
python -m unittest discover -s tests -p "*_tests.py" -v

# To also run test_*.py files
python -m pytest tests/ -v
```

### Lint

```bash
# Critical errors only (CI-blocking)
flake8 . --select=E9,F63,F7,F82 --show-source

# Full style check (non-blocking)
flake8 . --max-complexity=10 --max-line-length=127
```

### Build documentation locally

```bash
pip install mkdocs-material pymdown-extensions
mkdocs serve   # http://127.0.0.1:8000
```

---

## Further Reading

- **CLI reference:** [pgatk-cli.md](pgatk-cli.md)
- **Use-case workflows:** [use-cases.md](use-cases.md)
- **File formats:** [formats.md](formats.md)
- **Graph engine design:** [plans/2026-03-01-pgatk-graph-engine-design.md](plans/2026-03-01-pgatk-graph-engine-design.md)
- **Accession design:** [plans/2026-03-03-protein-accession-design.md](plans/2026-03-03-protein-accession-design.md)
- **Published paper:** Umer et al., *Bioinformatics* 2022, 38(5), 1470–1472.
