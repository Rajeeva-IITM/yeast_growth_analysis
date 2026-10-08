# FINEMAP + CAVIAR fine-mapping: setup log

Detailed record of everything done to add FINEMAP and CAVIAR fine-mapping
alongside the existing SuSiE analysis. Date: **2026-07-28**.

**Goal:** a script equivalent to `src/SuSiE analysis/susie_conditionwise_analysis.R`
but driving FINEMAP and CAVIAR via the [`finemapr`](https://github.com/variani/finemapr)
package, with all third-party binaries installed into a gitignored folder and
the rest of the project left untouched.

**Status:** setup complete and verified. The analysis itself has *not* been run
(deliberately — only the steps preceding it were requested).

---

## 1. Environment

| Item | Value |
|---|---|
| Project root | `H:\Rajeeva\Project\Himanshu\subsets\R_analysis-yeast_growth_pred` |
| OS | Windows 11 Pro (26200) |
| R | 4.5.1 (`C:\Program Files\R\R-4.5.1\bin\Rscript.exe`), renv project library |
| WSL | Ubuntu 20.04, kernel 5.15.133.1-microsoft-standard-WSL2, gcc/g++ 9.4.0 |
| Project path inside WSL | `/mnt/h/Rajeeva/Project/Himanshu/subsets/R_analysis-yeast_growth_pred` |

**Why WSL is involved.** FINEMAP ships Linux/macOS builds only and CAVIAR is a
C++ program expecting a POSIX toolchain; neither runs natively on Windows. R,
however, is the native Windows build with the project's renv library. So the R
side stays on Windows and only the two binaries live in WSL, bridged by shims
(section 4). WSL was already installed and the `H:` drive is mounted, so the
project path maps 1:1 and no file copying is needed.

---

## 2. What was installed, and where

Everything lives under `tools/`, which was added to `.gitignore`.

```
tools/
├── bin/
│   ├── finemap            FINEMAP v1.1 (Linux x86_64 binary)
│   ├── finemap.cmd        Windows -> WSL shim
│   ├── CAVIAR             CAVIAR v2.2 (compiled here)
│   └── CAVIAR.cmd         Windows -> WSL shim
├── gsl/                   GSL 2.7.1, built from source (CAVIAR dependency)
├── lib/                   libblas.so / liblapack.so symlinks (see 3.2)
├── caviar-src/            github.com/fhormoz/caviar @ 135b58b
├── finemapr-src/          github.com/variani/finemapr @ 65ab4e8 (patched, see 3.3)
├── finemap_example/       FINEMAP's bundled example region (used by smoke tests)
├── downloads/             source tarballs + build logs (~138 MB)
└── work/                  scratch dirs for the tools' own I/O at run time
```

| Component | Version | Source |
|---|---|---|
| FINEMAP | **v1.1** | `http://www.christianbenner.com/finemap_v1.1_x86_64.tgz` |
| CAVIAR | v2.2 (10/Apr/2018) | `https://github.com/fhormoz/caviar` @ `135b58b`, compiled locally |
| GSL | 2.7.1 | `https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz`, built from source |
| finemapr | 0.1.0 | `https://github.com/variani/finemapr` @ `65ab4e8`, patched |

Total footprint ~185 MB, none of it versioned.

### Why FINEMAP v1.1 specifically — not the current v1.4.x

This is the single most important pin in the setup. `finemapr::run_finemap()`
writes a **v1.1-format** master file:

```
z;ld;snp;config;log;n-ind
region.z;region.ld;region.snp;region.config;region.log;<n>
```

and a 2-column, headerless `.z` file (`<snp> <zscore>`). FINEMAP ≥ 1.2 changed
both: the master gained a `cred` column and renamed `n-ind` → `n_samples`, and
the `.z` file must now be 7 columns with a header (`rsid chromosome position
allele1 allele2 maf beta se`). Installing the current v1.4.2 would make every
`run_finemap()` call fail. v1.1 is also exactly what finemapr's own install
notes specify.

---

## 3. Chronological record

### 3.1 Locating the install instructions

`README.md` only links to the tools' homepages. The actual instructions are in
**`misc/install-finemaping-tools.md`** in the finemapr repo, which pins
`finemap_v1.1_x86_64` and gives the CAVIAR `git clone` + `make` recipe. The
repo installs to `~/.local/apps/<tool>/`; per the requirement that binaries be
kept in a gitignored project folder, this was redirected to `tools/`.

Also read `R/finemap.R`, `R/caviar.R` and `R/zzz.R` to pin down the exact
contract: option names `finemapr_finemap` / `finemapr_caviar`, the expected
binary names (`finemap`, `CAVIAR`), the file formats above, and CAVIAR's
invocation `-z <zscores> -l <ld> -o <prefix>` producing `<prefix>_post`,
`<prefix>_set`, `<prefix>.log`.

### 3.2 Build prerequisites — GSL missing, no root

`make` in `caviar-src/CAVIAR-C++` failed immediately:

```
Util.cpp:10:10: fatal error: gsl/gsl_linalg.h: No such file or directory
```

Inspection showed the machine has only the *runtime* libraries (`libblas3`,
`liblapack3`) — no `-dev` packages and no GSL at all — and `sudo` is
password-protected, so `apt-get install libgsl-dev` could not be run
non-interactively.

Resolved without root, keeping everything inside `tools/`:

1. **GSL 2.7.1 built from source** into `tools/gsl` (`--enable-static
   --disable-shared`). Configure/build/install logs are in `tools/downloads/`.
2. CAVIAR's Makefile links `-llapack -lblas`, which needs the unversioned
   `.so` names that only the `-dev` packages provide. Symlinks
   `tools/lib/libblas.so` → `/usr/lib/x86_64-linux-gnu/libblas.so.3` and
   likewise for lapack were created, and the link step pointed at them.
3. CAVIAR is compiled by invoking `g++` directly rather than via `make`,
   because the shipped Makefile hardcodes its include/link flags and cannot
   pick up the local GSL. The source list and flags mirror `Makefile:17`
   (bundled header-only armadillo, `-DARMA_DONT_USE_WRAPPER`).

On a machine with root, all of this collapses to
`sudo apt-get install -y libgsl-dev liblapack-dev libblas-dev` — noted in the
install script.

### 3.3 finemapr patch for R ≥ 4.0

finemapr was last updated in 2018 and `run_caviar()` guards its inputs with:

```r
stopifnot(class(ld) == "matrix")
stopifnot(class(dir_run) == "character")
```

Since **R 4.0.0**, `class(<matrix>)` returns `c("matrix", "array")`, so the
first comparison evaluates to `c(TRUE, FALSE)` and `stopifnot()` aborts —
`run_caviar()` is simply unusable on modern R. This was confirmed empirically
(`Error: class(ld) == "matrix" are not all TRUE`) before patching.

The install script rewrites both guards to `is.matrix()` / `is.character()`,
preserving the intended check. `run_finemap()` already uses `is.matrix()` and
needed no patch. The patch is applied by `sed` in
`install_finemap_caviar.sh` (idempotent) and is re-applied on any fresh clone.

Other compatibility notes:
- `tibble::as_data_frame`, `tibble::data_frame` and `dplyr::select_` are all
  imported by finemapr's NAMESPACE. All three are deprecated but **still
  exported** in the installed versions (tibble 3.3.0, dplyr 1.1.4), so the
  package loads and runs; it emits deprecation warnings on use. Harmless, but
  it is what would break first on a future tibble/dplyr release.
- Install emits `replacing previous import 'Matrix::head' by 'utils::head'` —
  benign, from finemapr's own NAMESPACE.
- All other finemapr `Imports` (tibble, dplyr, readr, magrittr, ggplot2,
  cowplot, data.table, Matrix) were already in the renv library.

### 3.4 `finemapr::extract_credible_set()` is unusable on this path

`extract_credible_set()` only supports the newer `finemapr()` S3 pipeline: it
expects `x$snp` to be a *list* of tables carrying a `snp_prob_cumsum` column
and reads `x$prop_credible`. The deprecated `run_finemap()`/`run_caviar()`
functions return `x$snp` as a single tibble with neither field, so `lapply()`
iterates over columns and the call dies with

```
no applicable method for 'filter' applied to an object of class "c('double','numeric')"
```

Rather than fight this, credible sets are computed in our own script
(`credible_set_from_pip()`: sort by PIP descending, accumulate to the coverage
threshold). This has the side benefit of defining credible sets identically for
FINEMAP, CAVIAR and SuSiE, so the three are directly comparable.

### 3.5 Package installation

`finemapr` was installed from the patched local clone into the **renv project
library** with `R CMD INSTALL --no-multiarch --no-staged-install
tools/finemapr-src`. `renv.lock` was **not** touched (no `renv::snapshot()`),
and `renv/` is already gitignored, so this leaves no trace in version control.

---

## 4. The Windows → WSL bridge

`tools/bin/finemap.cmd` and `tools/bin/CAVIAR.cmd` are generated by the install
script. Each resolves the WSL path of the current directory and of its sibling
Linux binary, then runs the binary inside WSL with that working directory:

```bat
for /f "usebackq delims=" %%i in (`wsl -d Ubuntu -e wslpath -a "%CD%"`) do set "WSLCWD=%%i"
for /f "usebackq delims=" %%i in (`wsl -d Ubuntu -e wslpath -a "%~dp0finemap"`) do set "WSLBIN=%%i"
wsl -d Ubuntu -e bash -c "cd '!WSLCWD!' && '!WSLBIN!' %*"
```

This works because `finemapr` `setwd()`s into the run directory and then passes
only **relative** filenames (`region.z`, `region.ld`, `region.master`). So the
shim only has to translate the working directory — no argument rewriting is
needed at all. `%~dp0` keeps the shims relocatable rather than hardcoding an
absolute path.

`finemapr_setup.R` picks the `.cmd` shim on Windows and the bare binary
elsewhere, so the analysis script is portable to a Linux box unchanged.

---

## 5. Verification performed

| Check | Result |
|---|---|
| FINEMAP smoke test (bundled `region1`, direct in WSL) | PASS — 51-line `region1.snp` |
| CAVIAR smoke test (bundled `region1`, direct in WSL) | PASS — causal set `rs15 rs47`, matching the expected output in finemapr's install notes |
| `finemap.cmd` / `CAVIAR.cmd` invoked from Windows PowerShell | PASS — exit 0, identical outputs |
| `finemapr::run_finemap()` from Windows R via shim | PASS — rs15/rs47 at PIP 1.0 |
| `finemapr::run_caviar()` from Windows R via shim (after patch) | PASS — rs15/rs47 at PIP 1.0 |
| `parse()` of the analysis script | PASS — 36 top-level expressions |
| Synthetic end-to-end smoke test (see below) | PASS |

**Synthetic smoke test.** 250 strains, 60 markers across two mock chromosomes,
every third marker duplicated to create perfect-LD partners, with `YAL010W` as
the single causal marker. Result:

- LD tag collapsing reduced each 30-marker region to 20 tags, exactly as
  designed.
- **Both tools recovered the true causal marker**: FINEMAP PIP 0.998, CAVIAR
  PIP 0.99999, both in the credible set.
- The null region (`chr2`) produced diffuse PIPs and a correspondingly large
  credible set — the expected behaviour when there is no signal.

The smoke test used synthetic data only; no real dataset was fine-mapped.

---

## 6. Design decisions in the analysis script

The SuSiE script fine-maps all ~2,700–3,000 varying markers of a condition in
one `susie_rss()` call. **That approach does not transfer to FINEMAP/CAVIAR**,
for two independent reasons:

1. **n < m.** Each condition has ~600 strains but ~2,700 markers, so the
   in-sample LD matrix has rank ≤ 600 and is badly singular. CAVIAR factorises
   the LD matrix directly and degrades sharply when it is rank-deficient.
2. **CAVIAR enumerates causal configurations** (O(m^c)). With m = 2,700 and
   c = 2 that is ~3.6 M configurations per condition, each involving linear
   algebra on the full matrix — not tractable. CAVIAR's own manual uses ~100
   SNPs.

Both are locus fine-mappers, so the script splits markers into **regions**,
fine-maps each independently, and stitches the results back into one
genome-wide table per condition.

| Decision | Default | Rationale |
|---|---|---|
| Region definition | whole chromosomes (`region_by = "chromosome"`) | In a two-parent cross, linkage blocks span chromosome arms, so a chromosome is a natural near-independent unit. A sliding-window mode (`"window"`, ordered by genomic position) is also provided. |
| Collapse near-perfect LD | `prune_r2 = 0.99` | Neighbouring markers in a cross are frequently perfectly correlated, which makes CAVIAR singular. Markers are greedily collapsed to a tag (highest \|z\| wins); collapsed partners are carried in a `tagged_markers` column so credible sets can still be read as the full set of implicated genes. |
| Region size cap | `max_markers = 400` | Bounds CAVIAR's runtime. If exceeded, top-\|z\| markers are kept and the number dropped is recorded in the meta table rather than silently discarded. |
| Max causal variants | `n_causal = 3` | FINEMAP `--n-causal-max`, CAVIAR `-c`. |
| Credible-set coverage | `0.95` | Matches the SuSiE script. |
| z-scores | univariate logistic regression via `fastglm` | Byte-identical helper to the SuSiE script's `uni_logreg_z()`, so all three methods consume the same summary statistics. Deliberately duplicated rather than shared, so either script runs standalone — **keep the two in sync**. |

`uni_logreg_z()` and `marginal_stats_for_genes()` are duplicated from the SuSiE
script for this reason; that is the main maintenance liability introduced here.

---

## 7. Files created and modified

**Created** (all new, nothing overwritten):

| Path | Purpose |
|---|---|
| `src/finemapping/install_finemap_caviar.sh` | Fetches/builds FINEMAP, GSL, CAVIAR; patches finemapr; writes shims; runs smoke tests. Idempotent. |
| `src/finemapping/finemapr_setup.R` | Installs finemapr if absent, sets `options(finemapr_*)`, validates the binaries exist. |
| `src/finemapping/finemap_caviar_conditionwise_analysis.R` | The analysis script (the deliverable). |
| `src/finemapping/SETUP_LOG.md` | This file. |
| `tools/**` | All third-party binaries and sources (gitignored). |

**Modified** — exactly one file:

| Path | Change |
|---|---|
| `.gitignore` | Appended a `tools/` entry (with a comment). Nothing removed. |

**Not touched:** no existing script, dataset, result or `renv.lock` was
modified, and **no files were deleted**. The pre-existing uncommitted changes in
`src/contingency/`, `src/result_model_evaluation.qmd` and `src/SuSiE analysis/`
are untouched and unrelated. The synthetic smoke-test output left under
`tools/work/*/SmokeTest/` was intentionally left in place; it is gitignored and
safe to delete at any time.

---

## 8. How to run

One-off setup (already done on this machine; idempotent, safe to re-run):

```powershell
wsl -d Ubuntu -- bash /mnt/h/Rajeeva/Project/Himanshu/subsets/R_analysis-yeast_growth_pred/src/finemapping/install_finemap_caviar.sh
```

Then, from the **project root** (required — renv activates from there):

```r
source("src/finemapping/finemap_caviar_conditionwise_analysis.R")
setup_finemapr()

# recommended first real run: one dataset, two conditions
res <- run_conditionwise_finemapping("bloom2013", "Bloom2013",
                                     conditions = c("4NQO", "maltose"))
cmp <- compare_finemapping_shap(res, "Bloom2013", "finemap")
qtl <- compare_finemapping_qtl_bloom2013(res, "finemap")
mth <- compare_methods(res, "Bloom2013")     # FINEMAP vs CAVIAR vs SuSiE

# FINEMAP only - substantially faster than CAVIAR
res <- run_conditionwise_finemapping("bloom2013", "Bloom2013",
                                     methods = "finemap")

# everything (all datasets, all conditions, both tools) - slow
results <- main_conditionwise_finemapping()
```

Outputs mirror the SuSiE layout, under
`results/result_finemapping/0.5_condition-wise/<method>/<Dataset>/`:
per-condition `*_credible_sets.csv`, combined `*_conditionwise_pip.csv` and
`*_conditionwise_meta.csv`, a `comparison/` folder (Spearman table + bar chart,
per-condition PIP-vs-SHAP scatters, Venns, QTL overlap) and a `triangulation/`
folder. Cross-method results go to
`results/result_finemapping/0.5_condition-wise/method_comparison/<Dataset>/`.

---

## 9. Caveats and things to watch

- **Runtime is unmeasured.** No real condition has been fine-mapped, so the
  cost of a full run is unknown. CAVIAR is the bottleneck; start with a
  two-condition subset and check `n_markers_used` in the meta table before
  committing to all 39 conditions × 3 datasets.
- **`max_markers = 400` may still be slow for CAVIAR.** If so, lower it or
  switch to `region_by = "window"` with a smaller `window_size`. Whatever is
  dropped is always reported in the meta table.
- **Region choice is a real modelling assumption.** Chromosome-wide regions
  assume signals do not straddle chromosomes (safe) and that within-chromosome
  LD is handled by the tools (reasonable for a cross, but worth sanity-checking
  against the SuSiE credible sets via `compare_methods()`).
- **FINEMAP is pinned to v1.1 by finemapr's file format.** Do not "upgrade" the
  binary without also rewriting the master/`.z` writers.
- **finemapr is unmaintained** (last release 2018) and relies on three
  deprecated tidyverse functions. A future tibble/dplyr release that removes
  `data_frame()`, `as_data_frame()` or `select_()` will break it; the fix would
  be to vendor the two `run_*` functions locally.
- The `tools/downloads/` folder holds ~138 MB of tarballs and build logs. Safe
  to delete once the binaries are built.
