# cyftools

**A command-line toolkit for wrangling single-cell tables from multiplex imaging (CyCIF / mIF).**

cyftools is to multiplex-imaging single-cell data what **samtools** is to sequencing
alignments. After you segment cells in a high-plex image and measure per-marker
intensities, you get a big table — one row per cell, columns for coordinates, marker
intensities, and derived annotations. cyftools stores that table in a compact,
self-describing binary format (**CYF**) and gives you ~50 composable subcommands that
read a stream of cells, transform it, and write a stream of cells — so commands chain
together through Unix pipes.

```
            samtools  :  SAM/BAM  ::  cyftools  :  CYF
```

| Genomics (SAM/BAM) | cyftools (CYF) |
|---|---|
| read / alignment record | a **cell** record |
| reference position (POS) | `x`, `y` (slide coordinates) |
| FLAG bitfield | `cflag` (structure) + `pflag` (phenotype) |
| `@SQ` sequence dictionary | `@MA` marker columns |
| `@RG` read group | `@SA` sample |
| `@PG` program chain | `@PG` provenance (every command stamps one) |
| `samtools view \| samtools sort \| ...` | `cyftools convert \| cyftools filter \| cyftools view` |

The on-disk format is specified in [`docs/CYF_FORMAT.md`](docs/CYF_FORMAT.md).

---

## Install

The core (the CYF format plus the streaming and spatial commands) builds with just a
C++17 compiler and header-only libraries (`nanoflann` and `delaunator` are vendored in
`src/`; `cereal` and `knncolle` are git submodules). The only remaining system
dependencies are for image and plot output.

```sh
# clone with submodules
git clone --recursive https://github.com/walaj/cyftools
cd cyftools

# system deps (OSX) — only needed for image/plot commands
brew install libjpeg libtiff libomp cairo

# build (CMake — recommended)
cmake -S . -B build
cmake --build build -j
ctest --test-dir build          # run the format self-test
```

The build is configurable; turn features off if you don't need them (a build with
all of these off needs no system libraries at all):

```sh
cmake -S . -B build -DCYFTOOLS_WITH_TIFF=OFF -DCYFTOOLS_WITH_CAIRO=OFF
```

| Option | Default | Feature |
|---|---|---|
| `CYFTOOLS_WITH_OPENMP` | ON | OpenMP parallelism |
| `CYFTOOLS_WITH_TIFF`   | ON | TIFF image I/O (`convolve`, `png`) |
| `CYFTOOLS_WITH_CAIRO`  | ON | Cairo PDF/PNG plotting |
| `CYFTOOLS_WITH_LDA`    | OFF | LDA topic modeling (needs supervised-lda + Eigen) |

(The legacy hand-written `src/Makefile` still works — `cd src && make all` — but CMake
is the maintained path.) See the [dependency-cleanup notes](docs/DEPENDENCY_CLEANUP.md)
for what was slimmed and why.

Quick self-test of the format codec (no system deps needed):

```sh
cd src && make cyf-test
```

---

## Quickstart

cyftools reads and writes the binary CYF format. Use a filename, or `-` to stream
from stdin / to stdout so commands can be piped.

**Import a CSV** (must have `CellID`, `X`/`X_centroid`, `Y`/`Y_centroid` columns, then
one column per marker):

```sh
cyftools convert cells.csv cells.bcyf -s 1       # -s = sample id
```

cyftools has two file forms, mirroring SAM/BAM: **`.cyf`** is the human-readable
tab-delimited text form (the SAM analog), and **`.bcyf`** is the BGZF-compressed
binary form (the BAM analog — compact, and readable by `zcat`/`bgzip`). Reads
auto-detect the form, so you don't have to tell the tools which is which.

**Look at it:**

```sh
cyftools view -H cells.bcyf                        # header only (markers, provenance)
cyftools view -R cells.bcyf | head                 # first rows as CSV
cyftools count cells.bcyf                           # number of cells
cyftools summary cells.bcyf                         # per-column min/max/mean
cyftools view -l cells.bcyf                         # list marker names
```

**Pipe commands together** (the samtools idiom — `-` is stdin/stdout):

```sh
# import, keep cells whose phenotype flag has bit 256 set, view as CSV
# (-s is a phenotype-flag OR mask; which bit means which marker depends on your panel)
cyftools convert cells.csv - -s 1 \
  | cyftools filter - - -s 256 \
  | cyftools view -R -
```

**Subset columns and cells:**

```sh
cyftools cut cells.bcyf slim.cyf -x CD4,CD8,PD-1   # keep only these columns
cyftools head -n 1000 cells.bcyf first1k.cyf       # first 1000 cells
cyftools crop cells.bcyf box.cyf ...               # spatial crop to a bounding box
```

**Combine samples:**

```sh
cyftools cat sampleA.cyf sampleB.cyf > cohort.cyf
```

**Transform values / call phenotypes:**

```sh
cyftools log10 cells.bcyf logged.cyf               # log10 marker intensities
cyftools pheno cells.bcyf phenotyped.cyf ...        # gate markers into phenotype flags
```

**Spatial analysis:**

```sh
cyftools radialdens cells.bcyf dens.cyf ...         # neighbor density in radial rings
cyftools margin cells.bcyf margin.cyf ...           # tumor/stroma boundary cells
cyftools filter cells.bcyf - -M | cyftools dbscan - clusters.cyf   # cluster marked cells
cyftools tls cells.bcyf tls.cyf                      # call tertiary lymphoid structures
cyftools delaunay cells.bcyf graph.cyf ...           # Delaunay neighbor graph
```

Run any command with no arguments to see its full help.

> **Writing CYF vs. the legacy format.** During the format migration, output defaults to
> the legacy (cereal) encoding; set `export CYFTOOLS_FORMAT=cyf` to write the new CYF
> format. Reads auto-detect either format. See [`docs/CYF_INTEGRATION.md`](docs/CYF_INTEGRATION.md).

---

## The command set

Every command reads a CYF stream and (mostly) writes one, so they compose. Grouped:

**I/O & inspection** — `convert` (CSV→CYF), `view`, `cat`, `head`, `reheader`, `info`,
`check`, `count`, `cellcount`, `summary`, `synth`, `debug`

**Row / column selection** — `cut`, `clean`, `trim`, `filter`, `sampleselect`,
`subsample`, `crop`, `roi`

**Value transforms** — `log10`, `divide`, `rescale`, `magnify`, `offset`, `flip`,
`scramble`, `hallucinate`, `flagset`, `pheno`, `mean`

**Spatial & graph** — `radialdens`, `delaunay`, `dbscan`, `tls`, `island`, `margin`,
`annotate`, `dist`

**Numeric / stats** — `pearson`, `jaccard`, `histogram`

**Latent Dirichlet allocation** — `ldacreate`, `ldarun` *(optional build)*

**Imaging & plots** — `png`, `plot`, `scatter`, `convolve`

---

## The CYF format in one paragraph

A CYF file is a header followed by a stream of cell records. The header is a set of
SAM-style tagged lines: `@VN` (version), `@SA` (sample), `@MA` (marker columns),
`@CA` (metadata columns), `@GA` (graph columns), `@FL` (flag-bit meanings), and `@PG`
(the provenance chain). Each record carries an `id`, `x`/`y` coordinates, two 64-bit
flag registers (`cflag`, `pflag`), and one typed value per data column. There are two
on-disk forms (the SAM/BAM two-tier model): **`.cyf`** text and **`.bcyf`** binary, the
latter being the explicit little-endian record stream wrapped in **BGZF** (blocked gzip,
so it's `zcat`-readable and carries the virtual offsets needed for future indexing). The
layout is readable from the spec by any language (a Python reference reader lives in
[`docs/reference/cyf.py`](docs/reference/cyf.py)). Full details: [`docs/CYF_FORMAT.md`](docs/CYF_FORMAT.md).

---

## Documentation

- [`docs/CYF_FORMAT.md`](docs/CYF_FORMAT.md) — the on-disk format specification (text + binary).
- [`docs/ARCHITECTURE_REVIEW.md`](docs/ARCHITECTURE_REVIEW.md) — codebase architecture and design notes.
- [`docs/CYF_INTEGRATION.md`](docs/CYF_INTEGRATION.md) — how the CYF reader/writer is wired in (and the path off cereal).
- [`docs/DEPENDENCY_CLEANUP.md`](docs/DEPENDENCY_CLEANUP.md) — the dependency-slimming plan and status.
- [`docs/reference/cyf.py`](docs/reference/cyf.py) — a pure-Python reference codec.
