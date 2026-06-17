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
cyftools convert cells.csv cells.byf -s 1       # -s = sample id
```

cyftools picks the output form from the file extension, mirroring SAM/BAM:
**`.cyf`** is the human-readable tab-delimited text form (the SAM analog) and
**`.byf`** is the BGZF-compressed binary form (the BAM analog — compact, and
readable by `zcat`/`bgzip`). Reads auto-detect the form, so you don't have to
tell the tools which is which. A third extension, **`.ocyf`**, denotes the legacy
serialized form: it is read-only (for migrating old files) and is never written.

**Look at it:**

```sh
cyftools view -H cells.byf                        # header only (markers, provenance)
cyftools view -R cells.byf | head                 # first rows as CSV
cyftools count cells.byf                           # number of cells
cyftools summary cells.byf                         # per-column min/max/mean
cyftools view -l cells.byf                         # list marker names
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
cyftools cut cells.byf slim.cyf -x CD4,CD8,PD-1   # keep only these columns
cyftools head -n 1000 cells.byf first1k.cyf       # first 1000 cells
cyftools crop cells.byf box.cyf ...               # spatial crop to a bounding box
```

**Combine samples:**

```sh
cyftools cat sampleA.cyf sampleB.cyf > cohort.cyf
```

**Transform values / call phenotypes:**

```sh
cyftools log10 cells.byf logged.cyf               # log10 marker intensities
cyftools pheno cells.byf phenotyped.cyf ...        # gate markers into phenotype flags
```

**Spatial analysis:**

```sh
cyftools radialdens cells.byf dens.cyf ...         # neighbor density in radial rings
cyftools margin cells.byf margin.cyf ...           # tumor/stroma boundary cells
cyftools filter cells.byf - -M | cyftools dbscan - clusters.cyf   # cluster marked cells
cyftools tls cells.byf tls.cyf                      # call tertiary lymphoid structures
cyftools delaunay cells.byf graph.cyf ...           # Delaunay neighbor graph
```

Run any command with no arguments to see its full help.

> **Output form follows the extension.** Write `out.cyf` for text or `out.byf` for BGZF
> binary; stdout (`-`) defaults to binary. The legacy `.ocyf` form is read-only — reads
> auto-detect it, but writing `.ocyf` is refused. To migrate an old file, just convert it:
> `cyftools clean old.ocyf new.byf`. See [`docs/CYF_INTEGRATION.md`](docs/CYF_INTEGRATION.md).

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
SAM-style tagged lines: `@HD` (header/version), `@SA` (sample), `@MA` (marker columns),
`@CA` (metadata columns), `@GA` (graph columns), `@FL` (flag-bit meanings), and `@PG`
(the provenance chain). Each record carries an `id`, `x`/`y` coordinates, two 64-bit
flag registers (`cflag`, `pflag`), and one typed value per data column. There are two
on-disk forms (the SAM/BAM two-tier model): **`.cyf`** text and **`.byf`** binary, the
latter being the explicit little-endian record stream wrapped in **BGZF** (blocked gzip,
so it's `zcat`-readable and carries the virtual offsets needed for future indexing). The
layout is readable from the spec by any language (a Python reference reader lives in
[`docs/reference/cyf.py`](docs/reference/cyf.py)). Full details: [`docs/CYF_FORMAT.md`](docs/CYF_FORMAT.md).

---

## Viewing in the browser (cyfview)

`viewer/cyfview.html` is a browser GUI for CYF tables; `viewer/view.py` serves it and runs
`cyftools` to convert whatever you load. The **same `view.py` + `cyftools`** run three ways —
a quick local launch, a local Docker container, or Cloud Run with a GCS bucket for multi-GB
uploads. The deep cloud config (IAM, bucket CORS, sizing) lives in
[`viewer/DEPLOY.md`](viewer/DEPLOY.md).

### A. Local, no Docker (fastest to iterate)

Needs `build/cyftools` and Python 3 — nothing else:

```sh
python3 viewer/view.py                 # serve; open the printed URL, drag a .byf onto the page
python3 viewer/view.py cells.byf       # ...or auto-open cells.byf, preloaded
```

Edit `viewer/cyfview.html` and just **refresh** the page; edit `viewer/view.py` and
**restart** the command. No build step on this path.

### B. Local, in Docker (the exact image that deploys to the cloud)

```sh
git submodule update --init --recursive          # once
docker build -t cyftools-viewer .                # build the image
docker run --rm -p 8080:8080 cyftools-viewer     # open http://localhost:8080/
```

### Rebuild the image after a code change

The container holds a **baked-in copy** of the code, so after editing anything you must
rebuild the image (path A above does not — it reads the files live):

```sh
docker build -t cyftools-viewer .                # rebuild
docker run --rm -p 8080:8080 cyftools-viewer     # run the new image
```

- Changed only `viewer/` (HTML or Python)? The rebuild is **seconds** — the C++ compile
  layer is cached.
- Changed `src/` or `CMakeLists.txt`? `cyftools` **recompiles** (a few minutes).

### C. Deploy to Cloud Run + GCS (multi-GB uploads)

`gcloud run deploy --source .` rebuilds the image with Cloud Build and deploys it, so
shipping a code change is just re-running that one command.

```sh
# 1. a bucket for uploads + render packs
gcloud storage buckets create gs://MY-CYF-BUCKET --location=us-central1

# 2. let the Cloud Run service account read/write objects AND sign upload URLs
SA="cyftools-viewer@MY-PROJECT.iam.gserviceaccount.com"
gcloud storage buckets add-iam-policy-binding gs://MY-CYF-BUCKET \
  --member="serviceAccount:$SA" --role="roles/storage.objectAdmin"
gcloud iam service-accounts add-iam-policy-binding "$SA" \
  --member="serviceAccount:$SA" --role="roles/iam.serviceAccountTokenCreator"

# 3. build + deploy (re-run this to ship any later code change)
gcloud run deploy cyftools-viewer --source . --region us-central1 \
  --service-account "$SA" --allow-unauthenticated \
  --set-env-vars CYFVIEW_BUCKET=MY-CYF-BUCKET,CYFVIEW_MAX_CELLS=4000000 \
  --memory 4Gi --cpu 2 --timeout 900
```

After the first deploy, set the bucket's CORS to allow the printed `run.app` origin to PUT/GET
it (one command in [`viewer/DEPLOY.md`](viewer/DEPLOY.md)). Without a bucket the container still
works — uploads just stay local; the bucket is what lifts the request-size limit so multi-GB
`.byf` files can be uploaded.

---

## Documentation

- [`docs/CYF_FORMAT.md`](docs/CYF_FORMAT.md) — the on-disk format specification (text + binary).
- [`docs/ARCHITECTURE_REVIEW.md`](docs/ARCHITECTURE_REVIEW.md) — codebase architecture and design notes.
- [`docs/CYF_INTEGRATION.md`](docs/CYF_INTEGRATION.md) — how the CYF reader/writer is wired in (and the path off cereal).
- [`docs/DEPENDENCY_CLEANUP.md`](docs/DEPENDENCY_CLEANUP.md) — the dependency-slimming plan and status.
- [`docs/reference/cyf.py`](docs/reference/cyf.py) — a pure-Python reference codec.
- [`viewer/DEPLOY.md`](viewer/DEPLOY.md) — cyfview cloud deploy: Docker, Cloud Run, GCS bucket setup.
