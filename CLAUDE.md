# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

cyftools is a single C++17 command-line binary for wrangling single-cell tables from
multiplex imaging (CyCIF / Orion / mIF). The mental model is **samtools for cells**:
SAM/BAM is to sequencing alignments as the **CYF** format is to per-cell measurement
tables. ~50 subcommands each read a stream of cell records, transform it, and write a
stream — so commands chain through Unix pipes. `-` means stdin/stdout.

## Build, test, run

The build is CMake-first. The format core plus streaming/spatial commands need only a
C++17 compiler and header-only libs; image/plot/LDA features are optional and gated.

```sh
git submodule update --init                 # cereal + knncolle live in external/
cmake -S . -B build && cmake --build build -j
ctest --test-dir build                       # runs the CYF format self-test (cyf_test)
./build/cyftools <command>                   # run the tool; no args to a command prints its help
```

Run a single test (there is currently one CTest target, the format codec round-trip):

```sh
ctest --test-dir build -R cyf_test -V        # via ctest
./build/cyf_test                             # or run the test binary directly
```

Minimal core build with no system libraries (turn the optional features off):

```sh
cmake -S . -B build -DCYFTOOLS_WITH_TIFF=OFF -DCYFTOOLS_WITH_CAIRO=OFF
```

Feature options (each maps to a `HAVE_*` compile definition the sources gate on
internally — every source file is always compiled, optional code just disappears when
its macro is off): `CYFTOOLS_WITH_OPENMP` (ON), `CYFTOOLS_WITH_TIFF` (ON, `convolve`/`png`),
`CYFTOOLS_WITH_CAIRO` (ON, PDF/PNG plotting), `CYFTOOLS_WITH_LDA` (OFF, needs supervised-lda + Eigen),
`CYFTOOLS_64BIT` (ON, 64-bit flag registers → `USE_64_BIT`).

A legacy hand-written `src/Makefile` still works (`cd src && make all`, `make cyf-test`)
but CMake is the maintained path. On macOS, OpenMP needs `brew install libomp` (the build
auto-detects Homebrew's libomp); image/plot features need `brew install libjpeg libtiff cairo`.

## Architecture

The codebase has three layers; understanding the boundary between them is the key to
being productive.

**1. Command dispatch — `src/cyftools.cpp` (the big ~128KB file).**
`main()` parses the first argv as a module name and dispatches to a per-command
`<name>func(argc, argv)` function via a long `if (opt::module == "...")` chain. Each
`*func` parses that command's own options, then does one of two things:

- *Streaming commands* construct a `CellProcessor` subclass and call
  `table.StreamTable(processor, infile)`. Cells flow through one at a time — low memory,
  composable. This is the common path.
- *Whole-table commands* (spatial/graph algorithms that need all cells at once) call a
  `CellTable` method like `Distances()`, `CallTLS()`, `Delaunay()` and then `OutputTable()`.

When adding a command, you touch `cyftools.cpp` (dispatch + option parsing + usage
string) and add a processor or a `CellTable` method.

**2. Streaming processors — `src/cell_processor.{h,cpp}`.**
~40 small classes deriving from `CellProcessor` (or `LineProcessor`), each implementing
`ProcessHeader(CellHeader&)` and `ProcessLine(Cell&)`. This is the per-cell transform
interface that `StreamTable` drives. A processor that only reads (e.g. counts) is a
reader; one that emits writes cells through the table's output writer. New
streaming transforms are usually a new processor here.

**3. The table + heavy algorithms — `src/cell_table*.cpp`.**
`CellTable` (`cell_table.cpp`) is the in-memory representation and owns I/O setup. The
algorithm-heavy code is split out by concern: `cell_table_cluster.cpp` (DBSCAN on
nanoflann), `cell_table_delaunay.cpp` (Delaunay graph), `cell_table_graph.cpp` (neighbor
graphs), `cell_table_lda.cpp` (topic modeling, optional).

### The CYF format and I/O

The on-disk format is specified in `docs/CYF_FORMAT.md` (authoritative — read it before
touching serialization). A language-agnostic Python reference reader is in
`docs/reference/cyf.py`. Shape: a SAM-style tagged **header** (`@HD` version, `@SA`
sample, `@MA` markers, `@CA` metadata, `@GA` graph, `@FL` flag-bit meanings, `@PG`
provenance chain) followed by a stream of cell **records**. The header declares the
column schema and the meaning of every flag bit, so a file is self-describing. Each
record has `id`, `x`/`y`, two 64-bit flag registers (`cflag` = structure, `pflag` =
phenotype), plus one typed value per data column.

Three on-disk forms, picked by extension (mirroring SAM/BAM):

- `.cyf` — TAB-delimited text (the SAM analog).
- `.byf` — binary record stream wrapped in **BGZF** (the BAM analog; `zcat`-readable). `-`/stdout defaults to this.
- `.ocyf` — legacy cereal-serialized form, **read-only**; reads auto-detect it, writing it is refused. Migrate with `cyftools clean old.ocyf new.byf`.

I/O implementation lives in `cyf_io.cpp`, `cell_header.cpp`, and `bgzf.cpp`; format
field/flag definitions in `cyf_field.h`, `cyf_flags.h`, `cell_flag.cpp`. **cereal is still
a core dependency** during the migration off the legacy serialized format — see
`docs/CYF_INTEGRATION.md` and `docs/DEPENDENCY_CLEANUP.md` for the migration status and
direction. The `@PG` provenance chain is stamped by every command (`cmd_input` in
`cyftools.cpp`), so changes to dispatch should preserve that.

### Vendored vs submoduled dependencies

Header-only libs vendored directly in `src/`: `nanoflann.hpp` (KD-tree / nearest
neighbor), `delaunator.hpp` (Delaunay), `csv.h` (CSV reader — uses a background
`std::thread`, hence the `-pthread` requirement at compile *and* link). Git submodules
in `external/`: `cereal` (serialization), `knncolle` (KNN, used by `island`/`dist`/`tls`,
built with the kmknn/HNSW/Annoy backends disabled). zlib is the only always-required
system lib (BGZF).

## Browser viewer & cloud deploy (`viewer/`)

`viewer/cyfview.html` is a single-file deck.gl GUI for a CYF table; `viewer/view.py` is
the server that feeds it. **Invariant: the browser never reimplements the format or any
algorithm — every conversion shells out to the real `cyftools` binary.** The GUI renders a
`.cyfv` "pack" (`cyftools export` → magic `CYFV`, columnar id/x/y/cflag/pflag + markers,
see `docs/VIEWER_PACK.md`); `.cyfv` is an internal intermediary, never user-facing (users
only ever name `.byf`).

`view.py` runs two ways: a **local launch** (`python3 viewer/view.py [cells.byf]` →
`127.0.0.1`, auto port, opens a browser) or **server mode** when `$PORT` is set (binds
`0.0.0.0:$PORT`, threaded, no browser — for Docker/Cloud Run). Loading a `.byf` uses a
**staged flow** so multi-GB files aren't bound by any request-size limit: `POST
/api/upload-url` → browser **PUTs the file straight to an object store** → `POST /api/pack`
(server runs `cyftools count`, then `export`, or `subsample -r | export` when over
`CYFVIEW_MAX_CELLS` so the browser payload stays bounded) → browser fetches the pack. The
object store is pluggable: `LocalStore` (cache dir; PUT to `/blob/...`) or `GcsStore`
(signed URLs to a GCS bucket). A raw `POST /convert` (whole body) is the small-file
fallback. Browser-side magic-byte sniff routes `CYFV`→render-direct, BGZF→server.

A second, **standalone** browser tool lives in `docs/` and is served statically via
**GitHub Pages** (no server, no Docker): `docs/cyf_cohort.html` is the cohort phenotyper —
it loads a `cohort.json` (from `cyftools cohort`), groups/compares samples, and exports
publication-ready R/ggplot2 figures standardized by an optional plot config
(`docs/cohort.config.example.json`, authored in-page via "⚙ Build config…"). It is pure
client-side (drag a file in, or `?file=`/`?config=` URL params) and only consumes
`cohort.json` — it never re-derives metrics, honoring the same "never reimplement" spirit.
`docs/.nojekyll` + `docs/index.html` make Pages serve it byte-for-byte. Enable it under
repo Settings → Pages → *Deploy from branch* `main` `/docs`.

Deploy: `Dockerfile` is a two-stage build (minimal-core `cyftools`, i.e. TIFF/Cairo/LDA/
OpenMP off → runtime needs only `python3`+`zlib`+`google-cloud-storage`). **After any code
change you must rebuild the image** (container has a baked-in copy); viewer-only edits
rebuild in seconds (compile cached), `src/`/`CMakeLists.txt` edits recompile. The build
stage copies only `CMakeLists.txt`/`src/`/`external/`, so viewer edits don't bust the
compile cache. Full local+Cloud Run+GCS steps: `README.md` ("Viewing in the browser") and
`viewer/DEPLOY.md`. Note: fixed-width int headers must `#include <cstdint>` (macOS/clang
hides a missing include that gcc/libstdc++ surfaces — keep the tree Linux-buildable).

## Conventions

- Output form follows the **filename extension**, never a flag. Don't add format-select flags.
- Commands compose through pipes; prefer a streaming `CellProcessor` over an in-memory
  whole-table pass unless the algorithm genuinely needs all cells.
- Optional/heavy code is guarded by `HAVE_*` macros so the minimal core keeps building
  with zero system deps — keep new optional features behind a `CYFTOOLS_WITH_*` option
  and its `HAVE_*` define rather than making them unconditional.
- The version string has a single source of truth: `project(... VERSION ...)` in
  `CMakeLists.txt`, surfaced as `CYFTOOLS_VERSION`. Don't hard-code it elsewhere.

## Key docs

`docs/CYF_FORMAT.md` (format spec), `docs/ARCHITECTURE_REVIEW.md` (design notes),
`docs/CYF_INTEGRATION.md` (reader/writer wiring + path off cereal),
`docs/DEPENDENCY_CLEANUP.md` (what was slimmed and why),
`docs/VIEWER_PACK.md` (`.cyfv` pack byte layout), `viewer/DEPLOY.md` (cyfview Docker/Cloud
Run/GCS deploy).
