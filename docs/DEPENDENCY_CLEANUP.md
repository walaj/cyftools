# Dependency cleanup — lowest-hanging fruit

Goal: shrink cyftools' build dependencies toward an htslib-style "compiler + a couple
of headers" baseline. This maps each dependency to the feature(s) that need it,
grounded in the actual source, and ranks removals by value ÷ effort.

> **Status: Tier 0 + the mlpack→nanoflann swap are DONE** (verified by compiling the
> formerly-heavy files and a runtime DBSCAN test — see below). Removed from the build:
> **mlpack, armadillo, CGAL, boost, gmp, hdf5, libbdplus**, and the `supervised-lda`
> hard requirement. Spatial search now uses the vendored header-only `nanoflann`
> (`src/nanoflann.hpp` + `src/cyf_kdtree.h`); DBSCAN was reimplemented on it. The
> remaining Tier 2 item (gating Cairo + libtiff/jpeg behind flags) is not yet done.

## The current dependency surface

README asks for: `libjpeg hdf5 eigen libomp cgal boost mlpack libbdplus armadillo`
plus submodules `cereal umappp CppIrlba knncolle CppKmeans` plus a hand-built
`supervised-lda` at a hardcoded `$HOME/git` path.

What actually uses what (from the code):

| Dependency | Used by | Reachable? | Weight |
|---|---|---|---|
| **libbdplus** | nothing (`grep` finds zero references) | dead | — |
| **HDF5** | `CellTable::HDF5Write` | **dead** — no command calls it | heavy (libhdf5) |
| **mlpack GMM** | `CellTable::GMM_EM` | **dead** — there is no `gmm` command | (part of mlpack) |
| **CGAL** | Voronoi-diagram PDF only (`cell_table_delaunay.cpp`, behind `#ifdef HAVE_CGAL`, already commented off in the brew/Mac build) | optional | very heavy (pulls Boost, GMP, MPFR) |
| **Boost** | one call: `boost::hash_combine` ×2 in `cell_table_delaunay.cpp` (otherwise only transitive via CGAL/mlpack) | trivial | heavy |
| **mlpack + Armadillo** | `mlpack::RangeSearch`/`Range` spatial radius search in `cell_table_graph.cpp` (5 sites: radial density, KNN, neighbors) and `DBSCAN` in `cell_table_mlpack.cpp` (the `dbscan`/`tls` commands) | live, core | **very heavy** (mlpack + Armadillo + ensmallen + BLAS/LAPACK) |
| **Cairo** | all PDF/PNG plotting (`png`, `plot`, delaunay/umap renders) | live | heavy |
| **libtiff + libjpeg** | `convolve` and TIFF I/O (`tiff_*.cpp`) | live | medium |
| **Eigen** | LDA + `umappp`/`CppIrlba` (header-only consumers) | live | light (header-only) |
| **OpenMP** | parallel loops | live | light (compiler) |
| **supervised-lda** | `ldacreate`/`ldarun` (`cell_table_lda.cpp`) | live, niche | annoying (hardcoded path, manual build) |
| **cereal** | legacy serialization | being replaced by CYF | light (header-only) |
| umappp/knncolle/CppIrlba/CppKmeans | UMAP + KNN | live | header-only (need only Eigen) |

Two things already point the way out: the codebase **already vendors `delaunator.hpp`**
(header-only Delaunay, used for the actual triangulation — CGAL only draws the Voronoi
PDF), and **already uses header-only `knncolle`** for KNN alongside mlpack. The heavy
libraries are doing less than the README implies.

## Ranked plan

### Tier 0 — free wins, no real feature loss (do first)

1. **Delete `libbdplus` from the README.** It is referenced nowhere. Pure cruft.
2. **Delete `HDF5Write` and `GMM_EM`.** Both are dead (no caller / no command). Removing
   `HDF5Write` drops the **HDF5** dependency outright; removing `GMM_EM` shrinks the mlpack
   surface to just RangeSearch + DBSCAN.
3. **Drop CGAL.** It is already optional and only renders the Voronoi PDF; Delaunay itself
   runs on the vendored `delaunator`. Replace the two `boost::hash_combine` calls with a
   three-line local hash and **Boost disappears too** (it was otherwise only transitive).
   → removes **CGAL, Boost, GMP, MPFR** in one move. Cost: lose Voronoi-PDF output.
4. **Make `supervised-lda` opt-in.** It is already `#ifdef HAVE_LDAPLUSPLUS`; stop hardcoding
   `$HOME/git/supervised-lda` and default the build to LDA-off. `ldacreate`/`ldarun` become
   an optional add-on instead of a required prerequisite.

After Tier 0 the *required* heavy libraries drop from nine to roughly **mlpack, Armadillo,
Cairo, libtiff/jpeg** — with no loss beyond Voronoi-PDF and an unused HDF5 export.

### Tier 1 — the big one (moderate effort, high payoff)

5. **Replace mlpack with a single-header KD-tree (nanoflann).** mlpack is the most painful
   dependency and is used for exactly two things:
   - *radius/neighbor search* (`mlpack::RangeSearch`/`Range`) — nanoflann's `radiusSearch`
     and `knnSearch` are drop-in equivalents and header-only;
   - *DBSCAN* — ~50 lines on top of that KD-tree, or drop the `dbscan` command.

   This removes **mlpack + Armadillo + ensmallen + BLAS/LAPACK** — the single largest cut.
   Effort is concentrated in `cell_table_graph.cpp` (the 5 RangeSearch sites) plus the small
   `cell_table_mlpack.cpp`.

### Tier 2 — optional, feature trade-offs (later)

6. **Gate Cairo and libtiff/jpeg behind build flags** so a "core" build needs neither. This
   makes `png`/`plot`/`convolve` opt-in. The plotting could later move to a lighter backend
   (or be emitted as data for an external plotter).

## End state

A **core cyftools** that builds from a clean checkout with only a C++17 compiler, **Eigen
(header-only)** and **OpenMP** — everything else vendored or header-only (CYF I/O,
delaunator, knncolle/umappp, nanoflann). Plotting, TIFF, and LDA become optional modules.
That is the htslib-style portability target: the format + core analysis build anywhere, and
heavy/niche features are add-ons rather than prerequisites.

## Note on the sandbox

These libraries are third-party, not OS-provided. On the user's Mac they come from Homebrew;
in the Cowork sandbox nothing scientific is preinstalled and the package network is locked
down, so the full binary can't be linked here. The **core subset** (format + streaming:
`convert`, `view`, `cat`, `filter`, …) does build and run here — which is exactly the subset
this cleanup would turn into the default build.
