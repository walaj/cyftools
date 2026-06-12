# cyftools — Architecture Review & Modernization Roadmap

*A code-level review of the current codebase, with a roadmap toward (1) a rigorous SAM-inspired format specification, (2) htslib-style portability, and (3) a serialization layer that replaces cereal with an explicit, language-agnostic, BGZF-backed format.*

---

## 1. What cyftools is

cyftools is a command-line toolkit for wrangling **multiplex immunofluorescence (mIF / CyCIF) single-cell quantification tables** — the tabular output you get after segmenting cells in a high-plex image and measuring per-marker intensities. The design is explicitly modeled on **samtools/htslib**: a single binary exposing ~52 subcommands that each read a stream of records, transform them, and emit a stream of records, so commands compose through Unix pipes (`-` means stdin/stdout throughout).

The conceptual mapping to genomics is the heart of the project:

| Genomics (SAM/BAM) | cyftools | Meaning |
|---|---|---|
| A read / alignment record | A **`Cell`** record | One row = one cell |
| Reference position (`POS`) | `x`, `y` | Spatial location in the slide |
| FLAG bitfield | `cflag`, `pflag` | Cell-level and phenotype-level bit flags |
| `@HD VN:` | `VN` tag | Format version |
| `@SQ` (sequence dictionary) | `MA` tags | The measured markers (the "dimensions") |
| `@RG` (read group) | `SA` tags | Sample identity |
| `@PG` (program chain) | `PG` tags | Provenance — every command appends one |
| Optional `TAG:TYPE:VALUE` columns | `MA`/`CA` data columns | Per-cell marker & metadata values |

That mapping is already well thought out and is the right foundation. The work ahead is mostly about **making the format explicit, portable, and indexable** rather than reinventing the model.

---

## 2. Repository map

```
src/                    The C++ engine (only this is compiled into `cyftools`)
  cyftools.cpp          MAIN. Arg parsing + dispatch for all 52 subcommands (~3,700 lines)
  cysift.cpp            LEGACY. Old name for the same program. NOT in the Makefile. Dead. (~3,200 lines)
  cell_row.{h,cpp}      The Cell record (the on-disk row)
  cell_header.{h,cpp}   CellHeader + Tag (the SAM-style header)
  cell_column.h         In-memory columnar types (NumericColumn<T>); much is commented out
  cell_table*.cpp       The in-memory CellTable + analytics (graph, delaunay, LDA, mlpack)
  cell_processor.{h,cpp} The streaming processor framework (the "samtools view" engine)
  cell_selector.*       Flag/field selection logic for `filter`
  cell_flag.*           CellFlag bit helpers
  csv.h                 Vendored fast-cpp-csv-parser (CSV ingest)
  polygon.*             ROI polygons
  tiff_*.{h,cpp}        A hand-rolled TIFF reader/writer (image ops, convolve)
  color_map.*, lda_model.h, delaunator.hpp, ...
  Makefile              Single hand-written Makefile; builds `cyftools`
external/               Git submodules: cereal, umappp, CppIrlba, knncolle, CppKmeans
R/ python/ matlab/ omero/ html/ Rmd/   Bindings, viewers, and helper scripts
palettes/               Color palettes for plotting
test/                   Two CSVs. (Effectively no test suite.)
```

Two structural facts worth acting on early:

- **`cysift.cpp` is dead code.** It is the program's former name, is not referenced by the Makefile, and duplicates ~3,200 lines of `cyftools.cpp`. It should be deleted (or git-archived) so there is one source of truth.
- **There is no format specification anywhere in the repo** — no spec doc, no magic number, no version gate at the byte level. The format exists only as an implicit consequence of the C++ struct layout plus cereal. This is the central risk and the thing most worth fixing first.

---

## 3. The data model

### 3.1 The record — `Cell` (`cell_row.h`)

```cpp
class Cell {
  uint64_t id;            // cell id
  cy_uint  pflag;         // phenotype flags (bitfield)
  cy_uint  cflag;         // cell/structure flags (bitfield)
  float    x, y;          // spatial coordinates
  std::vector<float> cols; // ALL marker + metadata values, as float
  // serialize: ar(id, cflag, pflag, x, y, cols);
};
```

Every measured value — marker intensities *and* derived metadata — lives in a single `vector<float> cols`. The meaning of each slot is positional: `cols[i]` corresponds to the *i-th* data tag (`MA` or `CA`) in the header. There are no per-record column names; the header is the schema and the record is a bare value vector. This is exactly the BAM philosophy (records are compact and positional; the header carries the dictionary), and it's good for size and speed.

`cy_uint` is `uint64_t` when `USE_64_BIT` is defined (it is, by default in the Makefile), otherwise `uint32_t`. So the flag width is a **compile-time** decision — a file written by a 64-bit build and one written by a 32-bit build are not the same format, and nothing in the file records which one produced it. (See §6.)

A design observation: forcing everything to `float` means integer ids/counts and categorical metadata are stored lossily as floats. SAM/BAM instead carries typed optional fields (`i`, `f`, `Z`, `A`, `B`). If the metadata is ever more than "a number," the single-`float` model will constrain you. The commented-out `StringColumn` / `GraphColumn` / `FlagColumn` in `cell_column.h` suggest you already hit this wall once.

### 3.2 The header — `CellHeader` + `Tag` (`cell_header.h/.cpp`)

`CellHeader` is just a `std::vector<Tag>`. A `Tag` is a SAM-style two-letter record:

```cpp
class Tag {
  uint8_t type;       // MA, CA, GA, PG, VN, SA
  std::string id;     // e.g. "CD3"
  std::string data;   // remaining "FIELD:VALUE FIELD:VALUE" payload
  int i;              // sort/order index (column position for MA/CA; temporal order for PG)
};
```

Tag types:

- **`MA` (marker)** — a measured marker; defines a data column.
- **`CA` (meta / cell-annotation)** — a derived/metadata column (also a data column).
- **`GA` (graph)** — spatial-graph metadata.
- **`PG` (program)** — provenance; every command appends one (like `@PG`).
- **`VN` (version)** — format/program version.
- **`SA` (sample)** — sample identity (like `@RG`).

The text rendering is deliberately SAM-like: `@MA\tID:CD3\t...`, `@PG\t...`. The parser (`Tag::Tag(const std::string&)`) splits on whitespace, reads the two-letter type, pulls out `ID:`, and stuffs everything else into `data` as a raw string. `MA`/`CA` tags are "data tags" and their order (`i`) is what binds `cols[i]` to a column name.

This is a clean, faithful adaptation of the SAM header. The main gaps are that the grammar is informal (the `data` field is an unparsed blob), and validation is a stub (`CellHeader::validate()` just `return true;`).

### 3.3 The in-memory table — `CellTable` (`cell_table.h`)

For operations that need the whole table at once (spatial KNN, Delaunay/Voronoi, DBSCAN, Moran's I, UMAP, LDA, TIFF convolution), records are loaded into a **columnar** structure:

```cpp
std::unordered_map<std::string, FloatColPtr> m_table;  // named data columns
FloatColPtr m_x_ptr, m_y_ptr;                          // dedicated coordinate columns
IntColPtr   m_pflag_ptr, m_cflag_ptr;
IDColPtr    m_id_ptr;
```

So cyftools has **two representations**: a streaming row form (`Cell`, for the pipe-friendly transforms) and an in-memory column form (`CellTable`, for analytics). That's a sound split — it mirrors "samtools view" (streaming) vs. a stats tool that buffers — but right now the columnar side is entangled with the I/O code, and `cell_column.h` carries a lot of dead, commented-out class machinery.

---

## 4. The on-disk format, as it actually exists today

This is the precise current reality, because the roadmap depends on naming it exactly:

```
┌─────────────────────────────────────────────┐
│  cereal::PortableBinaryOutputArchive stream  │
│                                              │
│  [ CellHeader ]   ← one cereal object        │  = vector<Tag>, each Tag = {u8 type, string id, int i, string data}
│  [ Cell ]         ← cereal object            │  = {u64 id, cy_uint cflag, cy_uint pflag, f32 x, f32 y, vector<f32> cols}
│  [ Cell ]                                    │
│  [ Cell ]                                    │
│   ...                                        │
│  (stream ends; reader stops on cereal EOF    │
│   exception)                                 │
└─────────────────────────────────────────────┘
```

Concretely:

- **No magic number.** The file begins immediately with cereal's encoding of the header. There is no way to identify a cyf file by inspection, and no way to detect "this isn't a cyf file" except that deserialization throws.
- **No explicit version in the bytes.** A `VN` tag may sit *inside* the header, but the binary container/layout itself is unversioned. A change to the `Cell` struct silently breaks every existing file with no detectable signal.
- **No compression.** It's raw cereal portable-binary. SAM's binary sibling BAM is BGZF-compressed; cyf is not compressed at all.
- **No index.** Reading is strictly sequential, front to back. There's no equivalent of `.bai`/`.csi`, so "give me cells in this spatial region" always means a full scan.
- **EOF is signaled by exception.** The read loop (`CellTable::StreamTable`) calls `inputArchive(cell)` in a `try` and treats *any* `cereal::Exception` as end-of-stream (`break`). That conflates clean EOF with a truncated/corrupt file — a partial write reads as a successful short file.
- **The format is only readable by C++ + cereal**, and only when the reading program is compiled with a `Cell`/`Tag` struct identical to the writer's. cereal's "portable" binary fixes endianness and integer sizes, but it does **not** give you a documented byte layout that an independent Python/R/Rust reader can implement. The schema *is* the C++ code.

The write path is the mirror image. Each processor independently calls `SetupOutputStream()` and then `(*m_archive)(m_header)` before streaming cells. That header-write boilerplate is duplicated **24 times** across `cell_processor.cpp` — a strong signal it belongs in the `CellProcessor` base class.

---

## 5. The I/O & command architecture

### 5.1 The streaming engine

The pipe-friendly commands use a small framework in `cell_processor.h`:

```cpp
class CellProcessor {
  virtual int ProcessHeader(CellHeader&) = 0;   // called once
  virtual int ProcessLine(Cell&) = 0;           // called per record
  // return codes steer the driver: WRITE_CELL / SAVE_CELL / NO_WRITE_CELL / ...
};
```

`CellTable::StreamTable()` is the driver: open input (file or stdin) → deserialize header → `ProcessHeader` → loop deserializing `Cell`s → `ProcessLine` → optionally `OutputLine`. Each subcommand is a `CellProcessor` subclass (`ViewProcessor`, `FilterProcessor`, `CatProcessor`, `HeadProcessor`, `CleanProcessor`, …). This is precisely the htslib read/modify/write loop, and it's the strongest part of the codebase: cheap to add commands, naturally O(1) memory, naturally composable.

`CerealProcessor` is the one **ingest** path — it's a `LineProcessor` (operates on text lines, not `Cell`s) that reads a CSV with the vendored `csv.h`, builds the header, and writes the binary stream. That's the `convert` boundary between "human CSV" and "binary cyf."

### 5.2 The command surface (52 subcommands)

They cluster cleanly, which will matter when we split a library from the app:

- **I/O & format:** `convert`, `view`, `cat`, `head`, `reheader`, `info`, `check`, `count`, `cellcount`, `summary`, `synth`, `debug`
- **Row/column manipulation:** `cut`, `clean`, `trim`, `filter`, `sampleselect`, `subsample`, `crop`, `roi`, `select`
- **Transforms:** `log10`, `divide`, `rescale`, `magnify`, `offset`, `flip`, `scramble`, `hallucinate`, `flagset`, `pheno`, `mean`
- **Spatial / graph:** `radialdens`, `delaunay`, `dbscan`, `tls`, `island`, `margin`, `annotate`, `dist`, `umap`
- **Numeric / stats:** `pearson`, `jaccard`, `histogram`
- **ML:** `ldacreate`, `ldarun`
- **Imaging / output:** `png`, `plot`, `scatter`, `convolve`

The first two clusters are pure format/streaming work with essentially no heavy dependencies. The spatial/ML/imaging clusters are where the heavy dependencies (mlpack, CGAL, Eigen, Armadillo, supervised-lda, libtiff, Cairo, HDF5) come in. That dependency gradient is the natural seam for the htslib-style "thin format library vs. heavy analysis app" split (§7.2).

---

## 6. Strengths

The bones are good, and several decisions are genuinely right:

The **samtools mental model** — header + record stream + composable subcommands over pipes — is the correct abstraction for this domain, and it's already implemented coherently. The **SAM-style tagged header** (`MA/CA/GA/PG/VN/SA`) is a faithful, extensible adaptation, and the **PG provenance chain** (every command stamps a `PG` tag) is a real strength most bioinformatics file formats lack. The **streaming processor framework** keeps memory flat and makes new commands cheap. The **two-tier representation** (streaming rows for transforms, in-memory columns for analytics) is the right separation of concerns even if it needs tidying. And the conceptual groundwork for a spec is largely *already done* — this review's job is mostly to make explicit what the code already implies.

---

## 7. Weaknesses, risks, and the path forward

### 7.1 Goal 1 — A rigorous, SAM-inspired format specification

This is the highest-leverage next step, and it should come *before* any serialization rewrite, because the spec is what the rewrite implements.

**Problems to resolve in the spec:**

- **No magic / no version / no self-description.** Adopt BAM's approach: a magic string (e.g. `CYF\x01` for binary, a `@HD VN:` line for text) and an explicit format-version integer in the bytes, independent of the program version.
- **Compile-time format variance.** `cy_uint` being `uint32_t` *or* `uint64_t` depending on `USE_64_BIT` means flag width is not fixed by the format. The spec must fix field widths (pick `uint64` flags and commit), so a file is a file regardless of build.
- **Hardcoded biology.** `cysift.h` compiles in project-specific semantics: `#define PROSTATE_CD4 1024`, `ORION_PANCK`, Gleason grade groups, etc. Flag-bit *meanings* are domain data, not language constants — they belong in header tags (or a sidecar dictionary), not in the binary. The spec should define *how flag semantics are declared in the header*, and the `#define` soup should migrate there. This is also a portability blocker: any non-prostate, non-Orion dataset is currently second-class.
- **Informal tag grammar.** `Tag::data` is an unparsed string blob and `validate()` is a stub. The spec should define each tag's required/optional fields (e.g. `MA` requires `ID`; what else? units? excitation? threshold?), the `cols[i] ↔ i-th data tag` binding rule, ordering guarantees, and escaping.
- **Typing.** Decide whether records stay all-`float` or gain typed columns (BAM-style `i/f/Z/A/B`). At minimum the spec should reserve space for typed metadata so you're not boxed in.

**Recommended shape:** specify a **text format** (the SAM analog — human-readable, diff-able, the lingua franca) and a **binary format** (the BAM analog — compact, BGZF-framed, indexable), with **lossless 1:1 conversion** between them. That duality is exactly what made SAM/BAM succeed, and cyftools already has both halves in embryo (CSV ingest + cereal binary).

### 7.2 Goal 2 — htslib-style portability

htslib is portable because it is a **small, self-contained library with a stable, documented format and a narrow C API**, cleanly separated from the *samtools application* that uses it. cyftools today is the inverse: one binary that fuses the format with a very heavy analysis stack (mlpack, CGAL, Eigen, Armadillo, HDF5, Cairo, libtiff, plus `supervised-lda` **hardcoded to `$HOME/git/supervised-lda`**). You cannot read a cyf file without, in principle, dragging all of that along.

The portability program, in order of value:

1. **Split `libcyf` (format) from `cyftools` (app).** Pull the record, header, reader/writer, index, and the streaming `CellProcessor` framework into a dependency-light library — ideally depending on nothing heavier than a compression lib. Everything in the "I/O & format" and "row/column manipulation" command clusters (§5.2) can live on top of just `libcyf`. The spatial/ML/imaging commands link the heavy stack separately. This mirrors `libhts` vs `samtools`.
2. **Put it behind a stable C ABI.** A C API (not a C++-template API) is what lets Python (ctypes/cffi), R, Rust, and Julia bind without ABI pain — this is precisely how htslib reaches pysam, Rsamtools, rust-htslib. Your existing `R/`, `python/`, `matlab/`, `omero/` directories become thin wrappers over one ABI instead of bespoke reimplementations.
3. **Go on a dependency diet.** Make heavy deps optional and feature-gated (the Makefile already has `HAVE_*` flags — formalize them), and remove hardcoded home-directory paths. Consider CMake + proper `find_package`/vendoring so a clean checkout builds the *format* library with zero exotic prerequisites.
4. **Document the byte layout** (falls out of Goal 1) so a third party can write a conformant reader from the spec alone, without linking your code. That, more than anything, is what "portable like htslib" means.

### 7.3 Goal 3 — Off cereal, toward an explicit format + BGZF

Your instinct is right, with one important clarification: **cereal and BGZF solve different problems, and BAM uses an explicit version of both.**

- **cereal** is *serialization*: how a `Cell` struct becomes bytes. Its weakness here isn't speed — it's that the byte layout is **defined by C++ code, not by a document**. Only a C++ program compiled with the identical struct can read it. That is fundamentally incompatible with "portable like htslib."
- **BGZF** is the *container/compression*: blocked gzip with virtual file offsets, which is what enables both compression *and* random access + indexing. cyf has neither today.

BAM = **explicit little-endian record layout** (the serialization) **wrapped in BGZF** (the container). So the migration is two distinct moves, and they're independent:

1. **Replace cereal with an explicit, documented binary encoding.** Define exact field order, widths, and endianness for the header block and the record (`id:u64, cflag:u64, pflag:u64, x:f32, y:f32, ncol:u32, cols:f32[ncol]`, etc.), with a leading magic + version. Hand-rolled fixed-width read/write is straightforward, removes the cereal submodule, and — crucially — becomes implementable in any language directly from the spec. This is the move that actually delivers portability; do it even if you defer compression.
2. **Adopt BGZF as the container.** Either vendor `bgzf.c`/`bgzf.h` from htslib (small, permissively licensed) or link `libhts`. BGZF gives you gzip-compatible compression *and* the virtual-offset machinery to build a spatial or sample index later — the cyf analog of `.bai`/`.csi`. Random spatial queries (§4) become possible instead of full scans.

A safe sequencing: **(a)** freeze the explicit binary layout behind a magic+version and a thin reader/writer that replaces cereal, validating round-trip equivalence against current files; **(b)** wrap that byte stream in BGZF; **(c)** add an index format once virtual offsets exist. Steps are independently shippable, and step (a) alone unblocks cross-language readers.

While you're in there, fold the **24× duplicated header-write boilerplate** into the base class, fix **EOF-by-exception** (a sized/sentinel record stream distinguishes clean EOF from truncation), and delete the **dead `cysift.cpp`** and commented-out column classes.

---

## 8. Suggested near-term sequence

1. **Delete the legacy/dead code** (`cysift.cpp`, commented-out classes) so there's one source of truth. (Low risk, high clarity.)
2. **Write the format specification** (`docs/CYF_FORMAT.md`) — text + binary, magic, version, header grammar, record layout, flag-semantics-in-header. This is the anchor document everything else implements.
3. **Extract `libcyf`** — the format/record/header/streaming core, dependency-light, behind a C ABI.
4. **Replace cereal** with the explicit binary encoding defined by the spec; round-trip-test against existing files.
5. **Add BGZF** as the container; then an index.
6. **Build out a real test suite** (currently two CSVs) — golden files, round-trip, and cross-language reader conformance.

Steps 1–2 are pure design/cleanup and unblock everything else; 3–5 are the portability/serialization payoff; 6 should grow alongside all of them.

---

*Generated as an initial deep-dive review. Section 7 is the part to argue with — the sequencing and the "explicit layout vs. BGZF are separate moves" framing are the load-bearing claims.*
