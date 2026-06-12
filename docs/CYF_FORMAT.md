# CYF Format Specification

**Version:** 1.0 (draft)
**Status:** Proposed. Defines the on-disk format for cyftools, replacing the implicit cereal encoding.
**Scope of v1:** the text format (`.cyf`), the typed data model, the self-describing header, the binary byte layout, and its **BGZF compression** container (`.bcyf`) — see §10. Indexing on top of BGZF virtual offsets is the one remaining deferred piece.

A companion reference implementation in pure Python lives at [`docs/reference/cyf.py`](reference/cyf.py); every byte offset in the worked example (§11) is generated and verified by it.

---

## 1. Design goals and lineage

CYF is to multiplex-imaging single-cell tables what **SAM/BAM** is to sequencing alignments. The format is deliberately modeled on SAM and inherits its best ideas, while fixing the things that make the current cereal-based encoding non-portable.

Concretely, CYF commits to five principles:

1. **A file is self-describing.** The header is the schema. A reader needs nothing but the file itself — no matching C++ struct, no compiled `#define`s — to know what every column and every flag bit means. This is the property the current format lacks.
2. **The byte layout is explicit and language-agnostic.** Every field has a documented type, width, and byte order, so a conformant reader can be written in Python, R, Rust, or Julia directly from this document (as `cyf.py` demonstrates).
3. **Two encodings, one model.** A human-readable **text** form (the SAM analog) and a compact **binary** form (the BAM analog) describe the same logical content and convert **losslessly** in both directions.
4. **Streaming-first.** Like SAM/BAM, a file is a header followed by an unbounded stream of records, processable in O(1) memory and composable through Unix pipes.
5. **Provenance is first-class.** Every transformation appends a program (`@PG`) record, preserving an auditable chain of how a file was produced.

Where CYF departs from BAM, it does so because a quantification table has a **fixed schema across all rows**, whereas a BAM record carries a variable set of self-typed optional tags. CYF therefore declares column types **once** in the header and stores **bare positional values** in each record — closer to a database row or an Arrow schema than to BAM's per-record type codes. This is more compact and faster to parse, at the cost of requiring the header to decode the records (which is always present).

---

## 2. Conceptual model

```
            ┌──────────────────────────────────────────────┐
  HEADER →  │ ordered list of tag records (the schema)      │
            │  VN  format/version             (like @HD)    │
            │  SA  sample / specimen          (like @RG)    │
            │  MA  marker  ─┐                 (like @SQ)    │
            │  CA  metadata ├─ define data columns, in     │
            │  GA  graph   ─┘   header order = the schema   │
            │  FL  flag-bit semantics  (NEW; self-describes │
            │                           cflag/pflag bits)   │
            │  PG  program / provenance       (like @PG)    │
            └──────────────────────────────────────────────┘
  RECORDS → cell, cell, cell, …   (positional typed values per the schema)
```

The ordered set of `MA` + `CA` + `GA` tags **is** the column schema. The *k*-th data value in every record corresponds to the *k*-th data tag in header order, encoded according to that tag's declared type.

---

## 3. Type system

Every data column declares a single-character **type code**. Codes are chosen to echo BAM's optional-field types.

| Code | Meaning | Binary encoding | Text rendering |
|------|---------|-----------------|----------------|
| `i` | signed 32-bit int | 4 bytes, little-endian | decimal |
| `I` | unsigned 32-bit int | 4 bytes, LE | decimal |
| `l` | signed 64-bit int | 8 bytes, LE | decimal |
| `L` | unsigned 64-bit int | 8 bytes, LE | decimal |
| `f` | 32-bit float (IEEE-754) | 4 bytes, LE | shortest round-trip decimal |
| `d` | 64-bit float (IEEE-754) | 8 bytes, LE | shortest round-trip decimal |
| `Z` | UTF-8 string | `u32` length + that many bytes | raw text (tabs/newlines escaped) |
| `A` | categorical / enum | `u32` index into the column's `LV:` level list | the level string |
| `B` | typed array | `u8` subtype code + `u32` count + count×subtype values | `<sub>:v0,v1,…` (e.g. `L:2,3,4`) |

Notes:

- **Markers (`MA`) are conventionally `f`.** They may declare any numeric type, but `f` matches the measured-intensity model and is the default if `TY:` is omitted.
- **`A` (categorical)** stores a small integer index; the human-readable levels live in the header (`LV:tumor,stroma,immune`). This is how cell-type calls, region labels, and cluster ids should be stored — compact on disk, self-describing in the header.
- **`B` (array)** carries variable-length per-cell vectors. Its subtype must be one of the fixed-width numeric codes (`i I l L f d`). This is the natural home for spatial-graph data: e.g. a `GA` column `knn_ids` of type `B`/`L` (neighbor cell ids) paired with `knn_dist` of type `B`/`f` (distances).
- **Missing values.** For `f`/`d`, missing is IEEE NaN. For integer and `A` columns, a column may declare `NA:<value>` in its header tag to reserve a sentinel; absent that, all values are significant. (v1 keeps NA handling deliberately minimal; this is an area expected to firm up in 1.x.)

---

## 4. The header

### 4.1 Grammar

The header is a sequence of **tag records**, one per line in the text form. Each tag record is:

```
@<CLASS> <TAB> <KEY>:<VALUE> <TAB> <KEY>:<VALUE> …
```

- `<CLASS>` is a two-letter uppercase code: `VN`, `SA`, `MA`, `CA`, `GA`, `FL`, `PG`.
- Each subsequent field is a two-or-more-letter `KEY`, a colon, then the value. Keys are unique within a record. Field order within a record is preserved on write but not semantically significant (except that `ID`, when present, conventionally comes first).
- `VALUE` may contain spaces but not raw tabs or newlines; those are escaped as `\t`, `\n`, `\\` in text and stored decoded in the binary header block.

This is the SAM header grammar, generalized so that field semantics are defined per class below rather than left as an opaque blob.

### 4.2 Tag reference

#### `@VN` — version / file head (required, must be first)
The CYF analog of SAM's `@HD`. Exactly one, first line.

| Key | Req | Meaning |
|-----|-----|---------|
| `VN` | yes | format version string, e.g. `1.0` |
| `SO` | no | sort order: `unsorted` (default), `coordinate`, `id` |
| `CW` | no | flag register width in bits; v1 fixes this at `64` |

#### `@SA` — sample / specimen
The analog of `@RG`. Zero or more.

| Key | Req | Meaning |
|-----|-----|---------|
| `ID` | yes | unique sample id, referenced elsewhere |
| `SM` | no | human sample name |
| `SL` | no | slide id |
| `PA` | no | panel name (e.g. `Orion`, `prostate-12plex`) |
| `DS` | no | free-text description |

#### `@MA` — marker (defines a data column)
The analog of `@SQ`. One per measured marker; order defines column position.

| Key | Req | Meaning |
|-----|-----|---------|
| `ID` | yes | marker name (e.g. `CD3`, `DAPI`) |
| `TY` | no | type code (default `f`) |
| `CH` | no | imaging channel |
| `EX`/`EM` | no | excitation / emission wavelength |
| `AB` | no | antibody / clone |
| `UN` | no | measurement units |
| `DS` | no | description |

#### `@CA` — cell annotation / metadata (defines a data column)
Derived or imported per-cell values that are not raw marker intensities.

| Key | Req | Meaning |
|-----|-----|---------|
| `ID` | yes | column name (e.g. `nucleus_area`, `celltype`) |
| `TY` | yes | type code |
| `LV` | if `TY:A` | comma-separated category levels; record stores the index |
| `NA` | no | reserved missing-value sentinel (non-float types) |
| `PG` | no | id of the `@PG` program that created this column |
| `DS` | no | description |

#### `@GA` — graph column (defines a data column)
A per-cell spatial-graph payload, normally a `B` array.

| Key | Req | Meaning |
|-----|-----|---------|
| `ID` | yes | column name (e.g. `knn_ids`) |
| `TY` | yes | type code, normally `B` |
| `KD` | no | graph kind: `knn`, `delaunay`, `radial` |
| `K` | no | neighbor count / parameter |
| `DS` | no | description |

#### `@FL` — flag-bit semantics (NEW)
Declares the meaning of a single bit in the `cflag` or `pflag` register. This **replaces the hardcoded `#define`s** (`PROSTATE_CD4`, `ORION_PANCK`, Gleason groups, …) that currently compile dataset-specific biology into the binary. With `@FL`, a file is panel-agnostic and a reader can label flags without the source.

| Key | Req | Meaning |
|-----|-----|---------|
| `RG` | yes | register: `cflag` or `pflag` |
| `BI` | yes | bit index, 0-based |
| `ID` | yes | symbolic name (e.g. `TUMOR`, `CD3pos`) |
| `DS` | no | description |

Example: `@FL␉RG:pflag␉BI:11␉ID:CD3pos` means "bit 11 of `pflag` denotes CD3-positive."

#### `@PG` — program / provenance
The analog of `@PG`, including the `PP` back-link that forms an ordered chain.

| Key | Req | Meaning |
|-----|-----|---------|
| `ID` | yes | unique invocation id |
| `PN` | no | program name (`cyftools`) |
| `VN` | no | program version |
| `CL` | no | full command line |
| `PP` | no | `ID` of the previous `@PG` (chain link) |
| `DS` | no | description |

### 4.3 Ordering rules

- `@VN` is first.
- The relative order of `MA`/`CA`/`GA` records defines the data-column order and therefore the record value order. It must be stable across a round-trip.
- `@PG` records are emitted in chain order (oldest first).
- Other classes may appear in any position after `@VN`; writers SHOULD group by class for readability.

---

## 5. The record (logical)

A record (one cell) has five **mandatory fields** followed by one value per data column:

| Field | Type | Notes |
|-------|------|-------|
| `id` | `u64` | cell identifier; value equal to the EOF sentinel (§8.3) is illegal |
| `x` | `f32` | x coordinate (see precision note) |
| `y` | `f32` | y coordinate |
| `cflag` | `u64` | cell/structure flag register; bit meanings from `@FL RG:cflag` |
| `pflag` | `u64` | phenotype flag register; bit meanings from `@FL RG:pflag` |
| *data…* | per schema | one value per `MA`/`CA`/`GA` tag, in header order |

**Coordinate precision.** v1 stores `x`/`y` as `f32`, matching the current model and keeping records compact. `f32` resolves to ≤ 1 µm out to ~16 million µm and to ~0.01 µm out to ~160 mm, which covers typical slides. Datasets needing more can carry high-precision coordinates as `d`-typed `CA` columns; a future minor version may add an `@VN` flag selecting `f64` mandatory coordinates.

**Flag width.** `cflag`/`pflag` are fixed at 64 bits, removing the current compile-time `USE_64_BIT` ambiguity. A file is a file regardless of how the writer was built.

---

## 6. Text format (the SAM analog)

A text CYF file is the header (the `@`-lines from §4) followed by one record per line. Each record line is **tab-delimited**:

```
id <TAB> x <TAB> y <TAB> cflag <TAB> pflag <TAB> data0 <TAB> data1 <TAB> …
```

- `cflag`/`pflag` render as unsigned decimal integers (their bit meanings are recoverable from `@FL`). A writer MAY offer hex (`0x…`) output, but decimal is canonical.
- Float fields use the shortest decimal string that round-trips to the same IEEE value.
- `Z` strings are emitted with tabs/newlines/backslashes escaped.
- `A` categoricals MAY render either as the level string or the integer index; the level string is canonical for human readability, the index is accepted on read.
- `B` arrays render as `<subtype>:v0,v1,…` (e.g. `L:2,3,4`), or an empty field for a zero-length array.

The text format is the interchange and debugging format; it is line-oriented, greppable, and diff-able.

---

## 7. Binary format (the BAM analog), v1 — uncompressed

All multi-byte integers and floats are **little-endian**. The file is a fixed preamble, then a stream of records, then an EOF marker.

### 7.1 Preamble

| Offset | Field | Type | Value / meaning |
|--------|-------|------|-----------------|
| 0 | `magic` | `char[4]` | `43 59 46 01` = `"CYF\x01"` |
| 4 | `version` | `u16` | format major version (1 for this spec) |
| 6 | `reserved` | `u16` | 0 |
| 8 | `l_header` | `u32` | byte length of the header text block |
| 12 | `header_text` | `char[l_header]` | the full text header (the `@`-lines), UTF-8, newline-separated, **stored verbatim** |
| 12+`l_header` | `n_dcol` | `u32` | number of data columns |
| … | `dcol_types` | `u8[n_dcol]` | the type code of each data column, in header order |

Storing the header as the same UTF-8 text the text format uses (exactly as BAM stores the SAM header verbatim) is what makes text↔binary conversion lossless and trivial. The redundant `dcol_types` array lets a reader decode records by offset arithmetic without parsing the header text, and lets a reader **self-check** that its header parse agrees with the writer (`cyf.py` asserts this).

### 7.2 Record stream

Each record is written back-to-back with no separator:

| Field | Type |
|-------|------|
| `id` | `u64` |
| `x` | `f32` |
| `y` | `f32` |
| `cflag` | `u64` |
| `pflag` | `u64` |
| each data value | encoded per its `dcol_types[k]` (see §3) |

Fixed-width columns make a record fixed-size; `Z`/`A`/`B` columns make it variable-size, but always self-delimiting given the schema.

### 7.3 End-of-file marker and truncation detection

A well-formed file ends with the 8-byte marker:

```
EOF_MARKER = 43 59 46 01 45 4F 46 00   ("CYF\x01EOF\x00")
```

A reader sits at a record boundary (it has just finished decoding a complete record using the known schema) and reads the next 8 bytes:

- if they equal `EOF_MARKER` → **clean end of file**;
- otherwise → those 8 bytes are the next record's `id`; continue.

To make this unambiguous, the `u64` value of the marker — **19790406162471235** — is a **reserved, illegal cell id**. Writers must never emit a cell with that id; the reference writer raises an error if asked to.

This directly fixes the current "EOF by catching a cereal exception" problem: a file that ends **without** the marker is a **truncation/corruption** and a reader must report it rather than silently treating it as success. (`cyf.py` raises `EOFError` in that case; see the verified truncation test in `docs/reference/example.py`.)

---

## 8. Versioning and compatibility

- The binary `version` (`u16`) and the text `@VN VN:` carry the **format** version, independent of the cyftools **program** version (which lives in `@PG VN:`).
- **Major** version bumps signal an incompatible byte layout; a reader MUST refuse a major version it does not implement.
- **Minor** additions (new optional header keys, new type codes used only when present) must be backward-compatible: an older reader can skip unknown header keys and SHOULD error only if it encounters a data column whose **type code** it does not understand.

---

## 9. Lossless conversion

Because the binary preamble stores the header as the identical UTF-8 text used by the text format, and because both encodings carry the same typed values, conversion is lossless in both directions:

```
text  ──parse header, parse rows──▶  binary
binary ──emit verbatim header, format rows──▶  text
```

A conformant tool MUST guarantee `text → binary → text` and `binary → text → binary` are identity transforms (modulo the canonical float rendering of §6). The reference codec demonstrates both directions.

---

## 10. The two file forms, compression, and indexing

CYF mirrors SAM/BAM's two-tier model:

| Tier | SAM/BAM | CYF | Extension | Compressed |
|------|---------|-----|-----------|------------|
| Text | SAM | text (§6) | `.cyf` | no — human-readable, tab-delimited |
| Binary | BAM | binary (§7) **wrapped in BGZF** | `.bcyf` | yes |

So the "uncompressed binary" of §7 is not a stored file form on its own — exactly as BAM never stores bare binary, the §7 byte stream is the *payload inside the BGZF blocks*. You store either the **text** form (`.cyf`) or the **BGZF binary** form (`.bcyf`).

- **BGZF container (implemented).** The §7 byte stream is wrapped, unchanged, in BGZF — blocked gzip with a `BC` extra field per block and the standard 28-byte EOF block — exactly as SAM's bytes are wrapped to form BAM. A `.bcyf` file therefore begins with the gzip magic `1f 8b`; decompressing it yields the §7 stream beginning with the `CYF\x01` magic. Because it is ordinary gzip, `.bcyf` is readable by `zcat`/`bgzip`/`tabix`. Readers auto-detect compression from the gzip magic, so the tools accept `.cyf`, `.bcyf`, (and, transitionally, bare §7 binary and the legacy cereal stream) interchangeably. BGZF also carries the *virtual file offsets* (block offset + in-block offset) needed for random access.
- **Index (`.cyi`, deferred).** With BGZF virtual offsets in place, a coordinate index (binning the slide in x/y, à la `.bai`/`.csi`) or a sample/id index will enable "fetch the cells in this region/sample" without a full scan. The `@VN SO:` sort-order field is the hook that makes an index meaningful. This is the remaining piece.

---

## 11. Worked example (verified)

The following is produced and verified by `docs/reference/cyf.py`. The header (10 tag records) and two cells:

```
@VN	VN:1.0	SO:unsorted
@SA	ID:S1	SM:tonsil_01	PA:Orion
@MA	ID:DAPI	TY:f	CH:1
@MA	ID:CD3	TY:f	CH:5
@CA	ID:nucleus_area	TY:i
@CA	ID:celltype	TY:A	LV:tumor,stroma,immune
@GA	ID:knn_ids	TY:B	KD:knn	K:3
@FL	RG:cflag	BI:0	ID:TUMOR
@FL	RG:pflag	BI:1	ID:CD3pos
@PG	ID:pg1	PN:cyftools	VN:2.0	CL:convert in.csv
1	100.5	200.25	1	2	12.0	3.5	42	2	L:2,3,4
2	101.0	202.0	0	0	8.0	0.0	37	1	L:1,4,5
```

The binary encoding of the same content is 484 bytes. Annotated hexdump of the preamble and first record (offsets in decimal):

```
   0: 43 59 46 01                       magic 'CYF\x01'
   4: 01 00                             version u16 = 1
   6: 00 00                             reserved u16
   8: 2d 01 00 00                       l_header u32 = 301
  12: <301 bytes of UTF-8 header text>
 313: 05 00 00 00                       n_dcol u32 = 5
 317: 66 66 69 41 42                    dcol type codes = "ffiAB"

 --- first record ---
 322: 01 00 00 00 00 00 00 00           id u64 = 1
 330: 00 00 c9 42                       x f32 = 100.5
 334: 00 40 48 43                       y f32 = 200.25
 338: 01 00 00 00 00 00 00 00           cflag u64 = 1   (bit0 TUMOR set)
 346: 02 00 00 00 00 00 00 00           pflag u64 = 2   (bit1 CD3pos set)
 354: 00 00 40 41                       DAPI  f32 = 12.0
 358: 00 00 60 40                       CD3   f32 = 3.5
 362: 2a 00 00 00                       nucleus_area i32 = 42
 366: 02 00 00 00                       celltype A idx = 2  -> "immune"
 370: 4c                                knn_ids B subtype 'L' (u64)
 371: 03 00 00 00                       knn_ids B count = 3
 375: 02 00 00 00 00 00 00 00           knn_ids[0] = 2
      …                                 knn_ids[1], knn_ids[2], then record 2, then EOF_MARKER
```

Note how the categorical `celltype` stores just the index `2` on disk while the header's `LV:tumor,stroma,immune` supplies the label, and how the `knn_ids` graph column carries a self-describing typed array — both impossible to express in the current all-`float` model.

---

## 12. Migration from the current cereal format

The present format is `cereal::PortableBinary(CellHeader)` followed by `cereal::PortableBinary(Cell)*`, where `Cell = {u64 id, cy_uint cflag, cy_uint pflag, f32 x, f32 y, vector<f32> cols}`. The mapping to CYF v1 is mechanical:

- The cereal `Tag` vector maps onto `@`-line tag records; the existing `MA/CA/GA/PG/VN/SA` types carry over unchanged, with `FL` added to absorb the `#define` flag table.
- Each `Cell.cols[i]` becomes the *i*-th data column, declared `TY:f` (or a richer type where the column is known to be integer/categorical/graph).
- `cy_uint` flags are written as fixed `u64`, ending the `USE_64_BIT` ambiguity.
- A one-time `cyftools convert --from cereal` path can read the old stream (via the existing C++ structs) and emit CYF, stamping a `@PG` migration record.

A practical first implementation step is to add a CYF reader/writer alongside cereal and make `view`/`convert` able to round-trip both, then flip the default once parity is proven against the reference codec.

---

## Appendix A — constants

| Constant | Value |
|----------|-------|
| Binary magic | `43 59 46 01` (`"CYF\x01"`) |
| EOF marker (8 bytes) | `43 59 46 01 45 4F 46 00` (`"CYF\x01EOF\x00"`) |
| Reserved illegal cell id | `19790406162471235` (the EOF marker as LE `u64`) |
| Byte order | little-endian throughout |
| Mandatory record fields | `id:u64, x:f32, y:f32, cflag:u64, pflag:u64` |

## Appendix B — tag classes at a glance

| Class | SAM analog | Defines a data column? | Purpose |
|-------|-----------|------------------------|---------|
| `VN` | `@HD` | no | format version, sort order |
| `SA` | `@RG` | no | sample / specimen |
| `MA` | `@SQ` | yes | marker (intensity) column |
| `CA` | — | yes | derived/metadata column |
| `GA` | — | yes | spatial-graph column |
| `FL` | — | no | self-describing flag-bit semantics (replaces `#define`s) |
| `PG` | `@PG` | no | program / provenance chain |

---

*This is a v1.0 **draft** for review. The most consequential open questions are: final file extensions (§10), NA semantics for non-float columns (§3), and whether to promote `x`/`y` to `f64` (§5). Everything in §7 (the binary layout) is exercised end-to-end by the reference codec.*
