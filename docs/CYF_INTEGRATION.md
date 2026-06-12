# CYF integration (parallel-to-cereal phase)

This describes how the CYF reader/writer is wired into cyftools alongside cereal,
how to switch formats, and what remains before cereal can be removed. The format
itself is specified in [CYF_FORMAT.md](CYF_FORMAT.md); a language-agnostic
reference implementation is in [reference/cyf.py](reference/cyf.py).

## What was added

- **`src/cyf_field.h`** — the typed value model: `CyfValueType` (`i/I/l/L/f/d/Z/A/B`)
  and `CyfField` (int / float / string / category / array). Header-only, dependency-light.
- **`src/cyf_io.{h,cpp}`** — the CYF v1 codec (`cyf::CyfWriter`, `cyf::CyfReader`,
  `detectCyf`, header text render/parse, and the typed `CyfRow` path). Self-contained
  and cereal-free; this is the seed of a future portable `libcyf`.
- **`src/cell_archive.h`** — `OutArchive` / `InArchive`, thin wrappers that expose
  cereal's exact `operator()` interface but dispatch to either cereal **or** CYF.
  This is the single seam that lets both formats coexist without rewriting call
  sites. `InArchive` auto-detects the input format from the stream's magic bytes.
- **`src/cyf_test.cpp`** + **`make cyf-test`** — a dependency-light self-test
  (codec round-trip, header round-trip, truncation detection, and round-trips
  through both adapter backends). Builds with only cereal.

## What changed in existing code

The change is deliberately small and additive — cereal remains the default:

- `cell_processor.h` — `CellProcessor::m_archive` and `CerealProcessor::m_archive`
  are now `std::unique_ptr<OutArchive>`; `SetupOutputStream()` constructs an
  `OutArchive`. `OutputLine()` and the ~24 `(*m_archive)(m_header)` sites are
  unchanged (they go through `OutArchive::operator()`).
- `cell_processor.cpp` — `CerealProcessor::ProcessHeader` constructs an `OutArchive`.
- `cell_table.h` / `cell_table.cpp` — `CellTable::m_archive` is now `OutArchive`;
  `SetupOutputWriter`/`OutputTable` use it. `StreamTable` now reads through
  `InArchive`, so it transparently accepts cereal **or** CYF input and reports
  truncation instead of treating it as EOF.
- `Makefile` — adds `cyf_io.cpp` to the build and a `cyf-test` target.
- `cell_header.h` — adds `appendRawTag()` (order-preserving tag append used by the
  codec when rebuilding a header from on-disk text).

`cyftools.cpp::cyftools_cat` still uses cereal directly, but it is **dead code**
(the `cat` command dispatches to `catfunc` → `CatProcessor` → `StreamTable`, which
is wired). It can be deleted.

## How to use it

Reading is always automatic — `StreamTable`/`InArchive` detect the format per file,
so existing cereal files keep working and mixed inputs are fine.

Writing is selected by an environment variable (no CLI changes needed yet):

```sh
# default: writes the legacy cereal stream
cyftools convert in.csv - | cyftools filter ... - out.cyf

# opt in to CYF output for the whole pipeline
export CYFTOOLS_FORMAT=cyf
cyftools convert in.csv - | cyftools filter ... - out.cyf
```

`cyf::useCyfOutput()` reads `CYFTOOLS_FORMAT` once and is mutable, so a future
`--cyf` flag can set it directly from `cyftools.cpp` arg parsing.

Run the self-test anytime:

```sh
cd src && make cyf-test
```

## Typed columns

The format, header schema, and codec are now **fully typed**:

- `cyf_field.h` provides `CyfValueType` and `CyfField` (int / float / string /
  categorical / array).
- The header is a typed schema: `Tag::ValueType()` parses `TY:`, `Tag::CategoryLevels()`
  parses `LV:`, and `CellHeader::ColumnTypes()` returns the per-column types. Columns
  with no `TY:` default to float, so **legacy all-float headers are unchanged** (and
  still byte-identical to before — verified).
- `CyfWriter::writeRow` / `CyfReader::readRow` encode/decode a typed `CyfRow` driven by
  that schema. Typed files are byte-for-byte interoperable with the Python reference
  codec (verified for int/string/categorical/array), and each side reads the other's
  typed files.

The legacy `writeCell(Cell)` / `readCell(Cell)` path remains for the all-float model:
`writeCell` requires an all-float schema; `readCell` converts numeric typed columns to
float and only errors on a String/Array column (which a `float`-only `Cell` can't hold).

**What is still all-float: the in-memory app.** `Cell::cols` is still `vector<float>` and
`CellTable` still stores float columns, so the 52 commands continue to produce/consume
float data. Making the commands actually *emit* typed columns (e.g. CSV `convert`
inferring int/categorical, the graph ops writing `B` arrays) is the next phase:
propagate `CyfField` into `Cell::cols` (~70 sites) and the `CellTable` store (~97 sites),
behavior-preserving for all-float data. Until then, the typed machinery is exercised via
`writeRow`/`readRow` and the self-test, and the on-disk format already supports it end to
end.

**Self-describing flags (`@FL`) are now implemented.** The `FL` tag type is in the model
(`Tag::FL_TAG`, parsed/rendered/round-tripped), `CellHeader::GetFlagTags()` exposes them,
and `src/cyf_flags.h` declares the structural `cflag` vocabulary (Tumor=bit0 … BuildGraph=bit21)
as `@FL` tags that `convert` now stamps into every new file's header — so a reader can
interpret `cflag` without the source. The code still uses the `#define`s internally; the
remaining step is moving the panel-specific *phenotype* (`pflag`) bits (ORION_*/PROSTATE_*)
to per-dataset `@FL pflag` declarations. Graph-as-`B`-column still awaits the in-memory migration.

## Path to removing cereal

1. Flip the default: make `useCyfOutput()` return `true` (or wire a `--cyf` flag and
   default it on), regenerate any test fixtures as CYF.
2. Bake-in period: reads still auto-detect, so legacy cereal files remain readable.
3. Delete the cereal branches in `cell_archive.h` (the adapters collapse onto
   `CyfWriter`/`CyfReader`), drop the `<cereal/...>` includes and the `serialize()`
   methods on `Cell`/`CellHeader`/`Tag`, and remove the `external/cereal` submodule.
4. Then layer on the deferred extensions from the spec: typed columns, `@FL`, and —
   when ready — the BGZF container + index.
