# Cohort summary JSON

`cyftools cohort [options] <in1.byf> <in2.byf> ...` reads a set of CYF tables and
writes a single JSON document that drives a cohort-level viewer (densities,
plots, comparisons across samples). It is a derived, read-only summary: the bulk
per-cell data stays in the `.byf` files — only **precomputed compartment areas**
and a **joint flag histogram** land here, which together let a viewer answer
"cells/mm² of phenotype X in compartment Y" with a quick lookup, for *any* X.

```sh
cyftools cohort -o cohort.json *.byf       # write a file (progress -> stderr)
cyftools cohort *.byf > cohort.json        # or stream to stdout
cyftools cohort -t 8 -o cohort.json *.byf  # paint with 8 threads
cyftools cohort -n "Orion CRC" -o c.json *.byf   # label the cohort (top-level "name")
```

`-n`/`--name <str>` is an optional human-readable label echoed to the JSON's top-level
`"name"` field (empty string when omitted); the `docs/cyf_cohort.html` viewer shows it,
over the loaded file name, in its header after a load.

Each input is streamed independently (cells are never held across files), so
memory is bounded to one file at a time. Progress is logged to stderr.

## The model: two flag registers

Every cell carries two 64-bit flag registers — `cflag` and `pflag`. They are the
same *kind* of thing (sets of bits); the convention is only that **cflag bits
mark regions/compartments** (Tumor, Margin, TLS, …) and **pflag bits mark
phenotypes** (CD3, CD8, …). A *compartment* is one cflag bit's region, or the
synthetic **`All`** (every cell, ignoring cflag). A *cell of interest* is any
boolean predicate over pflag bits (e.g. `CD3 AND CD8 AND NOT FOXP3`).

Density is then just:

```
density(compartment, pflagPredicate) = count(compartment, pflagPredicate)
                                        ─────────────────────────────────
                                              area_mm2[compartment]
```

The numerator is computed in the viewer from the joint histogram; the
denominator is **precomputed here** by painting (below). Compartments are limited
to single cflag bits + `All`, because those are the regions whose areas were
painted — see "Limitations".

## The "painting" area estimate

A compartment's area is estimated by **painting**: every cell in it stamps a
filled disc (default radius **20 µm**, `-r`) onto a raster grid (default **1 µm**
cells, `-p`/`--grid`); the area is the count of touched grid cells × cell². Overlapping
discs paint shared cells once, so the result is the **union** area of the discs
— a compact approximation of the tissue those cells occupy.

- `-r`, `--radius` (µm): larger fills gaps more aggressively → larger area, lower density.
- `-p`, `--grid` (µm): the area-raster resolution — the internal rasterization cell
  size, **not** an image pixel size. It trades accuracy for speed/memory; area is
  insensitive to it (1.6 M cells: grid=1 → 198.8 mm² vs grid=4 → 200.2 mm², ≈0.7 %).
  `--grid 2` is ≈4× faster and ≈4× less memory. Painting is multithreaded (`-t`) and
  the result is independent of thread count.

## Coordinate units (per file, not per cohort)

All painting is done in **microns**. Each input is converted to microns using **its
own** `@HD` calibration — `UN:micron` is used as-is; `UN:pixel` is scaled by that
file's `MP` (microns per pixel). The calibration therefore lives in each file, never
at the cohort level, so a cohort can mix differently-scaled inputs safely. A file in
`UN:pixel` with a missing/invalid `MP` cannot be converted and is **skipped** (add it
with `cyftools addtag in.byf out.byf -t HD -f MP:<µm-per-px> -f UN:pixel`). A legacy
file with no `UN` tag is assumed to already be in microns.

## Computing a density in the viewer

`count(compartment, pflagPredicate)` is a masked sum over `flag_histogram` rows:

```js
function count(hist, compBit, pred) {        // compBit = -1 means "All"
  let n = 0;
  for (const [cflagStr, pflagStr, c] of hist) {
    const cflag = BigInt(cflagStr), pflag = BigInt(pflagStr);
    if (compBit >= 0 && !(cflag & (1n << BigInt(compBit)))) continue;
    if (pred(pflag)) n += c;
  }
  return n;
}
// e.g. CD3 AND CD8 in the Tumor compartment:
const cd3 = bitOf("CD3"), cd8 = bitOf("CD8"), tumor = bitOf("Tumor");
const n   = count(hist, tumor, pf => (pf & (1n<<BigInt(cd3))) && (pf & (1n<<BigInt(cd8))));
const dens = n / areaOf(tumor);              // cells / mm^2
```

`bitOf`/`areaOf` come from `pflag_bits`/`cflag_bits` and `compartments`. Single-
marker densities are just the one-bit case of the same sum, so no separate
per-marker table is stored.

## Stratifying by a categorical column (`catcols`)

A CYF file may carry categorical `@CA` columns declared `TY:A` — labels such as an
`IC` immune-cluster id that are *not* numbers to average (see `docs/CYF_FORMAT.md` §3;
make one with `cyftools settype in out --set IC:A`). For each such column, cohort emits
a `catcols[]` entry with the column's `levels` (the `LV:` labels, in level-index order)
and a **3-way histogram** `[cflag_value, pflag_value, level_index, count]` — the joint
flag table with the category as a third key.

This is the same masked-sum model as `flag_histogram`, plus a category filter, so a
consumer can count *any* phenotype in *any* compartment **split by category label**:

```js
// cells of IC label L, phenotype `pred`, in compartment compBit
function countByLevel(cat, compBit, pred, L) {   // cat = a catcols[] entry
  const li = cat.levels.indexOf(String(L));
  let n = 0;
  for (const [cflagStr, pflagStr, lvl, c] of cat.histogram) {
    if (lvl !== li) continue;
    const cflag = BigInt(cflagStr), pflag = BigInt(pflagStr);
    if (compBit >= 0 && !(cflag & (1n << BigInt(compBit)))) continue;
    if (pred(pflag)) n += c;
  }
  return n;
}
// e.g. CD8+ cells of immune cluster 17 in Tumor, as a density:
const dens = countByLevel(catOf("IC"), bitOf("Tumor"),
                          pf => pf & (1n<<BigInt(bitOf("CD8"))), 17) / areaOf("Tumor");
```

Summing a column's `histogram` over `level_index` (ignoring the third field) reproduces
`flag_histogram` exactly, so existing (non-stratified) consumers need no change. Note the
histogram is up to `n_levels`× larger than `flag_histogram`.

### Membership compartments (`<col>!=0` / `<col>=0`)

Per-level counts answer "how many cells of label L", but a *density* also needs a painted
**area**, and areas are only painted for cflag compartments (below). So for each categorical
column that carries a `"0"` baseline level — the "not-an-IC" convention — cohort also claims
two otherwise-unused **cflag bits** and sets them per cell: label `!= "0"` gets the `<col>!=0`
bit, label `== "0"` gets the `<col>=0` bit. These are set *before* painting and histogramming,
so they flow through the ordinary compartment path: each appears in `cflag_bits`, gets a
painted `area_mm2` in `compartments`, and lands in `flag_histogram`. A consumer therefore gets
`IC != 0` / `IC == 0` as first-class, gateable compartments with real density denominators —
no categorical-aware code needed, they behave exactly like `Tumor`/`TLS`.

The two bits partition every cell (each cell sets exactly one), so
`count(IC!=0) + count(IC==0) == n_cells`, and their `flag_histogram` marginals equal the
`catcols` level-0-vs-rest split. The bit *numbers* are chosen per sample (lowest free, after
reserving the `@FL`-declared bits); align them across samples by **name** (`IC!=0`), as with
any compartment. Columns with no `"0"` level, or with one partition empty, get no membership
bits.

## Schema

```jsonc
{
  "tool": "cyftools",
  "version": "0.1.0",
  "command": "cohort",
  "name": "Orion CRC",              // optional cohort label from -n/--name ("" if unset)
  "params": { "paint_radius_um": 20, "paint_pixel_um": 1 },
  "n_samples": 2,
  "samples": [
    {
      "file": "path/as-passed/LSP10353.byf",  // input path, verbatim
      "sample": "LSP10353",        // @SA sample name, else the filename stem
      "sample_id": 10353,          // numeric id = cell.id >> 32 ((sample_id<<32)|cell_id)
      "n_cells": 1620375,          // total cells streamed
      "bbox_um": [xmin, ymin, xmax, ymax],   // cell-coordinate bounds (µm)

      // legends: bit -> name for each register (only bits present in the data)
      "cflag_bits": [ {"name": "Tumor", "bit": 0}, {"name": "Margin", "bit": 2}, ... ],
      "pflag_bits": [ {"name": "CD3", "bit": 12}, {"name": "CD8", "bit": 14}, ... ],

      // precomputed painted areas — the density denominators. "All" (bit -1) first.
      "compartments": [
        { "cflag": "All",   "bit": -1, "n_cells": 1620375, "area_mm2": 198.81, "density_per_mm2": 8150.4 },
        { "cflag": "Tumor", "bit":  0, "n_cells":  942588, "area_mm2":  87.74, "density_per_mm2": 10743.2 },
        ...
      ],

      // contingency table of the two registers: [cflag_value, pflag_value, count],
      // one row per distinct (cflag,pflag) pattern. Flag VALUES are decimal STRINGS
      // (a 64-bit register can exceed JS's 2^53 safe int -> parse with BigInt).
      "flag_histogram": [
        ["0", "0", 41233],
        ["1", "131072", 8821],
        ["453", "521048", 7],
        ...
      ],

      // categorical stratifiers: one entry per @CA TY:A column (e.g. IC immune
      // clusters). "levels" is the LV: label list; "histogram" rows are
      // [cflag_value, pflag_value, level_index, count] — the joint table with the
      // category as a third key. Summing a column's histogram over level_index
      // reproduces flag_histogram exactly. Only present if the file has ≥1 A column.
      "catcols": [
        {
          "name": "IC",
          "type": "A",
          "levels": ["0", "1", "2", "17", "28", ...],   // real labels (LV: from the header)
          "histogram": [
            ["1", "131072", 0, 8102],    // cflag=1, pflag=131072, IC level 0 ("0")
            ["1", "131072", 3, 41],      // same (cflag,pflag), IC level 3 ("17")
            ...
          ]
        }
      ]
    }
    // ... one entry per input file that read successfully
  ]
}
```

Notes for consumers:

- `compartments[i].density_per_mm2` is the convenience "all cells in the
  compartment ÷ area" value; per-phenotype densities are computed from
  `flag_histogram` as above.
- Compartments overlap (cflag bits are independent), so their `n_cells`/areas
  need not sum to `All`.
- The `pflag_bits`/`cflag_bits` legends list only bits present in *that* sample;
  **different samples may expose different bits** — align across samples by name
  (or bit), not by array index.
- A compartment with no cells reports `area_mm2: 0` and `density_per_mm2: 0`
  (no division by zero). A file that fails to read is logged and skipped;
  `n_samples` counts only emitted samples.

## Limitations

- Areas (hence valid density denominators) exist only for **single cflag bits and
  `All`**. The histogram itself can count an *arbitrary* cflag predicate too (e.g.
  `Tumor AND Margin`, or `NOT Tumor`), but no painted area is stored for such a
  combination or complement, so it has no density (a consumer can still report
  counts / % — `docs/cyf_cohort.html` does, auto-switching its y-axis to `%`). If a
  combined-, complement-, or `NOT`-compartment area is needed, `cohort` would paint
  it as an extra compartment (e.g. the `IC=0` membership compartment *is* the painted
  complement of `IC!=0`).
- pflag/cflag names come from the header's `@FL` tags (`RG:cflag|pflag`,
  `BI:<bit>`). Files predating `@FL` fall back to the built-in cflag vocabulary
  (`Tumor`, `Margin`, `TLS`, …); any bit with no name is `cflag_bit_<n>` /
  `pflag_bit_<n>`.
