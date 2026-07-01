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
```

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

## Schema

```jsonc
{
  "tool": "cyftools",
  "version": "0.1.0",
  "command": "cohort",
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
  `Tumor AND Margin`), but no painted area is stored for such a combination, so it
  has no density. If a combined-compartment area is needed later, `cohort` would
  paint it as an extra compartment.
- pflag/cflag names come from the header's `@FL` tags (`RG:cflag|pflag`,
  `BI:<bit>`). Files predating `@FL` fall back to the built-in cflag vocabulary
  (`Tumor`, `Margin`, `TLS`, …); any bit with no name is `cflag_bit_<n>` /
  `pflag_bit_<n>`.
