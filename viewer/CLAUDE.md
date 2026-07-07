# viewer/CLAUDE.md — cyf_cohort.html notes & session progress

Scoped notes for `viewer/cyf_cohort.html` (the cohort comparison viewer, a.k.a.
"cycif_cohort"). Read the repo-root `CLAUDE.md` first for the overall project.

## Status (as of 2026-07-07)

**All work described here is committed** as `ca3bc62 "updated ROI processing"`.
Working tree is clean — nothing is lost on restart. The commit spans
`src/` (C++), `viewer/cyf_cohort.html`, and `docs/COHORT_JSON.md`.

⚠️ **Rebuild the binary before using the C++ changes** — `build/cyftools` on disk
predates them:
```sh
cmake --build build -j && ctest --test-dir build   # cyf_test should pass
```

## How to test the viewer (the harness used this session)

The page loads a `cohort.json` (produced by `cyftools cohort`). Two gotchas:
- `file://` is blocked; **serve over HTTP**: `python3 -m http.server` in a dir
  containing `cyf_cohort.html` + a `cohort.json`, then open
  `http://127.0.0.1:PORT/cyf_cohort.html?file=cohort.json` (the `?file=` triggers
  auto-load).
- The browser MCP `javascript_tool` runs in an **isolated world** — page globals
  (`state`, functions) are NOT visible; only the shared DOM is. Drive via real
  DOM (click `#btn-export-r-box`, read `#r-modal-text`). To capture large R
  scripts past the content filter, POST them from the page to a tiny local writer
  server and inspect the file on disk (that's how the R was verified).
- Generate a grouped `cohort.json`: build a few `.byf` with
  `cyftools convert … -u micron -M CD4`, tag groups via
  `cyftools addtag in.byf out.byf --group <Label>`, then `cyftools cohort -o cohort.json *.byf`.
- Inspect exported PDFs: `pdffonts f.pdf` (embedded/subset fonts),
  `pdftotext f.pdf -` (selectable text = editable, empty = outlined paths),
  `pdftotext -bbox` (per-word boxes / glyph heights), `pdfinfo` (page size).

## R export subsystem (the part with the most hard-won context)

Two export buttons → `buildBoxPlotR()` (`#btn-export-r-box`) and
`buildWaterfallR()` (`#btn-export-r-wf`); `buildRExportOriginal()` is **dead
code** (faceted, not wired to a button) but kept in sync. Shared helpers:
`rExportHeader()`, `rEditableText()`, `rSavePlot()`, plus the duplicated
`create_beeswarm_plot` R helper (two copies — edit both, usually via replace_all).

### Font / device design — DO NOT casually "simplify" this; each piece fixes a real failure

The goal: **editable, embedded Arial text at absolute 7pt** when opened in
Illustrator. The path there was painful; the rules:

1. **Device = `cairo_pdf`** (in `rSavePlot()`), never base `pdf()`.
   `base_family = "Arial"` on the base `pdf()` device throws
   `"invalid font type"` (Arial isn't a base-14 PDF font). cairo_pdf embeds/
   subsets Arial and keeps text editable. The script **saves itself** via
   `ggsave(..., device = cairo_pdf)` — users must **NOT** wrap it in
   `pdf(); …; dev.off()`.
2. **`if ("showtext" %in% loadedNamespaces()) showtext::showtext_auto(FALSE)`**
   at the top of every script (`rEditableText()`). showtext is a **sticky,
   session-global** switch: once any earlier run called `showtext_auto()`, it
   draws every glyph as vector **outlines** (not editable) on ALL devices for the
   rest of the session. This line neutralizes a poisoned session so text stays
   editable. (This was the root cause of "it STILL renders as paths".)
3. **`print(g)` is guarded to `if (interactive())`** — a headless `Rscript`
   default device is base `pdf()`, which would error on `base_family="Arial"`.
4. **Absolute 7pt, uniform, no scaling with plot area.** cairo_pdf renders point
   sizes as physical points → 7pt is 7pt at any width/height (verified: identical
   glyph height at 5×4 and 10×8). Every text element is pinned explicitly to
   `size = 7` (theme `text`/`axis.title`/`axis.text`, and `base_size = 7`), so no
   `rel(0.8)` scaling makes ticks 5.6pt. If the user says "size 7" they mean an
   absolute 7 everywhere.
5. **Text is full black** (`colour = "black"` on text/axis elements) — theme_bw's
   default grey30 was rejected.

### Known limitation — ggsignif / significance labels stay Helvetica

`ggpubr::stat_compare_means(label="p.signif", …)` draws its brackets via
**ggsignif**, whose label text is stubbornly **Helvetica** — it ignores the
theme's `base_family`, the `family=` arg, and even
`update_geom_defaults(...)`. So the `ns`/`*` labels are the ONE non-Arial bit in
the box PDF. Everything else is Arial. Options if it matters: retype those few
labels in Illustrator, or drop the comparisons from the export. Not yet resolved.

### SVG vs PDF (context if the "merged labels" complaint returns)

cairo_pdf groups adjacent same-font labels into one `BT…ET` text object →
Illustrator imports the group axis labels as ONE merged text box. `svglite`
(SVG) writes one `<text>` per label → separate editable objects, and keeps Arial.
We briefly switched to SVG for that reason, then the user chose **PDF** again
(accepting the merge). If they want separate label objects back, swap
`rSavePlot()` to `device = svglite::svglite` and `.svg` filenames — it's a
localized change. (svglite is installed; SVG font-family stays "Arial", editable.)

### Current box-plot export settings

- Save size **1.7 × 1.6 in** (`rSavePlot('cohort_boxplot.pdf', 1.7, 1.6)`);
  waterfall 6×5, markers 7×5.
- Box/whisker line width **0.25** (`geom_boxplot(linewidth = 0.25)`, half the 0.5
  default).
- Y-axis ticks **~2× denser**: `scale_y_continuous(breaks = scales::pretty_breaks(n = 10))`.
- Significance label `size = 2.46` (≈7pt; ggpubr size is in mm) `family = "Arial"`
  (harmless even though ggsignif ignores it — see limitation).

## Group ordering feature

Box-plot / waterfall / stats / R-export group order is user-controllable via
**◀ ▶ arrows** on the "Group order & colors" chip row. Core: `state.groupOrder`
(persisted in settings + optional `group_order` key in a config JSON) and
`orderedGroups(list)` — listed groups first, then the rest alphabetically; new
groups fall to the alphabetical tail; stale entries are skipped. `moveGroup()`
swaps neighbours. `getGroupColor()` intentionally still indexes colors
**alphabetically** so reordering doesn't recolor groups. Plotly box uses
`categoryorder:'array', categoryarray: groups` to honor it.

## Mann–Whitney p-value fix (statistics bug)

`mannWhitneyU()` could return **p > 1** (user saw 1.06) for near-tied groups: the
0.5 continuity correction made `|U1-meanU| - 0.5` negative → negative z →
`2*(1-normCdf(z)) > 1`. Fixed by clamping the corrected gap at 0 and p at 1:
```js
const z = Math.max(0, Math.abs(U1 - meanU) - 0.5) / sdU;
const p = Math.min(1, 2 * (1 - normCdf(z)));
```
Only the in-browser "MW p" column was affected; the R export uses R's own
`wilcox.test` (always ≤1). Real p-values unchanged; only impossible >1 → ≤1.

## Related C++ changes this session (in commit ca3bc62)

All in `src/`, documented at their usage/help strings:
- **`flagroi`** — arg is now the **decimal flag value** (a power of two, e.g. `8`
  = bit index 3), validated; added **`-O/--overwrite`** (clear the bit on all
  cells first, then set inside ROIs) vs default additive OR. See `flagroifunc`
  and `FlagRoiProcessor`.
- **`clearroi`** — new **`-i <id>`** to remove a single @RO by exact ID
  (disambiguates duplicate names). `RemoveRoiTags` gained an `id_filter` param.
- **`roiinfo`** — NEW module: lists each @RO's id / name / sample / vertices /
  area (shoelace, `Polygon::Area()` in `polygon.cpp`) plus a microns² column when
  `@HD UN/MP` allow. Header-only reader (`RoiInfoProcessor`). Dispatch + the
  "is-implemented" allowlist in `cyftools.cpp` both updated. IDs are shown by
  `cyftools view -H` (NOT `info`).
- **`cohort`** — now converts each file's coords to **microns per-file** using its
  own `@HD UN/MP` (was silently assuming microns); skips a pixel-unit file with no
  usable MP. Renamed the raster-resolution flag `-p` → **`--grid`** (`-p` still an
  alias) since it's the area-raster cell size, not an image pixel size. See
  `CohortProcessor::ProcessHeader/ProcessLine` and `docs/COHORT_JSON.md`.

## Immediate open items / next steps

- Decide the ggsignif Helvetica question (leave / drop significance brackets).
- `cyftools roiinfo` not yet documented in root `CLAUDE.md` command list or README.
- Consider whether cohort's `paint_pixel_um` JSON key should be renamed to
  `paint_grid_um` (kept as-is to avoid breaking consumers).
- Rebuild + `ctest` after restart (binary on disk is stale).
