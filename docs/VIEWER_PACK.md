# CYFV viewer pack

`cyftools export <in> <out.cyfv>` writes this format: a dependency-free,
column-major binary that a GPU front-end (deck.gl, regl, OpenGL, …) can load
directly into typed arrays — no Apache Arrow dependency on either side. It is a
**derived, read-only view** of a `.byf`/`.cyf` table for rendering; edits are
applied back to the `.byf` by cyftools (e.g. `flagroi`, `pheno`), never by
writing this file.

**You normally never create one by hand.** The launcher `viewer/view.py
cells.byf` treats the `.cyfv` as an *invisible render cache*: it runs `export`
into a temp directory (`$TMPDIR/cyfview-cache/`, overridable with
`CYFVIEW_CACHE`) keyed on the source file's path + mtime + size and the
cyftools binary's mtime. A cache hit is reused instantly; editing the `.byf`
(or rebuilding cyftools) changes the key and rebuilds the pack on the next
browser refresh. So the workflow is `.byf` in, `.byf` out — the `.cyfv` is a
throwaway cache, not a file you manage. (A `.cyfv` passed directly is served
as-is for debugging.)

All multi-byte values are **little-endian**.

## Layout

```
offset  field        type           notes
------  -----------  -------------  ----------------------------------------
0       magic        char[4]        "CYFV"
4       version      u16            5
6       reserved     u16            0
8       n_cells      u64            number of cells (N)
16      n_marker     u32            number of marker/meta columns (M)
20      names        M × (u16 len, char[len])   UTF-8 column names, in order
...     kinds        char[M]        per column: 'M' = marker (MA), 'C' = calculated (CA)
        <pad>        0x00 bytes     pad to the next multiple of 8 (align columns)
...     id           u64[N]         packed id: (sample_id<<32)|cell_id
...     x            f32[N]
...     y            f32[N]
...     cflag        u64[N]         cell flag register
...     pflag        u64[N]         phenotype flag register
...     marker_0     f32[N]         column-major: all N values of column 0
...     marker_1     f32[N]         then all N of column 1, …
...     marker_{M-1} f32[N]
...     n_roi        u32            number of @RO polygons (R)
...     roi          R × { u16 id_len, char[id_len],
                            u16 nm_len, char[nm_len],
                            u32 n_pts, f32[2*n_pts] (x,y interleaved) }
...     n_pflag      u32            number of @FL pflag (phenotype) bits (P)
...     pflag_map    P × { u16 name_len, char[name_len], u8 bit }   name -> pflag bit
...     n_cflag      u32            number of @FL cflag (structural) bits (C)
...     cflag_map    C × { u16 name_len, char[name_len], u8 bit }   name -> cflag bit
```

The per-cell `cflag`/`pflag` registers are columns above; these trailing maps
name the *bits*. `pflag_map` lets a viewer find the phenotype bit for a marker
(by name) and derive a gate from existing calls — the lowest marker value among
cells whose pflag bit is set. `cflag_map` names the structural flags for labeling.

The `<pad>` keeps the column block 8-byte aligned; with the column order above
(u64, then f32s, then u64s, then f32 markers) every array lands on a valid
typed-array boundary, so they can be viewed zero-copy. The trailing `n_roi`/ROI
block is read with a `DataView` (and the points copied), so it needs no padding.

Each array is contiguous, so after parsing the variable-length `names` block to
find the body offset, every column is a single typed-array view (and a ready GPU
attribute buffer). Columns appear in this fixed order: `id, x, y, cflag, pflag`,
then the `M` marker columns in `names` order.

## Minimal JS loader

```js
function loadCyfv(buf) {
  const dv = new DataView(buf);
  if (String.fromCharCode(dv.getUint8(0),dv.getUint8(1),dv.getUint8(2),dv.getUint8(3)) !== "CYFV")
    throw new Error("not a CYFV pack");
  const N = Number(dv.getBigUint64(8, true));   // n_cells  (bytes 8..15)
  const M = dv.getUint32(16, true);             // n_marker (bytes 16..19)
  let o = 20, names = [];
  for (let i = 0; i < M; i++) {
    const len = dv.getUint16(o, true); o += 2;
    names.push(new TextDecoder().decode(new Uint8Array(buf, o, len))); o += len;
  }
  o = (o + 7) & ~7;                              // skip the alignment pad
  const id    = new BigUint64Array(buf, o, N); o += N * 8;
  const x     = new Float32Array(buf, o, N);   o += N * 4;
  const y     = new Float32Array(buf, o, N);   o += N * 4;
  const cflag = new BigUint64Array(buf, o, N); o += N * 8;
  const pflag = new BigUint64Array(buf, o, N); o += N * 8;
  const markers = {};
  for (const name of names) { markers[name] = new Float32Array(buf, o, N); o += N * 4; }
  const nRoi = dv.getUint32(o, true); o += 4;    // @RO polygons
  const rois = [];
  for (let r = 0; r < nRoi; r++) {
    const idLen = dv.getUint16(o, true); o += 2;
    const rid = new TextDecoder().decode(new Uint8Array(buf, o, idLen)); o += idLen;
    const nmLen = dv.getUint16(o, true); o += 2;
    const name = new TextDecoder().decode(new Uint8Array(buf, o, nmLen)); o += nmLen;
    const nPts = dv.getUint32(o, true); o += 4;
    const pts = new Float32Array(buf.slice(o, o + nPts * 8)); o += nPts * 8;  // copy (alignment)
    rois.push({ id: rid, name, pts });
  }
  return { N, id, x, y, cflag, pflag, markers, rois };
}
```

## Round-trip with the source

The packed `id` is the same `(sample_id<<32)|cell_id` as in the `.byf` (see
[CYF_FORMAT.md §5](CYF_FORMAT.md)). A selection in the viewer is therefore either
a small set of `id`s or a polygon. To write it back, send the polygon to
`cyftools addroi` (as an `@RO`) and apply it with `cyftools flagroi`, or send a
gate definition to `cyftools pheno` — the bulk point data never leaves C++.
