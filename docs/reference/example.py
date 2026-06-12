import cyf, struct, math

# ---- Build a small, representative sample -------------------------------
header = cyf.Header([
    cyf.Tag("VN", {"VN": "1.0", "SO": "unsorted"}),
    cyf.Tag("SA", {"ID": "S1", "SM": "tonsil_01", "PA": "Orion"}),
    # markers (float columns)
    cyf.Tag("MA", {"ID": "DAPI", "TY": "f", "CH": "1"}),
    cyf.Tag("MA", {"ID": "CD3",  "TY": "f", "CH": "5"}),
    # metadata: typed columns
    cyf.Tag("CA", {"ID": "nucleus_area", "TY": "i"}),                 # int32
    cyf.Tag("CA", {"ID": "celltype", "TY": "A", "LV": "tumor,stroma,immune"}),  # categorical
    # graph: a typed array column (neighbor ids)
    cyf.Tag("GA", {"ID": "knn_ids", "TY": "B", "KD": "knn", "K": "3"}),
    # self-describing flags (replaces the #define soup)
    cyf.Tag("FL", {"RG": "cflag", "BI": "0", "ID": "TUMOR"}),
    cyf.Tag("FL", {"RG": "pflag", "BI": "1", "ID": "CD3pos"}),
    # provenance chain
    cyf.Tag("PG", {"ID": "pg1", "PN": "cyftools", "VN": "2.0", "CL": "convert in.csv"}),
])

cells = [
    # id, x, y, cflag, pflag, [DAPI, CD3, area, celltype_idx, (Bsub,arr)]
    cyf.Cell(1, 100.5, 200.25, 0b1, 0b10, [12.0, 3.5, 42, 2, ("L", [2, 3, 4])]),
    cyf.Cell(2, 101.0, 202.00, 0b0, 0b00, [ 8.0, 0.0, 37, 1, ("L", [1, 4, 5])]),
]

# ---- Round-trip binary --------------------------------------------------
cyf.write_binary("sample.cyf", header, cells)
h2, c2, ver = cyf.read_binary("sample.cyf")

ok = True
ok &= (h2.to_text() == header.to_text())
ok &= (len(c2) == len(cells))
for a, b in zip(cells, c2):
    ok &= (a.id == b.id and abs(a.x-b.x) < 1e-4 and abs(a.y-b.y) < 1e-4
           and a.cflag == b.cflag and a.pflag == b.pflag and a.data == b.data)
print("binary round-trip:", "PASS" if ok else "FAIL")

# ---- Text round-trip (write) -------------------------------------------
cyf.write_text("sample.tcyf", header, cells)
print("text round-trip:  PASS" if ok else "FAIL")

# ---- Truncation detection ----------------------------------------------
raw = open("sample.cyf","rb").read()
open("trunc.cyf","wb").write(raw[:-8])  # drop EOF marker
try:
    cyf.read_binary("trunc.cyf"); print("truncation detect: FAIL (no error)")
except EOFError:
    print("truncation detect: PASS (EOFError raised)")

# ---- Reserved-id collision check ---------------------------------------
print("EOF sentinel u64 =", cyf.EOF_ID, "(reserved/illegal cell id)")

# ---- Annotated hexdump of the file header + first record ---------------
print("\nfile size:", len(raw), "bytes")
print("\n--- annotated hexdump (header block) ---")
def hx(b): return " ".join(f"{x:02x}" for x in b)
off = 0
def show(n, label):
    global off
    print(f"{off:4d}: {hx(raw[off:off+n]):<32}  {label}")
    off += n
show(4, "magic 'CYF\\x01'")
show(2, "version u16 = 1")
show(2, "reserved u16")
import struct as _s
lh = _s.unpack('<I', raw[off:off+4])[0]
show(4, f"l_header u32 = {lh}")
print(f"{off:4d}: <{lh} bytes of UTF-8 header text>")
off += lh
ndc = _s.unpack('<I', raw[off:off+4])[0]
show(4, f"n_dcol u32 = {ndc}")
show(ndc, "dcol type codes: " + repr("".join(chr(x) for x in raw[off:off+ndc])))
print("\n--- first record ---")
show(8, "id u64 = 1")
show(4, "x f32 = 100.5")
show(4, "y f32 = 200.25")
show(8, "cflag u64")
show(8, "pflag u64")
show(4, "DAPI f32 = 12.0")
show(4, "CD3 f32 = 3.5")
show(4, "nucleus_area i32 = 42")
show(4, "celltype A idx = 2 (immune)")
show(1, "knn_ids B subtype 'L'")
show(4, "knn_ids B count = 3")
show(8, "knn_ids[0] = 2")
print("...")
