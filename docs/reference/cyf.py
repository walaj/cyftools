#!/usr/bin/env python3
"""
Reference codec for the CYF format (v1, uncompressed).
Pure-Python, no C++ / cereal. Proves the byte layout is unambiguous and
implementable in any language directly from the spec.
"""
import struct, io, math

MAGIC      = b"CYF\x01"          # 4 bytes
EOF_MARKER = b"CYF\x01EOF\x00"   # 8 bytes; its LE-u64 value is a reserved/illegal cell id
EOF_ID     = struct.unpack("<Q", EOF_MARKER)[0]
VERSION    = 1

# type code -> (struct fmt for fixed-width, or special)
FIXED = {
    "i": ("<i", 4), "I": ("<I", 4),
    "l": ("<q", 8), "L": ("<Q", 8),
    "f": ("<f", 4), "d": ("<d", 8),
}
# Z = string (u32 len + utf8), A = categorical (u32 index into header LV), B = typed array

# ---------------------------------------------------------------- header model
class Tag:
    def __init__(self, cls, fields):  # cls in {HD,SA,MA,CA,GA,FL,PG}; fields = dict (ordered)
        self.cls = cls
        self.fields = dict(fields)
    def to_text(self):
        parts = ["@" + self.cls]
        for k, v in self.fields.items():
            parts.append(f"{k}:{v}")
        return "\t".join(parts)
    @staticmethod
    def from_text(line):
        toks = line.rstrip("\n").split("\t")
        cls = toks[0][1:]
        fields = {}
        for t in toks[1:]:
            k, _, v = t.partition(":")
            fields[k] = v
        return Tag(cls, fields)

DATA_CLASSES = ("MA", "CA", "GA")  # classes that define a data column, in header order

class Header:
    def __init__(self, tags):
        self.tags = list(tags)
    def data_columns(self):
        """Ordered list of (Tag) that define data columns -> the record schema."""
        return [t for t in self.tags if t.cls in DATA_CLASSES]
    def type_codes(self):
        return [t.fields.get("TY", "f") for t in self.data_columns()]
    def levels(self, tag):  # categorical dictionary
        return tag.fields["LV"].split(",") if "LV" in tag.fields else []
    def to_text(self):
        return "".join(t.to_text() + "\n" for t in self.tags)
    @staticmethod
    def from_text(text):
        return Header([Tag.from_text(l) for l in text.splitlines() if l.startswith("@")])

# ---------------------------------------------------------------- record model
class Cell:
    __slots__ = ("id", "x", "y", "cflag", "pflag", "data")
    def __init__(self, id, x, y, cflag, pflag, data):
        self.id, self.x, self.y = id, x, y
        self.cflag, self.pflag = cflag, pflag
        self.data = list(data)  # one entry per data column, native python value

# ---------------------------------------------------------------- value codecs
def _write_value(buf, ty, val):
    if ty in FIXED:
        buf.write(struct.pack(FIXED[ty][0], val)); return
    if ty == "Z":
        b = val.encode("utf-8"); buf.write(struct.pack("<I", len(b))); buf.write(b); return
    if ty == "A":
        buf.write(struct.pack("<I", val)); return          # category index
    if ty == "B":
        sub, arr = val                                     # (subtype_code, list)
        buf.write(struct.pack("<B", ord(sub)))
        buf.write(struct.pack("<I", len(arr)))
        fmt = FIXED[sub][0]
        for e in arr:
            buf.write(struct.pack(fmt, e))
        return
    raise ValueError(f"bad type {ty}")

def _read_value(rd, ty):
    if ty in FIXED:
        fmt, n = FIXED[ty]; return struct.unpack(fmt, rd.read(n))[0]
    if ty == "Z":
        n = struct.unpack("<I", rd.read(4))[0]; return rd.read(n).decode("utf-8")
    if ty == "A":
        return struct.unpack("<I", rd.read(4))[0]
    if ty == "B":
        sub = chr(struct.unpack("<B", rd.read(1))[0])
        n = struct.unpack("<I", rd.read(4))[0]
        fmt, w = FIXED[sub]
        arr = [struct.unpack(fmt, rd.read(w))[0] for _ in range(n)]
        return (sub, arr)
    raise ValueError(f"bad type {ty}")

# ---------------------------------------------------------------- binary writer
def write_binary(path, header, cells):
    htext = header.to_text().encode("utf-8")
    tys = header.type_codes()
    with open(path, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<H", VERSION))
        f.write(struct.pack("<H", 0))                  # reserved
        f.write(struct.pack("<I", len(htext)))
        f.write(htext)
        f.write(struct.pack("<I", len(tys)))
        f.write(bytes(ord(t[0]) if t in ("i","I","l","L","f","d","Z","A","B") else ord(t) for t in tys))
        for c in cells:
            if c.id == EOF_ID:
                raise ValueError("cell id collides with reserved EOF sentinel")
            rec = io.BytesIO()
            rec.write(struct.pack("<Q", c.id))
            rec.write(struct.pack("<f", c.x))
            rec.write(struct.pack("<f", c.y))
            rec.write(struct.pack("<Q", c.cflag))
            rec.write(struct.pack("<Q", c.pflag))
            for ty, val in zip(tys, c.data):
                _write_value(rec, ty, val)
            f.write(rec.getvalue())
        f.write(EOF_MARKER)

# ---------------------------------------------------------------- binary reader
def read_binary(path):
    with open(path, "rb") as f:
        if f.read(4) != MAGIC:
            raise ValueError("not a CYF binary file (bad magic)")
        ver, = struct.unpack("<H", f.read(2))
        struct.unpack("<H", f.read(2))                 # reserved
        lh, = struct.unpack("<I", f.read(4))
        htext = f.read(lh).decode("utf-8")
        header = Header.from_text(htext)
        ndc, = struct.unpack("<I", f.read(4))
        tys = [chr(b) for b in f.read(ndc)]
        assert tys == header.type_codes(), "type-code array disagrees with header"
        cells = []
        saw_eof = False
        while True:
            peek = f.read(8)
            if len(peek) == 0:
                break                                  # truncated: no EOF marker
            if peek == EOF_MARKER:
                saw_eof = True
                break
            cid, = struct.unpack("<Q", peek)
            x, = struct.unpack("<f", f.read(4))
            y, = struct.unpack("<f", f.read(4))
            cflag, = struct.unpack("<Q", f.read(8))
            pflag, = struct.unpack("<Q", f.read(8))
            data = [_read_value(f, ty) for ty in tys]
            cells.append(Cell(cid, x, y, cflag, pflag, data))
        if not saw_eof:
            raise EOFError("stream ended without EOF marker (possible truncation)")
        return header, cells, ver

# ---------------------------------------------------------------- text writer/reader
def write_text(path, header, cells):
    tys = header.type_codes()
    with open(path, "w") as f:
        f.write(header.to_text())
        for c in cells:
            fields = [str(c.id), repr(c.x), repr(c.y), str(c.cflag), str(c.pflag)]
            for ty, val in zip(tys, c.data):
                if ty == "B":
                    sub, arr = val
                    fields.append(sub + ":" + ",".join(str(e) for e in arr))
                else:
                    fields.append(str(val))
            f.write("\t".join(fields) + "\n")
