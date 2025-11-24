#!/usr/bin/env python3
import sys
import csv
import argparse

def shift_all_points(all_points_str, dx, dy):
    """
    Shift all x,y pairs in the all_points string by dx, dy.
    all_points_str: "1368.83,4853.56 604.07,7358.43 ..."
    """
    if not all_points_str:
        return all_points_str
    tokens = all_points_str.strip().split()
    out = []
    for t in tokens:
        try:
            x_str, y_str = t.split(",")
            x = float(x_str) + dx
            y = float(y_str) + dy
            out.append(f"{x},{y}")
        except ValueError:
            # If anything is weird, keep token unchanged
            out.append(t)
    return " ".join(out)

def csv_quote(field, force=False):
    """
    Basic CSV quoting:
      - If force=True: always wrap in double quotes, escaping internal quotes.
      - If force=False: quote only if needed (comma, quote, or newline present).
    """
    if field is None:
        field = ""
    s = str(field)
    needs_quote = force or ("," in s or '"' in s or "\n" in s)
    if not needs_quote:
        return s
    return '"' + s.replace('"', '""') + '"'

def write_data_row(row, dx, dy, all_points_col_index):
    # Shift all_points if present
    if len(row) > all_points_col_index:
        row[all_points_col_index] = shift_all_points(row[all_points_col_index], dx, dy)

    out_fields = []

    # Column 1: no quotes
    if len(row) >= 1:
        out_fields.append(str(row[0]))
    else:
        out_fields.append("")

    # Columns 2-4: always quoted
    for idx in range(1, len(row)):
        colnum = idx + 1  # 1-based index
        val = row[idx]
        if colnum in (2, 3, 4):
            # force quotes
            out_fields.append(csv_quote(val, force=True))
        else:
            # quote only if needed
            out_fields.append(csv_quote(val, force=False))

    sys.stdout.write(",".join(out_fields) + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="Shift ROI polygon coordinates in CSV all_points column."
    )
    parser.add_argument("infile", help="Input CSV file")
    parser.add_argument("-x", type=float, default=0.0,
                        help="Shift in X direction (added to all x coords)")
    parser.add_argument("-y", type=float, default=0.0,
                        help="Shift in Y direction (added to all y coords)")
    args = parser.parse_args()

    dx, dy = args.x, args.y
    all_points_col_index = 4  # 0-based index: 5th column

    with open(args.infile, newline="") as f:
        reader = csv.reader(f)
        first_row = next(reader, None)
        if first_row is None:
            return

        # Detect header: check if first row looks like column names
        header_candidates = [c.lower() for c in first_row]
        has_header = any(c in ("id", "name", "text", "type", "all_points")
                         for c in header_candidates)

        if has_header:
            # Print header exactly as unquoted comma-joined text
            sys.stdout.write(",".join(first_row) + "\n")
        else:
            # No header: treat first_row as data
            write_data_row(first_row, dx, dy, all_points_col_index)

        # Process remaining rows as data
        for row in reader:
            write_data_row(row, dx, dy, all_points_col_index)

if __name__ == "__main__":
    main()
