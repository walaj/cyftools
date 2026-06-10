#!/usr/bin/env bash

# Usage:
# ./roi_shift.sh input.csv -x 1000 -y 1000 > shifted.csv
# or
# ./roi_shift.sh input.csv -x 1000 > shifted.csv   (y defaults to 0)

infile="$1"
shift

# defaults
shift_x=0
shift_y=0

# parse flags
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -x) shift_x="$2"; shift 2;;
    -y) shift_y="$2"; shift 2;;
    *) break;;   # ignore trailing args (no outfile expected)
  esac
done

awk -v sx="$shift_x" -v sy="$shift_y" -F, '
BEGIN { OFS="," }
NR==1 { print; next }
{
    split($5, pts, " ")
    newpts = ""
    for (i in pts) {
        split(pts[i], xy, ",")
        x = xy[1] + sx
        y = xy[2] + sy
        newpts = newpts (newpts=="" ? "" : " ") x "," y
    }
    $5 = newpts
    print
}' "$infile"
