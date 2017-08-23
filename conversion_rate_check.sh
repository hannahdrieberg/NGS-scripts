#!/bin/bash
set -eu
# Get sodium bisulfite conversion efficiency by calculating percent unconverted C's in Cp

echo "Conversion % in $1"

file=$1

echo "chrC"
grep "C"  ${file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
echo "chrPt"
grep "Pt"  ${file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
echo "chrMt"
grep "Mt"  ${file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
