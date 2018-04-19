#!/bin/bash
set -u
# Get sodium bisulfite conversion efficiency by calculating percent unconverted C's in Clp and Mt

echo "Conversion % in $1"

file=$1

echo "ChrC"
grep "ChrC"  ${file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
echo "ChrM"
grep "ChrM"  ${file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
# echo "Pt"
# grep "Pt"  ${file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
# echo "Mt"
# grep "Mt"  ${file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'

