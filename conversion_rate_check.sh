#!/bin/bash

echo "chrC"
grep "C"  *CHH.bed.bismark.cov | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
echo "chrPt"
grep "Pt"  *CHH.bed.bismark.cov | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
echo "chrMt"
grep "Mt"  *CHH.bed.bismark.cov | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
