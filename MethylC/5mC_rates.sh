#!/bin/bash

set -u

# Get methylation rates for all contexts across Chr1-5 as well as % CHH in chloroplast genome as an indication of sodium bisulfite conversion efficiency (unconverted CHH in Cp and Mt genome)

if [ "$#" -lt 1 ]; then
	echo "Missing required arguments!"
	echo "USAGE: methylation_rates.sh <sample ID>"
	echo "EXAMPLE: methylation_rates.sh col0-r1"
	exit 1
fi

file=$1
cg=${file}*CG*cov
chg=${file}*CHG*cov
chh=${file}*CHH*cov

echo "5mC % in $1"

echo "mCG Chr1-5: "$cg" "
grep -e "Chr1" -e "Chr2" -e "Chr3" -e "Chr4" -e "Chr5" $cg |  awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHG Chr1-5: "$chg" "
grep -e "Chr1" -e "Chr2" -e "Chr3" -e "Chr4" -e "Chr5" $chg |  awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHH Chr1-5: "$chh" "
grep -e "Chr1" -e "Chr2" -e "Chr3" -e "Chr4" -e "Chr5" $chh |  awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHH ChrC: "$chh" "
grep -e "ChrC"  $chh | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "mCHH ChrM: "$chh" "
grep -e "ChrM"  $chh | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print ((met / total))}'

echo "DONE"
