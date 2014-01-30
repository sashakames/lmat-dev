awk -v base=$3 '{if (($NF != "ReadTooShort") && ($(base+3) == $(base+4)) && (($(NF - 2) == 9606 && $(NF -1 ) ==  1) || ($(base+5) == 9606 && $(base+6) == 1))) sum += 1;} END {print sum " human matches"}' $1 >> $2.human_matches

