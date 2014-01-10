
awk -v base=$3 '{if ( ( $NF == "DirectMatch" || $NF == "MultiMatch" ) && $(NF - 2) != 9606  && !($(base+5) == 9606 && $(base+6) == 1 ) && $(NF-2) > 1 ) print ">" $1 " " $2 " " $(NF - 2) "\n" $(base)  ;}' $1 >> $2.non_hg.fa
