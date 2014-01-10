
awk -v base=$3 '{if ( $NF != "ReadTooShort" && ($NF == "NoDbHits" || ($(base+3) != $(base+4)) || ($(base+5) == 9606 && $(base+6) < 1) || ($(NF - 2) == 9606 && $(NF -1 ) < 1))) print ">" $1 " " $2 "\n" $(base) ;}' $1 >> $2.new.fa
