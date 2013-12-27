awk '{if ( $NF != "ReadTooShort" && ($NF == "NoDbHits" || ($6 != $7) || ($8 == 9606 && $9 < 1) || ($(NF - 2) == 9606 && $(NF -1 ) < 1))) print ">" $1 " " $2 "\n" $3 ;}' $1 > $1.new.fa
