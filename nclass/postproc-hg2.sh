awk '{if ( ( $NF == "DirectMatch" || $NF == "MultiMatch" ) && $(NF - 2) != 9606 ) print ">" $1 " " $2 " " $(NF - 2) "\n" $3  ;}' $1 > $1.mod
