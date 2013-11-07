single=`grep "Singleton" $1 | awk '{total+=$2;} END {print total;}'`
kmerct=`grep "kmer-count" $1 | awk '{total+=$2;} END {print total;}'`
total=`grep "Total-t" $1 | awk '{total+=$2;} END {print total;}'`


echo "$single $kmerct $total" 

echo "$single $kmerct $total $2" | awk '{print (((($2 * 8) + (($3 - $1) * $4) + (($2 - $1) * $4)) / 1000000000) + 2);}'

