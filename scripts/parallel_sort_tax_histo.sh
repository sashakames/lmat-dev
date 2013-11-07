for (( i = 0 ; i < 16  ; i = i + 1 ))
do
LOW=$(( $(( $i * 16)) ))
HI=$(( $LOW + 15 ))
echo "seq $LOW $HI | $HOME/parallel-20120122/src//parallel single_sort.sh {.}" > job_$i.sh

done

