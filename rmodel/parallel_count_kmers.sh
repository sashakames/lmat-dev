export LD_LIBRARY_PATH="$HOME/metag_repo/rel/LMAT-1.1/third-party/perm-je-0.9.3/lib"

for (( i = 0 ; i < 8 ; i = i + 1 ))
do
LOW=$(( $(( $i * 32)) ))
HI=$(( $LOW + 31 ))
echo "seq $LOW $HI | $HOME/parallel-20120122/src/parallel ./countTaxidFrequency -i /p/lscratchd/hysom/m5/m5.{.}.tax_histo.bin -f 32 -o /p/lscratchd/ames4/m5.20.{.}" > job_$i.sh
msub -A genome -q pbatch job_$i.sh -l walltime=24:00:00
done

