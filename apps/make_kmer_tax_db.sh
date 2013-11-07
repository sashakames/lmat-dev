#!/bin/sh

## Need to fix directory issues

db_file_lst=""
use_moab=0
usage="Generate freq file
Usage: $0 options 

option list:
   --db_file_lst=""
   --use_moab (default=$use_grid), use grid
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi

while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --db_file_lst=*)
      db_file_lst=$optarg;;
   --use_moab*)
      use_moab=1;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done


while read kmer_db_part ; do
   echo "proc: $kmer_db_part"
   metag_dir="/g/g21/allen99/repo/latest/metag_repo"
   bin="$metag_dir/dev/apps"
   tdir="$metag_dir/taxonomy/kpath"
   ofile="$kmer_db_part.db" 
   log_file="$kmer_db_part.log" 
   cmd="$bin/tax_histo_new.runme -d $kmer_db_part -o $ofile -t $tdir/kpath_taxonomy.dat -g $tdir/gid_to_kpath_tax_id.dat > $log_file"
   if [ $use_moab -eq 1 ] ; then 
      run_moab_job.sh --machine=aztec --bank=cbitic --walltime="12:00:00" --ttc=1 --job_str="$cmd"
   else 
      $bin/tax_histo_new.runme -d $kmer_db_part -o $ofile -t $tdir/kpath_taxonomy.dat -g $tdir/gid_to_kpath_tax_id.dat > $log_file
   fi
done < $db_file_lst
