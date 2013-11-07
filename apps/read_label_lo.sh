#!/bin/sh -xvf

## Need to fix directory issues



thresh=0.001
threads=6
dbfile="/local/ramfs/test"
query_file=""
prune_depth=0
use_moab=0
usage="Generate freq file
Usage: $0 options 

option list:
   --db_file=$dbfile
   --query_file=$query_file
   --filter=$thresh
   --threads=$threads

example usage:
$0 --db_file=/mnt/ramfs/microbe2-18 --query_file=HC1.fna --threads=1 

"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi

while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --db_file=*)
      dbfile=$optarg;;
   --query_file=*)
      query_fn=$optarg;;
   --filter=*)
      thresh=$optarg;;
   --threads=*)
      threads=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done


query_file="/siscnas/scratch/ames4/kmerdb-in/$query_fn"
metag_dir="$HOME/metag_repo"

bin="$HOME/test/1.0/src"
tdir="$metag_dir/dev/runtime_inputs"
gtmap="$tdir/gid_to_tid_all.txt"
depthf="$tdir/depth_for_kpath_taxonomy.dat"
dbname=`basename $dbfile`
ofile="~/scratch//$query_fn.$dbname.$thresh.lo.rl_output" 
cfiles="$tdir/cutoffs.kcnt.flst"
logfile="$ofile.log" 
#$bin/read_label_lo -n $cfiles -e $depthf -w -p -r -t $threads -i $query_file -d $dbfile -c $tdir/kpath_taxonomy.dat -o $ofile -v $thresh -g $gtmap -a 
/usr/bin/time -v $bin/read_label_lo -n $cfiles -e $depthf -w -p -r -t $threads -i $query_file -d $dbfile -c $tdir/kpath_taxonomy.dat -o $ofile -v $thresh -g $gtmap -a 