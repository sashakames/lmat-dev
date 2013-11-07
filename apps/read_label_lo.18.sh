#!/bin/sh

## Need to fix directory issues

#export LD_LIBRARY_PATH=/mnt/fusionio/lib
#export LD_LIBRARY_PATH=/media/KmerDB/lib
thresh=0.001
threads=80
dbfile=""
query_fn=""
prune_depth=0
use_moab=0
model_size=100000
cutoff=50
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
   --cutoff=*)
      cutoff=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

svn info

query_file="/siscnas/scratch/ames4/kmerdb-in/$query_fn"
metag_dir="$HOME/svn/metag_repo"
bin="$metag_dir/dev/apps"
tdir="$metag_dir/dev/runtime_inputs"
gtmap="$tdir/gid_to_tid_all.txt"
depthf="$tdir/depth_for_kpath_taxonomy.dat"
dbname=`basename $dbfile`
ofile="/mnt/fusionio/$query_fn.$dbname.$model_size.$cutoff.lo.rl_output" 
cfiles="$tdir/cutoffs.kcnt.$cutoff.$dbname.flst"
logfile="$ofile.log" 
#$bin/read_label_lo -n $cfiles -e $depthf -w -p -r -t $threads -i $query_file -d $dbfile -c $tdir/kpath_taxonomy.dat -o $ofile -v $thresh -g $gtmap -a 
$bin/read_label_lo -n $cfiles -e $depthf -w -p -r -t $threads -i $query_file -d $dbfile -c $tdir/kpath_taxonomy.dat -o $ofile -v $thresh -g $gtmap -a -h $cutoff -k 18
