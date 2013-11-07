#!/bin/sh -xvf

#/collab/usr/global/tools/mpi/install/bin/enable_cpus
#export GOMP_CPU_AFFINITY=0-79
## Need to fix directory issues
metag_dir=$METAG_DIR
bin="$metag_dir/dev/apps"
tdir="$metag_dir/dev/runtime_inputs"

export LD_LIBRARY_PATH=/usr/mic/post1/metagenomics/lib
#export LD_LIBRARY_PATH=/p/lscratchc/allen99/metag/metag_repo/rel/LMAT-1.1/third-party/perm-je-0.9.3/lib

#nullm=$tdir/cutoffs.kcnt.m5.100.flst
dbfile=/local/ramfs/genedb
depthf="$tdir/depth_for_ncbi_taxonomy.segment_20130204.dat"
taxtree="$tdir/ncbi_taxonomy.segment_20130204.dat.nohl"
g2t="$tdir/gene_name_tax_id.txt"
thresh=100

#nullm=$tdir/cutoffs.kcnt.50.m3_k=20_leaf_jan_3_2012.flst
#dbfile=/local/ramfs/m3_k=20_leaf_jan_3_2012
#depthf="$tdir/depth_for_ncbi_taxonomy.dat"
#taxtree="$tdir/ncbi_taxonomy.dat.nohl"
#thresh=50

verbose=0;
query_file=""
query_lst=""
prune_depth=0
use_moab=0
sdiff=1
odir=.
usage="Generate freq file
Usage: $0 options 

option list:
   --db_file=$dbfile
   --query_file=$query_file
   --query_lst=$query_lst
   --filter=$thresh
   --threads=$threads
   --nullm=$nullm
   --verbose=$verbose
   --odir=$odir
   --sdiff=$sdiff

example usage:
$0 --db_file=$dbfile --query_file=HC1.fna --threads=1 

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
      query_file=$optarg;;
   --query_lst=*)
      query_lst=$optarg;;
   --filter=*)
      thresh=$optarg;;
   --threads=*)
      threads=$optarg;;
   --nullm=*)
      nullm=$optarg;;
   --odir=*)
      odir=$optarg;;
   --sdiff=*)
      sdiff=$optarg;;
   --verbose)
      verbose=1;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done
cfiles=$nullm
dbname=`basename $dbfile`

qstr=""
if [ $query_file ] ; then
   query_file_name=`basename $query_file`
   qstr="-i $query_file -t $threads"
else
   query_file_name=`basename $query_lst`
   qstr="-l $query_lst"
fi
ofile="$odir/$query_file_name.$dbname.$thresh.lo.rl_output" 
vstr=""
if [ $verbose == 1 ] ; then
   vstr="-y"
fi
logfile="$ofile.log" 

/usr/bin/time -v $bin/gene_label $vstr -q $g2t -e $depthf -p -r $qstr -d $dbfile -c $taxtree -o $ofile >& $logfile 
