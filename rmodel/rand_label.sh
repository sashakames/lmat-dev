#!/bin/sh -xvf

## Need to fix directory issues
#export LD_LIBRARY_PATH=/usr/mic/post1/metagenomics/lib
/collab/usr/global/tools/mpi/install/bin/enable_cpus
export GOMP_CPU_AFFINITY=0-79


export LD_LIBRARY_PATH=/usr/mic/post1/metagenomics/lib
#export LD_LIBRARY_PATH=$METAG_DIR/third-party/perm-je-0.9.3/lib
#export LD_LIBRARY_PATH=$METAG_DIR/rel/LMAT-1.1/third-party/perm-je-0.9.3/lib
verbose=0;
thresh=50
threads=1
dbfile=""
query_file=""
usage="Run read_label on randomly generated sequences
Usage: $0 options 

option list:
   --db_file=$dbfile
   --query_file=$query_file
   --filter=$thresh
   --threads=$threads
   --verbose=$verbose

example usage:
$0 --db_file=/mnt/ramfs/microbe2-18 --query_file=HC1.fna --threads=1 

"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
odir=.
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --db_file=*)
      dbfile=$optarg;;
   --query_file=*)
      query_file=$optarg;;
   --filter=*)
      thresh=$optarg;;
   --threads=*)
      threads=$optarg;;
   --output_dir)
      odir=$optarg;;
   --verbose)
      verbose=1;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done



metag_dir=$METAG_DIR
bin="$metag_dir/dev/apps"
tdir="$metag_dir/dev/runtime_inputs"

depthf="$tdir/depth_for_ncbi_taxonomy.segment_20130204.dat"
tax_tree="$tdir/ncbi_taxonomy.segment_20130204.dat"

dbname=`basename $dbfile`
oname=$query_file.$dbname.$thresh.rl_output
ofile=$odir/$oname
logfile="$ofile.log" 
vstr=""
if [ $verbose == 1 ] ; then
   vstr="-y"
fi

$bin/read_label -l 0.0 -h $thresh $vstr -e $depthf -p -t $threads -i $query_file -d $dbfile -c $tax_tree -o $ofile >& $logfile

sfile=$odir/out.$oname
if [ -e $sfile ] ;  then
   rm -f $sfile
fi
for lfile in `find $odir -maxdepth 1 -name $oname\*.out` ; do
   echo "check: $lfile"
   cat $lfile >> $sfile
done
for lfile in `find $odir -maxdepth 1 -name $oname\*.out` ; do
   rm -f $lfile
done
echo "Done"
