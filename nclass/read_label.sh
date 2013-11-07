#!/bin/sh -xvf

#/collab/usr/global/tools/mpi/install/bin/enable_cpus
#export GOMP_CPU_AFFINITY=0-79
## Need to fix directory issues
metag_dir=$METAG_DIR
bin="$metag_dir/dev/apps"
tdir="$metag_dir/dev/runtime_inputs"

export LD_LIBRARY_PATH=/usr/mic/post1/metagenomics/lib
#export LD_LIBRARY_PATH=/p/lscratchc/allen99/metag/metag_repo/rel/LMAT-1.1/third-party/perm-je-0.9.3/lib

#nullm=$tdir/cutoffs.kcnt.m9.100.flst
#nullm=$tdir/cutoffs.gckcnt.m9.100.flst
nullm=$tdir/null_lst.txt
descript=
dbfile=/local/ramfs/m9.db
genedbfile=/local/ramfs/genedb
depthf="$tdir/depth_for_ncbi_taxonomy.segment_20130204.dat"
taxtree="$tdir/ncbi_taxonomy.segment_20130204.dat.nohl"
taxfile="$METAG_DIR/dev/runtime_inputs/ncbi_taxonomy_rank.segment_20130204.txt"
species_map="$tdir/species_map.ncbi.segment_2130204.dat"
thresh=500
fqfile=""
do_rl=1 # default is run everything, for debugging, may want to run parts
do_gl=1
sdiff=0.05
min_kmer=35
verbose=0
hbias=0.0
query_file=""
idir=`pwd`
odir=`pwd`
prune_depth=0
use_moab=0
odir=.
usage="Generate freq file
Usage: $0 options 

option list:
   --rl_off 
   --gl_off
   --db_file=$dbfile
   --genedb_file=$genedbfile
   --query_file=$query_file
   --fqfile=$fqfile
   --filter=$thresh
   --threads=$threads
   --nullm=$nullm
   --verbose=$verbose
   --odir=$odir
   --sdiff=$sdiff
   --descrip=$descrip
   --hbias=$hbias

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
   --genedb_file=*)
      genedbfile=$optarg;;
   --query_file=*)
      query_file=$optarg;;
   --filter=*)
      thresh=$optarg;;
   --threads=*)
      threads=$optarg;;
   --nullm=*)
      nullm=$optarg;;
   --fqfile=*)
      fqfile=$optarg;;
   --odir=*)
      odir=$optarg;;
   --sdiff=*)
      sdiff=$optarg;;
   --hbias=*)
      hbias=$optarg;;
   --verbose)
      verbose=1;;
   --gl_off)
      do_gl=0;;
   --rl_off)
      do_rl=0;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

dbname=`basename $dbfile`
query_file_name=`basename $query_file`
ofile="$odir/$query_file_name.$dbname.$thresh.lo.rl_output" 
cfiles=$nullm
logfile="$ofile.log" 

vstr=""
if [ $verbose == 1 ] ; then
   vstr="-y"
fi

min_score=0
if [ $do_rl == 1 ] ; then
   /usr/bin/time -v $bin/read_label_v3 -x $min_score -j $min_kmer -u $taxfile -w -r $species_map -l $hbias -b $sdiff -h $thresh $vstr -n $cfiles -e $depthf -p -t $threads -i $query_file -d $dbfile -c $taxtree -o $ofile >& $logfile 
fi
exit 0
lst=$ofile.flst
if [ -e $lst ] ; then
   rm -f $lst
fi
if [ $verbose == 1 ] ; then
   exit 0
fi
oname="$query_file_name.$dbname.$thresh.lo.rl_output" 
for lfile in `find . -maxdepth 1 -name $oname\*out` ; do
   echo $lfile >> $lst
done
mk=35
if [ $do_rl == 1 ] ; then
   losummary_fast_mc.sh --los=1 --file_lst=$lst --threads=$threads --min_read_len=$mk
   losummary_fast_mc.sh --los=0 --file_lst=$lst --threads=$threads --min_read_len=$mk
   losummary_fast_mc.sh --los=-0.5 --file_lst=$lst --threads=$threads --min_read_len=$mk
   losummary_fast_mc.sh --los=-1 --file_lst=$lst --threads=$threads --min_read_len=$mk
fi
exit 0

## assign gene names
if [ $genedbfile ] ; then
   if [ $do_gl == 1 ] ; then
      gene_label.sh --db_file=$genedbfile --query_lst=$lst --threads=80
   fi
fi
los=1
## identify SNPs with near neighbors
if [ -e $fqfile ] ; then
   lname=`basename $lst` 
   oname=`basename $ofile` 
   lmat2snp.sh --lmat_out_lst=$lname --lmatfsfile=$oname.out.$los.$mk.fastsummary --fqfile=$fqfile
else
   echo "fastq file must be given $fqfile not valid"
fi
