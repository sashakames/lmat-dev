#!/bin/sh -xvf
### default options
lmat_out_lst=""
min_kmers=25
min_seqlen=50
top_nn=100
los=0
threads=40
idfile=""
lmatfsfile=""
idir=`pwd`
## seems to need to be in another directory but specific for each sample
odir="" 
fqfile=""
maptype=long
qval=20
mapscore=30
usage="Generate SNP data 
Usage: $0 options


option list:
   --top_nn=$top_nn (default) : attempt to find SNPs for X most abundant near neighbors
   --min_seqlen=$min_seqlen
   --min_kmers=$min_kmers
   --los=$los (default)
   --lmat_out_lst=$lmat_out_lst (default)
   --lmatfsfile=$lmatfsfile (default)
   --threads=$threads (default)
   --idir=$idir (default)
   --odir=`pwd` + fastsummary file name (default)
   --fqfile=$fqfile (default)
   --maptype=$maptype (default)
   --qval=$qval (default)
   --mapscore=$mapscore (default)
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --top_nn=*)
      top_nn=$optarg;;
   --min_seqlen=*)
      min_seqlen=$optarg;;
   --min_kmers=*)
      min_kmers=$optarg;;
   --los=*)
      los=$optarg;;
   --lmat_out_lst=*)
      lmat_out_lst=$optarg;;
   --lmatfsfile=*)
      lmatfsfile=$optarg;;
   --threads=*)
      threads=$optarg;;
   --idir=*)
      idir=$optarg;;
   --odir=*)
      odir=$optarg;;
   --fqfile=*)
      fqfile=$optarg;;
   --mtype=*)
      mtype=$optarg;;
   --qval=*)
      qval=$optarg;;
   --mapscore=*)
      mapscore=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

if [ !$lmat_out_lst ]; then
   lmat_out_lst=`echo $lmatfsfile | sed -e 's/out\..*\..*\.fastsummary/flst/'`
fi
get_reads_auto.pl $lmatfsfile.content_call $lmatfsfile.nn_map $lmatfsfile.remap_id $lmat_out_lst $threads $min_kmers $los $min_seqlen $top_nn 
