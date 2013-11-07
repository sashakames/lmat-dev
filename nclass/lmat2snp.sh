#!/bin/sh -xvf
### default options
lmat_out_lst=""
min_seqlen=50
top_nn=30
mincov=3
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
   --mincov=*)
      mincov=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

odir=`pwd`-$lmatfsfile
echo "use for odir: $odir"
num_cpu=12
Variant_finder -r $idir -o $odir -f $fqfile -l $maptype -q $qval -m $mapscore -p $num_cpu -s $los -c $mincov

lname=`basename $lmatfsfile`
lst=$lmatfsfile.snp_lst.flst
if [ -e $lst ] ; then
   rm -f $lst
fi
for lfile in `find $odir -maxdepth 1 -name \*.$lname.content_call.xml` ; do
   echo $lfile >> $lst
done
