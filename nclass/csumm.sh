#!/bin/sh -xvf
### default options
rank_lst="plasmid species genus family"
rfile="$LMAT_DIR/ncbi_taxid_to_rank.pruned.txt"
pfile="$LMAT_DIR/low_numid_plasmids.txt"
usage="Generate bam file
Usage: $0 options

option list:
   --fsumm="" (default)
   --rank_lst=$rank_lst (default)
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --fsumm=*)
      fsumm=$optarg;;
   --rank_lst=*)
      rank_lst=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

bindir=""
if ! hash summary.py >& /dev/null; then
   #echo "Warning could not find $prog will try $LMAT_DIR/../bin"
   bindir="$LMAT_DIR/../bin/"
fi

summ="$fsumm.summ"
if [ ! -e $summ ] ; then
   echo "warning could not find [$fsumm.summ] file, cannot continue"
   exit 1
fi

ofile=$fsumm.ordered
${bindir}summary.py $summ $rfile $fsumm $pfile $ofile "$rank_lst"
