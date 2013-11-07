#!/bin/sh -xvf
### default options
file_lst=""
min_kmers=1
min_read_len=-1
prog=gene_summary.pl
los=0.1
num_threads=2
#taxfile="$LMAT_DIR/gene_info"
taxfile="$LMAT_DIR/gn_ref2.txt"
usage="Summarize gene call file
Usage: $0 options

option list:
   --kmers=$min_kmers (default)
   --los=$los (default)
   --file_lst=$file_lst (default)
   --taxfile=$taxfile (default)
   --min_read_len=$min_read_len (default)
   --threads=$num_threads (default)
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --kmers=*)
      min_kmers=$optarg;;
   --los=*)
      los=$optarg;;
   --file_lst=*)
      file_lst=$optarg;;
   --min_read_len=*)
      min_read_len=$optarg;;
   --taxfile=*)
      taxfile=$optarg;;
   --threads=*)
      num_threads=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

false=0
user=`whoami`
while read file ; do
   echo "$file | $prog $taxfile $los $min_kmers $min_read_len"
   echo $file | $prog $taxfile $los $min_kmers $min_read_len &

   doWait=1
   while [ $doWait -eq 1 ] ; do
      runcnt=`ps x | grep $prog |wc -l`
      if [ "$runcnt" -gt "$num_threads" ] ; then
         doWait=1
         sleep 1
      else
         doWait=0
      fi
   done
done < $file_lst
doWait=1
while [ $doWait -eq 1 ] ; do
   runcnt=`ps x | grep $prog |wc -l`
   ## there will always be 1 for the grep
   if [ "$runcnt" -gt 1 ] ; then
      doWait=1
      sleep 1
   else
      doWait=0
   fi
done

lst=""
save=""
while read file ; do
   lst="$file.$los.$min_read_len.genesummary $lst"
   save="$file.$los.$min_read_len.genesummary"
done < $file_lst
echo "file_lst $file_lst"

out=`echo $save | sed 's/output[0-9]*/output/g'`
echo "$lst" | combine_gs.pl | sort -k1nr,1nr > $out
