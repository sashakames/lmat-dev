#!/bin/sh  -xvf

#create 
db_file=/local/ramfs/m9.db
read_len=0
read_range=0
#num_reads=125000
num_reads=60000
tax_histo_cnt=../runtime_inputs/tcnt.m9.20.tax_histo
min_sample_size=100
filter=500
usage="Generate random null model 

Usage: $0 options 

option list:
   --db_file=$db_file
   --read_len=$read_len
   --read_range=$read_range (ex: 50:1000:20)
   --filter=$filter
   --num_reads=$num_reads (defualt)
   --tax_histo_cnt=$tax_histo_cnt
   --min_sample_size=$min_sample_size

example usage:
$0 --db_file=$db_file --read_len=80
WARNING MUST MANUALLY SET SCRIPT FOR CORRECT OUTPUT FORMAT (OLD OR NEW) current=$SCRPT

"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
debug=0
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --debug*)
      debug=1;;
   --db_file=*)
      db_file=$optarg;;
   --read_len=*)
      read_len=$optarg;;
   --read_range=*)
      read_range=$optarg;;
   --filter=*)
      filter=$optarg;;
   --num_reads=*)
      num_reads=$optarg;;
   --min_sample_size=*)
      min_sample_size=$optarg;;
   --tax_histo_cnt=*)
      tax_histo_cnt=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

if [ $read_len = 0 ] &&  [ $read_range = 0 ] ; then
   echo "${usage}"
   exit 1
elif [ $read_len -gt 0 ]; then
   beg=$read_len
   end=$read_len
   int=100
else
   beg=`echo $read_range | awk '{ split($0,a,":"); print a[1] }'` 
   end=`echo $read_range | awk '{ split($0,a,":"); print a[2] }'` 
   int=`echo $read_range | awk '{ split($0,a,":"); print a[3] }'` 
   echo "what $read_range"
fi
threads=80
tot_reads=$[num_reads*threads]
db_file_name=`basename $db_file`
echo "check: $beg $end $int $tot_reads"
while [ $beg -le $end ] ; do
   read_len=$beg
   #rand_label_v2.sh --db_file=$db_file --num_reads=$num_reads --read_len=$read_len --threads=$threads --filter=$filter
   oname="$db_file_name.$filter.$read_len.$tot_reads.rl_output.rand_lst"
   if [ -e $oname ] ; then
      merge_cnts.sh $oname $min_sample_size > log.$oname
      if [ $debug == 0 ] ; then
         gzip null.bin.10.$oname
      fi
   else 
      echo "warning no $oname found"
   fi
   let "beg += int"
done
exit 0
ls -1 null.bin.*rand_lst.gz | perl -ne 'if(/\d+\.(\d+)\.\d+\.rl_output/) { $t=$1-19; print "$t $_";}' | sort -k1n,1n > null_lst.txt
cp *rand_lst.gz ../runtime_inputs
cp null_lst.txt ../runtime_inputs
