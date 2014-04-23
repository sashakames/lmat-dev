#!/bin/sh 
### default options
file_lst=""
adapt_file=/usr/mic/bio/blastdb/vector/all-illumina-adapters.fa
useMcf=0
odir=.
num_threads=0
isPe=0
usage="Filtering steps for Illumina reads
Usage: $0 options

option list:
   --file_lst=$file_lst (default)  : list of fastq files 
   --threads=$num_threads (default) : default is to use all cores - customize the call to gnu parallel to limit resource use.
   --pe (default=off) : handle paired end reads
   --mcf (default=$useMcf) : use fastq-mcf
   --odir (default=$odir) : output directory
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi

## Some environments require explicit enabling of hyperthreading
## Other environments may already enable this
if [ -e /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus ] ; then
   /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus
   export GOMP_CPU_AFFINITY=0-79
fi

bindir=""
if ! hash trim_galore >& /dev/null; then
   echo "Warning could not find trim_galore will try $LMAT_DIR/../bin"
   bindir="$LMAT_DIR/../bin/"
fi

pbin=""
if ! hash parallel >& /dev/null; then
   echo "Warning could not find GNU parallel will try $LMAT_DIR/../bin"
   pbin="$LMAT_DIR/../bin/"
fi
   
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --file_lst=*)
      file_lst=$optarg;;
   --threads=*)
      num_threads=$optarg;;
   --pe*)
      isPe=1;;
   --mcf*)
      useMcf=1;;
   --odir=*)
      odir=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

jobstr=""
if [ $num_threads -gt 0 ] ; then
   jobstr="--jobs $num_threads"
fi

tfile="$file_lst.tmp.$$" 
if [ -e $tfile ] ; then
   rm -f $tfile
fi

if [ $isPe -eq 1 ] ; then
   echo "Run paired end protocol"
fi
t=1
if [ $t -eq 1 ]; then

if [ $isPe == 0 -a $useMcf == 1 ] ; then
   while read file ; do
      fname=`basename $file`
      echo "fastq-mcf $adapt_file $file -o $odir/fc.$fname" >> $tfile
   done < $file_lst
   ${pbin}parallel --load 150% --progress $jobstr -a $tfile 

   if [ -e $tfile ] ; then
      rm -f $tfile
   fi
   while read file ; do
      fname=`basename $file`
      rfile="$odir/fc.$fname"
      echo "$rfile" >> $tfile
   done < $file_lst
   ${pbin}parallel --load 150% --progress $jobstr -a $tfile /usr/gapps/kpath/seqtk-master/seqtk seq -A -q 10 -n N {} ">" {.}.fasta

   exit 0
fi
if [ -e $tfile ] ; then
   rm -f $tfile
fi
while read file ; do
   echo "fastq_illumina_filter -N -o fc.$file $file" >> $tfile
   if [ $isPe == 1 ] ; then
      f2=`echo "$file" | perl -ne 'chomp; s/R1/R2/g; print;'`
      echo "fastq_illumina_filter -N -o fc.$f2 $f2" >> $tfile
   fi
done < $file_lst

${pbin}parallel --load 150% --progress $jobstr -a $tfile 
   

if [ -e $tfile ] ; then
   rm -f $tfile
fi

while read file ; do
   file="fc.$file"
   if [ $isPe == 1 ] ; then
      f1=$file
      f2=`echo "$file" | perl -ne 'chomp; s/R1/R2/g; print;'`
      echo "trim_galore --retain_unpaired --phred33 --fastqc --paired $f1 $f2" >> $tfile
   else
      echo "$file" >> $tfile
   fi
done < $file_lst
if [ $isPe == 1 ] ; then
   ${pbin}parallel --load 150% --progress $jobstr -a $tfile 
else 
   ${pbin}parallel --load 150% --progress $jobstr -a $tfile trim_galore --phred33 --fastqc
fi

fi
if [ -e $tfile ] ; then
   rm -f $tfile
fi

if [ $isPe == 1 ] ; then
   while read file ; do
      file="fc.$file"
      op1=`echo "$file" | perl -ne 'chomp; s/\.fastq$/_val_1.fq/; print;'` 
      op2=`echo "$file" | perl -ne 'chomp; s/_R1_/_R2_/g; s/\.fastq$/_val_2.fq/; print;'` 
      pre_op3=`echo "$file" | perl -ne 'chomp; s/_R1_/_/g; s/\.fastq$/\.pre\.fastq/; print;'` 
      echo "$op1 $op2 $pre_op3"  >> $tfile
   done < $file_lst

   ${pbin}parallel --load 150% --progress $jobstr -a $tfile merge_fastq_reads_with_N_separator.pl

   if [ -e $tfile ] ; then
      rm -f $tfile
   fi
   while read file ; do
      file="fc.$file"
      up1=`echo "$file" | perl -ne 'chomp; s/\.fastq$/_unpaired_1.fq/; print;'`
      up2=`echo "$file" | perl -ne 'chomp; s/_R1_/_R2_/g; s/\.fastq$/_unpaired_2.fq/; print;'`
      pre_op3=`echo "$file" | perl -ne 'chomp; s/_R1_/_/g; s/\.fastq$/\.pre\.fastq/; print;'` 
      op3=`echo "$file" | perl -ne 'chomp; s/_R1_/_/g; s/\.fastq$/\.proc\.fastq/; print;'` 
      if [ -e $op3 ]; then
         rm -f $op3
      fi
      echo "cat $up1 $up2 $pre_op3 > $op3" >> $tfile
   done < $file_lst

   while read file ; do
      file="fc.$file"
      pre_op3=`echo "$file" | perl -ne 'chomp; s/_R1_/_/g; s/\.fastq$/\.pre\.fastq/; print;'` 
      #rm -f $pre_op3
   done < $file_lst

   ${pbin}parallel --load 150% --progress $jobstr -a $tfile 
fi

if [ -e $tfile ] ; then
   rm -f $tfile
fi

while read file ; do
   file="fc.$file"
   op3=`echo $file | perl -ne 'chomp; s/_R1_/_/g; s/\.fastq$/\.proc\.fastq/; print;'` 
   echo "$op3" >> $tfile
   fa=`echo $op3 | perl -ne 'chomp; s/_R1_/_/g; s/\.fastq$/\.fasta/; print;'` 
   rm -f $fa
done < $file_lst

${pbin}parallel --load 150% --progress $jobstr -a $tfile /usr/gapps/kpath/seqtk-master/seqtk seq -A -q 10 -n N {} ">" {.}.fasta

if [ -e $tfile ] ; then
   rm -f $tfile
fi
