#!/bin/sh  -xvf

#####################################
### SCRIPT TO RUN THE Gene labeling PIPELINE
###
### Steps:
###
### Call gene_label assigne gene names to each read and count the reads assigned to each label
###
#####################################
if [ -z "$LMAT_DIR" ] ; then
   echo "Please set LMAT_DIR environment variable to point to the directory where LMAT datafiles are stored"
   exit 1
fi

## Some environments require explicit enabling of hyperthreading
## Other environments may already enable this
if [ -e /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus ] ; then
   /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus 
   export GOMP_CPU_AFFINITY=0-79
fi

## Should improve this 
## Location of binaries
if hash read_label >& /dev/null ; then
   bin_dir=

elif [ -f read_label ] ; then
    bin_dir=./
elif [ `basename $PWD` == "nclass" ] ; then
    bin_dir="../apps/"
 else
   #echo "Could not find read_label in your path assume LMAT binaries/scripts are here: $bin_dir"
   bin_dir="$LMAT_DIR/../bin/"
fi
## Assume the perm-je library is here
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LMAT_DIR/../lib
pipe_cmd=""
overwrite=0
dbfile=""
markerdb=""
## Memory mapped gene database file
genedbfile=/local/ramfs/gene.20mer.db

## Used by gene_label to assign human readable names to genes
genefile="$LMAT_DIR/gn_ref2.txt.gz"

## ignore reads with less valid k-mers than this value
min_read_kmer=30

## Additional user input default settings
## Not yet supported
do_fasta=0 
# number of threads
threads=80
# set to 1 for debugging only (too much output for large runs)
verbose=0
# Must be specified by user on command line
query_file=
# specify directory to place output
odir=  

# in some cases you might not want to bother summarizing the human reads
# turn on if there are a huge # of human reads and you want to save time

# run content summary on marker
usage="Run LMAT pipeline 
Usage: $0 options 

option list:
   --fasta  : search a fasta file as input (not yet supported)
   --db_file=$dbfile : Memory mapped database file
   --marker_library=$markerdb : Memory mapped marker database file
   --genedb_file=$genedbfile : Memory mapped gene database file
   --query_file=$query_file : Metagenomic reads in fasta format
   --threads=$threads : Number of threads to use 
   --nullm=$nullm : File containing the list of null models
   --verbose=$verbose : Only used for debugging single read queries (too much output for larger datasets)
   --odir=$odir : Place output in this directory (defaults to current)
   --overwrite (default=$overwrite) : overwrite output file if it exists 
   --min_read_kmer (default=$min_read_kmer) : minimum number of valid k-mers present in read needed for analysis
   --prune_thresh : threshold of maximum taxonomy IDs allowed per k-mer. 
   --pipe_cmd=$pipe_cmd : provide a pipe cmd instead of query_file, example: cat file.fastq.gz | gunzip | seqtk -A |

example usage:
$0 --genedb_file=$genedbfile --query_file=HC1.fna --threads=$threads

"

if test $# = 0; then
   echo "${usage}"
   exit 1
fi

while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --odir=*)
      odir=$optarg;;
   --db_file=*)
      dbfile=$optarg;;
   --min_read_kmer=*)
      min_read_kmer=$optarg;;
   --marker_library=*)
      markerdb=$optarg;;
   --pipe_cmd=*)
      pipe_cmd=$optarg;;
   --genedb_file=*)
      genedbfile=$optarg;;
   --query_file=*)
      query_file=$optarg;;
   --prune_thresh=*)
      PTHRESH=$optarg;;
   --threads=*)
      threads=$optarg;;
   --nullm=*)
      nullm=$optarg;;
   --verbose)
      verbose=1;;
   --overwrite)
      overwrite=1;;
   --fasta)
      do_fasta=1;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

if [ $do_rl == 1 ] && [! -e $query_file ] ; then
   echo "Error $query_file not found"
   exit 0
fi
if [ ! -e $db_file ] && [ ! -e $markerdb ] ; then
   echo "Error need to supply a markery library or full database file"
   exit 0
fi
query_file_name=`basename $query_file`
min_gene_read=0
## minimum percentage of k-mer matches from a read found in a gene before a call can be made
gene_score=0.01
## minimum number of valid k-mers  present in read needed before considering the read for gene labeling
num_gene_kmers=20

vstr=""
if [ $verbose == 1 ] ; then
   vstr="-y"
fi

if [ ! $odir == '' ]; then
    odir="$odir/"
fi

## assign gene names
if [ $genedbfile ] ; then
      dbname=`basename $dbfile`
      if test -z $PTHRESH; then
	      rlofile="$query_file.$dbname.lo.rl_output" 
      else
	      rlofile="$query_file.$dbname.$PTHRESH.lo.rl_output"	 
      fi
      genedbname=`basename $genedbfile`
      rlofile_name=`basename $rlofile`
      lst=$rlofile_name.flst
      if [ $do_fasta == 1 ] ; then
         lstr=""
         qstr="-i $query_file"
      else
         if [ -e $lst ] ; then
             rm -f $lst
         fi
         counter=0
         while [ $counter -lt $threads ] ; do
            lfile=${rlofile}${counter}.out
            echo $lfile >> $lst
            let counter=counter+1
         done
         lstr="-l $lst"
         qstr=""
      fi

      vstr=""
      ## note need to fix verbose setting to get here (if needed)
      if [ $verbose == 1 ] ; then
         vstr="-y"
      fi
      genofile="${odir}$query_file_name.$genedbname.rl_output"
      logfile="${odir}$query_file_name.$genedbname.rl_output.log"
      res=$genofile.$gene_score.$num_gene_kmers.genesummary
      if [ ! -e $res ] || [ $overwrite == 1 ] ; then
         ${bin_dir}gene_label $vstr -q $num_gene_kmers -x $gene_score -p -r $lstr $qstr -d $genedbfile -o $genofile -g $genefile >& $logfile
         cat $res | sort -k1gr,1gr > tmp.$$
         mv tmp.$$ $res 
      fi
fi
