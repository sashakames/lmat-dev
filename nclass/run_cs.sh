#!/bin/sh 

#################################################################################################################
### SCRIPT TO RUN THE content summary PIPELINE
### 
### The primary purpose of this procedure is to generate k-mer distribution plots for each species and plasmid
### to better evaluate relative abundance of individual organism calls. k-mer distribution plots can 
### be produced for any named rank, but the default is for species, plasmid and genus calls
###
### k-mer distribution plots are tradditionally used used to estimate coverage for the single organism case
### and this becomes a fundamentally hard problem for the metagenomics case since
### the correct choice of k depends on two unknown variables
### relative abundance of the organism (presumably unknown) and the percentage of the organism's genome that is rank
### specific (this can potentially be calculated but is not currently done, and will likely continue to change
### as new organisms are sequenced).
###
### Therefore a range of k are computed for user evalutation. By default k-mer distributions are calculated for 
### k=8,10,12,14 and 17.
###
### A summary file is produced for the species level showing the top strain match (where available)
### and the "second peak" in the k-mer distribution plot, assumed to estimate coverage, or -1 when
### a second peak does not exist
###
### Script requires that the LMAT read label process was run 
###
### 2) Call content_summ to generate summary of organisms present using read_label output
###
#################################################################################################################
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
if hash content_summ.1.3 >& /dev/null ; then
   bin_dir=

elif [ -f content_summ.1.3 ] ; then
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
#####################################################
## DEFAULT PARAMETER SETTINGS FOR CLASSIFICATION
#####################################################
## Memory mapped taxonomy database file
dbfile=""
markerdb=""

## NCBI taxonomy tree given in LMAT format
## this version removes most of the intermediate human lineage nodes
## to save compute time
#taxtree="$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat.nohl"
noprune_taxtree=$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat

taxfile="$LMAT_DIR/ncbi_taxonomy_rank.segment.pruned.txt"

## Used by content summarization to map taxids to rank (can this be merged with $species_map?)
rankval="$LMAT_DIR/ncbi_taxid_to_rank.pruned.txt"

## For content summarization, and read counting, ignore reads with scores below this threshold
min_score=0
marker_min_score=0

## Working to deprecate this with improved null models, it increases the human tax scores by 1 standard deviation
#hbias=1.5
hbias=1.5


## ignore reads with less valid k-mers than this value
min_read_kmer=30

## Additional user input default settings
# 1== run content_summ, 0 = skip
do_cs=1 
# number of threads
threads=80
# set to 1 for debugging only (too much output for large runs)
verbose=0
# Must be specified by user on command line
query_file=
# specify directory to place output
odir=  
sdir=

# in some cases you might not want to bother summarizing the human reads
# turn on if there are a huge # of human reads and you want to save time
skipHumanSumm=0
skipMarkSumm=0

# run content summary on marker
usage="Run LMAT pipeline 
Usage: $0 options 

option list:
   --db_file=$dbfile : Memory mapped database file
   --marker_library=$markerdb : Memory mapped marker database file
   --query_file=$query_file : Metagenomic reads in fasta format
   --threads=$threads : Number of threads to use 
   --verbose=$verbose : Only used for debugging single read queries (too much output for larger datasets)
   --sdir=$sdir : Place output in this directory (defaults to current)
   --odir=$odir : Place output in this directory (defaults to current)
   --min_score=$min_score : minimum score assigned to read for it to be included in binning
   --marker_min_score=$marker_min_score : minimum score assigned to read for it to be included in binning for the marker library
   --skipHumanSumm (default=$skipHumanSumm) : In summarization step, there's an option to pass over human reads to speed up microbial content summarization
   --skipMarkSumm (default=$skipMarkSumm) : Turn off summarization of marker output (may be preferred when running the full database)
   --overwrite (default=$overwrite) : overwrite output file if it exists 
   --min_read_kmer (default=$min_read_kmer) : minimum number of valid k-mers present in read needed for analysis
   --prune_thresh : threshold of maximum taxonomy IDs allowed per k-mer. 
   --pipe_cmd=$pipe_cmd : provide a pipe cmd instead of query_file, example: cat file.fastq.gz | gunzip | seqtk -A |

example usage:
$0 --db_file=$dbfile --genedb_file=$genedbfile --query_file=HC1.fna --threads=$threads

"

if test $# = 0; then
   echo "${usage}"
   exit 1
fi

while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --min_score=*)
      min_score=$optarg;;
   --marker_min_score=*)
      marker_min_score=$optarg;;
   --skipHumanSumm*)
      skipHumanSumm=1;;
   --skipMarkSumm*)
      skipMarkSumm=1;;
   --hbias=*)
      hbias=$optarg;;
   --sdir=*)
      sdir=$optarg;;
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
   --verbose)
      verbose=1;;
   --overwrite)
      overwrite=1;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

if [ ! -e $db_file ] && [ ! -e $markerdb ] ; then
   echo "Error need to supply a markery library or full database file"
   exit 0
fi
query_file_name=`basename $query_file`
## plasmids that were not associated with a chromosome did not 
## get special taxids, this file identifies these plasmids
## so that during content summarization, the minimum coverage parameter for
## plasmids is correctly applied
xtra_plas_file="$LMAT_DIR/low_numid_plasmids.txt"
## likely should always be close to 0, this is only used in XML summary file now
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
if [ ! $sdir == '' ]; then
    sdir="$sdir/"
fi

dlst="$markerdb $dbfile"
for db in $dlst ; do 
      dbname=`basename $db`
      if test -z $PTHRESH; then
	      rlsfile="${sdir}$query_file_name.$dbname.lo.rl_output" 
	      rlofile="${odir}$query_file_name.$dbname.lo.rl_output" 
      else
	      rlsfile="${sdir}$query_file_name.$dbname.$PTHRESH.lo.rl_output"	 
	      rlofile="${odir}$query_file_name.$dbname.$PTHRESH.lo.rl_output"	 
      fi
      logfile="$rlofile.log" 
      ## File giving a list of null models - assumes this specific naming convention

      ## The higher the number the more conservative the read label call
      ## This value specifices how much higher (in standard deviation units) the score of the assigned label must be
      sdiff=1.0  
      if [ $db != "$markerdb" ] ; then
         ## use full library parameters
         use_min_score=$min_score
         echo "search full db: $db"
      else 
         ## use marker library parameters
         use_min_score=$marker_min_score
         echo "search marker library: $db"
      fi

      fastsum_file="$rlsfile.$use_min_score.$min_read_kmer.fastsummary"
      echo "check: $fastsum_file"
      if [ -e $fastsum_file ] ; then
         lst=$rlofile.flst
         if [ -e $lst ] ; then
             rm -f $lst
         fi
         counter=0
         while [ $counter -lt $threads ] ; do
            lfile=${rlsfile}${counter}.out
            echo $lfile >> $lst
            let counter=counter+1
         done
         fastsum_ofile="$rlofile.$use_min_score.$min_read_kmer.fastsummary"
         sumofile="$fastsum_ofile.summ"
         logfile="$sumofile.log"
         kcnt="8,10,12,14,17"
         hstr=""
         ## in rarer cases where large samples are dominanted by humans
         ## it may be faster to not attempt to summarize the human read contents
         if [ $skipHumanSumm == 1 ] ; then
            hstr="-s"
         fi
         if [ ! -e $sumofile ] || [ $overwrite == 1 ] ; then
            echo "Summary process $query_file [overwrite=$overwrite (1=yes, 0=no)] [outputfile=$sumofile]"
            /usr/bin/time -v ${bin_dir}content_summ.1.3 $hstr -p $xtra_plas_file -c $noprune_taxtree -l $fastsum_file -k $kcnt -f $lst -r $rankval -o $sumofile >& $logfile
            ofile="$fastsum_ofile.ordered"
            ${bin_dir}summary.py $sumofile $rankval $fastsum_ofile $xtra_plas_file $ofile "plasmid species genus"
         fi
     fi
done
