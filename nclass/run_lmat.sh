#!/bin/sh

#####################################
### SCRIPT TO RUN THE LMAT PIPELINE
###
### Steps:
###
### 1) Call read_label to taxonomically label all reads, and count reads assigned to each label
### 2) Call content_summ to generate summary of organisms present using read_label output
### 3) call tolineage.py generates human readable taxonomy lineage, can be input to Krona, which is then run if the binary is found, to produce an html file
### 4) Call gene_label assigne gene names to each read and count the reads assigned to each label
###
#####################################
if [ -z "$LMAT_DIR" ] ; then
   echo "Please set LMAT_DIR environment variable to point to the directory where LMAT datafiles are stored"
   exit 1
fi

seqtk="/usr/gapps/kpath/seqtk-master/seqtk seq -A -n N -q" 

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
else 
   #echo "Could not find read_label in your path assume LMAT binaries/scripts are here: $bin_dir"
   bin_dir="$LMAT_DIR/../bin/"
fi

## Assume the perm-je library is here
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LMAT_DIR/../lib

overwrite=0
#####################################################
## DEFAULT PARAMETER SETTINGS FOR CLASSIFICATION
#####################################################
## Memory mapped taxonomy database file
dbfile=""
markerdb=""
## Memory mapped gene database file
genedbfile=/local/ramfs/gene.20mer.db

## NCBI taxonomy tree given in LMAT format
## this version removes most of the intermediate human lineage nodes
## to save compute time
taxtree="$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat.nohl"
## content_caller_qsum.py  uses the non pruned version
noprune_taxtree=$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat
## Used by gene_label to assign human readable names to genes
genefile="$LMAT_DIR/gn_ref2.txt.gz"
## Stores tree depth of each taxnomy node
depthf="$LMAT_DIR/depth_for_ncbi_taxonomy.segment.pruned.dat"

## Stores human readable lineage information on each taxonomy node
taxfile="$LMAT_DIR/ncbi_taxonomy_rank.segment.pruned.txt"

## Used by content summarization to map taxids to rank (can this be merged with $species_map?)
rankval="$LMAT_DIR/ncbi_taxid_to_rank.pruned.txt"

## For content summarization, and read counting, ignore reads with scores below this threshold
min_score=0
marker_min_score=0

## Working to deprecate this with improved null models, it increases the human tax scores by 1 standard deviation
#hbias=1.5
hbias=0


## ignore reads with less valid k-mers than this value
min_read_kmer=30

## Additional user input default settings
# 1== run content_summ, 0 = skip
do_cs=1 
# 1 == run read_label, 0 = skip 
do_rl=1 
# 1 == run gene_label, 0 = skip
do_gl=1 
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
skipHumanSumm=0
skipMarkSumm=0

# run content summary on marker
usage="Run LMAT pipeline 
Usage: $0 options 

option list:
   --rl_off  : skip the read label step
   --gl_off  : skip the gene label step
   --cs_off  : skip the content summarization step
   --db_file=$dbfile : Memory mapped database file
   --marker_library=$markerdb : Memory mapped marker database file
   --genedb_file=$genedbfile : Memory mapped gene database file
   --query_file=$query_file : Metagenomic reads in fasta format
   --threads=$threads : Number of threads to use 
   --nullm=$nullm : File containing the list of null models
   --verbose=$verbose : Only used for debugging single read queries (too much output for larger datasets)
   --sdiff=$sdiff : Scoring differential. Tax ids with scores at or above the the maximum score - sdiff*STDEV are considered
   --hbias=$hbias : For human samples where human DNA concentration is high, human taxid score + hbias*STDEV, hbias > 0 may be conveniant
   --odir=$odir : Place output in this directory (defaults to current)
   --min_score=$min_score : minimum score assigned to read for it to be included in binning
   --marker_min_score=$marker_min_score : minimum score assigned to read for it to be included in binning for the marker library
   --skipHumanSumm (default=$skipHumanSumm) : In summarization step, there's an option to pass over human reads to speed up microbial content summarization
   --skipMarkSumm (default=$skipMarkSumm) : Turn off summarization of marker output (may be preferred when running the full database)
   --overwrite (default=$overwrite) : overwrite output file if it exists 
   --min_read_kmer (default=$min_read_kmer) : minimum number of valid k-mers present in read needed for analysis
   --prune_thresh : threshold of maximum taxonomy IDs allowed per k-mer. 

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
   --qual_filter=*)
      qual_filter=$optarg;;
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
   --odir=*)
      odir=$optarg;;
   --db_file=*)
      dbfile=$optarg;;
   --min_read_kmer=*)
      min_read_kmer=$optarg;;
   --marker_library=*)
      markerdb=$optarg;;
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
   --sdiff=*)
      sdiff=$optarg;;
   --verbose)
      verbose=1;;
   --overwrite)
      overwrite=1;;
   --gl_off)
      do_gl=0;;
   --rl_off)
      do_rl=0;;
   --cs_off)
      do_cs=0;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

if [ ! -e $query_file ] ; then
   echo "Error $query_file not found"
   exit 0
fi
if [ ! -e $db_file ] && [ ! -e $markerdb ] ; then
   echo "Error need to supply a markery library or full database file"
   exit 0
fi
query_file_name=`basename $query_file`
## Content summarization parameters
##
## contains a list of genomes that may have contamination
## and should show a higher abundance before be including 
## in the summary table of called genomes
suspect_genomes=$LMAT_DIR/supsect_m9_genomes.txt
min_contam_filt=0.05
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

dlst="$markerdb $dbfile"
for db in $dlst ; do 
      dbname=`basename $db`
      if test -z $PTHRESH; then
	  rlofile="${odir}$query_file_name.$dbname.lo.rl_output" 
      else
	  rlofile="${odir}$query_file_name.$dbname.$PTHRESH.lo.rl_output"	 
      fi
      logfile="$rlofile.log" 
      tidmap="$LMAT_DIR/m9.32To16.map"
      ## File giving a list of null models - assumes this specific naming convention
      nullm=$LMAT_DIR/$dbname.null_lst.txt
      ## The higher the number the more conservative the read label call
      ## This value specifices how much higher (in standard deviation units) the score of the assigned label must be
      sdiff=0.5  
      fstr="-f $tidmap"
      if [ ! -z $PTHRESH ]; then
	  pstr="-g $PTHRESH -m $LMAT_DIR/numeric_ranks"
      fi
      rprog=${bin_dir}read_label
      if [ $db != "$markerdb" ] ; then
         ## use full library parameters
         use_min_score=$min_score
         echo "search full db: $db"
      else 
         ## use marker library parameters
         use_min_score=$marker_min_score

         sdiff=1.0
         echo "search marker library: $db"
      fi
      if [ -e $nullm ] ; then
         nullmstr="-n $nullm"
      else 
         nullmstr=""
      fi
      fastsum_file="$rlofile.$use_min_score.$min_read_kmer.fastsummary"
      if [ -f $db ] && [ $do_rl == 1 ] ; then
         if [ ! -e $fastsum_file ] || [ $overwrite == 1 ] ; then
            echo "Process $query_file [overwrite=$overwrite (1=yes, 0=no)] [outputfile=$fastsum_file]"

	    if [ -n $qual_filter ] ; then
		echo Quality filter Q$qual_filter on.
		gunzip -c $query_file | $seqtk $qual_filter - | $rprog $min_kmer_str $fstr $pstr -u $taxfile -x $use_min_score -j $min_read_kmer -l $hbias -b $sdiff $vstr $nullmstr -e $depthf -p -t $threads -i - -d $db -c $taxtree -o $rlofile >& $logfile
	    else
		/usr/bin/time -v $rprog $min_kmer_str $fstr $pstr -u $taxfile -x $use_min_score -j $min_read_kmer -l $hbias -b $sdiff $vstr $nullmstr -e $depthf -p -t $threads -i $query_file -d $db -c $taxtree -o $rlofile >& $logfile
	    fi

            min_reads=1
            if [ ! -e $fastsum_file ] ; then
               echo "Error, did not create a fastsummary file [$fastsum_file]"
               exit 0
            fi
            min_num_reads=10 ## 
            ${bin_dir}tolineage.py $taxfile $fastsum_file $fastsum_file.lineage $min_num_reads all

            if hash ktImportText > /dev/null 2>&1 ; then
               ktImportText $fastsum_file.lineage -o $fastsum_file.lineage.html
            fi
            if [ $verbose == 1 ] ; then
               echo "Verbose setting is used for debuging one program at a time only"
               exit 0
            fi
         fi
      else
	  echo "Did not find an LMAT db at $db."
      fi
      if [ $do_cs == 1 ] && [ -e $fastsum_file ] ; then
         lst=$rlofile.flst
         if [ -e $lst ] ; then
             rm -f $lst
         fi
         counter=0
         while [ $counter -lt $threads ] ; do
            lfile=${rlofile}${counter}.out
            echo $lfile >> $lst
            let counter=counter+1
         done
         sumofile="$fastsum_file.summ"
         logfile="$sumofile.log"
         ###########################################################################
         ## Below are parameters that may need to be tweaked differently depending
         ## on whether a marker library or full database is being used
         ###########################################################################
         ## minimum abundance of reads required before making a summarization call
         ## note you can still check read counts with no abundance filtering
         #virus
         min_virus_read=0.001
         # bacteria/archaea
         min_prok_read=0.0001
         # eukaryotes
         min_euk_read=0.0001
         # plasmids
         min_plas=0.01
         # minimum weight irrespective of call type
         # Can be used to filter out an organism with lots of low scoring reads
         min_wrdc=2
         # once one strain is called, this parameter specifies the minimum number of reads
         # required to specify a second strain
         #a standard value for collapsing strains
         #same_strain=2.0
         #default is to set this high and not collapse strain data
         #The issue now is that when re-assigning higher ranked reads, it is possible for the called strain
         #to not have all of the observed variants. Unless another strain remains as a viable candidate
         #there's a higher chance that the variant could be mis-assigned to another species.
         same_strain=0.0
         # Can be used to filter out low scoring singleton read calls if desired
         min_avg_wght=0.85
         ###########################################################################
         if [ $db == "$markerdb" ] ; then
            if [ $skipMarkSumm == 1 ] ; then
               continue
            fi
            ###########################################################################
            ## must change these settings for marker library
            ###########################################################################
            min_wrdc=0.07
            min_avg_wght=0.03
            ## Stores the k-mer counts associated with each taxid
            kmercnt="$LMAT_DIR/tcnt.kML18"
            ## I'm noticing slight over specificity in the marker library,so consider lower this value
            sdiff=0.5  
         else 
            ## Stores the k-mer counts associated with each taxid
            kmercnt="$LMAT_DIR/tcnt.m9.20.tax_histo"
         fi
         min_thresh="$min_euk_read $min_prok_read $min_virus_read $min_contam_filt $min_plas $min_wrdc $same_strain $min_avg_wght"
         kcnt=12
         hstr=""
         ## in rarer cases where large samples are dominanted by humans
         ## it may be faster to not attempt to summarize the human read contents
         if [ $skipHumanSumm == 1 ] ; then
            hstr="-s"
         fi
         if [ ! -e $sumofile ] || [ $overwrite == 1 ] ; then
            echo "Summary process $query_file [overwrite=$overwrite (1=yes, 0=no)] [outputfile=$sumofile]"
            /usr/bin/time -v ${bin_dir}content_summ $hstr -p $xtra_plas_file -b $suspect_genomes -t "$min_thresh" -c $noprune_taxtree -l $fastsum_file -k $kcnt -f $lst -r $rankval -m $kmercnt -o $sumofile >& $logfile
         fi
     fi

done

## assign gene names
if [ $genedbfile ] ; then
   if [ $do_gl == 1 ] ; then
      genedbname=`basename $genedbfile`
      lst=$rlofile.ras.flst
      if [ -e $lst ] ; then
         rm -f $lst
      fi
      counter=0
      while [ $counter -lt $threads ] ; do
         lfile=${rlofile}${counter}.out.ras.1
         echo $lfile >> $lst
         let counter=counter+1
      done
      vstr=""
      ## note need to fix verbose setting to get here (if needed)
      if [ $verbose == 1 ] ; then
         vstr="-y"
      fi
      genofile="$rlofile.ras.flst.$genedbname.rl_output"
      logfile="$rlofile.ras.flst.$genedbname.rl_output.log"
      res=$genofile.$gene_score.$num_gene_kmers.genesummary
      if [ ! -e $res ] || [ $overwrite == 1 ] ; then
         ${bin_dir}gene_label $vstr -q $num_gene_kmers -x $gene_score -p -r -l $lst -d $genedbfile -o $genofile -g $genefile >& $logfile
         cat $res | sort -k1gr,1gr > tmp.$$
         mv tmp.$$ $res 
      fi
   fi
fi
