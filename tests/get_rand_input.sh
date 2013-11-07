#!/bin/sh

tlabel=""
rl=100
max_seq=5000000

usage="Generate test data
Usage: $0 options 

option list:
   --test_name=NAME (no default) 
   --read_len=$rl (default)
   --max_seq=$max_seq (default) Maximum amount of genome to break into reads (takes first X reads)
   --help :  print usage
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`
   case $opt in
   --test_name=*)
      tlabel=$optarg;;
   --read_len=*)
      rl=$optarg;;
   --max_seq=*)
      max_seq=$optarg;;
   --help)
      echo "${usage}"
      exit 0;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

if [ ! "$tlabel" ] ; then
   echo "need to provide a for test set"
   exit 0
fi

echo "input parameters: test_name=$tlabel read_len=$read_len max_seq=$max_seq"

rm -f *.rout *.reads *.query
rm -f $tlabel.random_vf.flst $tlabel.random_vnf.flst $tlabel.random_ff.flst $tlabel.random_bf.flst $tlabel.random_bnf.flst
find /usr/mic/bio/shea/2011Dec20/ViralFamilies_clustered -name \*.fasta | sel_rand_fn.pl 5 > $tlabel.random_vf.flst
find /usr/mic/bio/shea/2011Dec20/Viral_noFamily_clustered -name \*.fasta | sel_rand_fn.pl 2 > $tlabel.random_vnf.flst
find /usr/mic/bio/shea/2011Dec20/FungiFamilies -name \*.fasta | sel_rand_fn.pl 3 > $tlabel.random_ff.flst 
find /usr/mic/bio/shea/2011Dec20/BacteriaFamilies -name \*.fasta | sel_rand_fn.pl 5 > $tlabel.random_bf.flst 
find /usr/mic/bio/shea/2011Dec20/Bacteria_noFamily -name \*.fasta | sel_rand_fn.pl 2 > $tlabel.random_bnf.flst 

function gen_dat {
   ilst=$1
   rl=$2
   max_seq=$3
   while read raw_list_fasta ; do
      list_fasta=`basename $raw_list_fasta`
      rm -f $list_fasta.rout
      get_fa_rand.pl $raw_list_fasta $list_fasta.rout 10
      ## save the first entry (random choice) to use as a query
      get_fa_idx.pl $list_fasta.rout 0 $list_fasta.rout.query 
      make_simple_read.pl $list_fasta.rout.query $list_fasta.rout.query.reads $rl $max_seq
   done < $ilst
}

gen_dat $tlabel.random_vf.flst 10 $rl $max_seq
gen_dat $tlabel.random_vnf.flst 5 $rl $max_seq
gen_dat $tlabel.random_bf.flst 10 $rl $max_seq
gen_dat $tlabel.random_bnf.flst 5 $rl $max_seq
gen_dat $tlabel.random_ff.flst 1 $rl $max_seq

cat *.rout > $tlabel.ref.fasta
cat *.rout.query.reads > $tlabel.query.reads.all.fasta
