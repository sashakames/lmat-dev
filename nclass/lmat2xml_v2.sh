#!/bin/sh 

qr=$1
tdb=$2
filt=$3
gdb=$4 #genedb
call_los=$5
read_los=$6
mrl=$7
hdr=$8
min_gene_read=$9
gene_score=${10}
min_kmer=${11}
echo "gs: $gene_score, $min_kmer"
cidmap=$LMAT_DIR/microbe3.20130207.headers.mapping.ncbionly.multi.added_human.custom_ids
taxcall=$qr.$tdb.$filt.lo.rl_output.$call_los.$mrl.fastsummary.summ.lin
gencall=$qr.$tdb.$filt.lo.rl_output.ras.flst.$gdb.lo.rl_output.$gene_score.$min_kmer.genesummary
rlfile=$qr.$tdb.$filt.lo.rl_output.ras.flst
ofile=$qr.lmat.xml
lmat2xml_v2.pl $taxcall $gencall $rlfile $ofile $hdr $cidmap $min_gene_read $read_los
