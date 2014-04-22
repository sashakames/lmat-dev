dn=$2
bn=`basename $1 .species`

python genusspecies2html.py $dn/$1 $dn/$bn.genus $LMAT_DIR/ncbi_taxonomy_rank.segment.pruned.txt > $3/$1.html

