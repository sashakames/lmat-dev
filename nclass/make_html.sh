if [ $# -lt 3 ] ;
then
    echo "html creation wrapper (convenient for use with gnu parallel)"
    echo "Usage: make_html.sh <species-file> <path-to-file> <output-dir>"
    exit
fi 

dn=$2
bn=`basename $1 .species`

python genusspecies2html.py $dn/$1 $dn/$bn.genus $LMAT_DIR/ncbi_taxonomy_rank.segment.pruned.txt > $3/$1.html

