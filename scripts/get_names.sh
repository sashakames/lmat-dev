PATH_TO_INPUT_DATA=/p/lscratche/ames4/snglorg/orig
for f in $PATH_TO_INPUT_DATA/*fasta; do
    n=${f%.fasta}
    n=${n##*/}

#    fN=`wc -l $f|sed 's/ .*//'`
    echo $n
    done