#  Note that the number below depends on k, eg. 19 for 20-mers


if [ $# -lt 3 ] ; then

    echo $#
    echo "Usage $0 <path-to-files> <db-name>  <red-len-field>"
    echo "Where the field is separated bys .'s"
    exit 0
fi

cd $1

ls null.bin.10.$2.*.gz | sort -t '.' -k $3 -n > tmp
awk -v FS='.' -v rlf=$3 '{print $(rlf) - 19}' tmp > coltmp
paste -d ' ' coltmp tmp > $2.null_lst.txt
mv null.bin.10.$2.*.gz $2.null_lst.txt $LMAT_DIR
