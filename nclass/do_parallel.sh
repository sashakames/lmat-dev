
inf=`basename $2`
dirbase=`basename $1`



length=`head -n10 $1/$inf*0.out  | grep "ReadTooShort"| grep length | head -n1 | awk '{print $3;}'`

echo $length

if [ ! -e $length ] ; then
    base=4
else
    base=3
fi

echo base $base


for n in 1 2 3 ; do

    ls $1/$inf*.out | $HOME/bin/parallel sh postproc-hg$n.sh  {.}.out $1/$dirbase $base
done





