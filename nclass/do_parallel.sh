ls $1/*.out | $HOME/bin/parallel sh postproc-hg1.sh  {.}.out 
ls $1/*.out | $HOME/bin/parallel sh postproc-hg2.sh  {.}.out 