#!/bin/sh
### default options
file_lst=""
num_jobs=200
usage="Simple general utility to manage moab job submissions size > 200 (to adhere to LC's 200 job submission limit)
Usage: $0 options

option list:
   --num_jobs=$num_jobs (default)
   --file_lst=$file_lst (default)
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --num_jobs=*)
      num_jobs=$optarg;;
   --file_lst=*)
      file_lst=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done
user=`whoami`
while read file ; do
   msub $file
   doWait=1
   while [ $doWait -eq 1 ] ; do
      runcnt=`showq -w user=allen99 | grep allen99 | wc -l`
      if [ "$runcnt" -ge "$num_jobs" ] ; then
         doWait=1
         sleep 5
      else
         doWait=0
      fi
   done
done < $file_lst
