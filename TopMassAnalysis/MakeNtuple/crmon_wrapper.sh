#!/bin/bash

date="20140223"

crabdir[0]="crab_TopMassNtuple_data_${date}"
#crabdir[1]="crab_TopMassSkim_mc_${date}"

count=0
end_job=""
end_time=$(( $(date +%s) + 60*60*24 ))
while [ $(date +%s) -le $end_time ] && [ "$end_job" != "yes" ]
do

   for dir in ${crabdir[*]}
   do
      if [ $count == 0 ] ; then
         multicrab -status -c ${dir}
         multicrab -getoutput -c ${dir}
      fi
      for subdir in `ls ${dir}`
      do
         if [ -d ${dir}/${subdir} ] ; then
            crmon -o --resub-any --max-threads=1 ${dir}/${subdir}
         fi
      done
   done

   echo -e "\ncurrent time = $(date)"
   echo -e "will resume in 5 minutes..."
   echo -e "end job? (yes/no) \c"
   read -t 300 end_job

   (( count++ ))
done
