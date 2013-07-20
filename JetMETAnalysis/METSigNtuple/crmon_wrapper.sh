#!/bin/bash

date="20130626"

crabdir[0]="crab_ZmumuNtuples_data_${date}"
crabdir[1]="crab_ZmumuNtuples_mc_${date}"
crabdir[2]="crab_WenuNtuples_data_${date}"
crabdir[3]="crab_WenuNtuples_mc_${date}"
crabdir[4]="crab_WenuNtuples_data_split_20130626"
crabdir[5]="crab_TtbarNtuples_mc_20130630"
crabdir[6]="crab_WenuNtuples_mc_qcd_20130626"

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
