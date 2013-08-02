#!/bin/bash

date="20130723"

crabdir[0]="crab_Zmumu_data_20130726"
crabdir[1]="crab_Zmumu_mc_${date}"
crabdir[2]="crab_Wenu_data_${date}"
crabdir[3]="crab_Wenu_mc_${date}"
crabdir[4]="crab_Wenu_loose_data_${date}"
crabdir[5]="crab_Wenu_loose_mc_${date}"
crabdir[6]="crab_Dijet_data_20130726"
crabdir[7]="crab_Dijet_mc_${date}"
crabdir[8]="crab_Ttbar0lept_data_${date}"
crabdir[9]="crab_Ttbar0lept_mc_${date}"
crabdir[10]="crab_Ttbar1lept_data_${date}"
crabdir[11]="crab_Ttbar1lept_mc_${date}"
crabdir[12]="crab_Zmumu_data_20130728"
crabdir[13]="crab_Dijet_data_20130728"
crabdir[14]="crab_Zmumu_mc_20130728"

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

   #echo -e "\ncurrent time = $(date)"
   #echo -e "will resume in 5 minutes..."
   #echo -e "end job? (yes/no) \c"
   #read -t 300 end_job

   (( count++ ))
done
