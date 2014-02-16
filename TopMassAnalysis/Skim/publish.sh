#!/bin/bash

date="20130708"

#crabdir[0]="crab_TopMassSkim_data_${date}"
crabdir[0]="crab_TopMassSkim_mc_${date}"

maxprocs=4
for dir in ${crabdir[*]}
do
   for subdir in `ls ${dir}`
   do

      # limit number of simultaneous processes
      while [ $(jobs -p | wc -l) -ge $max_procs ]; do
         sleep 1
      done

      (
      if [ -d ${dir}/${subdir} ] ; then
         echo "Publishing ${dir}/${subdir}:"
         crab -publish -USER.dbs_url_for_publication=https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_01_writer/servlet/DBSServlet -c ${dir}/${subdir}
      fi
      ) &

   done
done

