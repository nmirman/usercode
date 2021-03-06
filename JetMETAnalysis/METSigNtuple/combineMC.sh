#!/bin/bash

crabout="/eos/uscms/store/user/nmirman/CRABjobs/crab_ZmumuNtuples_mc_20130427"

for dir in `ls ${crabout}`
do

   echo "Opening ${crabout}/${dir}"

   filelist=""
   numfiles=`ls -l ${crabout}/${dir} | wc -l`

   for i in `seq $numfiles`
   do
      files=(`ls ${crabout}/${dir}/ntuple_${i}_*`)
      if [ "${#files[@]}" > "0" ]
      then
         filelist=${filelist}"${files[0]} "
      fi
   done

   hadd /uscmst1b_scratch/lpc1/3DayLifetime/nmirman/${dir}.root ${filelist}
   mv /uscmst1b_scratch/lpc1/3DayLifetime/nmirman/${dir}.root /eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130427/

done
