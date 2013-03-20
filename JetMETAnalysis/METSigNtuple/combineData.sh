#!/bin/bash

crabout="/eos/uscms/store/user/nmirman/CRABjobs/crab_ZmumuNtuples_data_20130315"

filelist=""
for dir in `ls ${crabout}`
do

   echo "Opening ${crabout}/${dir}"

   numfiles=`ls -l ${crabout}/${dir} | wc -l`

   for i in `seq $numfiles`
   do
      files=`ls ${crabout}/${dir}/ntuple_${i}_*`
      if [ "${#files[@]}" > "0" ]
      then
         filelist=${filelist}"${files[0]} "
      fi
   done

done

hadd ~/nobackup/Data.root ${filelist}
mv ~/nobackup/Data.root /eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130315/

