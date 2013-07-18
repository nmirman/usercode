#!/bin/bash

crabout="/eos/uscms/store/user/xguo93/CRABjobs/crab_TopMassNtuple_mc_20130625"
ntuple="/eos/uscms/store/user/xguo93/Ntuples/20130625/"

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

   hadd ~/nobackup/${dir}.root ${filelist}
   mv ~/nobackup/${dir}.root ${ntuple}

done
