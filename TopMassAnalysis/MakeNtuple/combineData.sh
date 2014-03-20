#!/bin/bash

crabout="/eos/uscms/store/user/nmirman/CRABjobs/crab_TopMassNtuple_data_20140223/"
ntuple="/eos/uscms/store/user/nmirman/Ntuples/TopMass/20140223/"

filelist=""
for dir in `ls ${crabout}`
do

   echo "Opening ${crabout}/${dir}"

   numfiles=`ls -l ${crabout}/${dir} | grep root | wc -l`

   for i in `seq $numfiles`
   do
      files=(`ls ${crabout}/${dir}/ntuple_${i}_*`)
      if [ "${#files[@]}" > "0" ]
      then
         filelist=${filelist}"${files[0]} "
      fi
   done

done

hadd ~/nobackup/Data.root ${filelist}
mkdir ${ntuple}
mv ~/nobackup/Data.root ${ntuple}

