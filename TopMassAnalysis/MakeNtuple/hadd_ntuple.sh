#! /bin/bash

export WORKING_DIR=temp
export OUT_DIR=root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/Ntuples/TopMassNtuples/20150526/

cwd=$(pwd)

for sample in `ls $WORKING_DIR | grep txt`
do

   sample=${sample#"filelist_"}
   sample=${sample%".txt"}

   #if [[ $sample == *"TTJets"* ]]
   #then
   #   continue
   #fi
   #if [[ $sample != *"MGDecays"* ]]
   #then
   #   continue
   #fi

   echo 'Submitting '$sample
   fileout=ntuple_$sample.root
   #if [ -e $fileout ]; then
   #   rm $fileout
   #fi
   #hadd $fileout $WORKING_DIR/ntuple_${sample}_?_numEvent10000.root $WORKING_DIR/ntuple_${sample}_??_numEvent10000.root

   filelist=""
   for file in `ls $WORKING_DIR/ntuple_${sample}_[0-9]*_numEvent50000.root`
   do
      filelist=$filelist" "$file
   done
   #qsub -d . -e logs_hadd/err_${sample}.txt -o logs_hadd/out_${sample}.txt \
   #<<< "hostname; cd ${cwd}/../..; eval `scramv1 runtime -sh` cd ${cwd}; hadd temp_hadd/$fileout $filelist; xrdcp temp_hadd/$fileout $OUT_DIR; rm temp_hadd/$fileout"
   hadd temp_hadd/$fileout $filelist
   #xrdcp temp_hadd/$fileout $OUT_DIR
   #rm temp_hadd/$fileout
   echo 'deleting files'
   for file in `ls $WORKING_DIR/ntuple_${sample}_[0-9]*_numEvent50000.root`
   do
      rm $file
   done
done
