#! /bin/bash

export WORKING_DIR=temp

for sample in `ls $WORKING_DIR | grep txt`
do

   sample=${sample#"filelist_"}
   sample=${sample%".txt"}

   fileout=$WORKING_DIR/ntuple_$sample.root
   if [ -e $fileout ]; then
      rm $fileout
   fi
   hadd $fileout $WORKING_DIR/ntuple_${sample}_?_numEvent10000.root $WORKING_DIR/ntuple_${sample}_??_numEvent10000.root

done
