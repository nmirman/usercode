#! /bin/bash

xrd='root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/Ntuples/TopMassSkim/20141030'
skimdir='/mnt/xrootd/user/nmirman/Ntuples/TopMassSkim/20141030'

export WORKING_DIR=temp
cwd=$(pwd)

count=1
for sample in `ls $skimdir`
do

   #if [[ $count != 5 ]]
   #then
   #   let count=count+1
   #   continue
   #fi

   echo 'Submitting '$sample'...'
   filename=filelist_$sample.txt
   if [ -a ${WORKING_DIR}/$filename ] ; then
      rm ${WORKING_DIR}/$filename
   fi
   
   for file in `ls $skimdir/$sample`
   do
      echo $xrd/$sample/$file >> ${WORKING_DIR}/$filename
   done

   # run options
   optMC='True'
   optTtbar='False'
   optGT='START53_V7A'
   if [[ $sample == *TTJets* ]]
   then
      optTtbar='True'
   fi
   if [[ $sample == *Run2012* ]]
   then
      optMC='False'
      optGT='FT_R_53_V18'
      if [[ $sample == *Run2012D* ]]
      then
         optGT='FT_R_53_V21'
      fi
   fi

   for i in {0..9}
   do
      numevts=10000
      skip=$(($i*$numevts))
      qsub -d . -e logs/err_${sample}_$i.txt -o logs/out_${sample}_$i.txt \
      <<< "hostname; cd ${cwd}/../..; eval `scramv1 runtime -sh` cd ${cwd}; cmsRun makentuple_cfg.py runOnMC=$optMC runTtbar=$optTtbar globalTag=$optGT fileList=${WORKING_DIR}/$filename outputFile=${WORKING_DIR}/ntuple_${sample}_$i.root randSeed=$count maxEvents=$numevts skipEvents=$skip"
   done

   let count=count+1
done

