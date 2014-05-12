#! /bin/bash

xrd='root://cmseos:1094/'
skimdir='/eos/uscms/store/user/nmirman/CRABjobs/Skim20140322'
outdir='/eos/uscms/store/user/nmirman/Ntuples/TopMass/20140507'

if [[ $1 == 1 ]]
then
   cd /uscms/home/nmirman/nobackup/CMSSW_5_3_9_patch3/src/
   eval `scramv1 runtime -sh`
   cd TopMassAnalysis/MakeNtuple/
   export WORKING_DIR=$_CONDOR_SCRATCH_DIR
fi

count=-1
for sample in `ls $skimdir`
do

   count=$[count+1]

   if [[ $2 != $count && $1 == 1 ]]
   then
      continue
   fi

   echo 'Processing '$sample'...'
   filename=filelist_$sample.txt
   if [ -a ${WORKING_DIR}$filename ] ; then
      rm ${WORKING_DIR}$filename
   fi
   
   for file in `tree $skimdir/$sample -i -f | grep .root`
   do
      echo $xrd$file >> ${WORKING_DIR}$filename
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

   cmsRun makentuple_cfg.py runOnMC=$optMC runTtbar=$optTtbar globalTag=$optGT fileList=${WORKING_DIR}$filename outputFile=${WORKING_DIR}ntuple.root
   mv ${WORKING_DIR}ntuple.root $outdir/ntuple_$sample.root
   rm ${WORKING_DIR}$filename

done

