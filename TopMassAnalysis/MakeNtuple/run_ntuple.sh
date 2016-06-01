#! /bin/bash

indate='20150510'
outdate='20160520'
xrd='root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/Ntuples/TopMassSkim/'$indate
skimdir='/mnt/xrootd/user/nmirman/Ntuples/TopMassSkim/'$indate
outdir='/mnt/xrootd/user/nmirman/Ntuples/TopMassNtuples/'$outdate
outxrd='root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/Ntuples/TopMassNtuples/'$outdate

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
   #if [[ $sample != "TTJets_FullLeptMGDecays_TuneP11" ]]
   #if [[ $sample != "TTJets_msdecay" ]]
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
      if [[ $file == "topmassSkim"* ]]
      then
         echo $xrd/$sample/$file >> ${WORKING_DIR}/$filename
      fi
   done

   # make output directories
   if [ ! -d "$outdir" ]; then
      mkdir $outdir
   fi
   if [ ! -d "$outdir/$sample" ]; then
      mkdir $outdir/$sample
   fi

   # run options
   optMC='True'
   optTtbar='False'
   optGT='START53_V27'
   if [[ $sample == *TTJets* ]]
   then
      optTtbar='True'
   fi
   if [[ $sample == *Run2012* ]]
   then
      optMC='False'
      optGT='FT_53_V21_AN6'
      #if [[ $sample == *Run2012D* ]]
      #then
      #   optGT='FT_53_V21_AN6'
      #fi
   fi

   njobs=10
   if [[ $sample == *"msdecay"* ]] || [[ $sample == *"FullLeptMGDecays"* ]]
   then
      njobs=100
      echo 'msdecays'
   fi
   i=0
   while [ $i -lt $njobs ]
   #for i in {0..10}
   do
      numevts=50000
      skip=$(($i*$numevts))
      qsub -d . -e logs/err_${sample}_$i.txt -o logs/out_${sample}_$i.txt \
      <<< "hostname; cd ${cwd}/../..; eval `scramv1 runtime -sh` cd ${cwd}; cmsRun makentuple_cfg.py runOnMC=$optMC runTtbar=$optTtbar globalTag=$optGT fileList=${WORKING_DIR}/$filename outputFile=${WORKING_DIR}/ntuple_${sample}_$i.root randSeed=$count maxEvents=$numevts skipEvents=$skip; xrdcp ${WORKING_DIR}/ntuple_${sample}_${i}_numEvent$numevts.root ${outxrd}/${sample}/ntuple_${sample}_${i}.root; rm ${WORKING_DIR}/ntuple_${sample}_${i}_numEvent$numevts.root"
      #cmsRun makentuple_cfg.py runOnMC=$optMC runTtbar=$optTtbar globalTag=$optGT fileList=${WORKING_DIR}/$filename outputFile=${WORKING_DIR}/ntuple_${sample}_$i.root randSeed=$count maxEvents=$numevts skipEvents=$skip
      let i=i+1
      #break
   done

   let count=count+1
done

