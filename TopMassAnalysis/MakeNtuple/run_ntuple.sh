#! /bin/bash

#indate='20160607'
indate='20160523'
outdate='20160714'
#xrd='root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/Ntuples/TopMassSkim/'$indate
xrd='root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/'
#skimdir='/mnt/xrootd/user/nmirman/Ntuples/TopMassSkim/'$indate
skimdir='/mnt/xrootd/user/nmirman/'
outdir='/mnt/xrootd/user/nmirman/Ntuples/TopMassNtuples/'$outdate
outxrd='root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/Ntuples/TopMassNtuples/'$outdate

export WORKING_DIR=temp
cwd=$(pwd)

count=1
for sample in `ls $skimdir`
do

   if [[ $sample == "Ntuples" ]] || [[ $sample == "corrfiles" ]]
   then
      continue
   fi

#   if [[ $sample != "TTJets_MSDecays_mass166_5_TuneZ2star_8TeV-madgraph-tauola" ]] && [[ $sample != "TTJets_MSDecays_mass175_5_TuneZ2star_8TeV-madgraph-tauola" ]]
   #if [[ $sample != "DoubleElectron" ]] && [[ $sample != "DoubleMu" ]] && [[ $sample != "DoubleMuParked" ]] && [[ $sample != "MuEG" ]]
   #if [[ $sample != "DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball" ]]
      #&& [[ $sample != "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball" ]] && [[ $sample != "WW_TuneZ2star_8TeV_pythia6_tauola" ]] && [[ $sample != "WZ_TuneZ2star_8TeV_pythia6_tauola" ]] && [[ $sample != "ZZ_TuneZ2star_8TeV_pythia6_tauola" ]]
   #if [[ $sample != "TTJets_MSDecays_mass169_5_TuneZ2star_8TeV-madgraph-tauola" ]]
   #if [[ $sample != "DoubleElectron" ]] 
   #then
   #   let count=count+1
   #   continue
   #fi

   #if [[ $count != 5 ]]
   #then
   #   let count=count+1
   #   continue
   #fi
   #if [[ $sample != "TTJets_FullLeptMGDecays_TuneP11" ]]
   #if [[ $sample != "TTJets_MSDecays_mass171_5_TuneZ2star_8TeV-madgraph-tauola" ]]
   if [[ $sample != "TT_CT10_TuneZ2star_8TeV-powheg-tauola" ]]
   then
      let count=count+1
      continue
   fi

   echo 'Submitting '$sample'...'
   filename=filelist_$sample.txt
   if [ -a ${WORKING_DIR}/$filename ] ; then
      rm ${WORKING_DIR}/$filename
   fi
   
   datetmp=0
   echo 'Dataset dates:'
   for subdir in `ls $skimdir/$sample`
   do
      # get latest date
      i=$((${#subdir}-8))
      date=${subdir:$i:8}
      echo ' * '$date
      if [[ $date -gt $datetmp ]] ; then
         datetmp=$date
      fi
   done
   echo 'Use date '$datetmp'...'

   for subdir in `ls $skimdir/$sample`
   do
      #if [[ $subdir == *$indate ]] ; then
      if [[ $subdir == *$datetmp ]] ; then
         for subsubdir in `ls $skimdir/$sample/$subdir`
         do
            for file in `ls $skimdir/$sample/$subdir/$subsubdir/`
            do
               if [[ $file == "topmassSkim"* ]]
               then
                  echo $xrd/$sample/$subdir/$subsubdir/$file >> ${WORKING_DIR}/$filename
               fi
            done
         done
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
   #if [[ $sample == *Run2012* ]]
   if [[ $sample == "DoubleElectron" ]] || [[ $sample == "DoubleMu" ]] || [[ $sample == "DoubleMuParked" ]] || [[ $sample == "MuEG" ]]
   then
      optMC='False'
      optGT='FT_53_V21_AN6'
      #if [[ $sample == *Run2012D* ]]
      #then
      #   optGT='FT_53_V21_AN6'
      #fi
   fi

   njobs=10
   if [[ $sample == "TBarToDilepton_tW-channel-DR_mass178_5_8TeV-powheg-tauola" ]]
   then
      njobs=20
   fi
   if [[ $sample == *"msdecay"* ]] || [[ $sample == *"FullLeptMGDecays"* ]] || [[ $sample == *"MSDecays"* ]] || [[ $sample == "TT_CT10_TuneZ2star_8TeV-powheg-tauola" ]]
   then
      njobs=100
      echo 'msdecays'
   fi
   i=0
   while [ $i -lt $njobs ]
   #for i in {0..10}
   do
      numevts=20000
      skip=$(($i*$numevts))
      qsub -d . -e logs/err_${sample}_$i.txt -o logs/out_${sample}_$i.txt \
      <<< "hostname; cd ${cwd}/../..; eval `scramv1 runtime -sh` cd ${cwd}; cmsRun makentuple_cfg.py runOnMC=$optMC runTtbar=$optTtbar globalTag=$optGT fileList=${WORKING_DIR}/$filename outputFile=${WORKING_DIR}/ntuple_${sample}_$i.root randSeed=$count maxEvents=$numevts skipEvents=$skip; xrdcp ${WORKING_DIR}/ntuple_${sample}_${i}_numEvent$numevts.root ${outxrd}/${sample}/ntuple_${sample}_${i}.root; rm ${WORKING_DIR}/ntuple_${sample}_${i}_numEvent$numevts.root"
      #cmsRun makentuple_cfg.py runOnMC=$optMC runTtbar=$optTtbar globalTag=$optGT fileList=${WORKING_DIR}/$filename outputFile=${WORKING_DIR}/ntuple_${sample}_$i.root randSeed=$count maxEvents=$numevts skipEvents=$skip
      #break
      let i=i+1
   done

   let count=count+1
done

