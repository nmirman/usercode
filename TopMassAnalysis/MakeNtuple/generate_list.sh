#! /bin/bash

dir='/eos/uscms/store/user/nmirman/CRABjobs/Skim20140322/TTJets_allevts/nmirman/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/TopMass_TTJets_allevts_20140322/3a67ca30416f008729590be490327b85/'

rm filelist.txt

for file in `ls ${dir}`
do
      echo 'file:'${dir}${file} >> filelist.txt
done
