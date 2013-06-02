#!/bin/bash

dir_data="crab_WenuNtuples_data_20130504"
dir_mc="crab_WenuNtuples_mc_20130504"

dir_top="/eos/uscms/store/user/nmirman/CRABjobs/"

path=${dir_top}${dir_data}
mkdir ${path}
mkdir ${path}/Run2012A-13Jul2012
mkdir ${path}/Run2012A-recover-06Aug2012
mkdir ${path}/Run2012B-13Jul2012
mkdir ${path}/Run2012C-24Aug2012
mkdir ${path}/Run2012C-PromptReco-v2
mkdir ${path}/Run2012C-EcalRecover-11Dec2012
mkdir ${path}/Run2012D-PromptReco-v1
#mkdir ${path}/Run2012D-16Jan2013

path=${dir_top}${dir_mc}
mkdir ${path}
mkdir ${path}/WJetsToLNu
mkdir ${path}/DYJetsToLL_M-50
mkdir ${path}/DYJetsToLL_M-10To50
mkdir ${path}/WW
mkdir ${path}/WZ
mkdir ${path}/ZZ
mkdir ${path}/TTJets
mkdir ${path}/Tbar_tW-channel
mkdir ${path}/T_tW-channel
mkdir ${path}/QCD_EMEnriched_20_30
mkdir ${path}/QCD_EMEnriched_30_80
mkdir ${path}/QCD_EMEnriched_80_170
mkdir ${path}/QCD_BCtoE_20_30
mkdir ${path}/QCD_BCtoE_30_80
mkdir ${path}/QCD_BCtoE_80_170
mkdir ${path}/Gamma_0_15
mkdir ${path}/Gamma_15_30
mkdir ${path}/Gamma_30_50
mkdir ${path}/Gamma_50_80
mkdir ${path}/Gamma_80_120
mkdir ${path}/Gamma_120_170
