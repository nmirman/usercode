#!/bin/bash

dir_data="crab_ZmumuNtuples_data_20130501"
dir_mc="crab_ZmumuNtuples_mc_20130501"

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
mkdir ${path}/Run2012D-16Jan2013

path=${dir_top}${dir_mc}
mkdir ${path}
mkdir ${path}/DYJetsToLL
mkdir ${path}/TTJets
mkdir ${path}/QCD
mkdir ${path}/WW
mkdir ${path}/WZ
mkdir ${path}/ZZ
mkdir ${path}/Tbar_tW-channel
mkdir ${path}/T_tW-channel
mkdir ${path}/WJetsToLNu
