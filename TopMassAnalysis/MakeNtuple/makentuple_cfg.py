import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_9_0_FyU.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_99_1_S17.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_98_1_BX0.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_97_1_ufP.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_96_1_POI.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_95_1_lC9.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_94_1_72q.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_93_1_bFP.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_92_1_8ZC.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_91_1_Lu1.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_126_0_rvv.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_121_0_mTW.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_12_0_26g.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_13_0_nLt.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_16_0_mNf.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_17_0_Mks.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_125_0_Hds.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_127_0_lSE.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_1_1_MkE.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_20_0_iNg.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_14_0_iyR.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_15_0_sBX.root',
       #'/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_18_0_Gxi.root',
       '/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_24_0_5To.root',
       '/store/user/neggert/TTJets_TuneZ2_7TeV-madgraph-tauola/TopDil_Skim_Jan2012_Signal/99d3aece7d55e016e7d52ecacf0253ac/BsmMassesSkimJan2012_25_0_l45.root'
    )
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START42_V13::All'

process.load("BsmMasses.MakeNtuple.makentuple_cfi")
process.makentuple.outFileName = 'ntuple.root'
process.makentuple.runOnMC = True
process.makentuple.negTagCut = 0.244

process.p = cms.Path(process.makentuple)
