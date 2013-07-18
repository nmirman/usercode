import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.setDefault( 'maxEvents', -1)

options.setDefault( 'outputFile',
      'ntuple.root' )

options.register( 'runOnMC',
      False,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      'True for MC')

options.register( 'globalTag',
      'START53_V7A',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      'CMS Global Tag')

options.register( 'wantSummary',
      False,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      "Print summary at end of job")

options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) #-1

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       #'file:BsmMassesSkim_Summer11_Sync_15_1_ZHz.root'
       '/store/user/nmirman/DoubleElectron/TopMass_DoubleElectron_Run2012A-22Jan2013/e31ff7515f7d356ae1a23ccf4a2e42cf/BsmMassesSkim_Summer11_Sync_217_1_LOI.root',
       '/store/user/nmirman/DoubleElectron/TopMass_DoubleElectron_Run2012A-22Jan2013/e31ff7515f7d356ae1a23ccf4a2e42cf/BsmMassesSkim_Summer11_Sync_101_1_XTb.root',
       '/store/user/nmirman/DoubleElectron/TopMass_DoubleElectron_Run2012A-22Jan2013/e31ff7515f7d356ae1a23ccf4a2e42cf/BsmMassesSkim_Summer11_Sync_103_1_ICO.root'
    )
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = ( options.globalTag+'::All' )

process.load("TopMassAnalysis.MakeNtuple.makentuple_cfi")
process.makentuple.outFileName = options.outputFile
process.makentuple.runOnMC = options.runOnMC
process.makentuple.negTagCut = 0.244

process.p = cms.Path(process.makentuple)
