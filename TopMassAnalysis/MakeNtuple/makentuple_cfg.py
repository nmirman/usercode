import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

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

options.register( 'fileList',
      'filelist.txt',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      'File list.')

options.register( 'runTtbar',
      False,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      'True for Ttbar MC')

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

options.setDefault( 'inputFiles',
      'root://cmseos:1094//eos/uscms/store/user/nmirman/CRABjobs/Skim20140322/TTJets/nmirman/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/TopMass_TTJets_20140322/6bc4e5e44b306d7fc6bf3112f474f486/topmassSkim_1_1_QfX.root' )

options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

fileList = FileUtils.loadListFromFile( options.fileList )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
       #'file:BsmMassesSkim_Summer11_Sync_15_1_ZHz.root'
       #'/store/user/nmirman/DoubleElectron/TopMass_DoubleElectron_Run2012A-22Jan2013/e31ff7515f7d356ae1a23ccf4a2e42cf/BsmMassesSkim_Summer11_Sync_217_1_LOI.root',
       #'/store/user/nmirman/DoubleElectron/TopMass_DoubleElectron_Run2012A-22Jan2013/e31ff7515f7d356ae1a23ccf4a2e42cf/BsmMassesSkim_Summer11_Sync_101_1_XTb.root',
       #'/store/user/nmirman/DoubleElectron/TopMass_DoubleElectron_Run2012A-22Jan2013/e31ff7515f7d356ae1a23ccf4a2e42cf/BsmMassesSkim_Summer11_Sync_103_1_ICO.root'
    #)
    fileNames = cms.untracked.vstring( *fileList )
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = ( options.globalTag+'::All' )

process.load("TopMassAnalysis.MakeNtuple.makentuple_cfi")
process.makentuple.outFileName = options.outputFile
process.makentuple.runOnMC = options.runOnMC
process.makentuple.runTtbar = cms.bool( options.runTtbar )
process.makentuple.negTagCut = 0.244
#process.makentuple.jetScale = cms.double(0.85)

process.p = cms.Path(process.makentuple)
