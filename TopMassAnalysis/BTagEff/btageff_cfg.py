import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.setDefault( 'maxEvents', -1)

options.setDefault( 'outputFile',
      'btageff.root' )

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

options.register( 'randSeed',
      0,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.int,
      'Random Seed for MC smearing')

options.register( 'skipEvents',
      0,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.int,
      'Skip first n events.')

options.register( 'wantSummary',
      False,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      "Print summary at end of job")

options.setDefault( 'inputFiles',
      'root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/nmirman/Ntuples/TopMassNtuples/20160601/DoubleElectron/ntuple_DoubleElectron_0.root' 
      )

options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

fileList = FileUtils.loadListFromFile( options.fileList )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( *fileList ),
    skipEvents = cms.untracked.uint32( options.skipEvents )
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = ( options.globalTag+'::All' )

process.load("TopMassAnalysis.BTagEff.btageff_cfi")
process.demo.outFileName = options.outputFile
process.demo.runOnMC = options.runOnMC
process.demo.runTtbar = cms.bool( options.runTtbar )
process.demo.negTagCut = 0.244
process.demo.randSeed = options.randSeed

process.TFileService = cms.Service('TFileService',
         fileName = cms.string(options.outputFile)
         )

process.p = cms.Path(process.demo)

