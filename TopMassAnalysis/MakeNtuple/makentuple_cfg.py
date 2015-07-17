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
      'root://cmseos:1094//eos/uscms/store/user/nmirman/CRABjobs/Skim20140322/TTJets/nmirman/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/TopMass_TTJets_20140322/6bc4e5e44b306d7fc6bf3112f474f486/topmassSkim_1_1_QfX.root' )

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

process.load("TopMassAnalysis.MakeNtuple.makentuple_cfi")
process.makentuple.outFileName = options.outputFile
process.makentuple.runOnMC = options.runOnMC
process.makentuple.runTtbar = cms.bool( options.runTtbar )
process.makentuple.negTagCut = 0.244
process.makentuple.randSeed = options.randSeed

# Produce PDF weights (maximum is 3)
#process.pdfWeights = cms.EDProducer("PdfWeightProducer",
#      # Fix POWHEG if buggy (this PDF set will also appear on output,
#      # so only two more PDF sets can be added in PdfSetNames if not "")
#      #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
#      #GenTag = cms.untracked.InputTag("genParticles"),
#      PdfInfoTag = cms.untracked.InputTag("generator"),
#      PdfSetNames = cms.untracked.vstring(
#         "cteq66.LHgrid"
#         , "MRST2006nnlo.LHgrid"
#         , "NNPDF10_100.LHgrid"
#         )
#      )

process.load('TopAnalysis.TopUtils.EventWeightBJES_cfi')

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("histo.root"),
      #closeFileFast = cms.untracked.bool(True)
      )

process.p = cms.Path(process.makentuple)

