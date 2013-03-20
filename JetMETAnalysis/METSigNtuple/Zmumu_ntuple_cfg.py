import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.setDefault( 'outputFile',
      'ntuple.root'
      )

options.register( 'globalTag',
      'START53_V7A',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "CMS Global Tag"
      )

options.register( 'runOnMC',
      False,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      "True for MC"
      )

options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(
         '/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root'
         #'/store/data/Run2012A/DoubleMu/AOD/13Jul2012-v1/00000/0048B245-B9D2-E111-A8DC-0018F3D096BE.root'
         #'/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v2/10000/00A4899E-666B-E211-A2AC-E0CB4E29C50D.root'
         )
      )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = ( options.globalTag+'::All' )
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

process.load("CommonTools.ParticleFlow.PF2PAT_cff")

process.pfAllMuons.src="particleFlow"
process.pfMuonsFromVertex.dzCut=9999
process.pfSelectedMuons.cut = cms.string('pt>20 && abs(eta)<2.1')
process.pfNoMuon.bottomCollection   = "particleFlow"
process.pfNoMuon.topCollection      = "pfSelectedMuons"
process.pfJets.doAreaFastjet        = True
process.pfJets.jetPtMin             = 0
process.pfJets.src                  = "pfNoMuon"

process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = cms.bool(True)

process.mypf2pat = cms.Sequence(
      process.pfAllMuons * 
      process.pfMuonsFromVertex *
      process.pfSelectedMuons *
      process.pfNoMuon *
      process.pfJets
      )

# met corrections and filters
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

if options.runOnMC :
   process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc
else :
   process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data

process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
      cms.InputTag('pfJetMETcorr', 'type1') ,
      cms.InputTag('pfMEtSysShiftCorr')  
      )
process.pfType1p2CorrectedMet.srcType1Corrections = cms.VInputTag(
      cms.InputTag('pfJetMETcorr', 'type1') ,
      cms.InputTag('pfMEtSysShiftCorr')       
      )

process.mymet = cms.Sequence(
      process.pfMET *
      process.pfMEtSysShiftCorrSequence *
      process.producePFMETCorrections
      )

metList = []
metList.append(cms.untracked.InputTag("pfMet", "", ""))

trigger_paths = ["HLT_Mu17_Mu8_v"]
trigger_pattern = [path+"*" for path in trigger_paths]

process.demo = cms.EDAnalyzer('METSigNtuple',
      runOnMC              = cms.untracked.bool(options.runOnMC),
      output_file          = cms.untracked.string(options.outputFile),

      pfjetsTag            = cms.untracked.InputTag('pfJets'),
      pfjetCorrectorL1     = cms.untracked.string('ak5PFL1Fastjet'),
      pfjetCorrectorL123   = cms.untracked.string('ak5PFL1FastL2L3'),
      jetResAlgo           = cms.string('AK5PF'),
      jetResEra            = cms.string('Spring10'),

      muonTag              = cms.untracked.InputTag("muons"),

      genparticlesTag      = cms.untracked.InputTag("genParticles"),
      pfcandidatesTag      = cms.untracked.InputTag("particleFlow"),

      genjetsTag           = cms.untracked.InputTag('ak5GenJets'),

      metsTag              = cms.untracked.VInputTag(metList),
      genmetTag            = cms.untracked.InputTag('genMetCalo'),

      verticesTag          = cms.untracked.InputTag('offlinePrimaryVertices'),
      pileupTag            = cms.untracked.InputTag('addPileupInfo')
      )
if not options.runOnMC:
   process.demo.pfjetCorrectorL123 = 'ak5PFL1FastL2L3Residual'

# MET filters for ICHEP 2012
# obtained from
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#A_Central_Filter_Package_RecoMET

#from PhysicsTools.PatAlgos.patTemplate_cfg import *
## The good primary vertex filter ____________________________________________||

process.primaryVertexFilter = cms.EDFilter(
      "VertexSelector",
      src = cms.InputTag("offlinePrimaryVertices"),
      cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
      filter = cms.bool(True)
      )

## The beam scraping filter __________________________________________________||
process.noscraping = cms.EDFilter(
      "FilterOutScraping",
      applyfilter = cms.untracked.bool(True),
      debugOn = cms.untracked.bool(False),
      numtrack = cms.untracked.uint32(10),
      thresh = cms.untracked.double(0.25)
      )

## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

## The CSC beam halo tight filter ____________________________________________||
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")

## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

## The ECAL laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
#process.MessageLogger.suppressError = cms.untracked.vstring ('ecalLaserCorrFilter')

## The Good vertices collection needed by the tracking failure filter ________||
process.goodVertices = cms.EDFilter(
      "VertexSelector",
      filter = cms.bool(False),
      src = cms.InputTag("offlinePrimaryVertices"),
      cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
      )

## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

## The tracking POG filters __________________________________________________||
process.load('RecoMET.METFilters.trackingPOGFilters_cff')

process.filtersSeq = cms.Sequence(
      process.primaryVertexFilter *
      process.noscraping *
      process.HBHENoiseFilter *
      process.CSCTightHaloFilter *
      process.hcalLaserEventFilter *
      process.EcalDeadCellTriggerPrimitiveFilter *
      process.goodVertices * process.trackingFailureFilter *
      process.eeBadScFilter *
      process.ecalLaserCorrFilter *
      process.trkPOGFilters
      )

# trigger filter                
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::HLT",
      HLTPaths = trigger_pattern,
      throw=False
      )

process.p = cms.Path(
      process.filtersSeq *
      process.triggerSelection *
      process.mypf2pat *
      process.mymet *
      process.demo
      )
if options.runOnMC :
   process.p.remove( process.filtersSeq )
