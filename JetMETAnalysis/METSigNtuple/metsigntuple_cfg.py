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

options.register( 'channel',
      'Zmumu',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "Physics channel: Zmumu or Wenu."
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
         #'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root'
         #'/store/data/Run2012A/DoubleMu/AOD/13Jul2012-v1/00000/0048B245-B9D2-E111-A8DC-0018F3D096BE.root'
         #'/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v2/10000/00A4899E-666B-E211-A2AC-E0CB4E29C50D.root'
         #'/store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v2/0000/000869F4-59EE-E111-9BF6-003048D47752.root'
         #'/store/data/Run2012A/SingleElectron/AOD/13Jul2012-v1/0000/001A2EB8-47D4-E111-B527-003048679070.root'
         #'/store/mc/Summer12_DR53X/TTJets_HadronicMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00000/002A756C-FA15-E211-9FA6-485B39800B75.root'
         '/store/data/Run2012A/DoubleMu/AOD/22Jan2013-v1/20000/001AE30A-BA81-E211-BBE7-003048FFD770.root'
         )
      )

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = ( options.globalTag+'::All' )
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

process.load("CommonTools.ParticleFlow.PF2PAT_cff")
process.pfPileUp.Enable = False

process.pfAllMuons.src="particleFlow"
process.pfMuonsFromVertex.dzCut=9999
process.pfNoMuon.bottomCollection   = "particleFlow"
if options.channel == 'Zmumu':
   process.pfNoMuon.topCollection      = "pfSelectedMuons"
process.pfJets.doAreaFastjet        = True
process.pfJets.jetPtMin             = 0
process.pfJets.src                  = "pfNoElectron"

process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = cms.bool(True)

process.mypf2pat_zmumu = cms.Sequence(
      process.pfNoPileUpSequence * # pfPileUp enable is false
      process.pfParticleSelectionSequence *
      process.pfAllMuons * 
      process.pfMuonsFromVertex *
      process.pfSelectedMuons *
      process.pfNoMuon *
      process.pfElectronSequence *
      process.pfNoElectron *
      process.pfJets
      )
process.mypf2pat_wenu = cms.Sequence(
      process.pfNoPileUpSequence *
      process.pfParticleSelectionSequence *
      process.pfMuonSequence *
      process.pfNoMuon *
      process.pfElectronSequence *
      process.pfNoElectron *
      process.pfJets
      )
#
## met corrections and filters
#

# type-I corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
if not options.runOnMC:
      process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

# type-0 corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")

process.pfType1CorrectedMetType0 = process.pfType1CorrectedMet.clone(
      applyType0Corrections = cms.bool(False),
      srcType1Corrections = cms.VInputTag(
         cms.InputTag('pfMETcorrType0'),
         cms.InputTag('pfJetMETcorr', 'type1')        
         )
      )
process.producePFMETCorrectionsType0 = cms.Sequence( process.producePFMETCorrections )
process.producePFMETCorrectionsType0.replace(
      process.pfType1CorrectedMet,
      process.pfType1CorrectedMetType0
      )

# x/y shift corrections
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
if not options.runOnMC:
   process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data
if options.runOnMC:
   process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc

process.pfType1CorrectedMetXYshift = process.pfType1CorrectedMet.clone(
      srcType1Corrections = cms.VInputTag(
         cms.InputTag('pfJetMETcorr', 'type1'),
         cms.InputTag('pfMEtSysShiftCorr')
         )
      )
process.producePFMETCorrectionsXYshift = cms.Sequence( process.producePFMETCorrections )
process.producePFMETCorrectionsXYshift.replace(
      process.pfType1CorrectedMet,
      process.pfType1CorrectedMetXYshift
      )

# type-0 and x/y shift corrections
process.pfType1CorrectedMetType0XYshift = process.pfType1CorrectedMet.clone(
      applyType0Corrections = cms.bool(False),
      srcType1Corrections = cms.VInputTag(
         cms.InputTag('pfMETcorrType0'),
         cms.InputTag('pfJetMETcorr', 'type1'),
         cms.InputTag('pfMEtSysShiftCorr')
         )
      )
process.producePFMETCorrectionsType0XYshift = cms.Sequence( process.producePFMETCorrections )
process.producePFMETCorrectionsType0XYshift.replace(
      process.pfType1CorrectedMet,
      process.pfType1CorrectedMetType0XYshift
      )

process.mymet = cms.Sequence(
      process.pfMET *
      process.producePFMETCorrections *
      process.type0PFMEtCorrection *
      process.producePFMETCorrectionsType0 *
      process.pfMEtSysShiftCorrSequence *
      process.producePFMETCorrectionsXYshift *
      process.producePFMETCorrectionsType0XYshift
      )

metList = []
metList.append(cms.untracked.InputTag("pfMet", "", ""))
metList.append(cms.untracked.InputTag("pfType1CorrectedMet", "", ""))
metList.append(cms.untracked.InputTag("pfType1CorrectedMetType0", "", ""))
metList.append(cms.untracked.InputTag("pfType1CorrectedMetXYshift", "", ""))
metList.append(cms.untracked.InputTag("pfType1CorrectedMetType0XYshift", "", ""))

# jet pileup id
from CMGTools.External.pujetidsequence_cff import puJetId, puJetMva

process.recoPuJetId = puJetId.clone(
      jets = cms.InputTag("pfJets"),
      applyJec = cms.bool(True),
      inputIsCorrected = cms.bool(False),                
      )

process.recoPuJetMva = puJetMva.clone(
      jets = cms.InputTag("pfJets"),
      jetids = cms.InputTag("recoPuJetId"),
      applyJec = cms.bool(True),
      inputIsCorrected = cms.bool(False),                
      )

process.recoPuJetIdSequence = cms.Sequence(process.recoPuJetId * process.recoPuJetMva )

# rho value for isolation
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# particle flow isolation
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

if options.channel == 'Zmumu':
   trigger_paths = ["HLT_Mu17_Mu8_v"]
if options.channel == 'Wenu':
   trigger_paths = ["HLT_Ele27_WP80_v"]
if options.channel == 'Ttbar':
   trigger_paths = ["HLT_SixJet45_v"]
if options.channel == 'Dijet':
   trigger_paths = ["HLT_PFJet320_v"]
trigger_pattern = [path+"*" for path in trigger_paths]

process.demo = cms.EDAnalyzer('METSigNtuple',
      runOnMC              = cms.untracked.bool(options.runOnMC),
      output_file          = cms.untracked.string(options.outputFile),

      saveJetInfo          = cms.untracked.bool(True),
      saveBTags          = cms.untracked.bool(False),

      selectionChannel     = cms.untracked.string('Zmumu'),

      pfjetsTag            = cms.untracked.InputTag('pfJets'),
      pfjetCorrectorL1     = cms.untracked.string('ak5PFL1Fastjet'),
      pfjetCorrectorL123   = cms.untracked.string('ak5PFL1FastL2L3'),
      jetResAlgo           = cms.string('AK5PF'),
      jetResEra            = cms.string('Spring10'),

      muonTag              = cms.untracked.InputTag("pfSelectedMuons"),
      electronTag          = cms.untracked.InputTag("pfIsolatedElectrons"),

      conversionsInputTag     = cms.InputTag("allConversions"),
      rhoIsoInputTag          = cms.InputTag("kt6PFJetsForIsolation", "rho"),
      isoValInputTags         = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
         cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
         cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),

      genparticlesTag      = cms.untracked.InputTag("genParticles"),
      pfcandidatesTag      = cms.untracked.InputTag("particleFlow"),

      genjetsTag           = cms.untracked.InputTag('ak5GenJets'),

      metsTag              = cms.untracked.VInputTag(metList),
      genmetTag            = cms.untracked.InputTag('genMetTrue'),

      verticesTag          = cms.untracked.InputTag('offlinePrimaryVertices'),
      pileupTag            = cms.untracked.InputTag('addPileupInfo'),

      metSig               = cms.untracked.InputTag('pfMetSig','METSignificance'),
      metSigMatrix00       = cms.untracked.InputTag('pfMetSig','CovarianceMatrix00'),
      metSigMatrix01       = cms.untracked.InputTag('pfMetSig','CovarianceMatrix01'),
      metSigMatrix10       = cms.untracked.InputTag('pfMetSig','CovarianceMatrix10'),
      metSigMatrix11       = cms.untracked.InputTag('pfMetSig','CovarianceMatrix11')
      )
if not options.runOnMC:
   process.demo.pfjetCorrectorL123 = 'ak5PFL1FastL2L3Residual'
if options.channel == 'Zmumu':
   process.demo.saveJetInfo = cms.untracked.bool(True)
   process.demo.selectionChannel = cms.untracked.string('Zmumu')
   process.demo.muonTag = cms.untracked.InputTag("pfSelectedMuons")
if options.channel == 'Wenu':
   process.demo.saveJetInfo = cms.untracked.bool(False)
   process.demo.selectionChannel = cms.untracked.string('Wenu')
   process.demo.muonTag = cms.untracked.InputTag("pfIsolatedMuons")
if options.channel == 'Ttbar':
   process.demo.saveJetInfo = cms.untracked.bool(True)
   process.demo.selectionChannel = cms.untracked.string('Ttbar')
   process.demo.muonTag = cms.untracked.InputTag("pfIsolatedMuons")
   process.demo.saveBTags = cms.untracked.bool(True)
if options.channel == 'Dijet':
   process.demo.saveJetInfo = cms.untracked.bool(True)
   process.demo.selectionChannel = cms.untracked.string('Dijet')
   process.demo.muonTag = cms.untracked.InputTag("pfIsolatedMuons")
   process.demo.saveBTags = cms.untracked.bool(True)

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

# MET Significance producer
process.load("JetMETAnalysis.METSignificance.metsignificance_cfi")
process.pfMetSig.runOnMC = options.runOnMC
if not options.runOnMC:
      process.pfMetSig.pfjetCorrectorL123 = 'ak5PFL1FastL2L3Residual'

process.mypf2pat = cms.Sequence()
if options.channel == 'Zmumu':
   process.mypf2pat *= process.mypf2pat_zmumu
if options.channel == 'Wenu':
   process.mypf2pat *= process.mypf2pat_wenu
if options.channel == 'Ttbar':
   process.mypf2pat *= process.mypf2pat_wenu
if options.channel == 'Dijet':
   process.mypf2pat *= process.mypf2pat_wenu

process.p = cms.Path(
      process.triggerSelection *
      process.filtersSeq *
      process.mypf2pat *
      process.mymet *
      process.recoPuJetIdSequence *
      process.kt6PFJetsForIsolation *
      process.eleIsoSequence *
      process.pfiso *
      process.pfMetSig *
      process.demo
      )
if options.runOnMC :
   process.p.remove( process.filtersSeq )

#process.out = cms.OutputModule( "PoolOutputModule"
#      , fileName = cms.untracked.string( "patTuple.root" ),
#      SelectEvents = cms.untracked.PSet(
#         SelectEvents = cms.vstring('p')
#         )
#      )
#process.outpath = cms.EndPath( process.out )
