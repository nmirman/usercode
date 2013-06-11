import sys

import FWCore.ParameterSet.Config as cms

# parse command line arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')

options.setDefault( 'maxEvents', -1)

options.setDefault( 'outputFile',
                  "BsmMassesSkim_Summer11_Sync.root" )
                  
options.setDefault( 'inputFiles',
                    '/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v2/10000/00A4899E-666B-E211-A2AC-E0CB4E29C50D.root')

options.register( 'inputType',
                  'signalMC',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'Type of input file. Options are data, signalMC, tauMC, topNonDilMC, backgroundMC.')

options.register( 'wantSummary',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Print summary at end of job")

options.parseArguments()

process = cms.Process('PAT')

if (options.inputType == 'signalMC') :
    runOnMC = True
    signalMC = True
    TauMC = False
    TopNonDilMC = False
elif (options.inputType == 'tauMC') :
    runOnMC = True
    signalMC = False
    TauMC = True
    TopNonDilMC = False
elif (options.inputType == 'topNonDilMC') :
    runOnMC = True
    signalMC = False
    TauMC = False
    TopNonDilMC = True
elif (options.inputType == 'backgroundMC') :
    runOnMC = True
    signalMC = False
    TauMC = False
    TopNonDilMC = False
elif (options.inputType == 'data') :
    runOnMC = False
    signalMC = False
    TauMC = False
    TopNonDilMC = False
else :
    sys.exit("inputType value was invalid. Options are data, signalMC, tauMC, topNonDilMC, backgroundMC.")

postfix = 'PFlow'

inputfiles = options.inputFiles

globalTagData = 'FT_P_V43E_AN3'
globalTagMC = 'START53_V7A'

fwkReportEvery = 1000

outCommands = ['drop *']

muonCutPF = "isGlobalMuon & isTrackerMuon & pt>20 & abs(eta)<2.4 & numberOfValidHits>10. & muonID('GlobalMuonPromptTight')"
electronCutPF = "pt>20 & abs(eta)<2.5 & test_bit( electronID('eidTightMC'), 0)"
jetCutPF = 'pt>30 & abs(eta)<2.5'
jetCutPF += ' && numberOfDaughters > 1'                                  # PF jet ID:
jetCutPF += ' && chargedEmEnergyFraction < 0.99'                         # PF jet ID:
jetCutPF += ' && neutralHadronEnergyFraction < 0.99'                     # PF jet ID:
jetCutPF += ' && neutralEmEnergyFraction < 0.99'                         # PF jet ID:
jetCutPF += ' && (chargedHadronEnergyFraction > 0. || abs(eta) >= 2.4)'  # PF jet ID:
jetCutPF += ' && (chargedMultiplicity > 0 || abs(eta) >= 2.4)'           # PF jet ID:


####################################################
######### end user configurables ###################
####################################################


# basics
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool( options.wantSummary )
)
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if runOnMC:
  process.GlobalTag.globaltag = globalTagMC   + '::All'
else:
  process.GlobalTag.globaltag = globalTagData + '::All'


# input stuff
process.source = cms.Source( "PoolSource"
, noEventSort        = cms.untracked.bool( True )
, duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
, fileNames          = cms.untracked.vstring( inputfiles )
)
# maximum number of events
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32( options.maxEvents )
)


# output stuff
process.out = cms.OutputModule( "PoolOutputModule"
, fileName       = cms.untracked.string( options.outputFile )
, SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring( 'p' ) )
, outputCommands = cms.untracked.vstring( outCommands)
, dropMetaData   = cms.untracked.string( 'ALL' )
)
process.outpath = cms.EndPath( process.out )

# event selection is configured later
process.out.SelectEvents.SelectEvents = []


# pick decay channel at gen level
if signalMC:
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")
    process.ttDilFilter = process.ttDecaySelection.clone()
    process.ttDilFilter.allowedTopDecays.decayBranchA.electron = True
    process.ttDilFilter.allowedTopDecays.decayBranchA.muon     = True
    process.ttDilFilter.allowedTopDecays.decayBranchA.tau      = False

    process.ttDilFilter.allowedTopDecays.decayBranchB.electron = True
    process.ttDilFilter.allowedTopDecays.decayBranchB.muon     = True
    process.ttDilFilter.allowedTopDecays.decayBranchB.tau     = False
    process.ttDilFilter.restrictTauDecays.leptonic = cms.bool(True)
elif TauMC:
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")
    process.ttDilFilter = process.ttDecaySelection.clone()
    process.ttDilFilter.allowedTopDecays.decayBranchA.electron = False
    process.ttDilFilter.allowedTopDecays.decayBranchA.muon     = False
    process.ttDilFilter.allowedTopDecays.decayBranchA.tau      = True

    process.ttDilFilter.allowedTopDecays.decayBranchB.electron = True
    process.ttDilFilter.allowedTopDecays.decayBranchB.muon     = True
    process.ttDilFilter.allowedTopDecays.decayBranchB.tau     = True
    process.ttDilFilter.restrictTauDecays.leptonic = cms.bool(True)

elif TopNonDilMC:
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")
    process.ttDilFilter = process.ttDecaySelection.clone()
    process.ttDilFilter.allowedTopDecays.decayBranchA.electron = True
    process.ttDilFilter.allowedTopDecays.decayBranchA.muon     = True
    process.ttDilFilter.allowedTopDecays.decayBranchA.tau      = True

    process.ttDilFilter.allowedTopDecays.decayBranchB.electron = True
    process.ttDilFilter.allowedTopDecays.decayBranchB.muon     = True
    process.ttDilFilter.allowedTopDecays.decayBranchB.tau     = True
    process.ttDilFilter.restrictTauDecays.leptonic = cms.bool(True)
    process.ttDilFilter.invert = True

# Gen Kinematics Filter for events with Pythia bug
process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")

# event cleaning
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *
process.HBHENoiseFilter = HBHENoiseFilter
# s. https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1196.html
HBHENoiseFilter.minIsolatedNoiseSumE        = 999999.
HBHENoiseFilter.minNumIsolatedNoiseChannels = 999999
HBHENoiseFilter.minIsolatedNoiseSumEt       = 999999.

process.scrapingFilter = cms.EDFilter(
  "FilterOutScraping"
, applyfilter = cms.untracked.bool( True )
, debugOn     = cms.untracked.bool( False )
, numtrack    = cms.untracked.uint32( 10 )
, thresh      = cms.untracked.double( 0.25 )
)

process.eventCleaning = cms.Sequence(
  process.HBHENoiseFilter
+ process.scrapingFilter
)

# primary vertex selection
pvSelection = cms.PSet(
  minNdof = cms.double( 4. )
, maxZ    = cms.double( 24. )
, maxRho  = cms.double( 2. )
)

process.goodOfflinePrimaryVertices = cms.EDFilter(
  "PrimaryVertexObjectFilter" # checks for fake PVs automatically
, filterParams = pvSelection
, filter       = cms.bool( False ) # use only as producer
, src          = cms.InputTag( 'offlinePrimaryVertices' )
)

process.goodOfflinePrimaryVertexFilter = cms.EDFilter(
  "PrimaryVertexFilter" # checks for fake PVs automatically
, pvSelection
, NPV = cms.int32( 1 )
, pvSrc = cms.InputTag( 'goodOfflinePrimaryVertices' )
)

process.vertexing = cms.Sequence(process.goodOfflinePrimaryVertices*process.goodOfflinePrimaryVertexFilter)


# PAT configuration
process.load( "PhysicsTools.PatAlgos.patSequences_cff" )
from PhysicsTools.PatAlgos.tools.coreTools import *


# JEC Levels
jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
if not runOnMC :
    jecLevels.append('L2L3Residual')


# PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT( process,
           runOnMC = runOnMC,
           jetAlgo = 'AK5',
           postfix = postfix,
           jetCorrections = ( 'AK5PFchs', jecLevels )
         )


# Top projection
applyPostfix( process, 'pfNoPileUp'  , postfix ).enable = True
applyPostfix( process, 'pfNoMuon'    , postfix ).enable = True
applyPostfix( process, 'pfNoElectron', postfix ).enable = True
applyPostfix( process, 'pfNoJet'     , postfix ).enable = True
applyPostfix( process, 'pfNoTau'     , postfix ).enable = False
applyPostfix( process, 'pfPileUp', postfix ).Vertices = cms.InputTag( 'goodOfflinePrimaryVertices' )
applyPostfix( process, 'pfPileUp', postfix ).checkClosestZVertex = False

# jet configuration
applyPostfix( process, 'pfJets', postfix ).doAreaFastjet = True
applyPostfix( process, 'pfJets', postfix ).doRhoFastjet  = False
applyPostfix( process, 'selectedPatJets', postfix ).cut  = jetCutPF



# muon configuration
applyPostfix( process, 'pfMuonsFromVertex'    , postfix ).vertices = cms.InputTag( 'goodOfflinePrimaryVertices' )
applyPostfix( process, 'pfMuonsFromVertex'    , postfix ).d0Cut    = 0.02 # this is wrt the primary vertex
#applyPostfix( process, 'isoValMuonWithCharged', postfix ).deposits[0].deltaR = 0.3
#applyPostfix( process, 'isoValMuonWithNeutral', postfix ).deposits[0].deltaR = 0.3
#applyPostfix( process, 'isoValMuonWithPhotons', postfix ).deposits[0].deltaR = 0.3
#applyPostfix( process, 'pfIsolatedMuons'      , postfix ).combinedIsolationCut = 0.2

applyPostfix( process, 'patMuons', postfix ).usePV      = True
applyPostfix( process, 'patMuons', postfix ).embedTrack = True
applyPostfix( process, 'selectedPatMuons', postfix ).cut = muonCutPF


# electron configuration
applyPostfix( process, 'pfElectronsFromVertex'    , postfix ).vertices = cms.InputTag( 'goodOfflinePrimaryVertices' )
applyPostfix( process, 'pfElectronsFromVertex'    , postfix ).d0Cut    = 0.04
#applyPostfix( process, 'isoValElectronWithCharged', postfix ).deposits[0].deltaR = 0.3
#applyPostfix( process, 'isoValElectronWithNeutral', postfix ).deposits[0].deltaR = 0.3
#applyPostfix( process, 'isoValElectronWithPhotons', postfix ).deposits[0].deltaR = 0.3
#applyPostfix( process, 'pfIsolatedElectrons'    , postfix ).combinedIsolationCut = 0.17

applyPostfix( process, 'patElectrons', postfix ).usePV      = False
applyPostfix( process, 'patElectrons', postfix ).embedTrack = True
applyPostfix( process, 'selectedPatElectrons', postfix ).cut = electronCutPF
applyPostfix( process, 'patElectrons', postfix ).electronIDSources = electronIDSources = cms.PSet( eidTight            = cms.InputTag( 'eidTight' )
                                                                                                 , eidLoose            = cms.InputTag( 'eidLoose' )
                                                                                                 , eidRobustTight      = cms.InputTag( 'eidRobustTight' )
                                                                                                 , eidRobustHighEnergy = cms.InputTag( 'eidRobustHighEnergy' )
                                                                                                 , eidRobustLoose      = cms.InputTag( 'eidRobustLoose' )
                                                                                                 , eidVeryLooseMC      = cms.InputTag( 'eidVeryLooseMC' )
                                                                                                 , eidLooseMC          = cms.InputTag( 'eidLooseMC' )
                                                                                                 , eidMediumMC         = cms.InputTag( 'eidMediumMC' )
                                                                                                 , eidTightMC          = cms.InputTag( 'eidTightMC' )
                                                                                                 , eidSuperTightMC     = cms.InputTag( 'eidSuperTightMC' )
                                                                                                 , eidHyperTight1MC    = cms.InputTag( 'eidHyperTight1MC' )
                                                                                                 )


# PAT requires 2 leptons
process.countPatLeptons.minNumber = 2
process.countPatLeptons.maxNumber = 2


# remove MC Matching
if not runOnMC :
    runOnData( process,
               names = [ 'PFAll' ],
               postfix = postfix
             )


process.out.outputCommands += [ 'keep edmTriggerResults_*_*_*'
                              , 'keep *_hltTriggerSummaryAOD_*_*'
                              # vertices and beam spot
                              , 'keep *_offlineBeamSpot_*_*'
                              , 'keep *_offlinePrimaryVertices*_*_*'
                              , 'keep *_goodOfflinePrimaryVertices*_*_*'
                              ]
if runOnMC:
    process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                  , 'keep recoGenParticles_*_*_*'
                                  , 'keep *_addPileupInfo_*_*'
                                  ]


# stuff for running L1FastJet
from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets
kt6PFJetsPFChs = kt6PFJets.clone( src = cms.InputTag( 'pfNoElectron'+postfix )
                                , doAreaFastjet = cms.bool( True )
                                , doRhoFastjet  = cms.bool( True ) )
setattr( process, 'kt6PFJetsChs' + postfix, kt6PFJetsPFChs )
getattr( process, 'patPF2PATSequence' + postfix).replace( getattr( process, 'patJetCorrFactors' + postfix )
                                                        , getattr( process, 'kt6PFJetsChs' + postfix ) * getattr( process, 'patJetCorrFactors' + postfix )
                                                        )
applyPostfix( process, 'patJetCorrFactors', postfix ).rho = cms.InputTag( 'kt6PFJetsChs' + postfix, 'rho' )
process.out.outputCommands.append( 'keep double_*' + postfix + '*_*_' + process.name_() )


# ElectronID
process.load( "RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi" )
process.eidCiCSequence = cms.Sequence(
  process.eidVeryLooseMC
+ process.eidLooseMC
+ process.eidMediumMC
+ process.eidTightMC
+ process.eidSuperTightMC
+ process.eidHyperTight1MC
)


# pair selection
process.eePair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatElectrons"+postfix+"@- selectedPatElectrons"+postfix+"@+"),
    cut = cms.string("mass > 12"),
)
process.eePairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("eePair"),
    minNumber = cms.uint32(1),
)
process.step1ee = cms.Sequence(process.eePair*process.eePairFilter)

process.emuPair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatElectrons"+postfix+"@- selectedPatMuons"+postfix+"@+"),
    cut = cms.string("mass > 12"),
)
process.emuPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("emuPair"),
    minNumber = cms.uint32(1),
)
process.step1emu = cms.Sequence(process.emuPair*process.emuPairFilter)

process.mumuPair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatMuons"+postfix+"@- selectedPatMuons"+postfix+"@+"),
    cut = cms.string("mass > 12"),
)
process.mumuPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("mumuPair"),
    minNumber = cms.uint32(1),
)
process.step1mumu = cms.Sequence(process.mumuPair*process.mumuPairFilter)

# Z mass filter
process.eePairNoZ = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("eePair"),
    cut = cms.string("!(76<mass<106)")
)
process.eePairNoZFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("eePairNoZ"),
    minNumber = cms.uint32(1)
)
process.step2ee = cms.Sequence( process.eePairNoZ*process.eePairNoZFilter)

process.mumuPairNoZ = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("mumuPair"),
    cut = cms.string("!(76<mass<106)")
)
process.mumuPairNoZFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("mumuPairNoZ"),
    minNumber = cms.uint32(1)
)
process.step2mumu = cms.Sequence( process.mumuPairNoZ*process.mumuPairNoZFilter)


# Require at least 2 jets
process.step3 = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("selectedPatJets"+postfix),
    minNumber = cms.uint32(2)
)


# MET filter
process.highMET = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("patMETs"+postfix),
    cut = cms.string("pt>30")
)

process.highMETFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("highMET"),
    minNumber = cms.uint32(1)
)

process.highMETemu = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("patMETs"+postfix),
    cut = cms.string("pt>20")
)

process.highMETFilteremu = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("highMETemu"),
    minNumber = cms.uint32(1)
)


process.step4 = cms.Sequence(process.highMET*process.highMETFilter)
process.step4emu = cms.Sequence(process.highMETemu*process.highMETFilteremu)

# trigger filter
trigger_paths = ["HLT_Mu17_Mu8_v"," HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"," HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7","HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"]
trigger_pattern = [path+"*" for path in trigger_paths]

from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::HLT",
      HLTPaths = trigger_pattern,
      throw=False
      )

# diagnostic stuff
# process.TFileService = cms.Service("TFileService",
#     fileName = cms.string('analyzeTopMuon.root')
# )
# process.load("TopQuarkAnalysis.Examples.TopJetAnalyzer_cfi")
# process.analyzeJet.input = "selectedPatJets"+postfix
# process.load("TopQuarkAnalysis.Examples.TopMuonAnalyzer_cfi")
# process.analyzeMuon.input = "selectedPatMuons"+postfix
# 
# process.load("BsmMasses.PrintMET.printmet_cfi")

# Scheduling
process.p_ee = cms.Path()
process.p_emu = cms.Path()
process.p_mumu = cms.Path()

process.p_common = cms.Sequence()

if signalMC or TauMC or TopNonDilMC:
    process.p_ee += process.makeGenEvt
    process.p_emu += process.makeGenEvt
    process.p_mumu += process.makeGenEvt
    #process.p_ee += process.ttDilFilter
    #process.p_emu += process.ttDilFilter
    #process.p_mumu += process.ttDilFilter

if runOnMC :
    process.p_ee += process.totalKinematicsFilter
    process.p_emu += process.totalKinematicsFilter
    process.p_mumu += process.totalKinematicsFilter
    
process.p_ee += process.eventCleaning
process.p_emu += process.eventCleaning
process.p_mumu += process.eventCleaning

process.p_ee += process.vertexing
process.p_emu += process.vertexing
process.p_mumu += process.vertexing

process.p_ee += process.eidCiCSequence
process.p_emu += process.eidCiCSequence
process.p_mumu += process.eidCiCSequence # this is probably still needed to make pat not crash


process.p_ee += getattr( process, 'patPF2PATSequence' + postfix )
process.p_emu += getattr( process, 'patPF2PATSequence' + postfix )
process.p_mumu += getattr( process, 'patPF2PATSequence' + postfix )

process.p_ee += process.step1ee
process.p_emu += process.step1emu
process.p_mumu += process.step1mumu

process.p_ee += process.step2ee
process.p_mumu += process.step2mumu

process.p_ee += process.step3
process.p_emu += process.step3
process.p_mumu += process.step3

# process.p_mumu += process.analyzeMuon
# process.p_mumu += process.analyzeJet
# process.p_mumu += process.printMET

process.p_ee += process.step4
process.p_mumu += process.step4
process.p_emu += process.step4emu

process.p_ee += process.triggerSelection
process.p_mumu += process.triggerSelection
process.p_emu += process.triggerSelection

process.out.SelectEvents.SelectEvents.append( 'p_mumu' )
process.out.SelectEvents.SelectEvents.append( 'p_emu' )
process.out.SelectEvents.SelectEvents.append( 'p_ee' )

# process.out.outputCommands.append( 'keep *' )

