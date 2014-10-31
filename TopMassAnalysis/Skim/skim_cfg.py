import sys

import FWCore.ParameterSet.Config as cms

# parse command line arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')

options.setDefault( 'maxEvents', -1)

options.setDefault( 'outputFile',
                  "topmassSkim.root" )
                  
options.setDefault( 'inputFiles',
                    #'/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v2/10000/00A4899E-666B-E211-A2AC-E0CB4E29C50D.root')
                    #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'
                    #'/store/caf/user/tjkim/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00000/000C5D15-AB1A-E211-8BDE-00215E22053A.root',
                    #'/store/caf/user/tjkim/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00000/0010005A-421A-E211-9E3C-E41F13181DA4.root'
                    'root://xrootd.unl.edu//store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/001C868B-B2E1-E111-9BE3-003048D4DCD8.root'
                    )

options.register( 'inputType',
                  'data',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'Type of input file. Options are Data or MC.'
                  )

options.register( 'globalTag',
                  'START53_V7A',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'globalTag')

options.register( 'dataType',
                  'DoubleMuData',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'Type of data file.'
                  )

options.register( 'keepAllEvts',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  'Flag to remove event filters.'
                  )

options.register( 'wantSummary',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Print summary at end of job")

options.parseArguments()

process = cms.Process('PAT')

if ( 'MC' in options.dataType ) :
   runOnMC = True
   trigger_paths = ["HLT_Mu17_Mu8_v17","HLT_Mu17_TkMu8_v10","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17","HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7","HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"]
   options.globalTag = 'START53_V7A'

elif ( 'Data' in options.dataType ) :
   runOnMC = False
   options.globalTag = 'FT_53_V21_AN4'
   if ( 'DoubleMu' in options.dataType ) :
      trigger_paths = ["HLT_Mu17_Mu8_v","HLT_Mu17_TkMu8_v"]
   elif ( 'DoubleElectron' in options.dataType ) :
      trigger_paths = ["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"]
   elif ( 'MuEG' in options.dataType ) :
      trigger_paths = ["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"]
   else : 
      sys.exit('dataType value was invalid.  The Data type does not exist.')


postfix = 'PFlow'

inputfiles = options.inputFiles

globalTag = options.globalTag

fwkReportEvery = 1000

outCommands = ['drop *']


####################################################
######### end user configurables ###################
####################################################


# basics
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool( options.wantSummary )
)
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = options.globalTag + '::All'

# input stuff
process.source = cms.Source( "PoolSource"
, noEventSort        = cms.untracked.bool( True )
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


# met filters
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
process.load('RecoMET.METFilters.metFilters_cff')

# scraping filter for data only
process.scrapingFilter = cms.EDFilter(
  "FilterOutScraping"
, applyfilter = cms.untracked.bool( True )
, debugOn     = cms.untracked.bool( False )
, numtrack    = cms.untracked.uint32( 10 )
, thresh      = cms.untracked.double( 0.25 )
)

process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
      filterParams = cms.PSet(
         minNdof = cms.double( 4. ),
         maxZ = cms.double( 24. ),
         maxRho = cms.double( 2. )
         ),
      filter = cms.bool( True),
      src = cms.InputTag( 'offlinePrimaryVertices' )
      )

process.vertexing = cms.Sequence(process.goodOfflinePrimaryVertices)


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

# object configurations from
# https://twiki.cern.ch/twiki/bin/view/CMS/TwikiTopRefHermeticTopProjections

# jet configuration
jetCutPF = 'pt>30 & abs(eta)<2.5'
jetCutPF += ' && numberOfDaughters > 1'                                  # PF jet ID:
jetCutPF += ' && neutralHadronEnergyFraction < 0.99'                     # PF jet ID:
jetCutPF += ' && neutralEmEnergyFraction < 0.99'                         # PF jet ID:
jetCutPF += ' && (chargedHadronEnergyFraction > 0. || abs(eta) >= 2.4)'  # PF jet ID:
jetCutPF += ' && (chargedMultiplicity > 0 || abs(eta) >= 2.4)'           # PF jet ID:
jetCutPF += ' && (chargedEmEnergyFraction < 0.99 || abs(eta) >= 2.4)'    # PF jet ID:

applyPostfix( process, 'pfJets', postfix ).doAreaFastjet = True
applyPostfix( process, 'pfJets', postfix ).doRhoFastjet  = False
applyPostfix( process, 'selectedPatJets', postfix ).cut  = jetCutPF


# muon configuration
applyPostfix( process, 'pfSelectedMuons', postfix ).cut = 'abs(eta)<2.5 && pt>10.'

applyPostfix( process, 'pfIsolatedMuons', postfix ).doDeltaBetaCorrection = True
applyPostfix( process, 'pfIsolatedMuons', postfix ).deltaBetaFactor = -0.5
applyPostfix( process, 'pfIsolatedMuons', postfix ).isolationCut = 0.20

applyPostfix( process, 'patMuons', postfix ).embedTrack = True
applyPostfix( process, 'patMuons', postfix ).usePV      = False
applyPostfix( process, 'selectedPatMuons', postfix ).cut = '''abs(eta)<2.4 && pt>20. &&
   (chargedHadronIso+max(0.,neutralHadronIso+photonIso-0.50*puChargedHadronIso))/pt < 0.20 &&
   (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''


# electron configuration
applyPostfix( process, 'pfElectronsFromVertex', postfix ).vertices = (
      cms.InputTag( 'goodOfflinePrimaryVertices' ) )
applyPostfix( process, 'pfElectronsFromVertex', postfix ).d0Cut    = 0.04
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.pfIdentifiedElectrons = cms.EDFilter("ElectronIDPFCandidateSelector", 
      recoGsfElectrons = cms.InputTag("gsfElectrons"),
      electronIdMap = cms.InputTag("mvaTrigV0"),
      electronIdCut = cms.double(0.0),
      src = cms.InputTag('pfElectronsFromVertex'+postfix))
process.load('EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi')
process.elPFIsoValueEA03.pfElectrons = cms.InputTag( 'pfSelectedElectrons'+postfix, '', '')

applyPostfix( process, 'pfSelectedElectrons', postfix ).src = 'pfIdentifiedElectrons'
applyPostfix( process, 'pfSelectedElectrons', postfix ).cut = '''abs(eta)<2.5 && pt>20. &&
   gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'''

applyPostfix( process, 'pfIsolatedElectrons', postfix ).isolationValueMapsCharged = (
      cms.VInputTag(cms.InputTag( 'elPFIsoValueCharged03PFId'+postfix )) )
applyPostfix( process, 'pfIsolatedElectrons', postfix ).isolationValueMapsNeutral = (
      cms.VInputTag(cms.InputTag( 'elPFIsoValueNeutral03PFId'+postfix ),
      cms.InputTag( 'elPFIsoValueGamma03PFId'+postfix )) )
applyPostfix( process, 'pfIsolatedElectrons', postfix ).deltaBetaIsolationValueMap = (
   'elPFIsoValueEA03' ) # new Effective Area edm::ValueMap
applyPostfix( process, 'pfIsolatedElectrons', postfix ).doDeltaBetaCorrection = (
   True ) # not really a 'deltaBeta' correction, but it serves
applyPostfix( process, 'pfIsolatedElectrons', postfix ).deltaBetaFactor = -1.0
applyPostfix( process, 'pfIsolatedElectrons', postfix ).isolationCut = 0.15

getattr( process, 'patPF2PATSequence'+postfix).replace(
      getattr( process, 'pfSelectedElectrons'+postfix ),
   process.mvaTrigV0 + process.pfIdentifiedElectrons +
   getattr( process, 'pfSelectedElectrons'+postfix ) + process.elPFIsoValueEA03 )

applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfPhotons = (
   'elPFIsoValueGamma03PFId'+postfix )
applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfNeutralHadrons = (
   'elPFIsoValueNeutral03PFId'+postfix )
applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfChargedHadrons = (
   'elPFIsoValueCharged03PFId'+postfix )
applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfChargedAll = (
   'elPFIsoValueChargedAll03PFId'+postfix )
applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfPUChargedHadrons = (
   'elPFIsoValuePU03PFId'+postfix )
applyPostfix( process, 'patElectrons', postfix ).isolationValues.user = (
   cms.VInputTag(cms.InputTag("elPFIsoValueEA03")) )
applyPostfix( process, 'patElectrons', postfix ).electronIDSources = (
   cms.PSet( mvaTrigV0 = cms.InputTag("mvaTrigV0")) )

applyPostfix( process, 'selectedPatElectrons', postfix ).cut = '''abs(eta)<2.5 && pt>20. && (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/pt < 0.15 && electronID("mvaTrigV0") > 0.5 && passConversionVeto'''

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


# pair selection
process.eePair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatElectrons"+postfix+"@- selectedPatElectrons"+postfix+"@+"),
    cut = cms.string("mass > 20"),
)
process.eePairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("eePair"),
    minNumber = cms.uint32(1),
)
process.step1ee = cms.Sequence(process.eePair*process.eePairFilter)

process.emuPair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatElectrons"+postfix+"@- selectedPatMuons"+postfix+"@+"),
    cut = cms.string("mass > 20"),
)
process.emuPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("emuPair"),
    minNumber = cms.uint32(1),
)
process.step1emu = cms.Sequence(process.emuPair*process.emuPairFilter)

process.mumuPair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatMuons"+postfix+"@- selectedPatMuons"+postfix+"@+"),
    cut = cms.string("mass > 20"),
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


## type-1 corrected MET -- run w/ process.producePFMETCorrections
#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#if not runOnMC:
#         process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

#process.pfType1CorrectedMetPFlow = process.pfType1CorrectedMet.clone(
#      src = cms.InputTag("pfMET"+postfix)
#      )

# load MET corrections

#process.load("JetMETCorrections.Type1MET.correctionTermsCaloMet_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
#process.load("JetMETCorrections.Type1MET.correctedMet_cff")
#
#if( runOnMC ):
#   process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
#   process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
#else:
#   process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
#   process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
#
#process.metcorr = cms.Path(
#      process.correctionTermsPfMetType1Type2 +
#      process.correctionTermsPfMetType0RecoTrack +
#      process.correctionTermsPfMetType0PFCandidate +
#      process.correctionTermsPfMetShiftXY +
#      process.correctionTermsCaloMet +
#      process.caloMetT1 + 
#      process.caloMetT1T2 + 
#      process.pfMetT0rt +
#      process.pfMetT0rtT1 +
#      process.pfMetT0pc +
#      process.pfMetT0pcT1 +
#      process.pfMetT0rtTxy +
#      process.pfMetT0rtT1Txy +
#      process.pfMetT0pcTxy +
#      process.pfMetT0pcT1Txy +
#      process.pfMetT1 +
#      process.pfMetT1Txy
#)


#
## MET Configuration
#

# type-I corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
if not runOnMC:
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
#process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data.clone(
#      numJetsMin = cms.int32(-1),
#      numJetsMax = cms.int32(-1),
#      px = cms.string("+4.83642e-02 + 2.48870e-01*Nvtx"),
#      py = cms.string("-1.50135e-01 - 8.27917e-02*Nvtx")
#      )
#process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc.clone(
#      numJetsMin = cms.int32(-1),
#      numJetsMax = cms.int32(-1),
#      px = cms.string("+1.62861e-01 - 2.38517e-02*Nvtx"),
#      py = cms.string("+3.60860e-01 - 1.30335e-01*Nvtx")
#      )

# most recent numbers for correction JEC '13
#process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data.px = cms.string("+4.83642e-02 + 2.48870e-01*Nvtx")
#process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data.py = cms.string("-1.50135e-01 - 8.27917e-02*Nvtx")
#process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data.px = cms.string("+1.62861e-01 - 2.38517e-02*Nvtx")
#process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data.py = cms.string("+3.60860e-01 - 1.30335e-01*Nvtx")

if not runOnMC:
   process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
if runOnMC:
   process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc

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

process.pfMET.jets = cms.InputTag('pfJetsPFlow')

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

# MET cut
process.highMET = cms.EDFilter("CandViewSelector",
    #src = cms.InputTag("patMETs"+postfix),
    #src = cms.InputTag("pfType1CorrectedMetPFLow"),
    #src = cms.InputTag("pfType1CorrectedMet"),
    src = cms.InputTag("pfType1CorrectedMetType0XYshift"),
    cut = cms.string("pt>40")
)

process.highMETFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("highMET"),
    minNumber = cms.uint32(1)
)

# add corrected MET to PAT
from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
process.patPfType1CorrectedMetType0XYshift = process.patMETsPFlow.clone(
      metSource = cms.InputTag('pfType1CorrectedMetType0XYshift'),
      addMuonCorrections = cms.bool(False),
      addGenMET    = cms.bool(False)
)

process.out.outputCommands += [ 'keep *_pfMet_*_*',
      'keep *_pfType1*_*_*',
      'keep *_ak5PFJets_*_*'
      ]

process.step4 = cms.Sequence(process.highMET*process.highMETFilter)
process.step4emu = cms.Sequence()

# trigger filter
trigger_pattern = [path+"*" for path in trigger_paths]

from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::HLT",
      HLTPaths = trigger_pattern,
      throw=False
      )

# Scheduling
process.p_ee = cms.Path()
process.p_emu = cms.Path()
process.p_mumu = cms.Path()

process.p_ee += process.vertexing
process.p_emu += process.vertexing
process.p_mumu += process.vertexing

process.p_ee += getattr( process, 'patPF2PATSequence' + postfix )
process.p_emu += getattr( process, 'patPF2PATSequence' + postfix )
process.p_mumu += getattr( process, 'patPF2PATSequence' + postfix )

if not runOnMC:
    process.p_ee += process.scrapingFilter
    process.p_emu += process.scrapingFilter
    process.p_mumu += process.scrapingFilter

process.p_ee += process.metFilters
process.p_emu += process.metFilters
process.p_mumu += process.metFilters

process.p_ee += process.pfIdentifiedElectrons
process.p_emu += process.pfIdentifiedElectrons
process.p_mumu += process.pfIdentifiedElectrons

process.p_ee += process.step1ee
process.p_emu += process.step1emu
process.p_mumu += process.step1mumu

process.p_ee += process.step2ee
process.p_mumu += process.step2mumu

process.p_ee += process.step3
process.p_emu += process.step3
process.p_mumu += process.step3

#process.p_ee += process.producePFMETCorrections
#process.p_emu += process.producePFMETCorrections
#process.p_mumu += process.producePFMETCorrections
process.p_ee += process.mymet
process.p_emu += process.mymet
process.p_mumu += process.mymet

process.p_ee += process.step4
process.p_mumu += process.step4
process.p_emu += process.step4emu

process.p_ee += process.triggerSelection
process.p_mumu += process.triggerSelection
process.p_emu += process.triggerSelection

process.out.SelectEvents.SelectEvents.append( 'p_mumu' )
process.out.SelectEvents.SelectEvents.append( 'p_emu' )
process.out.SelectEvents.SelectEvents.append( 'p_ee' )

if options.keepAllEvts:
   print 'KEEPING ALL EVENTS'
   process.out.SelectEvents.SelectEvents = []

#process.out.outputCommands.append( 'keep *' )

