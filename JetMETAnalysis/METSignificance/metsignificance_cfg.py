import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load('JetMETAnalysis.METSignificance.metsignificance_cfi')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       'root://xrootd.unl.edu//store/data/Run2012A/DoubleMu/AOD/22Jan2013-v1/20000/001AE30A-BA81-E211-BBE7-003048FFD770.root'
    )
)

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = ( 'FT_P_V43E_AN3::All' )

########################################
###### PF2PAT Configuration BEGIN ######
########################################

# set data or MC
isMC = False

process.pfMetSig.runOnMC = cms.untracked.bool(isMC)
if not isMC:
   process.pfMetSig.pfjetCorrectorL123 = 'ak5PFL1FastL2L3Residual'

process.load("CommonTools.ParticleFlow.PF2PAT_cff")
process.pfPileUp.Enable = False

process.pfAllMuons.src="particleFlowPtrs"
process.pfMuonsFromVertex.dzCut=9999
process.pfNoMuon.bottomCollection   = "particleFlowPtrs"

process.pfJets.srcPVs = cms.InputTag("offlinePrimaryVertices")
process.pfJets.doAreaFastjet        = True
process.pfJets.jetPtMin             = 0
process.pfJets.src                  = "pfNoElectron"

process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = cms.bool(True)

process.mypf2pat = cms.Sequence(
      process.particleFlowPtrs *
      process.pfNoPileUpSequence * # pfPileUp enable is false
      process.pfParticleSelectionSequence *
      process.pfMuonSequence *
      process.pfNoMuon *
      process.pfElectronSequence *
      process.pfNoElectron *
      process.pfJets
      )

#
## MET Configuration
#

# type-I corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
if not isMC:
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
if not isMC:
   process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
if isMC:
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

# set met tag
process.pfMetSig.metsTag = cms.untracked.VInputTag(metList)

######################################
###### PF2PAT Configuration END ######
######################################

process.p = cms.Path(
      process.mypf2pat *
      process.mymet *
      process.pfMetSig
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

process.e = cms.EndPath(process.out)
