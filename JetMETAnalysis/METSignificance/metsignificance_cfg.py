import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load('JetMETAnalysis.METSignificance.metsignificance_cfi')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v2/10000/00A4899E-666B-E211-A2AC-E0CB4E29C50D.root'
    )
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = ( 'FT_P_V43E_AN3::All' )

##################################
###### PF2PAT Configuration ######
##################################

process.load("CommonTools.ParticleFlow.PF2PAT_cff")
process.pfPileUp.Enable = False

process.pfAllMuons.src="particleFlow"
process.pfMuonsFromVertex.dzCut=9999
process.pfNoMuon.bottomCollection   = "particleFlow"
process.pfJets.doAreaFastjet        = True
process.pfJets.jetPtMin             = 0
process.pfJets.src                  = "pfNoElectron"

process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = cms.bool(True)

process.mypf2pat = cms.Sequence(
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

# met corrections and filters
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

process.mymet = cms.Sequence(
      process.pfMET *
      process.producePFMETCorrections
      )

process.pfMetSig.runOnMC = cms.untracked.bool(False)
if not process.pfMetSig.runOnMC:
   process.pfMetSig.pfjetCorrectorL123 = 'ak5PFL1FastL2L3Residual'

process.p = cms.Path(
      process.mypf2pat *
      process.mymet *
      process.pfMetSig
)

##################################
###### PF2PAT Configuration ######
##################################

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

process.e = cms.EndPath(process.out)
