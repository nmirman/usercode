import FWCore.ParameterSet.Config as cms

pfMetSig = cms.EDProducer('METSignificance',
      runOnMC = cms.untracked.bool(True),

      muonTag              = cms.untracked.InputTag('pfSelectedMuons'),
      electronTag          = cms.untracked.InputTag('pfSelectedElectrons'),

      pfjetsTag            = cms.untracked.InputTag('pfJets'),
      pfjetCorrectorL1     = cms.untracked.string('ak5PFL1Fastjet'),
      pfjetCorrectorL123   = cms.untracked.string('ak5PFL1FastL2L3'),

      jetResAlgo           = cms.string('AK5PF'),
      jetResEra            = cms.string('Spring10'),

      jetThreshold = cms.double(20)
)
