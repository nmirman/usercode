import FWCore.ParameterSet.Config as cms

pfMetSig = cms.EDProducer('METSignificance',
      runOnMC = cms.untracked.bool(False),

      muonTag              = cms.untracked.InputTag('pfSelectedMuons'),
      electronTag          = cms.untracked.InputTag('pfSelectedElectrons'),
      metTag              = cms.untracked.InputTag('pfType1CorrectedMet'),

      pfjetsTag            = cms.untracked.InputTag('pfJets'),
      pfjetCorrectorL1     = cms.untracked.string('ak5PFL1Fastjet'),
      pfjetCorrectorL123   = cms.untracked.string('ak5PFL1FastL2L3'),

      jetResAlgo           = cms.string('AK5PF'),
      jetResEra            = cms.string('Spring10'),

      jetThreshold = cms.double(20),
      parA1 = cms.double(1.39669),
      parA2 = cms.double(1.32037),
      parA3 = cms.double(1.32047),
      parA4 = cms.double(1.38161),
      parA5 = cms.double(1.51508),
      parN1 = cms.double(0.0),
      parS1 = cms.double(0.639158)
)
