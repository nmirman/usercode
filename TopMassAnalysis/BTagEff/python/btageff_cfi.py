import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('BTagEff',
      muonSrc = cms.InputTag('selectedPatMuonsPFlow'),
      electronSrc = cms.InputTag('selectedPatElectronsPFlow'),
      jetSrc = cms.InputTag('selectedPatJetsPFlow'),
      pmetSrc = cms.InputTag('patMETsPFlow'),
      metSrc = cms.InputTag('pfType1CorrectedMetType0XYshift'),
      jetScale = cms.double( 0. ),
      jetResScale = cms.double( 1.1 ), # nominal value is actually 1.1*value in db
      leptonScale = cms.double( 0. ),
      bTagger = cms.string('combinedSecondaryVertexBJetTags'),
      bTagCut = cms.double( 0.244 ), # CSVT
      negTagCut = cms.double( 0.05 ),
      doMETCut = cms.bool( False ),
      runOnMC = cms.bool( True ),
      outFileName = cms.string('outputCSVL.root'),
      genParticleSrc = cms.InputTag('genParticles'),
      randSeed = cms.int32( 0 ),
      PtNBins            = cms.int32(10),
      PtMin              = cms.double(0.),
      PtMax              = cms.double(1000.),
      EtaNBins           = cms.int32(4),
      EtaMin             = cms.double(0.),
      EtaMax             = cms.double(3.)
)