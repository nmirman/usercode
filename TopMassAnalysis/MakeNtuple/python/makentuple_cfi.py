import FWCore.ParameterSet.Config as cms

makentuple = cms.EDAnalyzer('MakeNtuple',
    muonSrc = cms.InputTag('selectedPatMuonsPFlow'),
    electronSrc = cms.InputTag('selectedPatElectronsPFlow'),
    jetSrc = cms.InputTag('selectedPatJetsPFlow'),
    metSrc = cms.InputTag('patMETsPFlow'),
    jetScale = cms.double( 0. ),
    jetResScale = cms.double( 1.1 ), # nominal value is actually 1.1*value in db
    leptonScale = cms.double( 0. ),
    bTagger = cms.string('combinedSecondaryVertexBJetTags'),
    bTagCut = cms.double( 0.244 ), # CSVT
    negTagCut = cms.double( 0.05 ),
    doMETCut = cms.bool( False ),
    runOnMC = cms.bool( True ),
    outFileName = cms.string('outputCSVL.root'),
    genParticleSrc = cms.InputTag('genParticles')
)
