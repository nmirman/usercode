import FWCore.ParameterSet.Config as cms

pfMetSig = cms.EDProducer('METSignificance',
    type			= cms.vstring("PFJet", "PFCandidate"),
    src				= cms.VInputTag("ak5PFJets", "particleFlow"),

    jetResolAlgo		= cms.string('AK5PF'),
    jetResolEra			= cms.string('Spring10'),

    smearMet			= cms.bool(False),
    scaleResolutions		= cms.bool(False),
    smearJetPtThreshold		= cms.double(15.0),
    resolJetPtThreshold		= cms.double(5.0),
    resolJetPhThreshold		= cms.double(10.0),
    jetParticleBoundaryPt	= cms.double(3.0),
    clusterAllParticlesAfterEta = cms.double(3.0), 
    scaleNeutralResols		= cms.double(1.8),
)

