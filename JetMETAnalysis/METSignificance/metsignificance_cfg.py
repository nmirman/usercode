import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	'file:/cu1/aleko/CMSSW442/Asymmetry/111216/WJetsToLNu_TuneZ5D2_7TeV-madgraph-tauola_Fall11.root'
	#'file:/cu1/aleko/files/FA6FEF7F-8E7E-E011-AC35-001A92971B78.root'
	#'file:/cu1/aleko/files/FE4A8ED9-0A67-E011-9F19-00215E21D61E.root'
    )
)

process.load("JetMETAnalysis.METSignificance.metsignificance_cfi")


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('/cu1/aleko/files/myOutputFile.root')
)

process.load("RecoMET.METProducers.PFMET_cfi")
process.pfMet1 = process.pfMet.clone(alias = "")

  
process.p = cms.Path(process.pfMet1 * process.pfMetSig)

process.e = cms.EndPath(process.out)
