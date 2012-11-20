ISMC = False

import FWCore.ParameterSet.Config as cms
process = cms.Process("NTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'file:../DY_44X_V9B.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = ('GR_R_44_V13::All','START44_V5D::All')[ISMC]
process.GlobalTag.globaltag = ('GR_R_44_V15::All','START44_V13::All')[ISMC]
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")


from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
)

process.load("CommonTools.ParticleFlow.PF2PAT_cff")

process.pfAllMuons.src="particleFlow"
process.pfMuonsFromVertex.dzCut=9999
process.pfSelectedMuons.cut = cms.string('pt>15 && abs(eta)<2.4')
process.pfNoMuon.bottomCollection		 = "particleFlow"
process.pfNoMuon.topCollection			 = "pfSelectedMuons"
process.pfJets.doAreaFastjet			 = True
process.pfJets.jetPtMin				 = 0
process.pfJets.src				 = "pfNoMuon"

process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = cms.bool(True)

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJets = kt4PFJets.clone(
    src = cms.InputTag('pfNoMuon'),
    rParam = cms.double(0.6),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax = cms.double(2.5)
)

process.mypf2pat = cms.Sequence(process.pfAllMuons * 
				process.pfMuonsFromVertex*
				process.pfSelectedMuons*
				process.pfNoMuon*
				process.kt6PFJets * 
				process.pfJets*
				process.pfMET)

####MET modules
process.load("JetMETAnalysis.METSignificance.metsignificance_cfi")
process.pfMetSig.src = cms.VInputTag("pfJets", "particleFlow") 
#these should be reproducible from the ntuples...
process.pfMetSigSmeared    =  process.pfMetSig.clone(scaleResolutions=True) 
process.pfMetSigSmeared15  = (process.pfMetSigSmeared.clone(smearJetPtThreshold=15), process.pfMetSigSmeared.clone(smearJetPtThreshold=15, smearMet=True))[ISMC]

metList = []
metList.append(cms.untracked.InputTag("pfMet", "", ""))
metList.append(cms.untracked.InputTag("pfMET", "", ""))
metList.append(cms.untracked.InputTag("pfMetSig", "", "NTUPLE")) #this is the central one
metList.append(cms.untracked.InputTag("pfMetSigSmeared15", "","NTUPLE"))

process.mymets = cms.Sequence(process.pfMetSig * process.pfMetSigSmeared15)

#trigger_paths = ["HLT_IsoMu15_v",
#		 "HLT_IsoMu24_v",
#		 "HLT_Mu24_v"
#]
trigger_paths = ["HLT_Mu17_Mu8_v"]
trigger_pattern = [path+"*" for path in trigger_paths]


process.ntuple = cms.EDAnalyzer('Mtuple',
	debug		    = cms.untracked.bool(False),
	isMC		    = cms.untracked.bool(ISMC),
	OutputFileName	    = cms.untracked.string('ntuple.root'),

	saveTriggerResults  = cms.untracked.bool(False),
	TriggerResultsTag   = cms.untracked.InputTag("TriggerResults","","HLT"),
	TriggerEventTag	    = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
	TriggerPath	    = cms.untracked.vstring(trigger_paths),

	saveMuons	    = cms.untracked.bool(True),
	muonTag		    = cms.untracked.InputTag("muons"),

	saveParticles	    = cms.untracked.bool(False),
	genparticlesTag	    = cms.untracked.InputTag("genParticles"),
	pfcandidatesTag     = cms.untracked.InputTag("particleFlow"),
	
	saveJets	    = cms.untracked.int32(10000),
	pfjetsTag	    = cms.untracked.InputTag("pfJets"),
	pfjetCorrector	    = cms.untracked.string("ak5PFL1Fastjet"),
	pfjetCorrector2	    = cms.untracked.string("ak5PFL1FastL2L3"),
	jetResolAlgo	    = cms.string('AK5PF'),
	jetResolEra	    = cms.string('Spring10'),

   genjetsTag = cms.untracked.InputTag("ak5GenJets"),

	saveMETs	    = cms.untracked.bool(True),
	mets		    = cms.untracked.VInputTag(metList),
	genMet		    = cms.untracked.InputTag("myGenMetCalo"),

    	saveVertices	    = cms.untracked.bool(True),
	verticesTag	    = cms.untracked.InputTag("offlinePrimaryVertices"),

	pileupTag	    = cms.untracked.InputTag("addPileupInfo")
)
if not ISMC:
    process.ntuple.pfjetCorrector2 = "ak5PFL1FastL2L3Residual"


for i in xrange(len(metList)): print "%s: %s" %(i,metList[i])

#trigger filter			       
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = trigger_pattern, throw=False)

#I will only this for the data
process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )


process.load("RecoMET.Configuration.GenMETParticles_cff")
process.load("RecoMET.METProducers.genMetCalo_cfi")
process.myGenMetCalo = process.genMetCalo.clone(src = "genParticlesForMETAllVisible")


process.p = cms.Path( process.noscraping *
		      process.triggerSelection * 
		      process.mypf2pat *
		      process.mymets *
		      process.genParticlesForMETAllVisible*
		      process.myGenMetCalo*
		      process.ntuple
		      )

if ISMC:
    process.p.remove(process.noscraping)
    #process.p.remove(process.triggerSelection)
else:
    process.p.remove(process.genParticlesForMETAllVisible)
    process.p.remove(process.myGenMetCalo)

#process.out = cms.OutputModule( "PoolOutputModule"
#      , fileName = cms.untracked.string( "patTuple.root" )
#      )
#process.outpath = cms.EndPath( process.out )
