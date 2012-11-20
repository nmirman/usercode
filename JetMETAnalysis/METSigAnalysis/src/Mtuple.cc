// -*- C++ -*-
//
// Package:    Mtuple
// Class:      Mtuple
// 
/**\class Mtuple Mtuple.cc EWK/Mtuple/src/Mtuple.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Aleko Khukhunaishvili,6 R-029,+41227678914,
//         Created:  Fri Jun 25 23:47:30 CEST 2010
// $Id: Mtuple.cc,v 1.2 2012/07/18 21:25:31 nmirman Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "JetMETAnalysis/METSigAnalysis/interface/MtupleFormats.h"

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

using namespace std;
using namespace edm;

bool comparePT(const reco::Muon &mu1, const reco::Muon &mu2){
   return (mu1.pt() > mu2.pt());
}


class Mtuple : public edm::EDAnalyzer {
   public:
      explicit Mtuple(const edm::ParameterSet&);
      ~Mtuple();


   private:
      virtual void beginJob() ;
      virtual void beginRun(const edm::Run &, const edm::EventSetup &);
      virtual void analyze(const edm::Event &, const edm::EventSetup &);
      void useOriginalPtrs(const edm::ProductID& productID);
      virtual void endJob() ;

      Bool_t		debug_;
      Bool_t		isMC_;
      std::string	OutputFileName_;

      Bool_t			    saveTriggerResults_;
      edm::InputTag		    TriggerResultsTag_;
      edm::InputTag		    TriggerEventTag_;  
      std::vector<std::string>	    TriggerPath_; 
      std::vector<edm::InputTag>    TriggerPathFilter_;    
      std::vector<std::string>	    TriggerPathVersioned_;
      Int_t  nTriggerPaths_;

      Bool_t		saveMuons_;
      edm::InputTag	muonTag_;

      Bool_t		saveMETs_;
      std::vector<edm::InputTag> metsTag_;
      Int_t		metsSize_;
      edm::InputTag	genMetTag_;

      Bool_t		saveParticles_;
      edm::InputTag	genparticlesTag_;
      edm::InputTag	pfcandidatesTag_;

      Int_t		saveJets_;
      edm::InputTag	pfjetsTag_;
      std::string	pfjetCorrector_;
      std::string	pfjetCorrector2_;

      edm::InputTag genjetsTag_;

      Bool_t		saveVertices_;
      edm::InputTag	verticesTag_;

      HLTConfigProvider hltConfig_;
      JetResolution *ptResol_;
      JetResolution *phiResol_;

      TFile		*OutFile__file;
      TTree		*results_tree;	 

      GenInfo		geninfo;
      GenW		W;

      RecoMuon		mu;
      Bool_t		hlt[20];
      float		prescale[20];
      METs		mets;
      float		genmet_et, genmet_phi, genmet_sumEt;
      Long64_t		run, event, lumi;
      PFCandidates	pfps;
      PFCandidates	pfes;
      PFJets		pfjs;
      GenJets     genjs;

      Vertices		vtxs;
      int		puN[3];
      int		puBC[3];
      float             puNMean;
      float		puRho;

      std::vector<reco::CandidatePtr> clusteredParticlePtrs_;
      std::vector<int>	jetIndex_;
};

Mtuple::Mtuple(const edm::ParameterSet& iConfig)
{
   using namespace std;
   debug_		= iConfig.getUntrackedParameter<Bool_t>("debug");    
   isMC_		= iConfig.getUntrackedParameter<Bool_t>("isMC");    
   OutputFileName_	= iConfig.getUntrackedParameter<std::string>("OutputFileName");

   //trigger
   saveTriggerResults_ = iConfig.getUntrackedParameter<Bool_t>("saveTriggerResults");
   TriggerResultsTag_	= iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultsTag");
   TriggerEventTag_    = iConfig.getUntrackedParameter<edm::InputTag>("TriggerEventTag");
   TriggerPath_	= iConfig.getUntrackedParameter<std::vector<string> >("TriggerPath");
   nTriggerPaths_	= TriggerPath_.size();

   //muons
   saveMuons_		= iConfig.getUntrackedParameter<Bool_t>("saveMuons");
   muonTag_		= iConfig.getUntrackedParameter<edm::InputTag>("muonTag");

   saveMETs_		= iConfig.getUntrackedParameter<Bool_t>("saveMETs");
   metsTag_		= iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("mets");
   genMetTag_		= iConfig.getUntrackedParameter<edm::InputTag>("genMet");
   metsSize_		= metsTag_.size();

   saveParticles_	= iConfig.getUntrackedParameter<Bool_t>("saveParticles");
   genparticlesTag_	= iConfig.getUntrackedParameter<edm::InputTag>("genparticlesTag");
   pfcandidatesTag_	= iConfig.getUntrackedParameter<edm::InputTag>("pfcandidatesTag");

   saveJets_		= iConfig.getUntrackedParameter<Int_t>("saveJets");
   pfjetsTag_		= iConfig.getUntrackedParameter<edm::InputTag>("pfjetsTag");
   pfjetCorrector_	= iConfig.getUntrackedParameter<std::string>("pfjetCorrector");
   pfjetCorrector2_	= iConfig.getUntrackedParameter<std::string>("pfjetCorrector2");

   genjetsTag_   = iConfig.getUntrackedParameter<edm::InputTag>("genjetsTag");

   saveVertices_	= iConfig.getUntrackedParameter<Bool_t>("saveVertices");
   verticesTag_	= iConfig.getUntrackedParameter<edm::InputTag>("verticesTag");

   using namespace std;
   string alg  = iConfig.getParameter<std::string>("jetResolAlgo");     
   string era  = iConfig.getParameter<std::string>("jetResolEra");     

   string path = "CondFormats/JetMETObjects/data";
   string ptFileName  = path + "/" + era + "_PtResolution_" +alg+".txt";
   string phiFileName = path + "/" + era + "_PhiResolution_"+alg+".txt";

   edm::FileInPath fpt(ptFileName);
   edm::FileInPath fphi(phiFileName);

   ptResol_  = new JetResolution(fpt.fullPath().c_str(),false);
   phiResol_ = new JetResolution(fphi.fullPath().c_str(),false);
}

Mtuple::~Mtuple(){}


   void
Mtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   run = iEvent.id().run();
   event = iEvent.id().event();
   lumi = iEvent.id().luminosityBlock();

   //some initialization to identify if there was such object or not.
   for(Int_t i=0; i<nTriggerPaths_; i++) {
      hlt[i]=false; 
      prescale[i]=0.;
   }

   if(isMC_){
      Handle<GenEventInfoProduct> gi;
      iEvent.getByLabel("generator", gi);

      edm::Handle< GenRunInfoProduct > genInfoProduct;
      iEvent.getRun().getByLabel("generator", genInfoProduct );

      geninfo.pid	= (int)gi->signalProcessID();
      geninfo.pthat	= gi->qScale();
      geninfo.weight	= gi->weight();
      geninfo.xsec	= genInfoProduct->crossSection();
      geninfo.eff	= genInfoProduct->filterEfficiency();
      geninfo.alphaQCD = gi->alphaQCD();
      geninfo.alphaQED = gi->alphaQED();

      if(gi->hasPDF()){
         geninfo.scalePDF = gi->pdf()->scalePDF;
         geninfo.id1	     = gi->pdf()->id.first;
         geninfo.id2	     = gi->pdf()->id.second;
         geninfo.x1	     = gi->pdf()->x.first;
         geninfo.x2	     = gi->pdf()->x.second;
         geninfo.xPDF1    = gi->pdf()->xPDF.first;
         geninfo.xPDF2    = gi->pdf()->xPDF.second;
      }

      Handle<std::vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByLabel("addPileupInfo", PupInfo);
      int pusize = PupInfo->size();
      if(pusize!=3) {
         std::cout << "WARNING: something wrong, check the pile-up info format." << std::endl;
         return;
      }
      puNMean = PupInfo->at(0).getTrueNumInteractions(); //same for each bunch-crossing
      for(int i=0; i<pusize; ++i){
         puN[i]  = PupInfo->at(i).getPU_NumInteractions();
         puBC[i] = PupInfo->at(i).getBunchCrossing();
      }
   }
   //get rho
   edm::Handle<double> rhoH;
   iEvent.getByLabel(edm::InputTag("kt6PFJets","rho","NTUPLE"),rhoH);
   puRho = *rhoH;

   edm::Handle<edm::View<reco::Candidate> > particles;
   //GenParticles
   if(isMC_){
      iEvent.getByLabel(genparticlesTag_, particles);

      bool wdone = false;
      bool edone = false;
      bool ndone = false;
      for( edm::View<reco::Candidate>::const_iterator iParticle = (particles.product())->begin(); iParticle != (particles.product())->end(); ++iParticle){
         int status = iParticle->status();
         int id = iParticle->pdgId();
         int aid = abs(id);
         if(status!=3) continue;

         if((aid==24 || aid==23 || aid==22) && !wdone){
            W.id	    = id;
            W.M		    = iParticle->mass();
            W.we	    = iParticle->energy();
            W.wpx	    = iParticle->px();
            W.wpy	    = iParticle->py();
            W.wpz	    = iParticle->pz();
            W.wpt	    = iParticle->pt();
            W.wphi	    = iParticle->phi();
            W.weta	    = iParticle->eta();
            W.wy	    = iParticle->rapidity();
            wdone = true;
         }

         else if(aid==13 && !edone){
            W.m	    = iParticle->mass();
            W.ee	    = iParticle->energy();
            W.epx	    = iParticle->px();
            W.epy	    = iParticle->py();
            W.epz	    = iParticle->pz();
            W.ept	    = iParticle->pt();
            W.ephi	    = iParticle->phi();
            W.eeta	    = iParticle->eta();
            W.eid	    = id;
            edone = true;
         }

         else if((aid==14 || aid==13) && !ndone){
            W.ne	    = iParticle->energy();
            W.npx	    = iParticle->px();
            W.npy	    = iParticle->py();
            W.npz	    = iParticle->pz();
            W.npt	    = iParticle->pt();
            W.nphi	    = iParticle->phi();
            W.neta	    = iParticle->eta();
            W.nid	    = id;
            ndone = true;
         }

         if(wdone && edone && ndone) break;
      } 
   }


   //HLT Trigger
   edm::Handle<TriggerResults> trh;
   if(saveTriggerResults_){
      iEvent.getByLabel(TriggerResultsTag_, trh);
      if(trh.isValid()){
         edm::TriggerNames  triggerNames_ = iEvent.triggerNames(*trh);
         if(debug_){
            for(int i=0; i<(int)triggerNames_.size(); ++i){
               TString eleTrigger = triggerNames_.triggerName(i);
               if(eleTrigger.Contains("HLT_Mu"))
                  std::cout << triggerNames_.triggerName(i) << "   " << trh->accept(i) << "   " << hltConfig_.prescaleValue(iEvent, iSetup, std::string(eleTrigger)) << std::endl; 
            }
         }
         for(Int_t i=0; i<nTriggerPaths_; ++i){
            Int_t index = triggerNames_.triggerIndex(TriggerPathVersioned_[i]); 
            if(index<(Int_t)trh->size()) {
               hlt[i] = trh->accept(index);
               prescale[i] = hltConfig_.prescaleValue(iEvent, iSetup, TriggerPathVersioned_[i]);
            }
         }
      }
   }


   edm::Handle<edm::View<reco::Vertex> > vertices;
   if(saveVertices_){
      //GenParticles
      for(Int_t i=0; i<NVertices; i++){
         vtxs.isFake[i] = true;
         vtxs.ndof[i]= -1.;
         vtxs.x[i] = 999.;
         vtxs.y[i] = 999.;
         vtxs.z[i] = 999.;
         vtxs.Rho[i] = -1.;
      }
      Int_t icand=0;
      iEvent.getByLabel(verticesTag_, vertices);
      for( edm::View<reco::Vertex>::const_iterator v = vertices->begin() ; v!= vertices->end() && icand<NVertices; v++ ){
         vtxs.isFake[icand]  = v->isFake();
         vtxs.ndof[icand] = v->ndof();
         vtxs.x[icand] = v->x();
         vtxs.y[icand] = v->y();
         vtxs.z[icand] = v->z();
         vtxs.Rho[icand] = v->position().Rho();
         icand++;
      } 
      vtxs.size = icand;
   }

   if(saveMuons_){
      //For HLT matching
      edm::Handle<trigger::TriggerEvent> pHLT;
      iEvent.getByLabel(TriggerEventTag_, pHLT);
      const trigger::TriggerObjectCollection& TOC(pHLT->getObjects());
      const Int_t nF(pHLT->sizeFilters());
      std::vector<Int_t> filterIndex;
      for(Int_t itrig=0; itrig<nTriggerPaths_; itrig++){
         Int_t index = pHLT->filterIndex(TriggerPathFilter_[itrig]);
         filterIndex.push_back(index);
      }

      edm::Handle<std::vector<reco::Muon> > muonsHandle;
      iEvent.getByLabel(muonTag_, muonsHandle);
      std::vector<reco::Muon> muons = *muonsHandle;
      std::sort(muons.begin(), muons.end(), comparePT);
      mu.size = muons.size();

      edm::Handle<reco::BeamSpot> beamSpotHandle;
      iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle);

      int i=0;
      for(std::vector<reco::Muon>::const_iterator it=muons.begin(); it!=muons.end(); ++it, ++i){
         reco::Muon muon = *it;

         mu.charge[i]    = muon.charge();
         mu.pt[i]	    = muon.pt();
         mu.p[i]	    = muon.p();
         mu.e[i]	    = muon.energy();
         mu.phi[i]	    = muon.phi();
         mu.eta[i]	    = muon.eta();
         mu.px[i]	    = muon.px();
         mu.py[i]	    = muon.py();
         mu.pz[i]	    = muon.pz();

         mu.isGlobal[i]  = muon.isGlobalMuon() ? 1 : 0; 
         mu.isTracker[i] = muon.isTrackerMuon() ? 1 : 0;

         mu.glpt[i]=-1;
         if(muon.isGlobalMuon()){
            reco::TrackRef	gm  = muon.globalTrack();
            mu.dxy[i]	    = gm->dxy(beamSpotHandle->position());
            mu.chi2[i]	    = gm->normalizedChi2();
            mu.trackerHits[i]   = gm->hitPattern().numberOfValidTrackerHits();
            mu.pixelHits[i]	    = gm->hitPattern().numberOfValidPixelHits();
            mu.muonHits[i]	    = gm->hitPattern().numberOfValidMuonHits();
            mu.nMatches[i]	    = muon.numberOfMatchedStations();
            mu.glpt[i]	    = gm->pt();
         }
         mu.trpt[i]=-1;
         if(muon.isTrackerMuon()){
            reco::TrackRef	it = muon.innerTrack();
            mu.numberOfValidHits[i]			= it->hitPattern().numberOfValidHits();
            mu.numberOfValidPixelHits[i]		= it->hitPattern().numberOfValidPixelHits();
            mu.numberOfValidStripHits[i]		= it->hitPattern().numberOfValidStripHits();
            mu.numberOfValidPixelBarrelHits[i]	= it->hitPattern().numberOfValidPixelBarrelHits();
            mu.pixelLayersWithMeasurement[i]	= it->hitPattern().pixelLayersWithMeasurement();
            mu.pixelBarrelLayersWithMeasurement[i]	= it->hitPattern().pixelBarrelLayersWithMeasurement();
            mu.pixelEndcapLayersWithMeasurement[i]	= it->hitPattern().pixelEndcapLayersWithMeasurement();
            mu.numberOfValidTrackerLayers[i]	= it->hitPattern().trackerLayersWithMeasurement();
            mu.trpt[i]				= it->pt();
         }
         // hlt matching
         // Save the smallest DeltaR among all trigger objects
         if(trh.isValid()){
            for(Int_t itrig=0; itrig<nTriggerPaths_; itrig++){
               mu.drTO[itrig][i] = 999.;
               if(filterIndex[itrig] == nF) continue;
               const trigger::Keys& KEYS(pHLT->filterKeys(filterIndex[itrig]));
               const Int_t nK(KEYS.size());
               for (Int_t iTrig = 0;iTrig<nK; ++iTrig ) {
                  const trigger::TriggerObject& TO(TOC[KEYS[iTrig]]);
                  Double_t dr = reco::deltaR(muon.eta(),muon.phi(),TO.eta(),TO.phi());
                  if(dr<mu.drTO[itrig][i]) {
                     mu.drTO[itrig][i] = dr;
                  }
               }
            }
         }

         mu.dr03TkSumPt[i]		= muon.isolationR03().sumPt;
         mu.dr03EcalRecHitSumEt[i]	= muon.isolationR03().emEt;
         mu.dr03HcalTowerSumEt[i]	= muon.isolationR03().hadEt;

         if(!isMC_) continue;

         //gen matching... probably will be needed to scale correction later...
         int k=-1;
         int index=0;
         float drmin=9999;
         for(edm::View<reco::Candidate>::const_iterator iParticle = particles->begin(); iParticle != particles->end(); ++iParticle, ++index){
            if(abs(iParticle->pdgId())!=13) continue;
            float dr = reco::deltaR(muon.eta(),muon.phi(),iParticle->eta(),iParticle->phi());
            if(dr<drmin){
               k=index;
               drmin=dr;
            }
         }
         mu.mid[i]=0;
         if(k>=0){
            const reco::Candidate &gmu = (*particles)[k];
            mu.gdr[i]     = drmin;
            mu.mid[i]     = gmu.pdgId();
            mu.mpx[i]     = gmu.px();
            mu.mpy[i]     = gmu.py();
            mu.mpz[i]     = gmu.pz();
            mu.men[i]     = gmu.energy();
         }
      }
   }

   //MET's
   if(saveMETs_){
      for(Int_t imet=0; imet<metsSize_ && imet<NMETs;imet++){
         edm::Handle<edm::View<reco::MET> > metHandle;
         iEvent.getByLabel(metsTag_[imet], metHandle);
         reco::MET metiter = (*metHandle)[0];
         mets.et[imet]	= metiter.pt();
         mets.phi[imet]	= metiter.phi();
         mets.sumEt[imet]	= metiter.sumEt();
         mets.dxx[imet]	= metiter.getSignificanceMatrix()(0,0); 
         mets.dxy[imet]	= metiter.getSignificanceMatrix()(0,1); 
         mets.dyy[imet]	= metiter.getSignificanceMatrix()(1,1); 
         mets.sig[imet]	= metiter.significance();
      }
      //add gen met here
      if(isMC_){
         edm::Handle<edm::View<reco::MET> > genmetHandle;
         iEvent.getByLabel(genMetTag_, genmetHandle);
         reco::MET genmet	= (*genmetHandle)[0];
         genmet_et		= genmet.pt();
         genmet_phi		= genmet.phi();
         genmet_sumEt	= genmet.sumEt();
      }
   }

   if(saveJets_>0){
      //PFJets
      for(Int_t i=0; i<saveJets_; i++){
         pfjs.nco[i]=0;
         pfjs.l1[i]=0;
         pfjs.l1l2l3[i]=0;
         pfjs.pt[i]=0;
         pfjs.phi[i]=0;
         pfjs.eta[i]=0;
         pfjs.energy[i]=0;
         pfjs.dpt[i]=0;
         pfjs.dphi[i]=0;

         pfjs.neutralHadronFraction[i]=0;
         pfjs.neutralEmFraction[i]=0;
         pfjs.chargedHadronFraction[i]=0;
         pfjs.chargedHadronMultiplicity[i]=0;
         pfjs.chargedEmFraction[i]=0;
      }

      Handle<reco::PFJetCollection> inputUncorJets;
      iEvent.getByLabel( pfjetsTag_, inputUncorJets );
      const JetCorrector* corrector  = JetCorrector::getJetCorrector (pfjetCorrector_, iSetup);
      const JetCorrector* corrector2 = JetCorrector::getJetCorrector (pfjetCorrector2_, iSetup);

      clusteredParticlePtrs_.clear();
      jetIndex_.clear();

      Int_t icand = 0;
      for(reco::PFJetCollection::const_iterator jet = inputUncorJets->begin(); jet != inputUncorJets->end() && icand<saveJets_; ++jet) {
         int index = jet - inputUncorJets->begin();
         edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::PFJetCollection>(inputUncorJets,index));
         pfjs.nco[icand]  = jet->nConstituents();
         pfjs.l1[icand]   = corrector->correction (*jet, iEvent, iSetup);
         pfjs.l1l2l3[icand]= corrector2->correction (*jet, iEvent, iSetup);
         pfjs.pt[icand]   = jet->pt();
         pfjs.phi[icand]  = jet->phi();
         pfjs.eta[icand]  = jet->eta();
         pfjs.energy[icand] = jet->energy();

         double jpt  = jet->pt();
         double jeta = jet->eta();

         pfjs.neutralHadronFraction[icand] = jet->neutralHadronEnergyFraction();
         pfjs.neutralEmFraction[icand] = jet->neutralEmEnergyFraction();
         pfjs.chargedHadronFraction[icand] = jet->chargedHadronEnergyFraction();
         pfjs.chargedHadronMultiplicity[icand] = jet->chargedHadronMultiplicity();
         pfjs.chargedEmFraction[icand] = jet->chargedEmEnergyFraction();

         TF1* fPtEta		= ptResol_ -> parameterEta("sigma",jeta);
         TF1* fPhiEta	= phiResol_-> parameterEta("sigma",jeta);
         pfjs.dpt[icand]     = jpt*fPtEta->Eval(jpt); 
         pfjs.dphi[icand]    = jpt*fPhiEta->Eval(jpt);
         delete fPtEta;
         delete fPhiEta;

         std::vector<reco::PFCandidatePtr> pfs = jet->getPFConstituents();
         for(std::vector<reco::PFCandidatePtr>::const_iterator it=pfs.begin(); it!=pfs.end(); ++it){
            reco::CandidatePtr ptr(*it);
            clusteredParticlePtrs_.push_back(ptr);
            jetIndex_.push_back(icand);
         }

         icand++;
      }
      pfjs.size = icand;


      // genjets
      if(isMC_){
         Handle<reco::GenJetCollection> genjets;
         iEvent.getByLabel( genjetsTag_, genjets );

         icand = 0;
         for(reco::GenJetCollection::const_iterator jet = genjets->begin(); jet != genjets->end(); ++jet) {

            genjs.pt[icand] = jet->pt();
            genjs.phi[icand] = jet->phi();
            genjs.eta[icand] = jet->eta();

            icand++;
         }
         genjs.size = icand;
      }
   }

   if(saveParticles_){
      edm::Handle<edm::View<reco::PFCandidate> > pfcandidates;
      iEvent.getByLabel(pfcandidatesTag_, pfcandidates);
      useOriginalPtrs(pfcandidates.id()); //this is needed with PF2PAT jet constituents, to find the original pf-candidate reference
      edm::Handle<edm::View<reco::Vertex> > vertices;
      iEvent.getByLabel(verticesTag_, vertices);

      for(unsigned int i=0; i<pfcandidates->size(); ++i){
         reco::PFCandidatePtr pf(pfcandidates,i);
         Int_t type	    = pf->particleId();
         pfps.type[i]    = type;
         pfps.pt[i]	    = pf->pt();
         pfps.et[i]	    = pf->et();
         pfps.phi[i]	    = pf->phi();
         pfps.eta[i]	    = pf->eta();

         reco::CandidatePtr ptr(pfcandidates,i);
         pfps.jetIndex[i]=-1;
         for(unsigned int j=0; j<clusteredParticlePtrs_.size(); ++j){
            if(ptr==clusteredParticlePtrs_.at(j)){
               pfps.jetIndex[i]=jetIndex_[j];
               break;
            }
         }
      }
      pfps.size = pfcandidates->size();
   }

   bool mu_selection = (mu.size == 2) and (mu.pt[0] > 20) and (mu.pt[1] > 20);
   //bool mu_tight = true;
   //bool mu_iso = true;
   //for(int i=0; i < mu.size; i++){
      //tight muon selection
   //}
   //FILL THE TREE
   if( mu_selection ){
      results_tree -> Fill();
   }

}


void 
Mtuple::useOriginalPtrs(const edm::ProductID& productID){
   std::vector<reco::CandidatePtr>::const_iterator it=clusteredParticlePtrs_.begin();

   reco::CandidatePtr ptr(*it);
   if(ptr.id()==productID) return; //If the first element is from the right product, return

   std::vector<reco::CandidatePtr> temp;
   for(; it!=clusteredParticlePtrs_.end(); ++it){
      reco::CandidatePtr ptr(*it);
      while(ptr.id()!=productID){
         ptr = ptr->sourceCandidatePtr(0);
         if(ptr.isNull()) return; //if it does not get to the correct product, return
      }
      temp.push_back(ptr);
   }
   clusteredParticlePtrs_.clear();
   clusteredParticlePtrs_ = temp;
}

   void 
Mtuple::beginJob()
{
   OutFile__file  = new TFile( OutputFileName_.c_str(), "RECREATE" );

   results_tree = new TTree("events", "events");
   results_tree -> Branch("run", &run, "run/I");
   results_tree -> Branch("lumi", &lumi, "lumi/I");
   results_tree -> Branch("event", &event, "event/I");
   if(isMC_){
      results_tree -> Branch("gi_pid", &geninfo.pid, "gi_pid/I");
      results_tree -> Branch("gi_pthat", &geninfo.pthat, "gi_pthat/F");
      results_tree -> Branch("gi_alphaQCD", &geninfo.alphaQCD, "gi_alphaQCD/F");
      results_tree -> Branch("gi_alphaQED", &geninfo.alphaQED, "gi_alphaQED/F");
      results_tree -> Branch("gi_scalePDF", &geninfo.scalePDF, "gi_scalePDF/F");
      results_tree -> Branch("gi_id1", &geninfo.id1, "gi_id1/I");
      results_tree -> Branch("gi_id2", &geninfo.id2, "gi_id2/I");
      results_tree -> Branch("gi_x1", &geninfo.x1, "gi_x1/F");
      results_tree -> Branch("gi_x2", &geninfo.x2, "gi_x2/F");
      results_tree -> Branch("gi_xPDF1", &geninfo.xPDF1, "gi_xPDF1/F");
      results_tree -> Branch("gi_xPDF2", &geninfo.xPDF2, "gi_xPDF2/F");
      results_tree -> Branch("gi_weight", &geninfo.weight, "gi_weight/F");
      results_tree -> Branch("gi_xsec", &geninfo.xsec, "gi_xsec/F");
      results_tree -> Branch("gi_eff", &geninfo.eff, "gi_eff/F");

      results_tree -> Branch("puN",   puN,   "puN[3]/I");
      results_tree -> Branch("puBC",   puBC,   "puBC[3]/I");
      results_tree -> Branch("puNMean",   &puNMean,   "puNMean/F");
   }

   results_tree -> Branch("puRho",   &puRho,   "puRho/F");
   results_tree -> Branch("nTriggerPaths_", &nTriggerPaths_, "nTriggerPaths_/I");
   results_tree -> Branch("hlt", &hlt, "hlt[nTriggerPaths_]/O");
   results_tree -> Branch("prescale", &prescale, "prescale[nTriggerPaths_]/F");


   if(saveMuons_){
      results_tree -> Branch("mu_size",   &mu.size, "mu_size/I");
      results_tree -> Branch("mu_charge",   mu.charge, "mu_charge[mu_size]/I");
      results_tree -> Branch("mu_isGlobal", mu.isGlobal, "mu_isGlobal[mu_size]/I");
      results_tree -> Branch("mu_isTracker", mu.isTracker, "mu_isTracker[mu_size]/I");
      results_tree -> Branch("mu_trackerHits", mu.trackerHits, "mu_trackerHits[mu_size]/I");
      results_tree -> Branch("mu_pixelHits", mu.pixelHits, "mu_pixelHits[mu_size]/I");
      results_tree -> Branch("mu_muonHits", mu.muonHits, "mu_muonHits[mu_size]/I");
      results_tree -> Branch("mu_nMatches", mu.nMatches, "mu_nMatches[mu_size]/I");

      results_tree -> Branch("mu_numberOfValidHits", mu.numberOfValidHits, "mu_numberOfValidHits[mu_size]/I");
      results_tree -> Branch("mu_numberOfValidPixelHits", mu.numberOfValidPixelHits, "mu_numberOfValidPixelHits[mu_size]/I");
      results_tree -> Branch("mu_numberOfValidStripHits", mu.numberOfValidStripHits, "mu_numberOfValidStripHits[mu_size]/I");
      results_tree -> Branch("mu_numberOfValidPixelBarrelHits", mu.numberOfValidPixelBarrelHits, "mu_numberOfValidPixelBarrelHits[mu_size]/I");
      results_tree -> Branch("mu_pixelLayersWithMeasurement", mu.pixelLayersWithMeasurement, "mu_pixelLayersWithMeasurement[mu_size]/I");
      results_tree -> Branch("mu_pixelBarrelLayersWithMeasurement", mu.pixelBarrelLayersWithMeasurement, "mu_pixelBarrelLayersWithMeasurement[mu_size]/I");
      results_tree -> Branch("mu_pixelEndcapLayersWithMeasurement", mu.pixelEndcapLayersWithMeasurement, "mu_pixelEndcapLayersWithMeasurement[mu_size]/I");
      results_tree -> Branch("mu_numberOfValidTrackerLayers", mu.numberOfValidTrackerLayers, "mu_numberOfValidTrackerLayers[mu_size]/I");

      results_tree -> Branch("mu_chi2", mu.chi2, "mu_chi2[mu_size]/F");
      results_tree -> Branch("mu_pt", mu.pt, "mu_pt[mu_size]/F");
      results_tree -> Branch("mu_glpt", mu.glpt, "mu_glpt[mu_size]/F");
      results_tree -> Branch("mu_trpt", mu.trpt, "mu_trpt[mu_size]/F");
      results_tree -> Branch("mu_p", mu.p, "mu_p[mu_size]/F");
      results_tree -> Branch("mu_e", mu.e, "mu_e[mu_size]/F");
      results_tree -> Branch("mu_phi", mu.phi, "mu_phi[mu_size]/F");
      results_tree -> Branch("mu_eta", mu.eta, "mu_eta[mu_size]/F");
      results_tree -> Branch("mu_px", mu.px, "mu_px[mu_size]/F");
      results_tree -> Branch("mu_py", mu.py, "mu_py[mu_size]/F");
      results_tree -> Branch("mu_pz", mu.pz, "mu_pz[mu_size]/F");
      results_tree -> Branch("mu_dxy", mu.dxy, "mu_dxy[mu_size]/F");
      results_tree -> Branch("mu_dr03TkSumPt", mu.dr03TkSumPt, "mu_dr03TkSumPt[mu_size]/F");
      results_tree -> Branch("mu_dr03EcalRecHitSumEt", mu.dr03EcalRecHitSumEt, "mu_dr03EcalRecHitSumEt[mu_size]/F");
      results_tree -> Branch("mu_dr03HcalTowerSumEt", mu.dr03HcalTowerSumEt, "mu_dr03HcalTowerSumEt[mu_size]/F");

      //trigger matching
      for(Int_t j=0; j<nTriggerPaths_; j++){
         results_tree -> Branch(Form("mu_DeltaRHLT_%d", j), mu.drTO[j], Form("mu_DeltaRHLT_%d[mu_size]/F", j));
      }

      results_tree -> Branch("mu_gdr", mu.gdr, "mu_gdr[mu_size]/F");
      results_tree -> Branch("mu_mid", mu.mid, "mu_mid[mu_size]/I");
      results_tree -> Branch("mu_mpx", mu.mpx, "mu_mpx[mu_size]/F");
      results_tree -> Branch("mu_mpy", mu.mpy, "mu_mpy[mu_size]/F");
      results_tree -> Branch("mu_mpz", mu.mpz, "mu_mpz[mu_size]/F");
      results_tree -> Branch("mu_men", mu.men, "mu_men[mu_size]/F");
   }


   if(saveMETs_){
      results_tree -> Branch("met_size", &metsSize_, "met_size/I");
      results_tree -> Branch("met_et", mets.et, "met_et[met_size]/F");
      results_tree -> Branch("met_phi", mets.phi, "met_phi[met_size]/F");
      results_tree -> Branch("met_sig", mets.sig, "met_sig[met_size]/F");
      results_tree -> Branch("met_sumEt", mets.sumEt, "met_sumEt[met_size]/F");
      results_tree -> Branch("met_dxx", mets.dxx, "met_dxx[met_size]/F");
      results_tree -> Branch("met_dxy", mets.dxy, "met_dxy[met_size]/F");
      results_tree -> Branch("met_dyy", mets.dyy, "met_dyy[met_size]/F");

      if(isMC_){
         results_tree -> Branch("genmet_et", &genmet_et, "genmet_et/F");
         results_tree -> Branch("genmet_phi", &genmet_phi, "genmet_phi/F");
         results_tree -> Branch("genmet_sumEt", &genmet_sumEt, "genmet_sumEt/F");
      }
   }

   if(saveParticles_){
      results_tree -> Branch("pf_size", &pfps.size, "pf_size/I");
      results_tree -> Branch("pf_type", pfps.type, "pf_type[pf_size]/I");
      results_tree -> Branch("pf_et",   pfps.et,   "pf_et[pf_size]/F");
      results_tree -> Branch("pf_pt",   pfps.pt,   "pf_pt[pf_size]/F");
      results_tree -> Branch("pf_phi",  pfps.phi,  "pf_phi[pf_size]/F");
      results_tree -> Branch("pf_eta",  pfps.eta,  "pf_eta[pf_size]/F");
      results_tree -> Branch("pf_jetindex", pfps.jetIndex,  "pf_jetindex[pf_size]/I");
   }

   if(isMC_){
      results_tree -> Branch("w_id",     &W.id,     "w_id/I");
      results_tree -> Branch("w_M",	   &W.M,      "w_M/D");
      results_tree -> Branch("w_we",	   &W.we,     "w_we/D");
      results_tree -> Branch("w_wpx",	   &W.wpx,    "w_wpx/D");
      results_tree -> Branch("w_wpy",    &W.wpy,    "w_wpy/D");
      results_tree -> Branch("w_wpz",    &W.wpz,    "w_wpz/D");
      results_tree -> Branch("w_wpt",    &W.wpt,    "w_wpt/D");
      results_tree -> Branch("w_weta",   &W.weta,   "w_weta/D");
      results_tree -> Branch("w_wphi",   &W.wphi,   "w_wphi/D");
      results_tree -> Branch("w_wy",     &W.wy,     "w_wy/D");

      results_tree -> Branch("w_m",	    &W.m,      "w_m/D");
      results_tree -> Branch("w_ee",	    &W.ee,     "w_ee/D");
      results_tree -> Branch("w_epx",	    &W.epx,    "w_epx/D");
      results_tree -> Branch("w_epy",	    &W.epy,    "w_epy/D");
      results_tree -> Branch("w_epz",	    &W.epz,    "w_epz/D");
      results_tree -> Branch("w_ept",	    &W.ept,    "w_ept/D");
      results_tree -> Branch("w_eeta",    &W.eeta,   "w_eeta/D");
      results_tree -> Branch("w_ephi",    &W.ephi,   "w_ephi/D");
      results_tree -> Branch("w_eid",     &W.eid,    "w_eid/I");

      results_tree -> Branch("w_ne",      &W.ne,     "w_ne/D");   
      results_tree -> Branch("w_npx",     &W.npx,    "w_npx/D");
      results_tree -> Branch("w_npy",     &W.npy,    "w_npy/D");
      results_tree -> Branch("w_npz",     &W.npz,    "w_npz/D");
      results_tree -> Branch("w_npt",     &W.npt,    "w_npt/D");
      results_tree -> Branch("w_neta",    &W.neta,   "w_neta/D");
      results_tree -> Branch("w_nphi",    &W.nphi,   "w_nphi/D");
      results_tree -> Branch("w_nid",     &W.nid,    "w_nid/I");
   }

   if(saveJets_){
      results_tree -> Branch("pfj_size",   &pfjs.size, "pfj_size/I");
      results_tree -> Branch("pfj_l1",     pfjs.l1,  "pfj_l1[pfj_size]/F");
      results_tree -> Branch("pfj_l1l2l3",   pfjs.l1l2l3,  "pfj_l1l2l3[pfj_size]/F");
      results_tree -> Branch("pfj_pt",     pfjs.pt,    "pfj_pt[pfj_size]/F");
      results_tree -> Branch("pfj_phi",    pfjs.phi,   "pfj_phi[pfj_size]/F");
      results_tree -> Branch("pfj_eta",    pfjs.eta,   "pfj_eta[pfj_size]/F");
      results_tree -> Branch("pfj_energy",    pfjs.energy,   "pfj_eta[pfj_size]/F");
      results_tree -> Branch("pfj_dpt",    pfjs.dpt,   "pfj_dpt[pfj_size]/F");
      results_tree -> Branch("pfj_dphi",   pfjs.dphi,  "pfj_dphi[pfj_size]/F");

      results_tree -> Branch("pfj_neutralHadronFraction", pfjs.neutralHadronFraction,
            "pfj_neutralHadronFraction[pfj_size]/F");
      results_tree -> Branch("pfj_neutralEmFraction", pfjs.neutralEmFraction,
            "pfj_neutralEmFraction[pfj_size]/F");
      results_tree -> Branch("pfj_chargedHadronFraction", pfjs.chargedHadronFraction,
            "pfj_chargedHadronFraction[pfj_size]/F");
      results_tree -> Branch("pfj_chargedHadronMultiplicity", pfjs.chargedHadronMultiplicity,
            "pfj_chargedHadronMultiplicity[pfj_size]/F");
      results_tree -> Branch("pfj_chargedEmFraction", pfjs.chargedEmFraction,
            "pfj_chargedEmFraction[pfj_size]/F");
      results_tree -> Branch("pfj_numConstituents", pfjs.nco, "pfj_numConstituents[pfj_size]/I");

      if(isMC_){
         results_tree -> Branch("genj_size",  &genjs.size, "genj_size/I");
         results_tree -> Branch("genj_pt",  genjs.pt, "genj_pt[genj_size]/F");
         results_tree -> Branch("genj_phi",  genjs.phi, "genj_phi[genj_size]/F");
         results_tree -> Branch("genj_eta",  genjs.eta, "genj_eta[genj_size]/F");
      }

   }
   if(saveVertices_){
      results_tree -> Branch("v_size",   &vtxs.size, "v_size/I");
      results_tree -> Branch("v_isFake",  vtxs.isFake,   "v_isFake[v_size]/O");
      results_tree -> Branch("v_ndof",    vtxs.ndof,    "v_ndof[v_size]/F");
      results_tree -> Branch("v_x",       vtxs.x,   "v_x[v_size]/F");
      results_tree -> Branch("v_y",       vtxs.y,   "v_y[v_size]/F");
      results_tree -> Branch("v_z",       vtxs.z,   "v_z[v_size]/F");
      results_tree -> Branch("v_Rho",     vtxs.Rho,   "v_Rho[v_size]/F");
   }
}


void
Mtuple::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup){
   bool changed(true);
   if (hltConfig_.init(iRun, iSetup, "HLT", changed)) {
      if (changed) {
         unsigned int prescales = hltConfig_.prescaleSize();
         std::cout << "New run: trigger menu has changed..."<< std::endl;
         std::vector<std::string> alltriggers=hltConfig_.triggerNames();

         TriggerPathVersioned_.clear();
         TriggerPathFilter_.clear();

         for(int i=0; i<nTriggerPaths_; ++i){
            bool triggerFound=false;
            for(std::vector<std::string>::const_iterator it=alltriggers.begin(); it!=alltriggers.end() && !triggerFound; ++it){
               TString fullname=*it;
               if(fullname.BeginsWith(TriggerPath_[i])){
                  std::cout << fullname << ":    ";
                  string stname = string(fullname);
                  TriggerPathVersioned_.push_back(stname);
                  for(unsigned int j=0; j<prescales; ++j) std::cout <<  hltConfig_.prescaleValue(j, stname) << " ";
                  std::cout << std::endl;

                  std::vector<std::string> trigObjects = hltConfig_.moduleLabels(stname);
                  int size = trigObjects.size();
                  edm::InputTag inputTag(trigObjects[size-2], "", "HLT");
                  TriggerPathFilter_.push_back(inputTag);
                  triggerFound=true;
               }
            }
            if(!triggerFound){
               TriggerPathVersioned_.push_back(TriggerPath_[i]);	
               TriggerPathFilter_.push_back(edm::InputTag("","",""));	
               std::cout << TriggerPath_[i] << "*   Was not found in this menu" << std::endl;
            }
         }
      }
   } else {
      LogError("Mtuple: ") << " HLT config extraction failure with process name HLT.";
   }
}

void 
Mtuple::endJob() {
   OutFile__file -> Write();
   OutFile__file -> Close();
}

DEFINE_FWK_MODULE(Mtuple);
