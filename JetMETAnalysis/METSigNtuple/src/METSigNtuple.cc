// -*- C++ -*-
//
// Package:    METSigNtuple
// Class:      METSigNtuple
// 
/**\class METSigNtuple METSigNtuple.cc JetMETAnalysis/METSigNtuple/src/METSigNtuple.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  nathan mirman
//         Created:  Wed Mar  6 16:05:43 CST 2013
// $Id: METSigNtuple.cc,v 1.1 2013/03/20 18:22:28 nmirman Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include "JetMETAnalysis/METSigNtuple/interface/NtupleFormats.h"

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

//
// class declaration
//

class METSigNtuple : public edm::EDAnalyzer {
   public:
      explicit METSigNtuple(const edm::ParameterSet&);
      ~METSigNtuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      static const double PU2012_MCf[60];
      static const double PU2012_Dataf[60];

      Bool_t      runOnMC_;
      std::string OutputFileName_;

      edm::InputTag  muonTag_;

      std::vector<edm::InputTag> metsTag_;
      Int_t    metsSize_;
      edm::InputTag  genMetTag_;

      edm::InputTag  genparticlesTag_;
      edm::InputTag  pfcandidatesTag_;

      Int_t    saveJets_;
      edm::InputTag  pfjetsTag_;
      std::string pfjetCorrectorL1_;
      std::string pfjetCorrectorL123_;

      edm::InputTag genjetsTag_;

      edm::InputTag  verticesTag_;

      JetResolution *ptRes_;
      JetResolution *phiRes_;

      TFile    *OutFile__file;
      TTree    *results_tree;  

      RecoMuon    mu;
      METs     mets;
      float    genmet_et, genmet_phi, genmet_sumEt;
      Long64_t    run, event, lumi;
      PFJets      pfjs;
      GenJets     genjs;

      Vertices    vtxs;

      float MyWeight;
      float T_nvertices;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
METSigNtuple::METSigNtuple(const edm::ParameterSet& iConfig)

{
   using namespace std;
   runOnMC_      = iConfig.getUntrackedParameter<Bool_t>("runOnMC");    
   OutputFileName_  = iConfig.getUntrackedParameter<std::string>("output_file");

   muonTag_    = iConfig.getUntrackedParameter<edm::InputTag>("muonTag");
   metsTag_    = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("metsTag");
   genMetTag_     = iConfig.getUntrackedParameter<edm::InputTag>("genmetTag");
   metsSize_      = metsTag_.size();

   genparticlesTag_  = iConfig.getUntrackedParameter<edm::InputTag>("genparticlesTag");
   pfcandidatesTag_  = iConfig.getUntrackedParameter<edm::InputTag>("pfcandidatesTag");

   pfjetsTag_    = iConfig.getUntrackedParameter<edm::InputTag>("pfjetsTag");
   pfjetCorrectorL1_  = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL1");
   pfjetCorrectorL123_ = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL123");

   genjetsTag_   = iConfig.getUntrackedParameter<edm::InputTag>("genjetsTag");

   verticesTag_  = iConfig.getUntrackedParameter<edm::InputTag>("verticesTag");

   string alg  = iConfig.getParameter<std::string>("jetResAlgo");     
   string era  = iConfig.getParameter<std::string>("jetResEra");     

   string path = "CondFormats/JetMETObjects/data";
   string ptFileName  = path + "/" + era + "_PtResolution_" +alg+".txt";
   string phiFileName = path + "/" + era + "_PhiResolution_"+alg+".txt";

   edm::FileInPath fpt(ptFileName);
   edm::FileInPath fphi(phiFileName);

   ptRes_  = new JetResolution(fpt.fullPath().c_str(),false);
   phiRes_ = new JetResolution(fphi.fullPath().c_str(),false);

   saveJets_ = 10000;
}


METSigNtuple::~METSigNtuple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

const double METSigNtuple::PU2012_Dataf[60] = {
   // 'true' distribution for 2012 dataset
   // obtained with pileupCalc.py (03.13.2013)
   12261.2,
   32854.9,
   90669,
   337108,
   619232,
   3.04778e+06,
   1.75106e+07,
   5.16043e+07,
   1.21968e+08,
   2.46908e+08,
   4.349e+08,
   6.77315e+08,
   8.77848e+08,
   9.97702e+08,
   1.07728e+09,
   1.13602e+09,
   1.16955e+09,
   1.17906e+09,
   1.17358e+09,
   1.15686e+09,
   1.13041e+09,
   1.09686e+09,
   1.05725e+09,
   1.00743e+09,
   9.40232e+08,
   8.52041e+08,
   7.46247e+08,
   6.30151e+08,
   5.11342e+08,
   3.97319e+08,
   2.95141e+08,
   2.09602e+08,
   1.42219e+08,
   9.19719e+07,
   5.64979e+07,
   3.28876e+07,
   1.81389e+07,
   9.50631e+06,
   4.76147e+06,
   2.29842e+06,
   1.08083e+06,
   501765,
   233650,
   111076,
   54812.8,
   28397.7,
   15488.6,
   8844.95,
   5236.19,
   3180.09,
   1964.04,
   1225.15,
   767.779,
   481.279,
   300.644,
   186.558,
   114.687,
   69.6938,
   41.7929,
   24.6979
};

const double METSigNtuple::PU2012_MCf[60] = {
   // pileup distribution for Summer2012 MC
   // obtained at https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
   2.560E-06,
   5.239E-06,
   1.420E-05,
   5.005E-05,
   1.001E-04,
   2.705E-04,
   1.999E-03,
   6.097E-03,
   1.046E-02,
   1.383E-02,
   1.685E-02,
   2.055E-02,
   2.572E-02,
   3.262E-02,
   4.121E-02,
   4.977E-02,
   5.539E-02,
   5.725E-02,
   5.607E-02,
   5.312E-02,
   5.008E-02,
   4.763E-02,
   4.558E-02,
   4.363E-02,
   4.159E-02,
   3.933E-02,
   3.681E-02,
   3.406E-02,
   3.116E-02,
   2.818E-02,
   2.519E-02,
   2.226E-02,
   1.946E-02,
   1.682E-02,
   1.437E-02,
   1.215E-02,
   1.016E-02,
   8.400E-03,
   6.873E-03,
   5.564E-03,
   4.457E-03,
   3.533E-03,
   2.772E-03,
   2.154E-03,
   1.656E-03,
   1.261E-03,
   9.513E-04,
   7.107E-04,
   5.259E-04,
   3.856E-04,
   2.801E-04,
   2.017E-04,
   1.439E-04,
   1.017E-04,
   7.126E-05,
   4.948E-05,
   3.405E-05,
   2.322E-05,
   1.570E-05,
   5.005E-06
};  

// ------------ method called for each event  ------------
   void
METSigNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   run = iEvent.id().run();
   event = iEvent.id().event();
   lumi = iEvent.id().luminosityBlock();

   // pileup reweighting
   MyWeight = 1.0;
   if( runOnMC_ ){

      Handle<std::vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByLabel("addPileupInfo", PupInfo);

      std::vector< float > PU2012_MC;
      std::vector< float > PU2012_Data;

      for( int i=0; i<60; i++) {
         PU2012_MC.push_back( PU2012_MCf[i] );
         PU2012_Data.push_back( PU2012_Dataf[i] );
      }
      edm::LumiReWeighting LumiWeights_( PU2012_MC, PU2012_Data, 0);

      std::vector<PileupSummaryInfo>::const_iterator PVI;

      float Tnvtx = -1.0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

         int BX = PVI->getBunchCrossing();

         if(BX == 0) { 
            Tnvtx = PVI->getTrueNumInteractions(); 
            continue;
         }

         MyWeight = LumiWeights_.weight( Tnvtx );
         T_nvertices = Tnvtx;
      }

   }

   edm::Handle<edm::View<reco::Vertex> > vertices;
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
      vtxs.chi2[icand] = v->chi2();
      vtxs.x[icand] = v->x();
      vtxs.y[icand] = v->y();
      vtxs.z[icand] = v->z();
      vtxs.Rho[icand] = v->position().Rho();
      icand++;
   } 
   vtxs.size = icand;
   reco::Vertex primaryVertex = vertices->at(0);

   // muons
   edm::Handle<std::vector<reco::Muon> > muonsHandle;
   iEvent.getByLabel(muonTag_, muonsHandle);
   std::vector<reco::Muon> muons = *muonsHandle;
   std::sort(muons.begin(), muons.end(), comparePT);
   mu.size = muons.size(); 

   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle);

   int i=0;
   for(std::vector<reco::Muon>::const_iterator it=muons.begin(); it!=muons.end(); it++, i++){
      reco::Muon muon = *it;

      mu.charge[i]   = muon.charge();
      mu.pt[i]       = muon.pt();
      mu.p[i]        = muon.p();
      mu.e[i]        = muon.energy();
      mu.phi[i]      = muon.phi();
      mu.eta[i]      = muon.eta();
      mu.px[i]       = muon.px();
      mu.py[i]       = muon.py();
      mu.pz[i]       = muon.pz();

      mu.isGlobal[i]  = muon.isGlobalMuon(); 
      mu.isTracker[i] = muon.isTrackerMuon();
      mu.isPF[i]      = muon.isPFMuon();

      mu.glpt[i]=-1;
      if(muon.isGlobalMuon()){
         reco::TrackRef gm    = muon.globalTrack();
         mu.chi2[i]           = gm->normalizedChi2();
         mu.trackerHits[i]    = gm->hitPattern().numberOfValidTrackerHits();
         mu.pixelHits[i]      = gm->hitPattern().numberOfValidPixelHits();
         mu.muonHits[i]       = gm->hitPattern().numberOfValidMuonHits();
         mu.nMatches[i]       = muon.numberOfMatchedStations();
         mu.glpt[i]           = gm->pt();
         mu.dxy[i]            = gm->dxy(primaryVertex.position());
         mu.dz[i]             = gm->dz(primaryVertex.position());
      }
      mu.trpt[i]=-1;
      if(muon.isTrackerMuon()){
         reco::TrackRef it = muon.innerTrack();
         mu.numberOfValidHits[i]       = it->hitPattern().numberOfValidHits();
         mu.numberOfValidPixelHits[i]     = it->hitPattern().numberOfValidPixelHits();
         mu.numberOfValidStripHits[i]     = it->hitPattern().numberOfValidStripHits();
         mu.numberOfValidPixelBarrelHits[i]  = it->hitPattern().numberOfValidPixelBarrelHits();
         mu.pixelLayersWithMeasurement[i] = it->hitPattern().pixelLayersWithMeasurement();
         mu.pixelBarrelLayersWithMeasurement[i] = it->hitPattern().pixelBarrelLayersWithMeasurement();
         mu.pixelEndcapLayersWithMeasurement[i] = it->hitPattern().pixelEndcapLayersWithMeasurement();
         mu.numberOfValidTrackerLayers[i] = it->hitPattern().trackerLayersWithMeasurement();
         mu.trpt[i]           = it->pt();
      }
      // if best track is available, redefine dxy and dz
      reco::TrackRef bt = muon.muonBestTrack();
      if( bt.isNonnull() and bt.isAvailable() ){
         mu.dxy[i]            = bt->dxy(primaryVertex.position());
         mu.dz[i]             = bt->dz(primaryVertex.position());
      }

      mu.dr03TkSumPt[i]      = muon.isolationR03().sumPt;
      mu.dr03EcalRecHitSumEt[i] = muon.isolationR03().emEt;
      mu.dr03HcalTowerSumEt[i]  = muon.isolationR03().hadEt;

      mu.dr04chHad[i] = muon.pfIsolationR04().sumChargedHadronPt;
      mu.dr04neutHad[i] = muon.pfIsolationR04().sumNeutralHadronEt;
      mu.dr04photons[i] = muon.pfIsolationR04().sumPhotonEt;

   }

   // met
   for(Int_t imet=0; imet<metsSize_ && imet<NMETs;imet++){
      edm::Handle<edm::View<reco::MET> > metHandle;
      iEvent.getByLabel(metsTag_[imet], metHandle);
      reco::MET metiter = (*metHandle)[0];
      mets.pt[imet]  = metiter.pt();
      mets.phi[imet] = metiter.phi();
      mets.px[imet]  = metiter.px();
      mets.py[imet]  = metiter.py();
      mets.pz[imet]  = metiter.pz();
      mets.sumEt[imet]  = metiter.sumEt();
   }
   // gen met
   if(runOnMC_){
      edm::Handle<edm::View<reco::MET> > genmetHandle;
      iEvent.getByLabel(genMetTag_, genmetHandle);
      reco::MET genmet  = (*genmetHandle)[0];
      genmet_et    = genmet.pt();
      genmet_phi    = genmet.phi();
      genmet_sumEt   = genmet.sumEt();
   } 

   // pfjets
   for(Int_t i=0; i<saveJets_; i++){
      pfjs.nco[i]=0;
      pfjs.l1[i]=0;
      pfjs.l1l2l3[i]=0;
      pfjs.pt[i]=0;
      pfjs.phi[i]=0;
      pfjs.eta[i]=0;
      pfjs.energy[i]=0;
      pfjs.sigmapt[i]=0;
      pfjs.sigmaphi[i]=0;

      pfjs.neutralHadronFraction[i]=0;
      pfjs.neutralEmFraction[i]=0;
      pfjs.chargedHadronFraction[i]=0;
      pfjs.chargedHadronMultiplicity[i]=0;
      pfjs.chargedEmFraction[i]=0;
   }

   Handle<reco::PFJetCollection> inputUncorJets;
   iEvent.getByLabel( pfjetsTag_, inputUncorJets );

   // pileup jet id
   Handle<ValueMap<float> > puJetIdMVA;
   iEvent.getByLabel(edm::InputTag("recoPuJetMva","fullDiscriminant"), puJetIdMVA);

   Handle<ValueMap<int> > puJetIdFlag;
   iEvent.getByLabel(edm::InputTag("recoPuJetMva","fullId"), puJetIdFlag);

   // jet energy corrections
   const JetCorrector* corrector  = JetCorrector::getJetCorrector (pfjetCorrectorL1_, iSetup);
   const JetCorrector* corrector2 = JetCorrector::getJetCorrector (pfjetCorrectorL123_, iSetup);

   icand = 0;
   for(reco::PFJetCollection::const_iterator jet = inputUncorJets->begin(); jet != inputUncorJets->end() && icand<saveJets_; ++jet) {
      int index = jet - inputUncorJets->begin();
      edm::RefToBase<reco::PFJet> jetRef(edm::Ref<reco::PFJetCollection>(inputUncorJets,index));
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

      TF1* fPtEta    = ptRes_ -> parameterEta("sigma",jeta);
      TF1* fPhiEta   = phiRes_-> parameterEta("sigma",jeta);
      pfjs.sigmapt[icand]     = fPtEta->Eval(jpt); 
      pfjs.sigmaphi[icand]    = fPhiEta->Eval(jpt);
      delete fPtEta;
      delete fPhiEta;

      pfjs.puid_mva[icand]  = (*puJetIdMVA)[jetRef];
      pfjs.puid_idflag[icand] = (*puJetIdFlag)[jetRef];
      pfjs.puid_passloose[icand] = PileupJetIdentifier::passJetId( pfjs.puid_idflag[icand],
            PileupJetIdentifier::kLoose );
      pfjs.puid_passmedium[icand] = PileupJetIdentifier::passJetId( pfjs.puid_idflag[icand],
            PileupJetIdentifier::kMedium );
      pfjs.puid_passtight[icand] = PileupJetIdentifier::passJetId( pfjs.puid_idflag[icand],
            PileupJetIdentifier::kTight );

      icand++;
   }
   pfjs.size = icand;
   // genjets
   if(runOnMC_){
      Handle<reco::GenJetCollection> genjets;
      iEvent.getByLabel( genjetsTag_, genjets );

      icand = 0;
      for(reco::GenJetCollection::const_iterator jet = genjets->begin(); jet != genjets->end(); ++jet) {

         genjs.pt[icand] = jet->pt();
         genjs.phi[icand] = jet->phi();
         genjs.eta[icand] = jet->eta();
         genjs.energy[icand] = jet->energy();

         genjs.emEnergy[icand] = jet->emEnergy();
         genjs.hadEnergy[icand] = jet->hadEnergy();
         genjs.invEnergy[icand] = jet->invisibleEnergy();
         genjs.auxEnergy[icand] = jet->auxiliaryEnergy();

         icand++;
      }
      genjs.size = icand;
   }


   // selection criteria
   bool mu_tight = true;
   bool mu_iso = true;
   bool mu_checketa = true;
   bool mu_checkpt = true;
   bool mu_zpeak = true;

   for(int i=0; i < mu.size; i++){
      // tight muon selection
      if( !(mu.isGlobal[i] and mu.isPF[i] and mu.chi2[i] < 10 and mu.muonHits[i] > 0
               and mu.nMatches[i] > 1 and mu.dxy[i] < 0.2 and mu.dz[i] < 0.5
               and mu.pixelHits[i] > 0 and mu.numberOfValidTrackerLayers[i] > 5) ){
         mu_tight = false;
      }
      // isolation
      if( !( (mu.dr04chHad[i]+mu.dr04neutHad[i]+mu.dr04photons[i])
               / mu.pt[i] < 0.12) ){
         mu_iso = false;
      }
      if( !(fabs(mu.eta[i]) < 2.4) ) mu_checketa = false;
      // pt cut
      if( !(mu.pt[i] > 20) ) mu_checkpt = false;
   }
   if( mu.size == 2 ){
      // Z-mass window
      TLorentzVector mu1temp( mu.px[0], mu.py[0], mu.pz[0], mu.e[0] );
      TLorentzVector mu2temp( mu.px[1], mu.py[1], mu.pz[1], mu.e[1] );
      mu_zpeak = (mu1temp+mu2temp).M() > 60 and (mu1temp+mu2temp).M() < 120;
   }
   bool Zmumu_selection = (mu.size == 2) and mu_tight and mu_iso and mu_checketa
      and mu_checkpt and mu_zpeak;

   // FILL THE TREE
   if( Zmumu_selection ){
      results_tree -> Fill();
   }
}


// ------------ method called once each job just before starting event loop  ------------
   void 
METSigNtuple::beginJob()
{
   OutFile__file  = new TFile( OutputFileName_.c_str(), "RECREATE" );

   results_tree = new TTree("events", "events");
   results_tree -> Branch("run", &run, "run/I");
   results_tree -> Branch("lumi", &lumi, "lumi/I");
   results_tree -> Branch("event", &event, "event/I");

   if(runOnMC_){
      results_tree -> Branch("puMyWeight", &MyWeight, "puMyWeight/F");
      results_tree -> Branch("puTnvtx", &T_nvertices, "puTnvtx/F");
   }
   results_tree -> Branch("mu_size",   &mu.size, "mu_size/I");
   results_tree -> Branch("mu_charge",   mu.charge, "mu_charge[mu_size]/I");
   results_tree -> Branch("mu_isGlobal", mu.isGlobal, "mu_isGlobal[mu_size]/O");
   results_tree -> Branch("mu_isTracker", mu.isTracker, "mu_isTracker[mu_size]/O");
   results_tree -> Branch("mu_isPF", mu.isPF, "mu_isPF[mu_size]/O");
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
   results_tree -> Branch("mu_dr04chHad", mu.dr04chHad, "mu_dr04chHad[mu_size]/F");
   results_tree -> Branch("mu_dr04neutHad", mu.dr04neutHad, "mu_dr04neutHad[mu_size]/F");
   results_tree -> Branch("mu_dr04photons", mu.dr04photons, "mu_dr04photons[mu_size]/F");

   results_tree -> Branch("met_size", &metsSize_, "met_size/I");
   results_tree -> Branch("met_pt", mets.pt, "met_pt[met_size]/F");
   results_tree -> Branch("met_px", mets.px, "met_px[met_size]/F");
   results_tree -> Branch("met_py", mets.py, "met_py[met_size]/F");
   results_tree -> Branch("met_pz", mets.pz, "met_pz[met_size]/F");
   results_tree -> Branch("met_phi", mets.phi, "met_phi[met_size]/F");
   results_tree -> Branch("met_sumEt", mets.sumEt, "met_sumEt[met_size]/F");

   if(runOnMC_){
      results_tree -> Branch("genmet_et", &genmet_et, "genmet_et/F");
      results_tree -> Branch("genmet_phi", &genmet_phi, "genmet_phi/F");
      results_tree -> Branch("genmet_sumEt", &genmet_sumEt, "genmet_sumEt/F");
   }

   results_tree -> Branch("pfj_size",   &pfjs.size, "pfj_size/I");
   results_tree -> Branch("pfj_l1",     pfjs.l1,  "pfj_l1[pfj_size]/F");
   results_tree -> Branch("pfj_l1l2l3",   pfjs.l1l2l3,  "pfj_l1l2l3[pfj_size]/F");
   results_tree -> Branch("pfj_pt",     pfjs.pt,    "pfj_pt[pfj_size]/F");
   results_tree -> Branch("pfj_phi",    pfjs.phi,   "pfj_phi[pfj_size]/F");
   results_tree -> Branch("pfj_eta",    pfjs.eta,   "pfj_eta[pfj_size]/F");
   results_tree -> Branch("pfj_energy",    pfjs.energy,   "pfj_eta[pfj_size]/F");
   results_tree -> Branch("pfj_sigmapt",    pfjs.sigmapt,   "pfj_sigmapt[pfj_size]/F");
   results_tree -> Branch("pfj_sigmaphi",   pfjs.sigmaphi,  "pfj_sigmaphi[pfj_size]/F");

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

   results_tree -> Branch("pfj_puid_mva", pfjs.puid_mva, "pfj_puid_mva[pfj_size]/F");
   results_tree -> Branch("pfj_puid_idflag", pfjs.puid_idflag, "pfj_puid_idflag[pfj_size]/I");
   results_tree -> Branch("pfj_puid_passloose", pfjs.puid_passloose,
         "pfj_puid_passloose[pfj_size]/O");
   results_tree -> Branch("pfj_puid_passmedium", pfjs.puid_passmedium,
         "pfj_puid_passmedium[pfj_size]/O");
   results_tree -> Branch("pfj_puid_passtight", pfjs.puid_passtight,
         "pfj_puid_passtight[pfj_size]/O");

   if(runOnMC_){
      results_tree -> Branch("genj_size",  &genjs.size, "genj_size/I");
      results_tree -> Branch("genj_pt",  genjs.pt, "genj_pt[genj_size]/F");
      results_tree -> Branch("genj_phi",  genjs.phi, "genj_phi[genj_size]/F");
      results_tree -> Branch("genj_eta",  genjs.eta, "genj_eta[genj_size]/F");
      results_tree -> Branch("genj_energy",  genjs.energy, "genj_energy[genj_size]/F");
      results_tree -> Branch("genj_emEnergy",  genjs.emEnergy, "genj_emEnergy[genj_size]/F");
      results_tree -> Branch("genj_hadEnergy",  genjs.hadEnergy, "genj_hadEnergy[genj_size]/F");
      results_tree -> Branch("genj_invEnergy",  genjs.invEnergy, "genj_invEnergy[genj_size]/F");
      results_tree -> Branch("genj_auxEnergy",  genjs.auxEnergy, "genj_auxEnergy[genj_size]/F");
   }

   results_tree -> Branch("v_size",   &vtxs.size, "v_size/I");
   results_tree -> Branch("v_isFake",  vtxs.isFake,   "v_isFake[v_size]/O");
   results_tree -> Branch("v_ndof",    vtxs.ndof,    "v_ndof[v_size]/F");
   results_tree -> Branch("v_chi2",    vtxs.chi2,    "v_chi2[v_size]/F");
   results_tree -> Branch("v_x",       vtxs.x,   "v_x[v_size]/F");
   results_tree -> Branch("v_y",       vtxs.y,   "v_y[v_size]/F");
   results_tree -> Branch("v_z",       vtxs.z,   "v_z[v_size]/F");
   results_tree -> Branch("v_Rho",     vtxs.Rho,   "v_Rho[v_size]/F");

}

// ------------ method called once each job just after ending the event loop  ------------
   void 
METSigNtuple::endJob() 
{
   OutFile__file -> Write();
   OutFile__file -> Close();
}

// ------------ method called when starting to processes a run  ------------
   void 
METSigNtuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
   void 
METSigNtuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
   void 
METSigNtuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
   void 
METSigNtuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METSigNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(METSigNtuple);
