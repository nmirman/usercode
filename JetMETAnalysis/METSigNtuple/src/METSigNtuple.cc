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
// $Id: METSigNtuple.cc,v 1.7 2013/06/30 22:04:13 nmirman Exp $
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

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "CommonTools/ParticleFlow/test/PFIsoReaderDemo.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

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

bool compareMuPT(const reco::PFCandidate &mu1, const reco::PFCandidate &mu2){
   return (mu1.pt() > mu2.pt());
}
bool compareElecPT(const reco::PFCandidate &elec1, const reco::PFCandidate &elec2){
   return (elec1.pt() > elec2.pt());
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

      Bool_t      saveJetInfo_;

      std::string selectionChannel_;

      edm::InputTag muonTag_;
      edm::InputTag electronTag_;

      edm::InputTag conversionsInputTag_;
      edm::InputTag rhoIsoInputTag;
      std::vector<edm::InputTag> isoValInputTags_;

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

      Bool_t    saveBTags_;

      edm::InputTag  verticesTag_;

      edm::InputTag metSigTag_;
      edm::InputTag metSigMatrix00_;
      edm::InputTag metSigMatrix01_;
      edm::InputTag metSigMatrix10_;
      edm::InputTag metSigMatrix11_;

      HLTConfigProvider hltConfig_;
      JetResolution *ptRes_;
      JetResolution *phiRes_;

      TFile    *OutFile__file;
      TTree    *results_tree;  

      RecoMuon    mu;
      RecoElectron    elec;
      METs     mets;
      float    genmet_et, genmet_phi, genmet_sumEt;
      Long64_t    run, event, lumi;
      PFJets      pfjs;
      GenJets     genjs;
      GenW        genW;
      GenNu       gennu;
      GenMu       genmu;
      GenInfo     geninfo;
      BTags       btags;

      Vertices    vtxs;

      float metsig;
      float metsigmatrix00;
      float metsigmatrix01;
      float metsigmatrix10;
      float metsigmatrix11;

      float MyWeight;
      float T_nvertices;

      std::vector<std::string>      TriggerPath_; 
      std::vector<edm::InputTag>    TriggerPathFilter_;    
      std::vector<std::string>      TriggerPathVersioned_;
      Int_t  nTriggerPaths_;

      int metcount;
};

//
// constants, enums and typedefs
//

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;

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

   selectionChannel_ = iConfig.getUntrackedParameter<std::string>("selectionChannel");

   saveJetInfo_      = iConfig.getUntrackedParameter<Bool_t>("saveJetInfo");
   saveBTags_      = iConfig.getUntrackedParameter<Bool_t>("saveBTags");

   muonTag_    = iConfig.getUntrackedParameter<edm::InputTag>("muonTag");
   electronTag_ = iConfig.getUntrackedParameter<edm::InputTag>("electronTag");
   metsTag_    = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("metsTag");
   genMetTag_     = iConfig.getUntrackedParameter<edm::InputTag>("genmetTag");
   metsSize_      = metsTag_.size();

   conversionsInputTag_    = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
   rhoIsoInputTag          = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
   isoValInputTags_        = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

   genparticlesTag_  = iConfig.getUntrackedParameter<edm::InputTag>("genparticlesTag");
   pfcandidatesTag_  = iConfig.getUntrackedParameter<edm::InputTag>("pfcandidatesTag");

   pfjetsTag_    = iConfig.getUntrackedParameter<edm::InputTag>("pfjetsTag");
   pfjetCorrectorL1_  = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL1");
   pfjetCorrectorL123_ = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL123");

   genjetsTag_   = iConfig.getUntrackedParameter<edm::InputTag>("genjetsTag");
   genparticlesTag_   = iConfig.getUntrackedParameter<edm::InputTag>("genparticlesTag");

   verticesTag_  = iConfig.getUntrackedParameter<edm::InputTag>("verticesTag");

   metSigTag_ = iConfig.getUntrackedParameter<edm::InputTag>("metSig");
   metSigMatrix00_ = iConfig.getUntrackedParameter<edm::InputTag>("metSigMatrix00");
   metSigMatrix01_ = iConfig.getUntrackedParameter<edm::InputTag>("metSigMatrix01");
   metSigMatrix10_ = iConfig.getUntrackedParameter<edm::InputTag>("metSigMatrix10");
   metSigMatrix11_ = iConfig.getUntrackedParameter<edm::InputTag>("metSigMatrix11");

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
   // obtained with pileupCalc.py (06.26.2013)
   12260.8,
   32850.4,
   92330.3,
   339464,
   618478,
   3.0497e+06,
   1.77215e+07,
   5.41421e+07,
   1.30521e+08,
   2.58981e+08,
   4.46344e+08,
   6.8564e+08,
   8.81642e+08,
   9.99085e+08,
   1.07862e+09,
   1.13797e+09,
   1.17211e+09,
   1.18207e+09,
   1.17701e+09,
   1.16108e+09,
   1.13609e+09,
   1.10481e+09,
   1.06807e+09,
   1.02107e+09,
   9.55582e+08,
   8.6706e+08,
   7.58729e+08,
   6.38851e+08,
   5.16436e+08,
   3.99862e+08,
   2.96257e+08,
   2.10055e+08,
   1.42404e+08,
   9.20546e+07,
   5.65387e+07,
   3.29089e+07,
   1.815e+07,
   9.51188e+06,
   4.76417e+06,
   2.29967e+06,
   1.08138e+06,
   501998,
   233744,
   111112,
   54826,
   28402.3,
   15490.1,
   8845.44,
   5236.34,
   3180.14,
   1964.06,
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

      Handle<GenEventInfoProduct> gi;
      iEvent.getByLabel("generator", gi);

      edm::Handle< GenRunInfoProduct > genInfoProduct;
      iEvent.getRun().getByLabel("generator", genInfoProduct );

      geninfo.pid   = (int)gi->signalProcessID();
      geninfo.pthat = gi->qScale();
      geninfo.weight   = gi->weight();
      geninfo.xsec  = genInfoProduct->crossSection();
      geninfo.eff   = genInfoProduct->filterEfficiency();
      geninfo.alphaQCD = gi->alphaQCD();
      geninfo.alphaQED = gi->alphaQED();

      if(gi->hasPDF()){
         geninfo.scalePDF = gi->pdf()->scalePDF;
         geninfo.id1        = gi->pdf()->id.first;
         geninfo.id2        = gi->pdf()->id.second;
         geninfo.x1      = gi->pdf()->x.first;
         geninfo.x2      = gi->pdf()->x.second;
         geninfo.xPDF1    = gi->pdf()->xPDF.first;
         geninfo.xPDF2    = gi->pdf()->xPDF.second;
      }

      Handle<std::vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByLabel("addPileupInfo", PupInfo);

      std::vector< float > PU2012_MC;
      std::vector< float > PU2012_Data;

      for( int i=0; i<60; i++) {
         PU2012_MC.push_back( PU2012_MCf[i] );
         PU2012_Data.push_back( PU2012_Dataf[i] );
      }
      edm::LumiReWeighting LumiWeights_( PU2012_MC, PU2012_Data);

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
   edm::Handle<std::vector<reco::PFCandidate> > muonsHandle;
   iEvent.getByLabel(muonTag_, muonsHandle);
   std::vector<reco::PFCandidate> muons = *muonsHandle;
   std::sort(muons.begin(), muons.end(), compareMuPT);
   mu.size = muons.size(); 

   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle);
   const reco::BeamSpot &beamSpot = *(beamSpotHandle.product());

   int i=0;
   for(std::vector<reco::PFCandidate>::const_iterator it=muons.begin(); it!=muons.end(); it++, i++){
      if( !(it->muonRef().isNonnull() and it->muonRef().isAvailable()) ) continue;
      reco::Muon muon = *(it->muonRef());
      reco::PFCandidate pfmuon = *it;

      mu.charge[i]   = pfmuon.charge();
      mu.pt[i]       = pfmuon.pt();
      mu.p[i]        = pfmuon.p();
      mu.e[i]        = pfmuon.energy();
      mu.phi[i]      = pfmuon.phi();
      mu.eta[i]      = pfmuon.eta();
      mu.px[i]       = pfmuon.px();
      mu.py[i]       = pfmuon.py();
      mu.pz[i]       = pfmuon.pz();

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

   // conversions
   edm::Handle<reco::ConversionCollection> conversions_h;
   iEvent.getByLabel(conversionsInputTag_, conversions_h);

   // iso deposits
   IsoDepositVals isoVals(isoValInputTags_.size());
   for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
      iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
   }

   // vertices
   edm::Handle<reco::VertexCollection> vtx_h;
   iEvent.getByLabel(verticesTag_, vtx_h);

   // rho for isolation
   edm::Handle<double> rhoIso_h;
   iEvent.getByLabel(rhoIsoInputTag, rhoIso_h);
   double rhoIso = *(rhoIso_h.product());

   //electrons
   edm::Handle<std::vector<reco::PFCandidate> > electronsHandle;
   iEvent.getByLabel(electronTag_, electronsHandle);
   std::vector<reco::PFCandidate> electrons = *electronsHandle;
   std::sort(electrons.begin(), electrons.end(), compareElecPT);
   elec.size = electrons.size(); 
   for(std::vector<reco::PFCandidate>::const_iterator it=electrons.begin(); it!=electrons.end(); it++, i++){
      //if( !(it->gsfElectronRef().isNonnull() and it->gsfElectronRef().isAvailable()) )
      //reco::GsfElectron gsfelectron = *(it->gsfElectronRef());

      reco::PFCandidate pfelectron = *it;

      elec.charge[i]   = pfelectron.charge();
      elec.pt[i]       = pfelectron.pt();
      elec.p[i]        = pfelectron.p();
      elec.e[i]        = pfelectron.energy();
      elec.phi[i]      = pfelectron.phi();
      elec.eta[i]      = pfelectron.eta();
      elec.px[i]       = pfelectron.px();
      elec.py[i]       = pfelectron.py();
      elec.pz[i]       = pfelectron.pz();

      reco::SuperCluster elec_sc = *(it->superClusterRef());
      if( it->superClusterRef().isNonnull() ){
         elec.supercluster_eta[i] = elec_sc.eta();
      }else{
         elec.supercluster_eta[i] = -99;
      }

      double iso_ch=-1, iso_em=-1, iso_nh=-1;
      bool veto=0, loose=0, medium=0, tight=0;

      if( it->gsfElectronRef().isNonnull() and it->gsfElectronRef().isAvailable() ){
         // get particle flow isolation
         iso_ch = (*(isoVals)[0])[it->gsfElectronRef()];
         iso_em = (*(isoVals)[1])[it->gsfElectronRef()];
         iso_nh = (*(isoVals)[2])[it->gsfElectronRef()];
         // working points
         veto = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, it->gsfElectronRef(),
               conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
         loose = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, it->gsfElectronRef(),
               conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
         medium = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, it->gsfElectronRef(),
               conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
         tight = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, it->gsfElectronRef(),
               conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
         // for 2011 WP70 trigger
         // bool trigwp70 = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERWP70, ele);
      }

      elec.IDveto[i] = veto;
      elec.IDloose[i] = loose;
      elec.IDmedium[i] = medium;
      elec.IDtight[i] = tight;

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
      mets.dxx[imet]   = metiter.getSignificanceMatrix()(0,0); 
      mets.dxy[imet]   = metiter.getSignificanceMatrix()(0,1); 
      mets.dyy[imet]   = metiter.getSignificanceMatrix()(1,1); 
      mets.sig[imet]   = metiter.significance();
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
   iEvent.getByLabel(edm::InputTag("recoPuJetMva","full53xDiscriminant"), puJetIdMVA);

   Handle<ValueMap<int> > puJetIdFlag;
   iEvent.getByLabel(edm::InputTag("recoPuJetMva","full53xId"), puJetIdFlag);

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

      double jpt  = ( jet->pt()*pfjs.l1l2l3[icand] > 10 ) ? jet->pt()*pfjs.l1l2l3[icand] : jet->pt();
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

      TF1 fPtResol = *(ptRes_ -> resolutionEtaPt(jeta,jpt));
      pfjs.jetres_par0[icand] = fPtResol.GetParameter(0);
      pfjs.jetres_par1[icand] = fPtResol.GetParameter(1);
      pfjs.jetres_par2[icand] = fPtResol.GetParameter(2);
      pfjs.jetres_par3[icand] = fPtResol.GetParameter(3);
      pfjs.jetres_par4[icand] = fPtResol.GetParameter(4);
      pfjs.jetres_par5[icand] = fPtResol.GetParameter(5);
      pfjs.jetres_par6[icand] = fPtResol.GetParameter(6);

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

         // get neutrino momentum sum
         double nu_px = 0;
         double nu_py = 0;
         std::vector<const reco::GenParticle*> genConstituents = jet->getGenConstituents();
         for( std::vector<const reco::GenParticle*>::const_iterator aCandidate
               = genConstituents.begin();
               aCandidate != genConstituents.end(); aCandidate++ ){
            int aId = abs( (*aCandidate)->pdgId() );
            if( aId == 11 or aId == 13 or aId == 15 ){
               nu_px += (*aCandidate)->px();
               nu_py += (*aCandidate)->py();
            }
         }
         double nu_pt = sqrt( nu_px*nu_px + nu_py*nu_py );

         genjs.nu_px[icand] = nu_px;
         genjs.nu_py[icand] = nu_py;
         genjs.nu_pt[icand] = nu_pt;

         icand++;
      }
      genjs.size = icand;

      int icand_nu = 0;
      int icand_mu = 0;
      Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByLabel(genparticlesTag_, genParticles);
      for(size_t i = 0; i < genParticles->size(); ++ i) {
         const reco::GenParticle & p = (*genParticles)[i];
         int id = p.pdgId();
         int st = p.status();  
         const reco::Candidate * mom = p.mother();
         double pt = p.pt(), eta = p.eta(), phi = p.phi(), energy = p.energy();
         double vx = p.vx(), vy = p.vy(), vz = p.vz();
         int charge = p.charge();
         int n = p.numberOfDaughters();
         if( abs(id) == 24 and st == 3 ){
            genW.id = id;
            genW.pt = pt;
            genW.eta = eta;
            genW.phi = phi;
            genW.energy = energy;
            for(int j = 0; j < n; ++ j) {
               const reco::Candidate * d = p.daughter( j );
               int dauId = d->pdgId();
               int aId = abs(d->pdgId());
               if( d->status() == 3 and (aId == 11 or aId == 13 or aId == 15) ){
                  genW.l_id = dauId;
                  genW.l_pt = d->pt();
                  genW.l_eta = d->eta();
                  genW.l_phi = d->phi();
                  genW.l_energy = d->energy();
               }
               if( d->status() == 3 and (aId == 12 or aId == 14 or aId == 16) ){
                  genW.nu_id = dauId;
                  genW.nu_pt = d->pt();
                  genW.nu_eta = d->eta();
                  genW.nu_phi = d->phi();
                  genW.nu_energy = d->energy();
               }
            }
         }

         // mising energy
         if( abs(id) == 12 or abs(id) == 14 or abs(id) == 16 ){
            gennu.id[icand_nu] = id;
            gennu.pt[icand_nu] = pt;
            gennu.eta[icand_nu] = eta;
            gennu.phi[icand_nu] = phi;
            gennu.energy[icand_nu] = energy;
            gennu.status[icand_nu] = st;
            icand_nu++;
         }
         // muons
         if( abs(id) == 13 ){
            genmu.pt[icand_mu] = pt;
            genmu.eta[icand_mu] = eta;
            genmu.phi[icand_mu] = phi;
            genmu.energy[icand_mu] = energy;
            genmu.status[icand_mu] = st;
            icand_mu++;
         }

      }
      gennu.size = icand_nu;
      genmu.size = icand_mu;

   }

   // met significance
   edm::Handle<double> metsigHandle;
   iEvent.getByLabel(metSigTag_, metsigHandle);
   metsig = *(metsigHandle.product());
   
   iEvent.getByLabel(metSigMatrix00_, metsigHandle);
   metsigmatrix00 = *(metsigHandle.product());
   iEvent.getByLabel(metSigMatrix01_, metsigHandle);
   metsigmatrix01 = *(metsigHandle.product());
   iEvent.getByLabel(metSigMatrix10_, metsigHandle);
   metsigmatrix10 = *(metsigHandle.product());
   iEvent.getByLabel(metSigMatrix11_, metsigHandle);
   metsigmatrix11 = *(metsigHandle.product());

   // b tagging
   for(Int_t i=0; i<saveBTags_; i++){
      btags.size = 0;
      btags.pt[i] = 0;
      btags.eta[i] = 0;
      btags.phi[i] = 0;
      btags.energy[i] = 0;
      btags.discriminator[i] = -100;
   }

   std::vector<double> goodbtag_eta;
   std::vector<double> goodbtag_phi;

   if( saveBTags_ ){
      edm::Handle<reco::JetTagCollection> bTagHandle;
      iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTagHandle);
      const reco::JetTagCollection & bTags = *(bTagHandle.product());

      btags.size = bTags.size();

      for (int i = 0; i < int(bTags.size()); ++i) {

         btags.pt[i] = bTags[i].first->pt();
         btags.eta[i] = bTags[i].first->eta();
         btags.phi[i] = bTags[i].first->phi();
         btags.energy[i] = bTags[i].first->energy();

         btags.discriminator[i] = bTags[i].second;

         if( btags.discriminator[i] > 0.898){
            goodbtag_eta.push_back( btags.eta[i] );
            goodbtag_phi.push_back( btags.phi[i] );
         }

      }
   }
   

   // selection criteria
   bool mu_tight = true;
   bool mu_iso = true;
   bool mu_checketa = true;
   bool mu_checkpt = true;
   bool mu_zpeak = true;

   bool elec_primary = false;
   bool elec_veto = false;

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
      if( !(fabs(mu.eta[i]) < 2.1) ) mu_checketa = false;
      // pt cut
      if( !(mu.pt[i] > 20) ) mu_checkpt = false;
   }
   if( mu.size == 2 ){
      // Z-mass window
      TLorentzVector mu1temp( mu.px[0], mu.py[0], mu.pz[0], mu.e[0] );
      TLorentzVector mu2temp( mu.px[1], mu.py[1], mu.pz[1], mu.e[1] );
      mu_zpeak = (mu1temp+mu2temp).M() > 60 and (mu1temp+mu2temp).M() < 120;
   }
   if( elec.size > 0 ){
      elec_primary = elec.IDmedium[0] and elec.pt[0] > 30 and fabs(elec.eta[0]) < 2.5;
   }
   for(int i=1; i < elec.size; i++){
      if( elec.IDloose[i] ) elec_veto = elec.IDloose[i] and elec.pt[i] > 20 and elec.eta[i] < 2.5;
   }
   // loose jetID
   int numjets_pt400 = 0;
   int numjets_pt200 = 0;
   int numjets_pt60 = 0;
   int numjets_pt50 = 0;
   int numjets_pt45 = 0;
   int numbtags = 0;
   for( int i=0; i < pfjs.size; i++ ){
      int nco = pfjs.nco[i];
      double nhf = pfjs.neutralHadronFraction[i];
      double nef = pfjs.neutralEmFraction[i];
      double chf = pfjs.chargedHadronFraction[i];
      double cef = pfjs.chargedEmFraction[i];
      int nch = pfjs.chargedHadronMultiplicity[i];

      bool passID = false;
      if( nco>1 and nhf<0.99 and nef<0.99 ){
         if( fabs(pfjs.eta[i]) < 2.4 ){
            if( cef<0.99 and chf>0 and nch>0 ) passID = true;
         }else{
            passID = true;
         }
      }

      double jetpt_corr = (pfjs.pt[i]*pfjs.l1l2l3[i] > 10) ? pfjs.pt[i]*pfjs.l1l2l3[i] : pfjs.pt[i];

      if( jetpt_corr > 400 ) numjets_pt400++;
      if( jetpt_corr > 200 ) numjets_pt200++;
      if( jetpt_corr > 60 ) numjets_pt60++;
      if( jetpt_corr > 50 ) numjets_pt50++;
      if( jetpt_corr > 45 ) numjets_pt45++;

      // match to b-tag
      for(int j=0; j < int(goodbtag_eta.size()); j++ ){
         double dphi = TVector2::Phi_mpi_pi( pfjs.phi[i] - goodbtag_phi[j] ); 
         double deta = pfjs.eta[i] - goodbtag_eta[j];
         double dRtemp = sqrt( deta*deta + dphi*dphi );
         if( dRtemp < 0.3 and passID and jetpt_corr > 45 ){
            numbtags++;
            continue;
         }
      }

   }

   bool Zmumu_selection = (mu.size == 2) and mu_tight and mu_iso and mu_checketa
      and mu_checkpt and mu_zpeak;
   bool Wenu_selection = elec_primary and !elec_veto;
   bool Ttbar_selection = (numjets_pt60 >= 4) and (numjets_pt50 >= 5) and (numjets_pt45 >= 6)
      and (numbtags > 1);
   bool Dijet_selection = (numjets_pt400 >= 1) and (numjets_pt200 >= 2);

   // FILL THE TREE
   if( selectionChannel_ == "Zmumu" ){
      if( Zmumu_selection ){
         results_tree -> Fill();
      }
   }
   else if( selectionChannel_ == "Wenu" ){
      if( Wenu_selection ){
         results_tree -> Fill();
      }
   }
   else if( selectionChannel_ == "Ttbar" ){
      if( Ttbar_selection ){
         results_tree->Fill();
      }
   }
   else if( selectionChannel_ == "Dijet" ){
      if( Dijet_selection ){
         results_tree->Fill();
      }
   }
   else{
      std::cout << "Error: Selection channel unknown." << std::endl;
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

      results_tree -> Branch("gi_pid", &geninfo.pid, "gi_pid/I");
      results_tree -> Branch("gi_pthat", &geninfo.pthat, "gi_pthat/F");
      results_tree -> Branch("gi_weight", &geninfo.weight, "gi_weight/F");
      results_tree -> Branch("gi_xsec", &geninfo.xsec, "gi_xsec/F");
      results_tree -> Branch("gi_eff", &geninfo.eff, "gi_eff/F");
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

   results_tree -> Branch("elec_size",   &elec.size, "elec_size/I");
   results_tree -> Branch("elec_charge",   elec.charge, "elec_charge[elec_size]/I");
   results_tree -> Branch("elec_pt", elec.pt, "elec_pt[elec_size]/F");
   results_tree -> Branch("elec_p", elec.p, "elec_p[elec_size]/F");
   results_tree -> Branch("elec_e", elec.e, "elec_e[elec_size]/F");
   results_tree -> Branch("elec_phi", elec.phi, "elec_phi[elec_size]/F");
   results_tree -> Branch("elec_eta", elec.eta, "elec_eta[elec_size]/F");
   results_tree -> Branch("elec_px", elec.px, "elec_px[elec_size]/F");
   results_tree -> Branch("elec_py", elec.py, "elec_py[elec_size]/F");
   results_tree -> Branch("elec_pz", elec.pz, "elec_pz[elec_size]/F");
   results_tree -> Branch("elec_supercluster_eta", elec.supercluster_eta, "elec_supercluster_eta[elec_size]/F");
   results_tree -> Branch("elec_IDveto", elec.IDveto, "elec_IDveto[elec_size]/O");
   results_tree -> Branch("elec_IDloose", elec.IDloose, "elec_IDloose[elec_size]/O");
   results_tree -> Branch("elec_IDmedium", elec.IDmedium, "elec_IDmedium[elec_size]/O");
   results_tree -> Branch("elec_IDtight", elec.IDtight, "elec_IDtight[elec_size]/O");

   if(runOnMC_){
      results_tree -> Branch("genW_id", &genW.id, "genW_id/I");
      results_tree -> Branch("genW_pt", &genW.pt, "genW_pt/F");
      results_tree -> Branch("genW_eta", &genW.eta, "genW_eta/F");
      results_tree -> Branch("genW_phi", &genW.phi, "genW_phi/F");
      results_tree -> Branch("genW_energy", &genW.energy, "genW_energy/F");
      results_tree -> Branch("genW_l_id", &genW.l_id, "genW_l_id/I");
      results_tree -> Branch("genW_l_pt", &genW.l_pt, "genW_l_pt/F");
      results_tree -> Branch("genW_l_eta", &genW.l_eta, "genW_l_eta/F");
      results_tree -> Branch("genW_l_phi", &genW.l_phi, "genW_l_phi/F");
      results_tree -> Branch("genW_l_energy", &genW.l_energy, "genW_l_energy/F");
      results_tree -> Branch("genW_nu_id", &genW.nu_id, "genW_nu_id/I");
      results_tree -> Branch("genW_nu_pt", &genW.nu_pt, "genW_nu_pt/F");
      results_tree -> Branch("genW_nu_eta", &genW.nu_eta, "genW_nu_eta/F");
      results_tree -> Branch("genW_nu_phi", &genW.nu_phi, "genW_nu_phi/F");
      results_tree -> Branch("genW_nu_energy", &genW.nu_energy, "genW_nu_energy/F");

      results_tree -> Branch("genNu_size", &gennu.size, "genNu_size/I");
      results_tree -> Branch("genNu_id", gennu.id, "genNu_id[genNu_size]/I");
      results_tree -> Branch("genNu_pt", gennu.pt, "genNu_pt[genNu_size]/F");
      results_tree -> Branch("genNu_eta", gennu.eta, "genNu_eta[genNu_size]/F");
      results_tree -> Branch("genNu_phi", gennu.phi, "genNu_phi[genNu_size]/F");
      results_tree -> Branch("genNu_energy", gennu.energy, "genNu_energy[genNu_size]/F");
      results_tree -> Branch("genNu_status", gennu.status, "genNu_status[genNu_size]/I");

      results_tree -> Branch("genMu_size", &genmu.size, "genMu_size/I");
      results_tree -> Branch("genMu_pt", genmu.pt, "genMu_pt[genMu_size]/F");
      results_tree -> Branch("genMu_eta", genmu.eta, "genMu_eta[genMu_size]/F");
      results_tree -> Branch("genMu_phi", genmu.phi, "genMu_phi[genMu_size]/F");
      results_tree -> Branch("genMu_energy", genmu.energy, "genMu_energy[genMu_size]/F");
      results_tree -> Branch("genMu_status", genmu.status, "genMu_status[genMu_size]/I");
   }

   results_tree -> Branch("met_size", &metsSize_, "met_size/I");
   results_tree -> Branch("met_pt", mets.pt, "met_pt[met_size]/F");
   results_tree -> Branch("met_px", mets.px, "met_px[met_size]/F");
   results_tree -> Branch("met_py", mets.py, "met_py[met_size]/F");
   results_tree -> Branch("met_pz", mets.pz, "met_pz[met_size]/F");
   results_tree -> Branch("met_phi", mets.phi, "met_phi[met_size]/F");
   results_tree -> Branch("met_sumEt", mets.sumEt, "met_sumEt[met_size]/F");
   results_tree -> Branch("met_dxx", mets.dxx, "met_dxx[met_size]/F");
   results_tree -> Branch("met_dxy", mets.dxy, "met_dxy[met_size]/F");
   results_tree -> Branch("met_dyy", mets.dyy, "met_dyy[met_size]/F");
   results_tree -> Branch("met_sig", mets.sig, "met_sig[met_size]/F");

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
   results_tree -> Branch("pfj_energy",    pfjs.energy,   "pfj_energy[pfj_size]/F");
   results_tree -> Branch("pfj_sigmapt",    pfjs.sigmapt,   "pfj_sigmapt[pfj_size]/F");
   results_tree -> Branch("pfj_sigmaphi",   pfjs.sigmaphi,  "pfj_sigmaphi[pfj_size]/F");
   results_tree -> Branch("pfj_jetres_par0",    pfjs.jetres_par0,   "pfj_jetres_par0[pfj_size]/F");
   results_tree -> Branch("pfj_jetres_par1",    pfjs.jetres_par1,   "pfj_jetres_par1[pfj_size]/F");
   results_tree -> Branch("pfj_jetres_par2",    pfjs.jetres_par2,   "pfj_jetres_par2[pfj_size]/F");
   results_tree -> Branch("pfj_jetres_par3",    pfjs.jetres_par3,   "pfj_jetres_par3[pfj_size]/F");
   results_tree -> Branch("pfj_jetres_par4",    pfjs.jetres_par4,   "pfj_jetres_par4[pfj_size]/F");
   results_tree -> Branch("pfj_jetres_par5",    pfjs.jetres_par5,   "pfj_jetres_par5[pfj_size]/F");
   results_tree -> Branch("pfj_jetres_par6",    pfjs.jetres_par6,   "pfj_jetres_par6[pfj_size]/F");

   if( saveJetInfo_ ){
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
   }

   if( saveBTags_ ){
      results_tree -> Branch("btags_size",   &btags.size,    "btags_size/I");
      results_tree -> Branch("btags_pt",     btags.pt,    "btags_pt[btags_size]/F");
      results_tree -> Branch("btags_phi",    btags.phi,   "btags_phi[btags_size]/F");
      results_tree -> Branch("btags_eta",    btags.eta,   "btags_eta[btags_size]/F");
      results_tree -> Branch("btags_energy",    btags.energy,   "btags_energy[btags_size]/F");
      results_tree -> Branch("btags_discriminator",    btags.discriminator,   "btags_discriminator[btags_size]/F");
   }

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
      results_tree -> Branch("genj_nu_px",  genjs.nu_px, "genj_nu_px[genj_size]/F");
      results_tree -> Branch("genj_nu_py",  genjs.nu_py, "genj_nu_py[genj_size]/F");
      results_tree -> Branch("genj_nu_pt",  genjs.nu_pt, "genj_nu_pt[genj_size]/F");
   }

   results_tree -> Branch("v_size",   &vtxs.size, "v_size/I");
   results_tree -> Branch("v_isFake",  vtxs.isFake,   "v_isFake[v_size]/O");
   results_tree -> Branch("v_ndof",    vtxs.ndof,    "v_ndof[v_size]/F");
   results_tree -> Branch("v_chi2",    vtxs.chi2,    "v_chi2[v_size]/F");
   results_tree -> Branch("v_x",       vtxs.x,   "v_x[v_size]/F");
   results_tree -> Branch("v_y",       vtxs.y,   "v_y[v_size]/F");
   results_tree -> Branch("v_z",       vtxs.z,   "v_z[v_size]/F");
   results_tree -> Branch("v_Rho",     vtxs.Rho,   "v_Rho[v_size]/F");
   
   results_tree -> Branch("metsig", &metsig, "metsig/F");
   results_tree -> Branch("metsigmatrix00", &metsigmatrix00, "metsigmatrix00/F");
   results_tree -> Branch("metsigmatrix01", &metsigmatrix01, "metsigmatrix01/F");
   results_tree -> Branch("metsigmatrix10", &metsigmatrix10, "metsigmatrix10/F");
   results_tree -> Branch("metsigmatrix11", &metsigmatrix11, "metsigmatrix11/F");

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
METSigNtuple::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
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
