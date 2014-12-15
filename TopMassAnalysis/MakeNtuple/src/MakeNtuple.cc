// -*- C++ -*-
//
// Package:    MakeNtuple
// Class:      MakeNtuple
// 
/**\class MakeNtuple MakeNtuple.cc BsmMasses/MakeNtuple/src/MakeNtuple.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Nicholas Eggert
//         Created:  Tue Sep 20 15:51:21 CDT 2011
// $Id: MakeNtuple.cc,v 1.2 2013/07/18 20:06:43 xguo Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "TLorentzVector.h"
#include "TMatrix.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

#include <string>

//
// class declaration
//

class JetSorter {
   public:
      JetSorter(std::string bTagName) {bTagName_ = bTagName; }
      bool operator() (pat::Jet jet1, pat::Jet jet2) {
         // return true if jet1 should come before jet2

         // if either of the jets has a b-tag, the one with the larger tag goes first
         // this means that if one has no tag (negative discriminator), the tagged one goes first
         if ((jet1.bDiscriminator(bTagName_) > 0.) || jet2.bDiscriminator(bTagName_) > 0.) {
            return (jet1.bDiscriminator(bTagName_) > jet2.bDiscriminator(bTagName_));
         }
         else { // if neither jet is tagged, the one with the higher pt goes first
            return (jet1.pt() > jet2.pt());
         }
      }
   private:
      std::string bTagName_;
};

struct event_id {
   // struct to hold event identifiers so we can check for uniqueness
   int run;
   int lumi;
   int event;
};

struct compare_event_id {
   bool operator() (event_id id1, event_id id2) {
      // sorting class for event_id. Needed to use the set object with these
      if (id1.run != id2.run) return (id1.run > id2.run);
      else if (id1.lumi != id2.lumi) return (id1.lumi > id2.lumi);
      else return (id1.event > id2.event);
   }
};




class MakeNtuple : public edm::EDAnalyzer {
   public:
      explicit MakeNtuple(const edm::ParameterSet&);
      ~MakeNtuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      bool passOfflineSelection( const TLorentzVector met,
               const TLorentzVector jet1, const TLorentzVector jet2,
               const TLorentzVector lep1, const TLorentzVector lep2 );
      void getSmearingFactors( const edm::Handle<edm::View<pat::Jet> >& jets,
            const TLorentzVector jet1, const TLorentzVector jet2,
            double *jet1_rE, double *jet2_rE, double *met_dx, double *met_dy, double *met_dz, double *met_dE );
      void smearJetsMET(
            const TLorentzVector jet1uns, const TLorentzVector jet2uns, const TLorentzVector metuns,
            TLorentzVector& jet1, TLorentzVector& jet2, TLorentzVector& met,
            double jet1_pt, double jet2_pt, double met_dx, double met_dy, double met_dz, double met_dE );
      void calculateSystematics( const edm::Event&, const TLorentzVector met, const TLorentzVector jet1, const TLorentzVector jet2,
            const TLorentzVector lep1, const TLorentzVector lep2);

      // ----------member data ---------------------------
      double jetScale_, leptonScale_, jetResScale_;
      std::string bTagger_;
      double bTagCut_, negTagCut_;
      bool doMETCut_, runOnMC_, runTtbar_;
      std::string outFileName_;
      JetSorter js;
      event_id this_event_id;
      std::set<event_id, compare_event_id> event_ids;

      edm::InputTag               muonTag_;
      edm::InputTag               electronTag_;
      edm::InputTag               jetTag_;
      edm::InputTag               metTag_;
      edm::InputTag               pmetTag_;
      edm::InputTag               genParticleTag_;

      int randSeed_;

      std::vector<pat::Jet> taggedJets_;
      std::vector<pat::Jet> negTaggedJets_;
      std::vector<const reco::RecoCandidate*> goodLeptons_;

      // This will be a pointer to the tree we want to use for this event
      std::map<std::string, TTree*> trees;
      TTree *tree;
      TFile *file;

      JetResolution *ptResol_;
      JetResolution *phiResol_;
      JetResolution *etaResol_;

      const reco::Candidate *t, *tb, *b, *bb, *Wp, *Wm, *lp, *lm, *n, *nb;

      TLorentzVector jet1FourVector, jet2FourVector, lep1FourVector, lep2FourVector;
      TLorentzVector uncorrectedJet1FourVector, uncorrectedJet2FourVector;
      TLorentzVector generatedJet1FourVector, generatedJet2FourVector, metFourVector, generatedMetFourVector;
      TLorentzVector bGEN, bbarGEN, lpGEN, lmGEN, nGEN, nbarGEN;
      int lpPdgIdGEN, lmPdgIdGEN;
      int nPdgIdGEN, nbPdgIdGEN;

      double jet1jesuncertainty, jet2jesuncertainty;
      double jet1PtResolution, jet1PhiResolution, jet1EtaResolution;
      double jet2PtResolution, jet2PhiResolution, jet2EtaResolution;
      double uncorrectedJet1PtResolution, uncorrectedJet1PhiResolution, uncorrectedJet1EtaResolution;
      double uncorrectedJet2PtResolution, uncorrectedJet2PhiResolution, uncorrectedJet2EtaResolution;
      double PDG1, PDG2;
      int jet1GenId, jet2GenId;
      int jet1ParentIdGEN, jet2ParentIdGEN;
      int vertices, numBJets;
      int run, lumi, event;
      int jetcount;
      double jet1Vz, jet2Vz;
      double jet1bdisc, jet2bdisc;

      static const double PU2012_MCf[60];
      static const double PU2012_Dataf[60];
      edm::LumiReWeighting LumiWeights_;
      reweight::PoissonMeanShifter PShiftUp_;
      reweight::PoissonMeanShifter PShiftDown_;
      float MyWeight;
      float TotalWeight_plus;
      float TotalWeight_minus;
      float T_nvertices;

      int geninfo_pid;
      float geninfo_pthat;
      float geninfo_weight;
      float geninfo_xsec;
      float geninfo_eff;
      float geninfo_alphaQCD;
      float geninfo_alphaQED;

      std::vector<std::string> systnames;
      std::vector<std::string> jsystnames;
      std::map<std::string, TLorentzVector> metsyst, jet1syst, jet2syst, lep1syst, lep2syst;

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
MakeNtuple::MakeNtuple(const edm::ParameterSet& iConfig) :
   jetScale_ (iConfig.getParameter<double>("jetScale")),
   leptonScale_ (iConfig.getParameter<double>("leptonScale")),
   jetResScale_ (iConfig.getParameter<double>("jetResScale")),
   bTagger_ (iConfig.getParameter<std::string>("bTagger")),
   bTagCut_ (iConfig.getParameter<double>("bTagCut")),
   negTagCut_ (iConfig.getParameter<double>("negTagCut")),
   doMETCut_ (iConfig.getParameter<bool>("doMETCut")),
   runOnMC_ (iConfig.getParameter<bool>("runOnMC")),
   runTtbar_ (iConfig.getParameter<bool>("runTtbar")),
   outFileName_ (iConfig.getParameter<std::string>("outFileName")),
   js (JetSorter(bTagger_)),
   // tree (NULL);
   muonTag_         (iConfig.getParameter<edm::InputTag>("muonSrc") ),
   electronTag_     (iConfig.getParameter<edm::InputTag>("electronSrc") ),
   jetTag_          (iConfig.getParameter<edm::InputTag>("jetSrc") ),
   metTag_          (iConfig.getParameter<edm::InputTag>("metSrc") ),
   pmetTag_          (iConfig.getParameter<edm::InputTag>("pmetSrc") ),
   genParticleTag_         (iConfig.getParameter<edm::InputTag>("genParticleSrc") ),
   randSeed_        (iConfig.getParameter<int>("randSeed"))



{
   //now do what ever initialization is needed
   std::string PATH("CondFormats/JetMETObjects/data/");
   edm::FileInPath ptFileName(PATH+"Spring10_PtResolution_AK5PF.txt");
   edm::FileInPath phiFileName(PATH+"Spring10_PhiResolution_AK5PF.txt");
   edm::FileInPath etaFileName(PATH+"Spring10_EtaResolution_AK5PF.txt");

   ptResol_  = new JetResolution(ptFileName.fullPath().c_str(),true);
   phiResol_ = new JetResolution(phiFileName.fullPath().c_str(),true);
   etaResol_ = new JetResolution(etaFileName.fullPath().c_str(),true);

}


MakeNtuple::~MakeNtuple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

const double MakeNtuple::PU2012_Dataf[60] = {
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

const double MakeNtuple::PU2012_MCf[60] = {
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
MakeNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   this_event_id.run   = iEvent.eventAuxiliary().run();
   this_event_id.lumi  = iEvent.eventAuxiliary().luminosityBlock();
   this_event_id.event = iEvent.eventAuxiliary().event();

   // insert returns a pair
   // the second element is false if the id is already in the set
   if (!event_ids.insert(this_event_id).second) return;

   taggedJets_.clear();
   negTaggedJets_.clear();
   goodLeptons_.clear();

   // pile up reweighting
   if( runOnMC_ ){

      edm::Handle<GenEventInfoProduct> gi;
      iEvent.getByLabel("generator", gi);

      //edm::Handle< GenRunInfoProduct > genInfoProduct;
      //iEvent.getRun().getByLabel("generator", genInfoProduct );

      geninfo_pid   = (int)gi->signalProcessID();
      geninfo_pthat = gi->qScale();
      geninfo_weight   = gi->weight();
      //geninfo_xsec  = genInfoProduct->crossSection();
      //geninfo_eff   = genInfoProduct->filterEfficiency();
      geninfo_xsec = -1;
      geninfo_eff = -1;
      geninfo_alphaQCD = gi->alphaQCD();
      geninfo_alphaQED = gi->alphaQED();

      std::vector<PileupSummaryInfo>::const_iterator PVI;
      edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByLabel("addPileupInfo", PupInfo);

      float Tnvtx = -1.0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

         int BX = PVI->getBunchCrossing();

         if(BX == 0) { 
            Tnvtx = PVI->getTrueNumInteractions(); 
            continue;
         }

         MyWeight = LumiWeights_.weight( Tnvtx );
         T_nvertices = Tnvtx;

         TotalWeight_plus = MyWeight*PShiftUp_.ShiftWeight( Tnvtx );
         TotalWeight_minus = MyWeight*PShiftDown_.ShiftWeight( Tnvtx );
      }

   }

   // get number of vertices
   double primary_vertex_z = 0;
   edm::Handle<std::vector<reco::Vertex> > goodVertices;
   iEvent.getByLabel("goodOfflinePrimaryVertices", "", goodVertices);
   vertices = goodVertices->size();
   reco::Vertex primary_vertex = goodVertices->at(0);
   primary_vertex_z = primary_vertex.z();


   // Get jets that pass the b tag cut
   edm::Handle<edm::View<pat::Jet> > jets;
   iEvent.getByLabel(jetTag_,jets);
   jetcount = 0;
   for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet!=jets->end();++jet){
      if (jet->bDiscriminator(bTagger_) > bTagCut_) {
         taggedJets_.push_back(*jet);
         jetcount++;
      }
   }

   // Get jets that pass the negative b tag cut
   for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet!=jets->end();++jet){
      if (jet->bDiscriminator(bTagger_) < negTagCut_) {
         negTaggedJets_.push_back(*jet);
         jetcount++;
      }
   }

   // Now sort the jets by b-tag, then by pt for jets without b-tags
   // use the JetSorter object js, declared earlier
   std::sort(taggedJets_.begin(), taggedJets_.end(), js);

   pat::Jet jet1, jet2;
   if( taggedJets_.size() >= 2 ){
      jet1 = taggedJets_[0];
      jet2 = taggedJets_[1];
   }
   else return;

   jet1Vz = jet1.vz()-primary_vertex_z;
   jet2Vz = jet2.vz()-primary_vertex_z;
   jet1bdisc = jet1.bDiscriminator(bTagger_);
   jet2bdisc = jet2.bDiscriminator(bTagger_);

   pat::Jet jet1Uncor = jet1.correctedJet(0); // 0 = Uncorrected
   pat::Jet jet2Uncor = jet2.correctedJet(0);
   const reco::GenJet *jet1Gen = NULL;
   const reco::GenJet *jet2Gen = NULL;
   jet1GenId = 0;
   jet2GenId = 0;
   jet1ParentIdGEN = 0;
   jet2ParentIdGEN = 0;
   double t_pt=0, tb_pt=0;
   if (runTtbar_) {
      if (jet1.genJet()){
         jet1GenId = jet1.partonFlavour();
         jet1Gen = jet1.genJet();
      }
      if(jet1.genParton()){
         jet1ParentIdGEN = (jet1.genParton())->mother()->pdgId();
      }
      if (jet2.genJet()) {
         jet2GenId = jet2.partonFlavour();
         jet2Gen = jet2.genJet();
      }
      if(jet2.genParton()){
         jet2ParentIdGEN = (jet2.genParton())->mother()->pdgId();
      }

      // get gen particles
      edm::Handle<edm::View<reco::GenParticle> > genParticles;
      iEvent.getByLabel(genParticleTag_,genParticles);
      for(edm::View<reco::GenParticle>::const_iterator p = genParticles->begin(); p!=genParticles->end();++p){
         if( fabs(p->pdgId()) == 6 and (p->daughter(0)) and (p->daughter(1)) ){
            if (p->pdgId()==6 && (p->daughter(0)->pdgId()==24 || p->daughter(1)->pdgId() == 24)) {
               t = &*p;
            }
            else if (p->pdgId()==-6 && (p->daughter(0)->pdgId()==-24 || p->daughter(1)->pdgId() == -24)){
               tb = &*p;
            }
         }
      }
      if (t->daughter(0)->pdgId() == 5 ) {
         b = (t->daughter(0));
         Wp = (t->daughter(1));
      }
      else {
         b = (t->daughter(1));
         Wp = (t->daughter(0));
      }
      if (tb->daughter(0)->pdgId() == -5) {
         bb = (tb->daughter(0));
         Wm = (tb->daughter(1));
      }
      else {
         bb = (tb->daughter(1));
         Wm = (tb->daughter(0));
      }

      if (Wp->daughter(0) and (Wp->daughter(0)->pdgId() == -13 || Wp->daughter(0)->pdgId() == -11) ) {
         lp = (Wp->daughter(0));
         n = (Wp->daughter(1));
      }
      else if (Wp->daughter(1) and (Wp->daughter(1)->pdgId() == -13 || Wp->daughter(1)->pdgId() == -11)) {
         lp = (Wp->daughter(1));
         n = (Wp->daughter(0));
      }
      else {
         lp = (Wp->daughter(1));
         n = (Wp->daughter(0));
      }
      if (Wm->daughter(0)->pdgId() == 13 || Wm->daughter(0)->pdgId() == 11) {
         lm = (Wm->daughter(0));
         nb = (Wm->daughter(1));
      }
      else if (Wm->daughter(1)->pdgId() == 13 || Wm->daughter(1)->pdgId() == 11) {
         lm = (Wm->daughter(1));
         nb = (Wm->daughter(0));
      }
      else {
         lm = (Wm->daughter(1));
         nb = (Wm->daughter(0));
      }

      bGEN.SetPxPyPzE(b->px(), b->py(), b->pz(), b->energy());
      bbarGEN.SetPxPyPzE(bb->px(), bb->py(), bb->pz(), bb->energy());
      lpGEN.SetPxPyPzE(lp->px(), lp->py(), lp->pz(), lp->energy());
      lmGEN.SetPxPyPzE(lm->px(), lm->py(), lm->pz(), lm->energy());
      nGEN.SetPxPyPzE(n->px(), n->py(), n->pz(), n->energy());
      nbarGEN.SetPxPyPzE(nb->px(), nb->py(), nb->pz(), nb->energy());
      lpPdgIdGEN = lp->pdgId();
      lmPdgIdGEN = lm->pdgId();
      nPdgIdGEN = n->pdgId();
      nbPdgIdGEN = nb->pdgId();

      t_pt = t->pt();
      tb_pt = tb->pt();

   } // if runTtbar_
   else
   {
      bGEN.SetPxPyPzE(0.0,0.0,0.0,0.0);
      bbarGEN.SetPxPyPzE(0.0,0.0,0.0,0.0);
      lpGEN.SetPxPyPzE(0.0,0.0,0.0,0.0);
      lmGEN.SetPxPyPzE(0.0,0.0,0.0,0.0);
      nGEN.SetPxPyPzE(0.0,0.0,0.0,0.0);
      nbarGEN.SetPxPyPzE(0.0,0.0,0.0,0.0);
      lpPdgIdGEN = 0;
      lmPdgIdGEN = 0;
      nPdgIdGEN = 0;
      nbPdgIdGEN = 0;
   }

   // get jet resolutions
   jet1PtResolution  = ptResol_ ->resolutionEtaPt(jet1.eta(), jet1.pt())->GetParameter(2)*jetResScale_;
   jet1PhiResolution = phiResol_->resolutionEtaPt(jet1.eta(), jet1.pt())->GetParameter(2);
   jet1EtaResolution = etaResol_->resolutionEtaPt(jet1.eta(), jet1.pt())->GetParameter(2);
   jet2PtResolution  = ptResol_ ->resolutionEtaPt(jet2.eta(), jet2.pt())->GetParameter(2)*jetResScale_;
   jet2PhiResolution = phiResol_->resolutionEtaPt(jet2.eta(), jet2.pt())->GetParameter(2);
   jet2EtaResolution = etaResol_->resolutionEtaPt(jet2.eta(), jet2.pt())->GetParameter(2);

   uncorrectedJet1PtResolution  = ptResol_ ->resolutionEtaPt(jet1Uncor.eta(), jet1Uncor.pt())->GetParameter(2);
   uncorrectedJet1PhiResolution = phiResol_->resolutionEtaPt(jet1Uncor.eta(), jet1Uncor.pt())->GetParameter(2);
   uncorrectedJet1EtaResolution = etaResol_->resolutionEtaPt(jet1Uncor.eta(), jet1Uncor.pt())->GetParameter(2);
   uncorrectedJet2PtResolution  = ptResol_ ->resolutionEtaPt(jet2Uncor.eta(), jet2Uncor.pt())->GetParameter(2);
   uncorrectedJet2PhiResolution = phiResol_->resolutionEtaPt(jet2Uncor.eta(), jet2Uncor.pt())->GetParameter(2);
   uncorrectedJet2EtaResolution = etaResol_->resolutionEtaPt(jet2Uncor.eta(), jet2Uncor.pt())->GetParameter(2);

   // Now get leptons
   edm::Handle<edm::View<pat::Electron> > electrons;
   iEvent.getByLabel(electronTag_, electrons);
   for( edm::View<pat::Electron>::const_iterator electron = electrons->begin(); electron != electrons->end(); ++electron) {
      const reco::RecoCandidate* lepton = &*electron;
      goodLeptons_.push_back(lepton);
   }
   edm::Handle<edm::View<pat::Muon> > muons;
   iEvent.getByLabel(muonTag_, muons);
   for( edm::View<pat::Muon>::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
      const reco::RecoCandidate* lepton = &*muon;
      goodLeptons_.push_back(lepton);
   }

   if (goodLeptons_.size() < 2) return;
   const reco::RecoCandidate* lep1 = goodLeptons_[0];
   const reco::RecoCandidate* lep2 = goodLeptons_[1];

   // get met
   edm::Handle<edm::View<reco::PFMET> > mets;
   iEvent.getByLabel(metTag_, mets);
   reco::PFMET met = mets->front();

   // get pat met for genmet object
   edm::Handle<edm::View<pat::MET> > pmets;
   iEvent.getByLabel(pmetTag_, pmets);
   pat::MET pmet = pmets->front();

   // Now set the variables we haven't set yet.
   jet1FourVector.SetPxPyPzE(jet1.px(), jet1.py(), jet1.pz(), jet1.energy());
   jet2FourVector.SetPxPyPzE(jet2.px(), jet2.py(), jet2.pz(), jet2.energy());
   uncorrectedJet1FourVector.SetPxPyPzE(jet1Uncor.px(), jet1Uncor.py(), jet1Uncor.pz(), jet1Uncor.energy());
   uncorrectedJet2FourVector.SetPxPyPzE(jet2Uncor.px(), jet2Uncor.py(), jet2Uncor.pz(), jet2Uncor.energy());
   if (jet1Gen!=NULL) generatedJet1FourVector.SetPxPyPzE(jet1Gen->px(), jet1Gen->py(), jet1Gen->pz(), jet1Gen->energy());
   else generatedJet1FourVector.SetPxPyPzE(0, 0, 0, 0);
   if (jet2Gen!=NULL) generatedJet2FourVector.SetPxPyPzE(jet2Gen->px(), jet2Gen->py(), jet2Gen->pz(), jet2Gen->energy());
   else generatedJet2FourVector.SetPxPyPzE(0, 0, 0, 0);

   // met
   metFourVector.SetPxPyPzE(met.px(), met.py(), met.pz(), met.energy());

   // jet energy scale uncertainty
   //std::string path = std::getenv("CMSSW_BASE");
   JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( "data/Summer13_V5_DATA_Uncertainty_AK5PFchs.txt" );
   double factor = 1;
   double corr = 0;
   jecUnc->setJetEta(jet1.eta());
   jecUnc->setJetPt(jet1.pt());
   jet1jesuncertainty = jecUnc->getUncertainty(true);
   if (jetScale_ != 0.) {
      // back out value of jec factor
      factor = jet1Uncor.pt() != 0 ? jet1.pt() / jet1Uncor.pt() : 1.0;
      jecUnc->setJetEta(jet1.eta());
      jecUnc->setJetPt(jet1.pt());
      corr = jetScale_*jecUnc->getUncertainty(true)/factor;
      // correct met first
      metFourVector -= jet1FourVector*corr;
      jet1FourVector *= 1+corr;
   }
   jecUnc->setJetEta(jet2.eta());
   jecUnc->setJetPt(jet2.pt());
   jet2jesuncertainty = jecUnc->getUncertainty(true);
   if (jetScale_ != 0.) {
      factor = jet2Uncor.pt() != 0 ? jet2.pt() / jet2Uncor.pt() : 1.0;
      jecUnc->setJetEta(jet2.eta());
      jecUnc->setJetPt(jet2.pt());
      corr = jetScale_*jecUnc->getUncertainty(true)/factor;
      // correct met first
      metFourVector -= jet2FourVector*corr;
      jet2FourVector *= 1+corr;
   }
   delete jecUnc;

   numBJets = taggedJets_.size();

   lep1FourVector.SetPxPyPzE(lep1->px(), lep1->py(), lep1->pz(), lep1->energy());
   lep2FourVector.SetPxPyPzE(lep2->px(), lep2->py(), lep2->pz(), lep2->energy());

   PDG1 = lep1->pdgId();
   PDG2 = lep2->pdgId();

   if (runTtbar_) {
      generatedMetFourVector.SetPxPyPzE(pmet.genMET()->px(), pmet.genMET()->py(), pmet.genMET()->pz(), pmet.genMET()->energy());
   }

   // save unsmeared jets and met
   TLorentzVector jet1Unsmeared = jet1FourVector;
   TLorentzVector jet2Unsmeared = jet2FourVector;
   TLorentzVector metUnsmeared = metFourVector;
   bool pass_selection = false;

   // jets+met smearing factors
   // with JER systematic variations
   double met_dx [3] = {0};
   double met_dy [3] = {0};
   double met_dz [3] = {0};
   double met_dE [3] = {0};
   double jet1_rE [3] = {1.0, 1.0, 1.0};
   double jet2_rE [3] = {1.0, 1.0, 1.0};
   if( runOnMC_ ) getSmearingFactors( jets, jet1FourVector, jet2FourVector, jet1_rE, jet2_rE,
         met_dx, met_dy, met_dz, met_dE );

   // 
   // use original objects to calculate JES systematics (should be accessed with const)
   //
   calculateSystematics( iEvent, metFourVector, jet1FourVector, jet2FourVector, lep1FourVector, lep2FourVector );

   //
   // fill tree for central sample (no systematic variations applied)
   //
   smearJetsMET( jet1Unsmeared, jet2Unsmeared, metUnsmeared,
         jet1FourVector, jet2FourVector, metFourVector,
         jet1_rE[0], jet2_rE[0], met_dx[0], met_dy[0], met_dz[0], met_dE[0] ); 
   pass_selection = passOfflineSelection( metFourVector,
         jet1FourVector, jet2FourVector, lep1FourVector, lep2FourVector );
   if( pass_selection ) trees["Central"]->Fill();

   // jet energy resolution systematics
   smearJetsMET( jet1Unsmeared, jet2Unsmeared, metUnsmeared,
         jet1FourVector, jet2FourVector, metFourVector,
         jet1_rE[1], jet2_rE[1], met_dx[1], met_dy[1], met_dz[1], met_dE[1] ); 
   pass_selection = passOfflineSelection( metFourVector,
         jet1FourVector, jet2FourVector, lep1FourVector, lep2FourVector );
   if( pass_selection ) trees["JetEnergyResolutionUP"]->Fill();

   smearJetsMET( jet1Unsmeared, jet2Unsmeared, metUnsmeared,
         jet1FourVector, jet2FourVector, metFourVector,
         jet1_rE[2], jet2_rE[2], met_dx[2], met_dy[2], met_dz[2], met_dE[2] ); 
   pass_selection = passOfflineSelection( metFourVector,
         jet1FourVector, jet2FourVector, lep1FourVector, lep2FourVector );
   if( pass_selection ) trees["JetEnergyResolutionDN"]->Fill();

   //
   // fill trees for systematics samples
   //
   double weight_temp = MyWeight;
   for( std::map<std::string, TLorentzVector>::iterator syst = metsyst.begin(); syst != metsyst.end(); syst++ ){
      std::string name = syst->first;

      MyWeight  = weight_temp;

      // get syst varied quantities
      metFourVector = metsyst[name];
      jet1FourVector = jet1syst[name];
      jet2FourVector = jet2syst[name];
      lep1FourVector = lep1syst[name];
      lep2FourVector = lep2syst[name];

      jet1Unsmeared = jet1FourVector;
      jet2Unsmeared = jet2FourVector;
      metUnsmeared = metFourVector;

      // pile up systematic
      if( name == "PileUpUP" ) MyWeight = TotalWeight_plus;
      if( name == "PileUpDN" ) MyWeight = TotalWeight_minus;

      // top pt reweighting
      double a = 0.148;
      double b = -0.00129;
      if( name.find("PtTopReweighting") != std::string::npos ) MyWeight *= sqrt( exp(a+b*t_pt)*exp(a+b*tb_pt) );

      // apply jet smearing
      smearJetsMET( jet1Unsmeared, jet2Unsmeared, metUnsmeared,
            jet1FourVector, jet2FourVector, metFourVector,
            jet1_rE[0], jet2_rE[0], met_dx[0], met_dy[0], met_dz[0], met_dE[0] ); 

      // fill syst tree
      pass_selection = passOfflineSelection( metFourVector,
            jet1FourVector, jet2FourVector, lep1FourVector, lep2FourVector );
      if( pass_selection ) trees[name]->Fill();
   }

}

   bool
MakeNtuple::passOfflineSelection( const TLorentzVector met,
      const TLorentzVector jet1, const TLorentzVector jet2,
      const TLorentzVector lep1, const TLorentzVector lep2 ){

   double met_pt = 40; // only for ee and mumu events
   double jet_pt = 30;
   double jet_eta = 2.5;
   double mu_pt = 20;
   double mu_eta = 2.4;
   double e_pt = 20;
   double e_eta = 2.5;
   // skipping dilepton mass criteria

   // ##### TEMPORARY SOLUTION ##### TODO
   
   met_pt = 45;
   jet_pt = 35;
   mu_pt = 21;
   e_pt = 21;
   
   // ##### TEMPORARY SOLUTION ##### TODO

   bool met_ok = (abs(PDG1) != abs(PDG2)) or (met.Pt() > met_pt);
   bool jet1_ok = jet1.Pt() > jet_pt and jet1.Eta() < jet_eta;
   bool jet2_ok = jet2.Pt() > jet_pt and jet2.Eta() < jet_eta;
   bool lep1_ok = abs(PDG1) == 11 ? (lep1.Pt() > e_pt and lep1.Eta() < e_eta)
                                 : (lep1.Pt() > mu_pt and lep1.Eta() < mu_eta);
   bool lep2_ok = abs(PDG1) == 11 ? (lep2.Pt() > e_pt and lep2.Eta() < e_eta)
                                 : (lep2.Pt() > mu_pt and lep2.Eta() < mu_eta);

   return (met_ok and jet1_ok and jet2_ok and lep1_ok and lep2_ok);
}

   void
MakeNtuple::getSmearingFactors( const edm::Handle<edm::View<pat::Jet> >& jets,
      const TLorentzVector jet1, const TLorentzVector jet2,
      double *jet1_rE, double *jet2_rE, double *met_dx, double *met_dy, double *met_dz, double *met_dE ){

   for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet!=jets->end();++jet){
      for(int i=0; i < 3; i++){
         
         // get 8 TeV scale factors from
         // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
         double c = 1.0;
         if( i==0 ){ // central
            if( fabs(jet->eta()) < 0.5 ) c = 1.079;
            else if( fabs(jet->eta()) < 1.1 ) c = 1.099;
            else if( fabs(jet->eta()) < 1.7 ) c = 1.121;
            else if( fabs(jet->eta()) < 2.3 ) c = 1.208;
            else if( fabs(jet->eta()) < 2.8 ) c = 1.254;
            else if( fabs(jet->eta()) < 3.2 ) c = 1.395;
            else if( fabs(jet->eta()) < 5.0 ) c = 1.056;
         }else if( i==1 ){ // up
            if( fabs(jet->eta()) < 0.5 ) c = 1.053;
            else if( fabs(jet->eta()) < 1.1 ) c = 1.071;
            else if( fabs(jet->eta()) < 1.7 ) c = 1.192;
            else if( fabs(jet->eta()) < 2.3 ) c = 1.162;
            else if( fabs(jet->eta()) < 2.8 ) c = 1.192;
            else if( fabs(jet->eta()) < 3.2 ) c = 1.332;
            else if( fabs(jet->eta()) < 5.0 ) c = 0.865;
         }else if( i==2 ){ // down
            if( fabs(jet->eta()) < 0.5 ) c = 1.105;
            else if( fabs(jet->eta()) < 1.1 ) c = 1.127;
            else if( fabs(jet->eta()) < 1.7 ) c = 1.150;
            else if( fabs(jet->eta()) < 2.3 ) c = 1.254;
            else if( fabs(jet->eta()) < 2.8 ) c = 1.316;
            else if( fabs(jet->eta()) < 3.2 ) c = 1.458;
            else if( fabs(jet->eta()) < 5.0 ) c = 1.247;
         }

         // get the smeared pt by scaling the RECO-GEN pt difference
         const reco::GenJet *genjet = NULL;
         double smearedJetE = jet->energy();
         if(jet->genJet()){
            genjet = jet->genJet();
            double dE = jet->energy() - genjet->energy();
            smearedJetE = genjet->energy() + c*dE;
         }
         double rE = (jet->pt() > 10) ? (smearedJetE / jet->energy()) : 1.0;

         // smear MET
         met_dx[i] -= (rE-1)*jet->px();
         met_dy[i] -= (rE-1)*jet->py();
         met_dz[i] -= (rE-1)*jet->pz();
         met_dE[i] -= (rE-1)*jet->energy();

         if( jet->pt() == jet1.Pt() and jet->phi() == jet1.Phi() ) jet1_rE[i] = rE;
         if( jet->pt() == jet2.Pt() and jet->phi() == jet2.Phi() ) jet2_rE[i] = rE;

      }
   }

   return;
}

   void
MakeNtuple::smearJetsMET(
      const TLorentzVector jet1uns, const TLorentzVector jet2uns, const TLorentzVector metuns,
      TLorentzVector& jet1, TLorentzVector& jet2, TLorentzVector &met,
      double jet1_rE, double jet2_rE, double met_dx, double met_dy, double met_dz, double met_dE
      ){

   jet1 *= jet1_rE;
   jet2 *= jet2_rE;

   met.SetPx( metuns.Px()+met_dx );
   met.SetPy( metuns.Py()+met_dy );
   met.SetPz( metuns.Pz()+met_dz );
   met.SetE( metuns.E()+met_dE );

   return;
}

   void
MakeNtuple::calculateSystematics( const edm::Event& iEvent, const TLorentzVector met,
      const TLorentzVector jet1, const TLorentzVector jet2,
      const TLorentzVector lep1, const TLorentzVector lep2)
{
   // Instantiate uncertainty sources
   std::vector<JetCorrectionUncertainty*> vsrc(jsystnames.size());
   std::string sgn [2] = {"DN","UP"};

   // jet energy scale
   for (unsigned int isrc = 0; isrc < jsystnames.size(); isrc++) {

      const char *name = jsystnames[isrc].c_str();
      JetCorrectorParameters *p = new JetCorrectorParameters("data/Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt", name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc[isrc] = unc;

      for(int i=0; i < 2; i++){
         metsyst[name+sgn[i]] = met;
         jet1syst[name+sgn[i]] = jet1;
         jet2syst[name+sgn[i]] = jet2;
         lep1syst[name+sgn[i]] = lep1;
         lep2syst[name+sgn[i]] = lep2;
      }

   } // for isrc

   // other systematics
   for (unsigned int i = 0; i < systnames.size(); i++) {
      const char *name = systnames[i].c_str();
      for(int i=0; i < 2; i++){
         metsyst[name+sgn[i]] = met;
         jet1syst[name+sgn[i]] = jet1;
         jet2syst[name+sgn[i]] = jet2;
         lep1syst[name+sgn[i]] = lep1;
         lep2syst[name+sgn[i]] = lep2;
      }
   }

   edm::Handle<edm::View<pat::Jet> > jets;
   iEvent.getByLabel(jetTag_,jets);

   // calculate unclustered energy contribution to the met
   TLorentzVector met_uncl = met + lep1 + lep2;
   for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet!=jets->end();++jet){
      TLorentzVector jetFourVector;
      jetFourVector.SetPxPyPzE(jet->px(), jet->py(), jet->pz(), jet->energy());
      met_uncl += jetFourVector;
   }

   // loop over jes uncertainty sources
   for (unsigned int isrc = 0; isrc < jsystnames.size(); isrc++) {
      JetCorrectionUncertainty *unc = vsrc[isrc];

      // up and down variations
      for(int i=0; i < 2; i++){
         std::string name = jsystnames[isrc]+sgn[i];

         // loop over all jets
         for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet!=jets->end();++jet){

            TLorentzVector jetFourVector;
            jetFourVector.SetPxPyPzE(jet->px(), jet->py(), jet->pz(), jet->energy());

            unc->setJetPt(jet->pt());
            unc->setJetEta(jet->eta());
            int sign = i == 0 ? -1 : 1;
            double var = sign*unc->getUncertainty(i); // 0 for down, 1 for up

            // get the uncorrected jet
            pat::Jet jetUncor = jet->correctedJet(0);

            // back out the correction factor
            double factor = jetUncor.pt() != 0 ? jet->pt() / jetUncor.pt() : 1.0;

            bool isJet1 = jetFourVector.Pt() == jet1syst[name].Pt();
            bool isJet2 = jetFourVector.Pt() == jet2syst[name].Pt();

            // correct met first
            metsyst[name] -= jetFourVector*(var/factor);
            jetFourVector *= 1+(var/factor);
            if( isJet1 ) jet1syst[name] = jetFourVector;
            if( isJet2 ) jet2syst[name] = jetFourVector;

         }
      } // up and down
   } // for isrc

   // loop over other systematics
   for(unsigned int isrc = 0; isrc < systnames.size(); isrc++){
      for(int i=0; i < 2; i++){
         std::string name = systnames[isrc]+sgn[i];
         int s = 2*i-1;

         if( name.find("MuonMomentumScale") != std::string::npos ){
            if( abs(PDG1) == 13 ) lep1syst[name] *= 1.000+s*0.002;
            if( abs(PDG2) == 13 ) lep2syst[name] *= 1.000+s*0.002;
         }
         if( name.find("ElectronEnergyScale") != std::string::npos ){
            double unc1=0, unc2=0;
            if( fabs(lep1syst[name].Eta()) < 1.48 ) unc1 = 0.006;
            else unc1 = 0.015;
            if( fabs(lep2syst[name].Eta()) < 1.48 ) unc2 = 0.006;
            else unc2 = 0.015;
            if( abs(PDG1) == 11 ) lep1syst[name] *= 1.000+s*unc1;
            if( abs(PDG2) == 11 ) lep2syst[name] *= 1.000+s*unc2;
         }

         // propagate lepton uncertainties to MET
         if( name.find("MuonMomentumScale") != std::string::npos
               or name.find("ElectronEnergyScale") != std::string::npos ){
            metsyst[name] -= lep1syst[name] - lep1;
            metsyst[name] -= lep2syst[name] - lep2;
         }

         // unclustered energy uncertainty
         if( name.find("METUnclustered") != std::string::npos ){
            metsyst[name] += s*0.1*met_uncl;
         }

      }

   }


   return;
}

// ------------ method called once each job just before starting event loop  ------------
   void 
MakeNtuple::beginJob()
{

   // pile up reweighting
   std::vector< float > PU2012_MC;
   std::vector< float > PU2012_Data;

   for( int i=0; i<60; i++) {
      PU2012_MC.push_back( PU2012_MCf[i] );
      PU2012_Data.push_back( PU2012_Dataf[i] );
   }
   LumiWeights_ = edm::LumiReWeighting( PU2012_MC, PU2012_Data);
   PShiftDown_ = reweight::PoissonMeanShifter(-0.05);
   PShiftUp_ = reweight::PoissonMeanShifter(0.05);


   file = new TFile(outFileName_.c_str(), "RECREATE");
   file->cd();
   trees["Central"] = new TTree("Central", "Central");

   trees["Central"]->Branch("jet1FourVector", &jet1FourVector);
   trees["Central"]->Branch("jet1PtResolution", &jet1PtResolution);
   trees["Central"]->Branch("jet1PhiResolution", &jet1PhiResolution);
   trees["Central"]->Branch("jet1EtaResolution", &jet1EtaResolution);

   trees["Central"]->Branch("jet2FourVector", &jet2FourVector);
   trees["Central"]->Branch("jet2PtResolution", &jet2PtResolution);
   trees["Central"]->Branch("jet2PhiResolution", &jet2PhiResolution);
   trees["Central"]->Branch("jet2EtaResolution", &jet2EtaResolution);

   trees["Central"]->Branch("lep1FourVector", &lep1FourVector);
   trees["Central"]->Branch("lep2FourVector", &lep2FourVector);

   trees["Central"]->Branch("uncorrectedJet1FourVector", &uncorrectedJet1FourVector);
   trees["Central"]->Branch("uncorrectedJet1PtResolution", &uncorrectedJet1PtResolution);
   trees["Central"]->Branch("uncorrectedJet1PhiResolution", &uncorrectedJet1PhiResolution);
   trees["Central"]->Branch("uncorrectedJet1EtaResolution", &uncorrectedJet1EtaResolution);

   trees["Central"]->Branch("uncorrectedJet2FourVector", &uncorrectedJet2FourVector);
   trees["Central"]->Branch("uncorrectedJet2PtResolution", &uncorrectedJet2PtResolution);
   trees["Central"]->Branch("uncorrectedJet2PhiResolution", &uncorrectedJet2PhiResolution);
   trees["Central"]->Branch("uncorrectedJet2EtaResolution", &uncorrectedJet2EtaResolution);

   trees["Central"]->Branch("jet1JESUncertainty", &jet1jesuncertainty);
   trees["Central"]->Branch("jet2JESUncertainty", &jet2jesuncertainty);

   trees["Central"]->Branch("metFourVector", &metFourVector);

   trees["Central"]->Branch("PDG1", &PDG1);
   trees["Central"]->Branch("PDG2", &PDG2);
   trees["Central"]->Branch("vertices", &vertices);
   trees["Central"]->Branch("numBJets", &numBJets);
   trees["Central"]->Branch("jetnumber",&jetcount);
   trees["Central"]->Branch("jet1Vz",&jet1Vz);
   trees["Central"]->Branch("jet2Vz",&jet2Vz);
   trees["Central"]->Branch("jet1bdisc",&jet1bdisc);
   trees["Central"]->Branch("jet2bdisc",&jet2bdisc);

   if (runOnMC_) {
      trees["Central"]->Branch("puMyWeight", &MyWeight);
      trees["Central"]->Branch("puTnvtx", &T_nvertices);
   }
   if (runTtbar_) {
      trees["Central"]->Branch("jet1GenId", &jet1GenId);
      trees["Central"]->Branch("jet2GenId", &jet2GenId);
      trees["Central"]->Branch("generatedMetFourVector", &generatedMetFourVector);
      trees["Central"]->Branch("generatedJet1FourVector", &generatedJet1FourVector);
      trees["Central"]->Branch("generatedJet2FourVector", &generatedJet2FourVector);
      trees["Central"]->Branch("bGEN", &bGEN);
      trees["Central"]->Branch("bbarGEN", &bbarGEN);
      trees["Central"]->Branch("lpGEN", &lpGEN);
      trees["Central"]->Branch("lmGEN", &lmGEN);
      trees["Central"]->Branch("nGEN", &nGEN);
      trees["Central"]->Branch("nbarGEN", &nbarGEN);
      trees["Central"]->Branch("lpPdgIdGEN", &lpPdgIdGEN);
      trees["Central"]->Branch("lmPdgIdGEN", &lmPdgIdGEN);
      trees["Central"]->Branch("nPdgIdGEN", &nPdgIdGEN);
      trees["Central"]->Branch("nbPdgIdGEN", &nbPdgIdGEN);
      trees["Central"]->Branch("jet1ParentIdGEN", &jet1ParentIdGEN);
      trees["Central"]->Branch("jet2ParentIdGEN", &jet2ParentIdGEN);
      trees["Central"]->Branch("geninfo_pid", &geninfo_pid);
      trees["Central"]->Branch("geninfo_pthat", &geninfo_pthat);
      trees["Central"]->Branch("geninfo_weight", &geninfo_weight);
      trees["Central"]->Branch("geninfo_xsec", &geninfo_xsec);
      trees["Central"]->Branch("geninfo_eff", &geninfo_eff);
      trees["Central"]->Branch("geninfo_alphaQCD", &geninfo_alphaQCD);
      trees["Central"]->Branch("geninfo_alphaQED", &geninfo_alphaQED);
   }

   trees["Central"]->Branch("run", &(this_event_id.run));
   trees["Central"]->Branch("lumi", &(this_event_id.lumi));
   trees["Central"]->Branch("event", &(this_event_id.event));

   systnames.push_back("JetEnergyResolution");
   systnames.push_back("METUnclustered");
   systnames.push_back("PileUp");
   systnames.push_back("ElectronEnergyScale");
   systnames.push_back("MuonMomentumScale");
   systnames.push_back("ElectronId");
   systnames.push_back("MuonId");
   systnames.push_back("BTaggingEff");
   systnames.push_back("PDF");
   systnames.push_back("PtTopReweighting");

   jsystnames.push_back("CorrelationGroupMPFInSitu");
   jsystnames.push_back("CorrelationGroupFlavor");
   jsystnames.push_back("CorrelationGroupIntercalibration");
   jsystnames.push_back("CorrelationGroupUncorrelated");
   jsystnames.push_back("CorrelationGroupbJES");
   jsystnames.push_back("Total");

   std::string sgn [2] = {"DN","UP"};
   for(std::vector<std::string>::iterator syst = jsystnames.begin(); syst != jsystnames.end(); syst++){
      for(int i=0; i < 2; i++){
         std::string name = *syst+sgn[i];
         trees[name] = trees["Central"]->CloneTree(0);
         trees[name]->SetName(name.c_str());
         trees[name]->SetTitle(name.c_str());
      }
   }
   for(std::vector<std::string>::iterator syst = systnames.begin(); syst != systnames.end(); syst++){
      for(int i=0; i < 2; i++){
         std::string name = *syst+sgn[i];
         trees[name] = trees["Central"]->CloneTree(0);
         trees[name]->SetName(name.c_str());
         trees[name]->SetTitle(name.c_str());
      }
   }

}


// ------------ method called once each job just after ending the event loop  ------------
   void 
MakeNtuple::endJob() 
{
   file->Write();
   file->Close();
}

// ------------ method called when starting to processes a run  ------------
   void 
MakeNtuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
   void 
MakeNtuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
   void 
MakeNtuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
   void 
MakeNtuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeNtuple);
