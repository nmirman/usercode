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
      edm::InputTag               genParticleTag_;

      std::vector<pat::Jet> taggedJets_;
      std::vector<pat::Jet> negTaggedJets_;
      std::vector<const reco::RecoCandidate*> goodLeptons_;

      // This will be a pointer to the tree we want to use for this event
      TTree *tree, *treeData, *treeBkg, *treeBkg2;
      TFile *file;

      JetResolution *ptResol_;
      JetResolution *phiResol_;
      JetResolution *etaResol_;

      const reco::Candidate *t, *tb, *b, *bb, *Wp, *Wm, *lp, *lm, *n, *nb;

      TLorentzVector jet1FourVector, jet2FourVector, lepton1FourVector, lepton2FourVector;
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
      TMatrixD metSignificanceMatrix;
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
      float MyWeight;
      float T_nvertices;

     int geninfo_pid;
     float geninfo_pthat;
     float geninfo_weight;
     float geninfo_xsec;
     float geninfo_eff;
     float geninfo_alphaQCD;
     float geninfo_alphaQED;

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
metSignificanceMatrix (TMatrixD(2, 2)),
// tree (NULL);
muonTag_         (iConfig.getParameter<edm::InputTag>("muonSrc") ),
electronTag_     (iConfig.getParameter<edm::InputTag>("electronSrc") ),
jetTag_          (iConfig.getParameter<edm::InputTag>("jetSrc") ),
metTag_          (iConfig.getParameter<edm::InputTag>("metSrc") ),
genParticleTag_         (iConfig.getParameter<edm::InputTag>("genParticleSrc") )



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
  if (taggedJets_.size() >= 2) { // signal
    tree = treeData;
    jet1 = taggedJets_[0];
    jet2 = taggedJets_[1];
  }
  else if (negTaggedJets_.size() >= 2) { // uu background
    tree = treeBkg2;
    jet1 = negTaggedJets_[0];
    jet2 = negTaggedJets_[1];
  }
  else if (taggedJets_.size() == 1 && negTaggedJets_.size() == 1) { // bu background
    tree = treeBkg;
    jet1 = taggedJets_[0];
    jet2 = negTaggedJets_[0];
  }
  else return;

  jet1Vz = jet1.vz()-primary_vertex_z;
  jet2Vz = jet2.vz()-primary_vertex_z;
  //jet1VzCovariance =jet1.vertexCovariance(2,2);
  //jet2VzCovariance =jet2.vertexCovariance(2,2);
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
        //std::cout << "Wp\n";
        //std::cout << Wp->daughter(1)->pdgId() << std::endl;
        //std::cout << Wp->daughter(0)->pdgId() << std::endl;
        //return; // not a dilepton event
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
        //std::cout << "Wm\n";
        //std::cout << Wm->daughter(1)->pdgId() << std::endl;
        //std::cout << Wm->daughter(0)->pdgId() << std::endl;
        //return; // not a dilepton event
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
  edm::Handle<edm::View<pat::MET> > mets;
  //edm::Handle<edm::View<reco::MET> > mets;
  iEvent.getByLabel(metTag_, mets);
  pat::MET met = mets->front();
  
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
  std::string path = std::getenv("CMSSW_BASE");
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( (path+"/src/PhysicsTools/PatUtils/data/Summer13_V4_DATA_Uncertainty_AK5PF.txt").c_str() );
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

    //std::cout << "JES SCALE TEST: " << std::endl;
    //std::cout << "jet1 pt, unc: " << jet1FourVector.Pt() << " " << jet1jesuncertainty << std::endl;
    //std::cout << "jet2 pt, unc: " << jet2FourVector.Pt() << " " << jet2jesuncertainty << std::endl;
    //std::cout << "met: " << metFourVector.Pt() << std::endl;
    //fflush(stdout);

  numBJets = taggedJets_.size();

  lepton1FourVector.SetPxPyPzE(lep1->px(), lep1->py(), lep1->pz(), lep1->energy());
  lepton2FourVector.SetPxPyPzE(lep2->px(), lep2->py(), lep2->pz(), lep2->energy());

  PDG1 = lep1->pdgId();
  PDG2 = lep2->pdgId();

	if (runTtbar_) {
	  generatedMetFourVector.SetPxPyPzE(met.genMET()->px(), met.genMET()->py(), met.genMET()->pz(), met.genMET()->energy());
	}
  metSignificanceMatrix = met.getSignificanceMatrix();

  tree->Fill();
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


  file = new TFile(outFileName_.c_str(), "RECREATE");
  file->cd();
  treeData = new TTree("RealData", "RealData");

  treeData->Branch("jet1FourVector", &jet1FourVector);
  treeData->Branch("jet1PtResolution", &jet1PtResolution);
  treeData->Branch("jet1PhiResolution", &jet1PhiResolution);
  treeData->Branch("jet1EtaResolution", &jet1EtaResolution);

  treeData->Branch("jet2FourVector", &jet2FourVector);
  treeData->Branch("jet2PtResolution", &jet2PtResolution);
  treeData->Branch("jet2PhiResolution", &jet2PhiResolution);
  treeData->Branch("jet2EtaResolution", &jet2EtaResolution);

  treeData->Branch("lepton1FourVector", &lepton1FourVector);
  treeData->Branch("lepton2FourVector", &lepton2FourVector);

  treeData->Branch("uncorrectedJet1FourVector", &uncorrectedJet1FourVector);
  treeData->Branch("uncorrectedJet1PtResolution", &uncorrectedJet1PtResolution);
  treeData->Branch("uncorrectedJet1PhiResolution", &uncorrectedJet1PhiResolution);
  treeData->Branch("uncorrectedJet1EtaResolution", &uncorrectedJet1EtaResolution);

  treeData->Branch("uncorrectedJet2FourVector", &uncorrectedJet2FourVector);
  treeData->Branch("uncorrectedJet2PtResolution", &uncorrectedJet2PtResolution);
  treeData->Branch("uncorrectedJet2PhiResolution", &uncorrectedJet2PhiResolution);
  treeData->Branch("uncorrectedJet2EtaResolution", &uncorrectedJet2EtaResolution);

  treeData->Branch("jet1JESUncertainty", &jet1jesuncertainty);
  treeData->Branch("jet2JESUncertainty", &jet2jesuncertainty);

  treeData->Branch("metFourVector", &metFourVector);
  treeData->Branch("metSignificanceMatrix", &metSignificanceMatrix);

  treeData->Branch("PDG1", &PDG1);
  treeData->Branch("PDG2", &PDG2);
  treeData->Branch("vertices", &vertices);
  treeData->Branch("numBJets", &numBJets);
  treeData->Branch("jetnumber",&jetcount);
  treeData->Branch("jet1Vz",&jet1Vz);
  treeData->Branch("jet2Vz",&jet2Vz);
  treeData->Branch("jet1bdisc",&jet1bdisc);
  treeData->Branch("jet2bdisc",&jet2bdisc);

  if (runOnMC_) {
     treeData->Branch("puMyWeight", &MyWeight);
     treeData->Branch("puTnvtx", &T_nvertices);
  }
  if (runTtbar_) {
     treeData->Branch("jet1GenId", &jet1GenId);
     treeData->Branch("jet2GenId", &jet2GenId);
     treeData->Branch("generatedMetFourVector", &generatedMetFourVector);
     treeData->Branch("generatedJet1FourVector", &generatedJet1FourVector);
     treeData->Branch("generatedJet2FourVector", &generatedJet2FourVector);
     treeData->Branch("bGEN", &bGEN);
     treeData->Branch("bbarGEN", &bbarGEN);
     treeData->Branch("lpGEN", &lpGEN);
     treeData->Branch("lmGEN", &lmGEN);
     treeData->Branch("nGEN", &nGEN);
     treeData->Branch("nbarGEN", &nbarGEN);
     treeData->Branch("lpPdgIdGEN", &lpPdgIdGEN);
     treeData->Branch("lmPdgIdGEN", &lmPdgIdGEN);
     treeData->Branch("nPdgIdGEN", &nPdgIdGEN);
     treeData->Branch("nbPdgIdGEN", &nbPdgIdGEN);
     treeData->Branch("jet1ParentIdGEN", &jet1ParentIdGEN);
     treeData->Branch("jet2ParentIdGEN", &jet2ParentIdGEN);
     treeData->Branch("geninfo_pid", &geninfo_pid);
     treeData->Branch("geninfo_pthat", &geninfo_pthat);
     treeData->Branch("geninfo_weight", &geninfo_weight);
     treeData->Branch("geninfo_xsec", &geninfo_xsec);
     treeData->Branch("geninfo_eff", &geninfo_eff);
     treeData->Branch("geninfo_alphaQCD", &geninfo_alphaQCD);
     treeData->Branch("geninfo_alphaQED", &geninfo_alphaQED);
  }

  treeData->Branch("run", &(this_event_id.run));
  treeData->Branch("lumi", &(this_event_id.lumi));
  treeData->Branch("event", &(this_event_id.event));

  treeBkg = treeData->CloneTree(0);
  treeBkg->SetName("buBkg");
  treeBkg->SetTitle("buBkg");

  treeBkg2 = treeData->CloneTree(0);
  treeBkg2->SetName("uuBkg");
  treeBkg2->SetTitle("uuBkg");

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
