// -*- C++ -*-
//
// Package:    BTagEff
// Class:      BTagEff
// 
/**\class BTagEff BTagEff.cc TopMassAnalysis/BTagEff/src/BTagEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed Jul 13 16:26:11 EDT 2016
// $Id$
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "LHAPDF/LHAPDF.h"

#include "TLorentzVector.h"
#include "TMatrix.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH2.h"

#include <string>

//
// class declaration
//

class BTagEff : public edm::EDAnalyzer {
   public:
      explicit BTagEff(const edm::ParameterSet&);
      ~BTagEff();

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
            const TLorentzVector lep1, const TLorentzVector lep2 );

      double jetScale_, leptonScale_, jetResScale_;
      std::string bTagger_;
      double bTagCut_, negTagCut_;
      bool doMETCut_, runOnMC_, runTtbar_;
      std::string outFileName_;

      edm::InputTag               muonTag_;
      edm::InputTag               electronTag_;
      edm::InputTag               jetTag_;
      edm::InputTag               metTag_;
      edm::InputTag               pmetTag_;
      edm::InputTag               genParticleTag_;

      int randSeed_;

      std::vector<pat::Jet> taggedJets_;
      std::vector<pat::Jet> allJets_;
      std::vector<pat::Jet> negTaggedJets_;
      std::vector<const reco::RecoCandidate*> goodLeptons_;
      std::vector<const reco::RecoCandidate*> goodMuons_;
      std::vector<const reco::RecoCandidate*> goodElectrons_;

      const reco::Candidate *t, *tb, *b, *bb, *Wp, *Wm, *lp, *lm, *n, *nb;

      TLorentzVector jet1FourVector, jet2FourVector, lep1FourVector, lep2FourVector;
      TLorentzVector uncorrectedJet1FourVector, uncorrectedJet2FourVector;
      TLorentzVector generatedJet1FourVector, generatedJet2FourVector, metFourVector, generatedMetFourVector;
      TLorentzVector bGEN, bbarGEN, lpGEN, lmGEN, nGEN, nbarGEN;
      TLorentzVector metUnclustered;
      int lpPdgIdGEN, lmPdgIdGEN;
      int nPdgIdGEN, nbPdgIdGEN;

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


      // ----------member data ---------------------------
      const int     ptNBins;
      const double  ptMin;
      const double  ptMax;
      const int     etaNBins;
      const double  etaMin;
      const double  etaMax;
      edm::Service<TFileService>  fs;
      TH2D  *h2_BTaggingEff_Denom_b;
      TH2D  *h2_BTaggingEff_Denom_c;
      TH2D  *h2_BTaggingEff_Denom_udsg;
      TH2D  *h2_BTaggingEff_Num_b;
      TH2D  *h2_BTaggingEff_Num_c;
      TH2D  *h2_BTaggingEff_Num_udsg;
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
BTagEff::BTagEff(const edm::ParameterSet& iConfig):
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
   // tree (NULL);
   muonTag_         (iConfig.getParameter<edm::InputTag>("muonSrc") ),
   electronTag_     (iConfig.getParameter<edm::InputTag>("electronSrc") ),
   jetTag_          (iConfig.getParameter<edm::InputTag>("jetSrc") ),
   metTag_          (iConfig.getParameter<edm::InputTag>("metSrc") ),
   pmetTag_          (iConfig.getParameter<edm::InputTag>("pmetSrc") ),
   genParticleTag_         (iConfig.getParameter<edm::InputTag>("genParticleSrc") ),
   randSeed_        (iConfig.getParameter<int>("randSeed")),
   ptNBins(iConfig.getParameter<int>("PtNBins")),
   ptMin(iConfig.getParameter<double>("PtMin")),
   ptMax(iConfig.getParameter<double>("PtMax")),
   etaNBins(iConfig.getParameter<int>("EtaNBins")),
   etaMin(iConfig.getParameter<double>("EtaMin")),
   etaMax(iConfig.getParameter<double>("EtaMax"))


{
   //now do what ever initialization is needed
   h2_BTaggingEff_Denom_b    = fs->make<TH2D>("h2_BTaggingEff_Denom_b", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Denom_c    = fs->make<TH2D>("h2_BTaggingEff_Denom_c", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Denom_udsg = fs->make<TH2D>("h2_BTaggingEff_Denom_udsg", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Num_b    = fs->make<TH2D>("h2_BTaggingEff_Num_b", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Num_c    = fs->make<TH2D>("h2_BTaggingEff_Num_c", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Num_udsg = fs->make<TH2D>("h2_BTaggingEff_Num_udsg", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);

}


BTagEff::~BTagEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BTagEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   taggedJets_.clear();
   allJets_.clear();
   negTaggedJets_.clear();
   goodLeptons_.clear();
   goodMuons_.clear();
   goodElectrons_.clear();

   // leptons
   edm::Handle<edm::View<pat::Electron> > electrons;
   iEvent.getByLabel(electronTag_, electrons);
   int countmu=0, counte=0;
   for( edm::View<pat::Electron>::const_iterator electron = electrons->begin(); electron != electrons->end(); ++electron) {
      const reco::RecoCandidate* lepton = &*electron;
      goodLeptons_.push_back(lepton);
      goodElectrons_.push_back(lepton);
      counte++;
   }
   edm::Handle<edm::View<pat::Muon> > muons;
   iEvent.getByLabel(muonTag_, muons);
   for( edm::View<pat::Muon>::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
      const reco::RecoCandidate* lepton = &*muon;
      goodLeptons_.push_back(lepton);
      goodMuons_.push_back(lepton);
      countmu++;
   }

   if (goodLeptons_.size() < 2) return;
   if (goodLeptons_.size() > 2){
      std::cout << "ERROR: MORE THAN 2 LEPTONS: ";
      for(unsigned int i=0; i < goodLeptons_.size(); i++) std::cout << goodLeptons_[i]->pt() << " ";
      std::cout << std::endl;
      return;
   }
   const reco::RecoCandidate* lep1 = goodLeptons_[0];
   const reco::RecoCandidate* lep2 = goodLeptons_[1];

   lep1FourVector.SetPxPyPzE(lep1->px(), lep1->py(), lep1->pz(), lep1->energy());
   lep2FourVector.SetPxPyPzE(lep2->px(), lep2->py(), lep2->pz(), lep2->energy());

   PDG1 = lep1->pdgId();
   PDG2 = lep2->pdgId();


   // met
   edm::Handle<edm::View<reco::PFMET> > mets;
   iEvent.getByLabel(metTag_, mets);
   reco::PFMET met = mets->front();

   edm::Handle<edm::View<pat::MET> > pmets;
   iEvent.getByLabel(pmetTag_, pmets);
   pat::MET pmet = pmets->front();

   metFourVector.SetPxPyPzE(met.px(), met.py(), met.pz(), met.energy());

   // jets
   // pass b tag cut
   edm::Handle<edm::View<pat::Jet> > jets;
   iEvent.getByLabel(jetTag_,jets);
   for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet!=jets->end();++jet){
      if( jet->pt() >= 30 and fabs(jet->eta()) < 2.5 ){
         if (jet->bDiscriminator(bTagger_) > bTagCut_) {
            taggedJets_.push_back(*jet);
         }
         allJets_.push_back(*jet);
      }
   }

   bool pass_selection = false;
   pass_selection = passOfflineSelection( metFourVector, lep1FourVector, lep2FourVector ) and allJets_.size() > 1;

   if( !pass_selection ) return;

   // loop over jets
   for(std::vector<pat::Jet>::iterator it = allJets_.begin(); it != allJets_.end(); ++it){

      int partonFlavor = it->partonFlavour();
      int jpt = it->pt();
      int jeta = fabs(it->eta());

      if( abs(partonFlavor)==5 )
      {
         h2_BTaggingEff_Denom_b->Fill(jpt, jeta);
         if( it->bDiscriminator(bTagger_) > bTagCut_ ) h2_BTaggingEff_Num_b->Fill(jpt, jeta);
      }
      else if( abs(partonFlavor)==4 )
      {
         h2_BTaggingEff_Denom_c->Fill(jpt, jeta);
         if( it->bDiscriminator(bTagger_) > bTagCut_ ) h2_BTaggingEff_Num_c->Fill(jpt, jeta);
      }
      else
      {
         h2_BTaggingEff_Denom_udsg->Fill(jpt, jeta);
         if( it->bDiscriminator(bTagger_) > bTagCut_ ) h2_BTaggingEff_Num_udsg->Fill(jpt, jeta);
      }

   }


}

bool
BTagEff::passOfflineSelection( const TLorentzVector met,
      const TLorentzVector lep1, const TLorentzVector lep2 ){

   double met_pt = 40; // only for ee and mumu events
   double mu_pt = 20;
   double mu_eta = 2.4;
   double e_pt = 20;
   double e_eta = 2.5;

   bool met_ok = (abs(PDG1) != abs(PDG2)) or (met.Pt() > met_pt);
   bool lep1_ok = abs(PDG1) == 11 ? (lep1.Pt() > e_pt and fabs(lep1.Eta()) < e_eta)
      : (lep1.Pt() > mu_pt and fabs(lep1.Eta()) < mu_eta);
   bool lep2_ok = abs(PDG2) == 11 ? (lep2.Pt() > e_pt and fabs(lep2.Eta()) < e_eta)
      : (lep2.Pt() > mu_pt and fabs(lep2.Eta()) < mu_eta);
   bool Zpeak_ok = abs(PDG1) != abs(PDG2) or ((lep1+lep2).M() < 60 or (lep1+lep2).M() > 120);

   return (met_ok and lep1_ok and lep2_ok and Zpeak_ok );
}



// ------------ method called once each job just before starting event loop  ------------
void 
BTagEff::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagEff::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
BTagEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BTagEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BTagEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BTagEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BTagEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagEff);
