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

#include "LHAPDF/LHAPDF.h"

#include "TLorentzVector.h"
#include "TMatrix.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"

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
      std::vector<const reco::RecoCandidate*> goodMuons_;
      std::vector<const reco::RecoCandidate*> goodElectrons_;

      // This will be a pointer to the tree we want to use for this event
      std::map<std::string, TTree*> trees;
      TTree *tree;
      TFile *file;
      TFile *fmueff;

      JetResolution *ptResol_;
      JetResolution *phiResol_;
      JetResolution *etaResol_;

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

      int nmuons, nelectrons;

      static const double PU2012_MCf[60];
      static const double PU2012_Dataf[60];
      static const double PU2012_DatafUP[60];
      static const double PU2012_DatafDN[60];
      edm::LumiReWeighting LumiWeights_;
      edm::LumiReWeighting LumiWeightsUP_;
      edm::LumiReWeighting LumiWeightsDN_;
      float weight_pu;
      float weight_pu_UP;
      float weight_pu_DN;
      float weight_toppt;
      float weight_btag;
      float weight_mu;
      float weight_elec;
      float T_nvertices;
      float weight_bjes_nuUP;
      float weight_bjes_nuDN;
      float weight_bjes_rbLEP;
      float weight_bfrag;

      int geninfo_pid;
      float geninfo_pthat;
      float geninfo_weight;
      float geninfo_xsec;
      float geninfo_eff;
      float geninfo_alphaQCD;
      float geninfo_alphaQED;

      float geninfo_scalePDF;
      float geninfo_id1;
      float geninfo_id2; 
      float geninfo_x1;
      float geninfo_x2; 
      float geninfo_xPDF1;  
      float geninfo_xPDF2;
      std::vector<float> pdf_weights;

      std::vector<std::string> systnames;
      std::vector<std::string> jsystnames;
      std::map<std::string, TLorentzVector> metsyst, jet1syst, jet2syst, lep1syst, lep2syst;

      std::vector<JetCorrectionUncertainty*> vsrc;

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
   // obtained with pileupCalc.py (03.17.2015)
   12238.5,
   32499.7,
   92004.3,
   336891,
   615848,
   2.99224e+06,
   1.74697e+07,
   5.34979e+07,
   1.29072e+08,
   2.56502e+08,
   4.42325e+08,
   6.80775e+08,
   8.77629e+08,
   9.95861e+08,
   1.07561e+09,
   1.13527e+09,
   1.16992e+09,
   1.18034e+09,
   1.1756e+09,
   1.15999e+09,
   1.13529e+09,
   1.10426e+09,
   1.06784e+09,
   1.02145e+09,
   9.56904e+08,
   8.69452e+08,
   7.62038e+08,
   6.4278e+08,
   5.20651e+08,
   4.04012e+08,
   3.00029e+08,
   2.13254e+08,
   1.44959e+08,
   9.39843e+07,
   5.79133e+07,
   3.38274e+07,
   1.87238e+07,
   9.84719e+06,
   4.94827e+06,
   2.39536e+06,
   1.12892e+06,
   524821,
   244475,
   116129,
   57201.2,
   29562.3,
   16083.7,
   9165.98,
   5418.61,
   3288.31,
   2030.33,
   1266.68,
   794.193,
   498.231,
   311.566,
   193.591,
   119.194,
   72.5585,
   43.5935,
   25.8146,
};

const double MakeNtuple::PU2012_DatafUP[60] = {
   // 'true' distribution for 2012 dataset
   // obtained with pileupCalc.py (03.17.2015)
   // luminosity increased by 5%
   11510.3,
   22652.3,
   82217.2,
   259033,
   549655,
   1.5763e+06,
   1.04678e+07,
   3.54953e+07,
   8.76743e+07,
   1.83325e+08,
   3.24051e+08,
   5.22025e+08,
   7.34096e+08,
   8.83328e+08,
   9.74281e+08,
   1.04131e+09,
   1.09053e+09,
   1.11736e+09,
   1.12425e+09,
   1.11893e+09,
   1.10432e+09,
   1.08189e+09,
   1.05396e+09,
   1.02151e+09,
   9.81233e+08,
   9.26384e+08,
   8.52021e+08,
   7.5908e+08,
   6.53504e+08,
   5.42733e+08,
   4.33866e+08,
   3.33348e+08,
   2.46117e+08,
   1.74608e+08,
   1.18866e+08,
   7.74316e+07,
   4.81206e+07,
   2.84783e+07,
   1.60549e+07,
   8.64648e+06,
   4.47159e+06,
   2.23688e+06,
   1.09245e+06,
   526812,
   254249,
   124670,
   63037.7,
   33241,
   18366.5,
   10606.1,
   6353.01,
   3911.3,
   2454.36,
   1559.17,
   997.472,
   639.906,
   410.197,
   261.934,
   166.172,
   104.5,
};

const double MakeNtuple::PU2012_DatafDN[60] = {
   // 'true' distribution for 2012 dataset
   // obtained with pileupCalc.py (03.17.2015)
   // luminosity decreased by 5%
   13063.1,
   44782.4,
   108214,
   432674,
   770689,
   5.79686e+06,
   2.8364e+07,
   8.15956e+07,
   1.89599e+08,
   3.58821e+08,
   6.02838e+08,
   8.52318e+08,
   1.01239e+09,
   1.11001e+09,
   1.18219e+09,
   1.22677e+09,
   1.24206e+09,
   1.23825e+09,
   1.22156e+09,
   1.19423e+09,
   1.15955e+09,
   1.11835e+09,
   1.06428e+09,
   9.87551e+08,
   8.84235e+08,
   7.60075e+08,
   6.25723e+08,
   4.91962e+08,
   3.68576e+08,
   2.6303e+08,
   1.78765e+08,
   1.15453e+08,
   7.05755e+07,
   4.06879e+07,
   2.20968e+07,
   1.13309e+07,
   5.519e+06,
   2.57704e+06,
   1.1677e+06,
   521332,
   233595,
   107219,
   51375.5,
   26006.1,
   13923.5,
   7821.75,
   4555.09,
   2718.26,
   1646.87,
   1005.84,
   615.858,
   376.265,
   228.484,
   137.442,
   81.6777,
   47.8496,
   27.5889,
   15.6368,
   8.7044,
   4.75596,
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
   goodMuons_.clear();
   goodElectrons_.clear();

   weight_pu=1;
   weight_toppt=1;
   weight_btag=1;
   weight_mu=1;
   weight_elec=1;
   weight_bfrag=1;

   // pdfs
   /*
   if( runOnMC_ ){

      edm::Handle<reco::PdfInfo> pdfs;
      iEvent.getByType( pdfstuff );
      cout << pdfs->id1 << " " << pdfstuff->x1 << endl;
      cout << pdfs->id2 << " " << pdfstuff->x2 << endl;
      cout << pdfs->scalePDF << endl;

   }
   */

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

      geninfo_scalePDF = 0;
      geninfo_id1 = 0;
      geninfo_id2 = 0; 
      geninfo_x1 = 0;
      geninfo_x2 = 0; 
      geninfo_xPDF1 = 0;  
      geninfo_xPDF2 = 0;

      pdf_weights.clear();

      if(gi->hasPDF()){
         geninfo_scalePDF = gi->pdf()->scalePDF;
         geninfo_id1        = gi->pdf()->id.first;
         geninfo_id2        = gi->pdf()->id.second;
         geninfo_x1      = gi->pdf()->x.first;
         geninfo_x2      = gi->pdf()->x.second;
         geninfo_xPDF1    = gi->pdf()->xPDF.first;
         geninfo_xPDF2    = gi->pdf()->xPDF.second;

         LHAPDF::usePDFMember(1,0);
         double xpdf1 = LHAPDF::xfx(1, geninfo_x1, geninfo_scalePDF, geninfo_id1);
         double xpdf2 = LHAPDF::xfx(1, geninfo_x2, geninfo_scalePDF, geninfo_id2);
         double w0 = xpdf1 * xpdf2;
         for(int i=1; i <=50; ++i){
            LHAPDF::usePDFMember(1,i);
            double xpdf1_new = LHAPDF::xfx(1, geninfo_x1, geninfo_scalePDF, geninfo_id1);
            double xpdf2_new = LHAPDF::xfx(1, geninfo_x2, geninfo_scalePDF, geninfo_id2);
            double weight = xpdf1_new * xpdf2_new / w0;
            pdf_weights.push_back(weight);
         }

      }

      /*
      edm::InputTag pdfWeightTag("pdfWeights:cteq66");
      edm::Handle<std::vector<double> > weightHandle;
      iEvent.getByLabel(pdfWeightTag, weightHandle);

      std::vector<double> weights = (*weightHandle);
      std::cout << "Event weight for central PDF:" << weights[0] << std::endl;
      unsigned int nmembers = weights.size();
      for (unsigned int j=1; j<nmembers; j+=2) {
         std::cout << "Event weight for PDF variation +" << (j+1)/2 << ": " << weights[j] << std::endl;
         std::cout << "Event weight for PDF variation -" << (j+1)/2 << ": " << weights[j+1] << std::endl;
      }
      */

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

      }

      weight_pu = LumiWeights_.weight( Tnvtx );
      weight_pu_UP = LumiWeightsUP_.weight( Tnvtx );
      weight_pu_DN = LumiWeightsDN_.weight( Tnvtx );
      T_nvertices = Tnvtx;

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

   nmuons = goodMuons_.size();
   nelectrons = goodElectrons_.size();

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

   numBJets = taggedJets_.size();

   lep1FourVector.SetPxPyPzE(lep1->px(), lep1->py(), lep1->pz(), lep1->energy());
   lep2FourVector.SetPxPyPzE(lep2->px(), lep2->py(), lep2->pz(), lep2->energy());

   PDG1 = lep1->pdgId();
   PDG2 = lep2->pdgId();

   if (runTtbar_) {
      generatedMetFourVector.SetPxPyPzE(pmet.genMET()->px(), pmet.genMET()->py(), pmet.genMET()->pz(), pmet.genMET()->energy());
   }

   // calculate unclustered energy contribution to the met
   TLorentzVector met_uncl = metFourVector + lep1FourVector + lep2FourVector;
   for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet!=jets->end();++jet){
      TLorentzVector jetFourVector;
      jetFourVector.SetPxPyPzE(jet->px(), jet->py(), jet->pz(), jet->energy());
      met_uncl += jetFourVector;
   }
   metUnclustered = met_uncl;

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

   // for b tagging efficiency weight
   double weight_bjets [] = {-1,-1};
   weight_bjets[0] = 0.997942*((1.+(0.00923753*jet1FourVector.Pt()))/(1.+(0.0096119*jet1FourVector.Pt())));
   weight_bjets[1] = 0.997942*((1.+(0.00923753*jet2FourVector.Pt()))/(1.+(0.0096119*jet2FourVector.Pt())));
   weight_btag = weight_bjets[0]*weight_bjets[1];

   // muon scale factors
   TGraph *gmueff = (TGraph*)fmueff->Get("DATA_over_MC_Loose_eta_pt20-500");
   double weight_muons [] = {-1,-1};
   for(unsigned int i=0; i < goodMuons_.size(); i++){
      for(int n=0; n < gmueff->GetN(); n++){
         double x, y;
         gmueff->GetPoint(n,x,y);
         double low_edge = x - gmueff->GetErrorX(n);
         double high_edge = x + gmueff->GetErrorX(n);
         if( goodMuons_[i]->eta() >= low_edge and goodMuons_[i]->eta() < high_edge ) weight_muons[i] = y;
      }
   }

   // electron scale factors
   double elec_factors [3][4] = {
      {0.904, 0.977, 0.978, 0.968},
      {0.980, 0.979, 0.980, 0.964},
      {0.933, 0.963, 0.976, 0.975}
   };
   double elec_factors_errUP [3][4] = {
      {0.005, 0.001, 0.001, 0.001},
      {0.010, 0.003, 0.000, 0.001},
      {0.121, 0.144, 0.001, 0.001}
   };
   double elec_factors_errDN [3][4] = {
      {0.005, 0.001, 0.001, 0.001},
      {0.008, 0.003, 0.000, 0.001},
      {0.129, 0.144, 0.001, 0.002}
   };
   double elec_etamin [] = {0.00, 0.80, 1.48};
   double elec_etamax [] = {0.80, 1.48, 2.50};
   double elec_ptmin [] = {20, 30, 40, 50};
   double elec_ptmax [] = {30, 40, 50, 150};

   // b fragmentation weights
   edm::Handle<double> bjesweight_nuUP;
   edm::Handle<double> bjesweight_nuDN;
   edm::Handle<double> bjesweight_rbLEP;
   iEvent.getByLabel("bjesweightNUUP", bjesweight_nuUP);
   iEvent.getByLabel("bjesweightNUDN", bjesweight_nuDN);
   iEvent.getByLabel("bjesweightRBLEP", bjesweight_rbLEP);

   weight_bjes_nuUP = *bjesweight_nuUP;
   weight_bjes_nuDN = *bjesweight_nuDN;
   weight_bjes_rbLEP = *bjesweight_rbLEP;
   
   //
   // fill trees for systematics samples
   //
   double weight_pu_temp = weight_pu;
   double weight_toppt_temp = weight_toppt;
   double weight_btag_temp = weight_btag;
   double weight_mu_temp = weight_mu;
   double weight_elec_temp = weight_elec;
   double weight_bfrag_temp = weight_bfrag;
   for( std::map<std::string, TLorentzVector>::iterator syst = metsyst.begin(); syst != metsyst.end(); syst++ ){
      std::string name = syst->first;

      weight_pu  = weight_pu_temp;
      weight_toppt = weight_toppt_temp;
      weight_btag = weight_btag_temp;
      weight_mu = weight_mu_temp;
      weight_elec = weight_elec_temp;
      weight_bfrag = weight_bfrag_temp;

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
      if( name == "PileUpUP" ) weight_pu = weight_pu_UP;
      if( name == "PileUpDN" ) weight_pu = weight_pu_DN;

      // top pt reweighting
      double a = 0.148;
      double b = -0.00129;
      if( name.find("PtTopReweighting") != std::string::npos ) weight_toppt = sqrt( exp(a+b*t_pt)*exp(a+b*tb_pt) );

      // b tag scale factors
      // Tagger: CSVL within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
      float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
      float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
      double SFb_error [] = { 0.033299, 0.0146768, 0.013803, 0.0170145, 0.0166976, 0.0137879, 0.0149072, 0.0153068, 0.0133077, 0.0123737, 0.0157152, 0.0175161, 0.0209241, 0.0278605, 0.0346928, 0.0350099 };
      int index_jet1 = -1;
      int index_jet2 = -1;
      for(int i=0; i < 16; i++){
         if( jet1FourVector.Pt() >= ptmin[i] and jet1FourVector.Pt() < ptmax[i] ) index_jet1 = i;
         if( jet2FourVector.Pt() >= ptmin[i] and jet2FourVector.Pt() < ptmax[i] ) index_jet2 = i;
      }
      if( name == "BTaggingUP" ){
         weight_btag = (weight_bjets[0]+SFb_error[index_jet1])*(weight_bjets[1]+SFb_error[index_jet2]);
      }
      else if( name == "BTaggingDN" ){
         weight_btag = (weight_bjets[0]-SFb_error[index_jet1])*(weight_bjets[1]-SFb_error[index_jet2]);
      }else{
         weight_btag = weight_bjets[0]*weight_bjets[1];
      }

      // muon ID scale factors
      int count_mu = 0;
      for(int i=0; i < 2; i++){
         if( weight_muons[i] != -1 ){
            weight_mu *= weight_muons[i];
            count_mu++;
         }
      }
      if( name == "MuonIdUP" ) weight_mu *= pow(1.005,count_mu);
      if( name == "MuonIdDN" ) weight_mu *= pow(0.995,count_mu);

      // electron scale factors
      for(unsigned int i=0; i < goodElectrons_.size(); i++){

         double eta = goodElectrons_[i]->eta();
         double pt = goodElectrons_[i]->pt();
         int ipt=-1, ieta=-1;
         for(int j=0; j < 3; j++){
            if( fabs(eta) >= elec_etamin[j] and fabs(eta) < elec_etamax[j] ) ieta = j;
         }
         for(int j=0; j < 4; j++){
            if( pt >= elec_ptmin[j] and pt < elec_ptmax[j] ) ipt = j;
         }
         if( ieta != -1 and ipt != -1 ){
            if( name == "ElectronIdUP" ){
               weight_elec *= elec_factors[ieta][ipt] + elec_factors_errUP[ieta][ipt];
            }
            else if( name == "ElectronIdDN" ){
               weight_elec *= elec_factors[ieta][ipt] - elec_factors_errDN[ieta][ipt];
            }else{
               weight_elec *= elec_factors[ieta][ipt];
            }
         }

      }

      if( name == "BFRAGnuUP" ) weight_bfrag *= weight_bjes_nuUP;
      if( name == "BFRAGnuDN" ) weight_bfrag *= weight_bjes_nuDN;
      if( name == "BFRAGrbLEPUP" ) weight_bfrag *= weight_bjes_rbLEP;

      // apply jet smearing
      int ismear = 0;
      if( name == "JetEnergyResolutionUP" ){
         ismear = 1;
      }else if( name == "JetEnergyResolutionDN" ){
         ismear = 2;
      }else{
         ismear = 0;
      }
      smearJetsMET( jet1Unsmeared, jet2Unsmeared, metUnsmeared,
            jet1FourVector, jet2FourVector, metFourVector,
            jet1_rE[ismear], jet2_rE[ismear], met_dx[ismear], met_dy[ismear], met_dz[ismear], met_dE[ismear] ); 

      // a hack for the Central variation (no systematics applied)
      if( name == "CentralUP" ) name = "Central";
      if( name == "CentralDN" ) continue;

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

   //double met_pt = 40; // only for ee and mumu events
   double met_pt = 30;
   //double jet_pt = 30;
   double jet_pt = 20;
   double jet_eta = 2.5;
   double mu_pt = 20;
   double mu_eta = 2.4;
   double e_pt = 20;
   double e_eta = 2.5;

   bool met_ok = (abs(PDG1) != abs(PDG2)) or (met.Pt() > met_pt);
   bool jet1_ok = jet1.Pt() > jet_pt and fabs(jet1.Eta()) < jet_eta;
   bool jet2_ok = jet2.Pt() > jet_pt and fabs(jet2.Eta()) < jet_eta;
   bool lep1_ok = abs(PDG1) == 11 ? (lep1.Pt() > e_pt and fabs(lep1.Eta()) < e_eta)
                                 : (lep1.Pt() > mu_pt and fabs(lep1.Eta()) < mu_eta);
   bool lep2_ok = abs(PDG2) == 11 ? (lep2.Pt() > e_pt and fabs(lep2.Eta()) < e_eta)
                                 : (lep2.Pt() > mu_pt and fabs(lep2.Eta()) < mu_eta);
   bool Zpeak_ok = abs(PDG1) != abs(PDG2) or ((lep1+lep2).M() < 60 or (lep1+lep2).M() > 120);

   return (met_ok and jet1_ok and jet2_ok and lep1_ok and lep2_ok and Zpeak_ok );
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
   //std::vector<JetCorrectionUncertainty*> vsrc(jsystnames.size());
   std::string sgn [2] = {"DN","UP"};

   // jet energy scale
   for (unsigned int isrc = 0; isrc < jsystnames.size(); isrc++) {

      const char *name = jsystnames[isrc].c_str();

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

            // separate Flavor uncertainties
            int flavor = jet->partonFlavour();
            if( jsystnames[isrc] == "JESFlavorPureGluon" ){
               if( fabs(flavor) != 21 ) var = 0.0;
            }
            if( jsystnames[isrc] == "JESFlavorPureQuark" ){
               if( fabs(flavor) != 1 and fabs(flavor) != 2 and fabs(flavor) != 3 ) var = 0.0;
            }
            if( jsystnames[isrc] == "JESFlavorPureCharm" ){
               if( fabs(flavor) != 4 ) var = 0.0;
            }
            if( jsystnames[isrc] == "JESFlavorPureBottom" ){
               if( fabs(flavor) != 5 ) var = 0.0;
            }

            bool isJet1 = jetFourVector.Pt() == jet1syst[name].Pt();
            bool isJet2 = jetFourVector.Pt() == jet2syst[name].Pt();

            // correct met first
            metsyst[name] -= jetFourVector*(var/factor);

            jetFourVector *= 1+(var/factor);
            if( isJet1 ) jet1syst[name] = jetFourVector;
            if( isJet2 ) jet2syst[name] = jetFourVector;

         } // jet loop

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

         if( name.find("Central") != std::string::npos or name.find("JetEnergyResolution") != std::string::npos ){
            // do nothing
         }

      }

   }

   return;
}

// ------------ method called once each job just before starting event loop  ------------
   void 
MakeNtuple::beginJob()
{

   // pdfs
   //LHAPDF::initPDFSet(1, "ct10nnlo.LHgrid");
   LHAPDF::initPDFSet(1, "CT10.LHgrid");

   // pile up reweighting
   std::vector< float > PU2012_MC;
   std::vector< float > PU2012_Data;
   std::vector< float > PU2012_DataUP;
   std::vector< float > PU2012_DataDN;

   for( int i=0; i<60; i++) {
      PU2012_MC.push_back( PU2012_MCf[i] );
      PU2012_Data.push_back( PU2012_Dataf[i] );
      PU2012_DataUP.push_back( PU2012_DatafUP[i] );
      PU2012_DataDN.push_back( PU2012_DatafDN[i] );
   }
   LumiWeights_ = edm::LumiReWeighting( PU2012_MC, PU2012_Data);
   LumiWeightsUP_ = edm::LumiReWeighting( PU2012_MC, PU2012_DataUP);
   LumiWeightsDN_ = edm::LumiReWeighting( PU2012_MC, PU2012_DataDN);

   // muon scale factors
   fmueff = new TFile("data/MuonEfficiencies_Run2012ReReco_53X.root");

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

   trees["Central"]->Branch("metFourVector", &metFourVector);
   trees["Central"]->Branch("metUnclustered", &metUnclustered);

   trees["Central"]->Branch("PDG1", &PDG1);
   trees["Central"]->Branch("PDG2", &PDG2);
   trees["Central"]->Branch("vertices", &vertices);
   trees["Central"]->Branch("numBJets", &numBJets);
   trees["Central"]->Branch("jetnumber",&jetcount);
   trees["Central"]->Branch("jet1Vz",&jet1Vz);
   trees["Central"]->Branch("jet2Vz",&jet2Vz);
   trees["Central"]->Branch("jet1bdisc",&jet1bdisc);
   trees["Central"]->Branch("jet2bdisc",&jet2bdisc);

   trees["Central"]->Branch("nmuons", &nmuons);
   trees["Central"]->Branch("nelectrons", &nelectrons);

   if (runOnMC_) {
      trees["Central"]->Branch("puTnvtx", &T_nvertices);
      trees["Central"]->Branch("weight_pu", &weight_pu);
      trees["Central"]->Branch("weight_toppt", &weight_toppt);
      trees["Central"]->Branch("weight_btag", &weight_btag);
      trees["Central"]->Branch("weight_mu", &weight_mu);
      trees["Central"]->Branch("weight_elec", &weight_elec);
      trees["Central"]->Branch("weight_bfrag", &weight_bfrag);

      trees["Central"]->Branch("geninfo_scalePDF", &geninfo_scalePDF);
      trees["Central"]->Branch("geninfo_id1", &geninfo_id1);
      trees["Central"]->Branch("geninfo_id2", &geninfo_id2); 
      trees["Central"]->Branch("geninfo_x1", &geninfo_x1);
      trees["Central"]->Branch("geninfo_x2", &geninfo_x2); 
      trees["Central"]->Branch("geninfo_xPDF1", &geninfo_xPDF1); 
      trees["Central"]->Branch("geninfo_xPDF2", &geninfo_xPDF2);

      trees["Central"]->Branch("pdf_weights", &pdf_weights);
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

   systnames.push_back("Central");
   systnames.push_back("JetEnergyResolution");
   systnames.push_back("METUnclustered");
   systnames.push_back("PileUp");
   systnames.push_back("ElectronEnergyScale");
   systnames.push_back("MuonMomentumScale");
   systnames.push_back("ElectronId");
   systnames.push_back("MuonId");
   systnames.push_back("BTagging");
   systnames.push_back("PDF");
   systnames.push_back("PtTopReweighting");
   systnames.push_back("BFRAGnu");
   systnames.push_back("BFRAGnu");
   systnames.push_back("BFRAGrbLEP");

   jsystnames.push_back("JESCorrelationGroupMPFInSitu");
   jsystnames.push_back("JESCorrelationGroupFlavor");
   jsystnames.push_back("JESCorrelationGroupIntercalibration");
   jsystnames.push_back("JESCorrelationGroupUncorrelated");
   jsystnames.push_back("JESCorrelationGroupbJES");

   jsystnames.push_back("JESAbsoluteStat");
   jsystnames.push_back("JESAbsoluteScale");
   jsystnames.push_back("JESAbsoluteMPFBias");
   jsystnames.push_back("JESAbsoluteFlavMap");
   jsystnames.push_back("JESFlavorPureGluon");
   jsystnames.push_back("JESFlavorPureQuark");
   jsystnames.push_back("JESFlavorPureCharm");
   jsystnames.push_back("JESFlavorPureBottom");
   jsystnames.push_back("JESRelativeFSR");
   jsystnames.push_back("JESTotal");

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
      if( *syst == "Central" ) continue;
      for(int i=0; i < 2; i++){
         std::string name = *syst+sgn[i];
         trees[name] = trees["Central"]->CloneTree(0);
         trees[name]->SetName(name.c_str());
         trees[name]->SetTitle(name.c_str());
      }
   }

   // uncertainty sources
   vsrc.resize( jsystnames.size() );
   for (unsigned int isrc = 0; isrc < jsystnames.size(); isrc++) {

      std::string nametemp = jsystnames[isrc];
      nametemp.erase( 0, 3 ); // erase 'JES'

      const char *nametempc = nametemp.c_str();
      JetCorrectorParameters *p = new JetCorrectorParameters("data/Winter14_V8_DATA_UncertaintySources_AK5PFchs.txt", nametempc);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc[isrc] = unc;

      delete p;
   }

}


// ------------ method called once each job just after ending the event loop  ------------
   void 
MakeNtuple::endJob() 
{
   file->Write();
   file->Close();
   fmueff->Close();
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
