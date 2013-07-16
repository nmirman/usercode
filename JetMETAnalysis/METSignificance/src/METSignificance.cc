// -*- C++ -*-
//
// Package:    METSignificance
// Class:      METSignificance
// 
/**\class METSignificance METSignificance.cc JetMETAnalysis/METSignificance/src/METSignificance.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  nathan mirman
//         Created:  Thu May 30 16:39:52 CDT 2013
// $Id: METSignificance.cc,v 1.4 2013/07/16 21:12:45 nmirman Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

//
// class declaration
//

class METSignificance : public edm::EDProducer {
   public:
      explicit METSignificance(const edm::ParameterSet&);
      ~METSignificance();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      
      edm::InputTag muonTag_;
      edm::InputTag electronTag_;
      edm::InputTag pfjetsTag_;

      std::string pfjetCorrectorL1_;
      std::string pfjetCorrectorL123_;

      bool runOnMC_;

      double jetThreshold_;
      double parA1;
      double parA2;
      double parA3;
      double parA4;
      double parA5;
      double parN1;
      double parS1;

      JetResolution *ptRes_;
      JetResolution *phiRes_;
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
METSignificance::METSignificance(const edm::ParameterSet& iConfig)
{
   muonTag_    = iConfig.getUntrackedParameter<edm::InputTag>("muonTag");
   electronTag_ = iConfig.getUntrackedParameter<edm::InputTag>("electronTag");

   pfjetsTag_    = iConfig.getUntrackedParameter<edm::InputTag>("pfjetsTag");
   pfjetCorrectorL1_  = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL1");
   pfjetCorrectorL123_ = iConfig.getUntrackedParameter<std::string>("pfjetCorrectorL123");

   runOnMC_ = iConfig.getUntrackedParameter<bool>("runOnMC");

   jetThreshold_ = iConfig.getParameter<double>("jetThreshold");
   if( runOnMC_ ){
      parA1 = 1.12660;
      parA2 = 1.09322;
      parA3 = 1.10951;
      parA4 = 1.17178;
      parA5 = 1.12164;
      parN1 = 0.0;
      parS1 = 0.585145;
   } else {
      parA1 = 1.39669;
      parA2 = 1.32037;
      parA3 = 1.32047;
      parA4 = 1.38161;
      parA5 = 1.51508;
      parN1 = 0.0;
      parS1 = 0.639158;
   }

   std::string alg  = iConfig.getParameter<std::string>("jetResAlgo");
   std::string era  = iConfig.getParameter<std::string>("jetResEra");

   std::string path = "CondFormats/JetMETObjects/data";
   std::string ptFileName  = path + "/" + era + "_PtResolution_" +alg+".txt";
   std::string phiFileName = path + "/" + era + "_PhiResolution_"+alg+".txt";

   edm::FileInPath fpt(ptFileName);
   edm::FileInPath fphi(phiFileName);

   ptRes_  = new JetResolution(fpt.fullPath().c_str(),false);
   phiRes_ = new JetResolution(fphi.fullPath().c_str(),false);

   produces<double>("METSignificance");
   produces<double>("CovarianceMatrix00");
   produces<double>("CovarianceMatrix01");
   produces<double>("CovarianceMatrix10");
   produces<double>("CovarianceMatrix11");
}


METSignificance::~METSignificance()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
METSignificance::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //
   // met and covariance
   //
   double met_px = 0;
   double met_py = 0;

   double cov_xx = 0;
   double cov_xy = 0;
   double cov_yy = 0;

   
   //
   // muons
   //
   Handle<std::vector<reco::PFCandidate> > muonsHandle;
   iEvent.getByLabel(muonTag_, muonsHandle);
   std::vector<reco::PFCandidate> muons = *muonsHandle;
   for(std::vector<reco::PFCandidate>::const_iterator it=muons.begin(); it!=muons.end(); it++){
      if( !(it->muonRef().isNonnull() and it->muonRef().isAvailable()) ) continue;
      reco::PFCandidate pfmuon = *it;
      met_px -= pfmuon.px();
      met_py -= pfmuon.py();
   }

   //
   // electrons
   //
   Handle<std::vector<reco::PFCandidate> > electronsHandle;
   iEvent.getByLabel(electronTag_, electronsHandle);
   std::vector<reco::PFCandidate> electrons = *electronsHandle;
   for(std::vector<reco::PFCandidate>::const_iterator it=electrons.begin(); it!=electrons.end(); it++){
      reco::PFCandidate pfelectron = *it;
      met_px -= pfelectron.px();
      met_py -= pfelectron.py();
   }

   //
   // jets
   //
   Handle<reco::PFJetCollection> inputUncorJets;
   iEvent.getByLabel( pfjetsTag_, inputUncorJets );

   // jet energy corrections
   const JetCorrector* corrector  = JetCorrector::getJetCorrector (pfjetCorrectorL1_, iSetup);
   const JetCorrector* corrector2 = JetCorrector::getJetCorrector (pfjetCorrectorL123_, iSetup);

   // pseudo-jet initialization
   double pjet_px = 0;
   double pjet_py = 0;
   double pjet_scalpt = 0;

   for(reco::PFJetCollection::const_iterator jet = inputUncorJets->begin(); jet != inputUncorJets->end(); ++jet) {
      double jpt  = jet->pt();
      double jeta = jet->eta();
      double feta = fabs(jeta);
      double c = cos(jet->phi());
      double s = sin(jet->phi());

      double jcorrl1 = corrector->correction (*jet, iEvent, iSetup);
      double jcorrl123 = corrector2->correction (*jet, iEvent, iSetup);

      // corrected jet pt's -- corrections apply for jets w/ corrected pt > 10 GeV
      double jptL123 = (jpt*jcorrl123 > 10) ? jpt*jcorrl123 : jpt;
      double jptT1 = (jpt*jcorrl123 > 10) ? jpt*(jcorrl123+1-jcorrl1) : jpt;

      // jet energy resolutions
      double jeta_res = (fabs(jeta) < 9.9) ? jeta : 9.89; // JetResolutions defined for |eta|<9.9
      TF1* fPtEta    = ptRes_ -> parameterEta("sigma",jeta_res);
      TF1* fPhiEta   = phiRes_-> parameterEta("sigma",jeta_res);
      double sigmapt = fPtEta->Eval(jptL123);
      double sigmaphi = fPhiEta->Eval(jptL123);
      delete fPtEta;
      delete fPhiEta;

      met_px -= c*jptT1;
      met_py -= s*jptT1;

      // split into high-pt and low-pt sector
      if( jptL123 > jetThreshold_ ){
         // high-pt jets enter into the covariance matrix via JER
         
         double scale = 0;
         if(feta<0.5) scale = parA1;
         else if(feta<1.1) scale = parA2;
         else if(feta<1.7) scale = parA3;
         else if(feta<2.3) scale = parA4;
         else scale = parA5;

         double dpt = scale*jptT1*sigmapt;
         double dph = jptT1*sigmaphi;

         cov_xx += dpt*dpt*c*c + dph*dph*s*s;
         cov_xy += (dpt*dpt-dph*dph)*c*s;
         cov_yy += dph*dph*c*c + dpt*dpt*s*s;

      }else{
         // low-pt jets are lumped into the pseudo-jet

         pjet_px += jptL123*c;
         pjet_py += jptL123*s;
         pjet_scalpt += jptL123;

      }

   }

   // contribution to covariance matrix from pseudo-jet
   cov_xx += parN1*parN1 + parS1*parS1*pjet_scalpt;
   cov_yy += parN1*parN1 + parS1*parS1*pjet_scalpt;

   //
   // compute significance
   //
   double det = cov_xx*cov_yy - cov_xy*cov_xy;

   double ncov_xx = cov_yy / det;
   double ncov_xy = -cov_xy / det;
   double ncov_yy = cov_xx / det;

   double sig = met_px*met_px*ncov_xx + 2*met_px*met_py*ncov_xy + met_py*met_py*ncov_yy;

   std::auto_ptr<double> significance (new double);
   (*significance) = sig;

   std::auto_ptr<double> sigmatrix_00 (new double);
   (*sigmatrix_00) = cov_xx;
   std::auto_ptr<double> sigmatrix_01 (new double);
   (*sigmatrix_01) = cov_xy;
   std::auto_ptr<double> sigmatrix_10 (new double);
   (*sigmatrix_10) = cov_xy;
   std::auto_ptr<double> sigmatrix_11 (new double);
   (*sigmatrix_11) = cov_yy;

   iEvent.put( significance, "METSignificance" );
   iEvent.put( sigmatrix_00, "CovarianceMatrix00" );
   iEvent.put( sigmatrix_01, "CovarianceMatrix01" );
   iEvent.put( sigmatrix_10, "CovarianceMatrix10" );
   iEvent.put( sigmatrix_11, "CovarianceMatrix11" );
}

// ------------ method called once each job just before starting event loop  ------------
void 
METSignificance::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METSignificance::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
METSignificance::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
METSignificance::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
METSignificance::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
METSignificance::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METSignificance::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(METSignificance);
