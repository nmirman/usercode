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
// Original Author:  Aleko Khukhunaishvili,6 R-029,+41227678914,
//         Created:  Sun Aug 14 15:59:00 CEST 2011
// $Id$
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"


#include "TMatrixD.h"
#include "TRandom3.h"
#include "TMath.h"


const double pi = TMath::Pi();

class METSignificance : public edm::EDProducer {
    public:
	explicit METSignificance(const edm::ParameterSet&);
	~METSignificance();

    private:
	virtual void produce(edm::Event&, const edm::EventSetup&);

	void addObject(double pt, double phi, double dpt, double dpphi);
	void addPFCandidate(reco::PFCandidatePtr pf);
	void addPFJet(const reco::PFJet *jet);
	void useOriginalPtrs(const edm::ProductID& productID);
	double getEnergyResolutionEm(double CorrectedEnergy, double eta);
	double getEnergyResolutionHad(double energyHCAL, double eta, double phi);

	double set_;
	double mex_;
	double mey_;
	double Cxx_;
	double Cxy_;
	double Cyy_;

	std::set<reco::CandidatePtr> clusteredParticlePtrs_;

	//PFEnergyResolution *pfresol_; //will get electron, photon and neutral hadron resolutions from particle flow tools
	JetResolution *ptResol_;
	JetResolution *phiResol_;
	TF1* fPtEtaBinned[10]; 
	TF1* fPhiEtaBinned[10];
	double alpha[10]; //σ(p_T)=α[i]*sqrt(pT) for low pT jets
	double beta[10]; //σ(p_φ)=β[ι]*sqrt(pT) for low pT jets

	std::vector<std::string> type_;
	std::vector<edm::InputTag> src_;

	bool smearMet_;
	bool scaleResolutions_;
	double smearJetPtThreshold_;
	double resolJetPtThreshold_;
	double resolJetPhThreshold_;
	double scaleNeutralResols_;
	double jetParticleBoundaryPt_;
	double clusterAllParticlesAfterEta_;

	TRandom3 *random;

	int nClustered;
	int nUnclustered;
};

METSignificance::METSignificance(const edm::ParameterSet& iConfig)
{
    using namespace std;

    type_ = iConfig.getParameter<std::vector<std::string> >("type");
    src_  = iConfig.getParameter<std::vector<edm::InputTag> >("src");

    smearMet_	         = iConfig.getParameter<bool>("smearMet");
    scaleResolutions_	 = iConfig.getParameter<bool>("scaleResolutions");
    smearJetPtThreshold_ = iConfig.getParameter<double>("smearJetPtThreshold");
    resolJetPtThreshold_ = iConfig.getParameter<double>("resolJetPtThreshold");
    resolJetPhThreshold_ = iConfig.getParameter<double>("resolJetPhThreshold");
    scaleNeutralResols_  = iConfig.getParameter<double>("scaleNeutralResols");
    jetParticleBoundaryPt_  = iConfig.getParameter<double>("jetParticleBoundaryPt");
    clusterAllParticlesAfterEta_  = iConfig.getParameter<double>("clusterAllParticlesAfterEta");

    //jet resolutions
    string alg  = iConfig.getParameter<std::string>("jetResolAlgo");     
    string era  = iConfig.getParameter<std::string>("jetResolEra");     
    string path = "CondFormats/JetMETObjects/data";
    string ptFileName  = path + "/" + era + "_PtResolution_" +alg+".txt";
    string phiFileName = path + "/" + era + "_PhiResolution_"+alg+".txt";
    edm::FileInPath fpt(ptFileName);
    edm::FileInPath fphi(phiFileName);
    ptResol_ = new JetResolution(fpt.fullPath().c_str(),false);
    phiResol_ = new JetResolution(fphi.fullPath().c_str(),false);

    //this is just to make it much faster. Evaluate each of the 9(10) functions at the beginning
    float jetapt[10] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75};
    for(int i=0; i<10; ++i){ 
	fPtEtaBinned[i] = ptResol_->parameterEta("sigma", jetapt[i]);
	alpha[i] =  sqrt(resolJetPtThreshold_)*fPtEtaBinned[i]->Eval(resolJetPtThreshold_);
    }
    float jetaphi[10]= {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75};
    for(int i=0; i<10; ++i){ 
	fPhiEtaBinned[i] =  phiResol_->parameterEta("sigma", jetaphi[i]);
	beta[i] = sqrt(resolJetPhThreshold_)*fPhiEtaBinned[i]->Eval(resolJetPhThreshold_);
    }

    random = new TRandom3(0); //random seed for now, don't forget to fix this...
    produces<std::vector<reco::MET> >("");
}

METSignificance::~METSignificance(){}

void
METSignificance::addObject(double pt, double phi, double dpt, double dpphi){
    
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    double Ctt = dpt*dpt;
    double Cff = dpphi*dpphi;

    Cxx_ += Ctt*cosphi*cosphi + Cff*sinphi*sinphi;
    Cxy_ += cosphi*sinphi*(Ctt-Cff);
    Cyy_ += Cff*cosphi*cosphi + Ctt*sinphi*sinphi;
}


void
METSignificance::addPFCandidate(reco::PFCandidatePtr pf){
    double pt = pf->pt(); ////note, different from the standard one where et's are summed up
    double phi = pf->phi();
    int id = pf->particleId();
    double dpt=999;
    double dpphi=999;

    set_ += pt;
    mex_ -= pt*cos(phi);
    mey_ -= pt*sin(phi);

    if(clusteredParticlePtrs_.find(pf) != clusteredParticlePtrs_.end()) return; //pf candidate already added from jet collection

    nUnclustered++;

    double eta = pf->eta();
    double energy = pf->energy();

    if(id==1){ //charged hadrons
	reco::TrackRef trackRef = pf->trackRef();
	double dpttrack=999999;
	if(!trackRef.isNull()){
	    dpphi = pt*trackRef->phiError();
	    dpttrack=trackRef->ptError();
	}
	double dpthcal =pt*getEnergyResolutionHad(energy, eta, phi); 
	dpt = dpttrack<dpthcal ? dpttrack : dpthcal; //this is to add the protection when track energy errors are huge (resulting in a crash when calling significance())
    }
    else if(id==2){ //electrons
	reco::GsfTrackRef gsfTrackRef = pf->gsfTrackRef();
	reco::TrackRef trackRef = pf->trackRef();
	if(!gsfTrackRef.isNull()) dpphi = pt*gsfTrackRef->phiError();
	else if(!trackRef.isNull()) dpphi = pt*trackRef->ptError();
	dpt = getEnergyResolutionEm(energy, eta)/cosh(eta);
    }
    else if(id==3){ //muons
	reco::TrackRef trackRef = pf->trackRef();
	if(!trackRef.isNull()){
	    dpphi = pt*trackRef->phiError();
	    dpt   = trackRef->ptError();
	}
    }
    else if(id==4){ //photons
	dpt = getEnergyResolutionEm(energy, eta)/cosh(eta); //ignoring eta resolutins...
	dpphi=0;		 //ignoring phi resolutions... 
    }
    else if(id==5){ // neutral hadrons
	dpt = pt*getEnergyResolutionHad(energy, eta, phi)*scaleNeutralResols_; //ignoring eta resolutions... scaling by some (to be tuned)factor 
	dpphi=0;		 //ignoring phi resolutions...
    }

    else{ //forward deposits
	double feta = fabs(eta);
	int ietapt = feta<4.5? feta/0.5 : 9;
	dpt = alpha[ietapt]*sqrt(pt);
	dpphi = 0;		 
    }

    addObject(pt, phi, dpt, dpphi);
}


void
METSignificance::addPFJet(const reco::PFJet *jet){

    std::vector<reco::PFCandidatePtr> pfs = jet->getPFConstituents();
    for(std::vector<reco::PFCandidatePtr>::const_iterator it=pfs.begin(); it!=pfs.end(); ++it){
	reco::CandidatePtr ptr(*it);
	clusteredParticlePtrs_.insert(ptr);
	nClustered++;
    }

    double jpt  = jet->pt();
    double jphi = jet->phi();
    double jeta = jet->eta();
    double feta = fabs(jeta);
    double jdeltapt = 999.;
    double jdeltapphi = 999.;
    double c = cos(jphi);
    double s = sin(jphi);
    int ieta =feta<4.5? feta/0.5 : 9;

    jdeltapt   = jpt>resolJetPtThreshold_ ? jpt*fPtEtaBinned[ieta]->Eval(jpt)   : alpha[ieta]*sqrt(jpt);
    jdeltapphi = jpt>resolJetPhThreshold_ ? jpt*fPhiEtaBinned[ieta]->Eval(jpt)  : beta[ieta]*sqrt(jpt);

    double factor=1.0;
    if(jpt>smearJetPtThreshold_){		    
	if(feta<0.5) factor=1.052;
	else if(feta<1.1) factor=1.057;
	else if(feta<1.7) factor=1.096;
	else if(feta<2.3) factor=1.134;
	else factor = 1.288;
    }

    if(smearMet_){ //false for data, true for MC
	double newSigma = sqrt(factor*factor-1)*jdeltapt;
	double newjpt = random->Gaus(jpt, newSigma);
	if(newjpt<0) newjpt=0;
	mex_ = mex_ + (jpt-newjpt)*c;
	mey_ = mey_ + (jpt-newjpt)*s;
    }

    if(scaleResolutions_) jdeltapt*=factor; //resolution is increased for both, MC and data

    addObject(jpt, jphi, jdeltapt, jdeltapphi);
}

void 
METSignificance::useOriginalPtrs(const edm::ProductID& productID){
    std::set<reco::CandidatePtr>::const_iterator it=clusteredParticlePtrs_.begin();
    reco::CandidatePtr ptr(*it);
    if(ptr.id()==productID) return; //If the first element is from the right product, return

    std::set<reco::CandidatePtr> temp;
    for(; it!=clusteredParticlePtrs_.end(); ++it){
	reco::CandidatePtr ptr(*it);
	while(ptr.id()!=productID){
	    ptr = ptr->sourceCandidatePtr(0);
	    if(ptr.isNull()) return; //if it does not get to the correct product, return
	}
	temp.insert(ptr);
    }
    clusteredParticlePtrs_.clear();
    clusteredParticlePtrs_ = temp;
}

double 
METSignificance::getEnergyResolutionEm(double CorrectedEnergy, double eta){
  double C;
  double S;
  double N;
  if(TMath::Abs(eta)<1.48){C=0.35/100; S=5.51/100; N=98./1000.;}
  else{C=0; S=12.8/100; N=440./1000.;} 
  double result = TMath::Sqrt(C*C*CorrectedEnergy*CorrectedEnergy + S*S*CorrectedEnergy + N*N);
  return result; 
}

double 
METSignificance::getEnergyResolutionHad(double energyHCAL, double eta, double phi){
  return 1.49356/sqrt(energyHCAL) + 6.62527e-03*sqrt(energyHCAL) - 6.33966e-02;
}


void
METSignificance::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    set_ = 0;
    mex_ = 0;
    mey_ = 0;
    Cxx_ = 0;
    Cxy_ = 0;
    Cyy_ = 0;
    clusteredParticlePtrs_.clear();

    nClustered=0;
    nUnclustered=0;

    for(unsigned int i=0; i<src_.size(); ++i){
	if(type_[i]=="PFJet"){
	    Handle<reco::PFJetCollection> jets;
	    iEvent.getByLabel(src_[i], jets);
	    for (reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it) {
		if(it->pt()<jetParticleBoundaryPt_ && fabs(it->eta())<clusterAllParticlesAfterEta_) continue;
		addPFJet(&(*it));
	    }
	}
	else if(type_[i]=="PFCandidate"){
	    edm::Handle<edm::View<reco::PFCandidate> > particles;
	    iEvent.getByLabel(src_[i], particles);

	    useOriginalPtrs(particles.id());

	    for (edm::View<reco::PFCandidate>::const_iterator it = particles->begin(); it != particles->end(); ++it) {
		reco::PFCandidatePtr dau(particles, it - particles->begin());
		if(dau.isNonnull () && dau.isAvailable()){
		    reco::PFCandidatePtr pf(dau.id(), &(*it), dau.key());
		    addPFCandidate(pf);
		}
	    }
	}
    }

    reco::MET::LorentzVector ptvec(mex_, mey_, 0, sqrt(mex_*mex_+mey_*mey_));
    reco::MET met(set_, ptvec, reco::MET::Point());

    TMatrixD RM(2,2);
    RM(0,0) = Cxx_;
    RM(0,1) = Cxy_;
    RM(1,0) = Cxy_;
    RM(1,1) = Cyy_;

    met.setSignificanceMatrix(RM);

    /*
    edm::Handle<edm::View<reco::PFCandidate> > particles;
    iEvent.getByLabel("particleFlow", particles);
    std::cout << particles->size() << "      " << nClustered << "      " << nUnclustered << "      " << (particles->size()-nClustered-nUnclustered) << std::endl;
    
    edm::Handle<edm::View<reco::PFMET> > pfmet;
    iEvent.getByLabel("pfMET", pfmet);
    double significance = (pfmet->front() ).significance();
    double sigmaX2= (pfmet->front() ).getSignificanceMatrix()(0,0);
    double sigmaXY= (pfmet->front() ).getSignificanceMatrix()(0,1);
    double sigmaY2= (pfmet->front() ).getSignificanceMatrix()(1,1);
    double metet = (pfmet->front() ).pt();
    std::cout << sigmaX2 << " " <<  Cxx_ << "     " <<
		 sigmaY2 << " " <<  Cyy_ << "     " << 	
		 sigmaXY << " " <<  Cxy_ << "     " << 	
		 significance << " " <<  met.significance() << "     " << 	
		 metet << " " <<  met.pt() << std::endl; 
    */

    std::auto_ptr<std::vector<reco::MET> > metp(new std::vector<reco::MET>());
    metp->push_back(met);
    iEvent.put(metp);
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(METSignificance);
