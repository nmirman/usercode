//maximum number of particles and jets saved.
const Int_t NMuons=1000;
const Int_t NParticles = 50000;
const Int_t NJets = 10000;  //for Wenu study it is more than enough
const Int_t NMETs = 100;
const Int_t NVertices = 100;


struct RecoMuon{
   int size;
   int charge[NMuons];
   bool isGlobal[NMuons];
   bool isTracker[NMuons];
   bool isPF[NMuons];
   int trackerHits[NMuons];
   int pixelHits[NMuons];
   int muonHits[NMuons];
   int nMatches[NMuons];
   int numberOfValidHits[NMuons];
   int numberOfValidPixelHits[NMuons];
   int numberOfValidStripHits[NMuons];
   int numberOfValidPixelBarrelHits[NMuons];
   int pixelLayersWithMeasurement[NMuons];
   int pixelBarrelLayersWithMeasurement[NMuons];
   int pixelEndcapLayersWithMeasurement[NMuons];
   int numberOfValidTrackerLayers[NMuons];
   float chi2[NMuons];
   float pt[NMuons];
   float glpt[NMuons];
   float trpt[NMuons];
   float p[NMuons];
   float e[NMuons];
   float phi[NMuons];
   float eta[NMuons];
   float px[NMuons];
   float py[NMuons];
   float pz[NMuons];
   float dxy[NMuons];
   float dz[NMuons];
   float dr03TkSumPt[NMuons];
   float dr03EcalRecHitSumEt[NMuons];
   float dr03HcalTowerSumEt[NMuons];
   float dr04chHad[NMuons];
   float dr04neutHad[NMuons];
   float dr04photons[NMuons];
};

struct Vertices{
   Int_t size;
   Bool_t isFake[NVertices];
   float ndof[NVertices];
   float chi2[NVertices];
   float x[NVertices];
   float y[NVertices];
   float z[NVertices];
   float Rho[NVertices];
};

struct METs{
   Int_t size;
   float pt[NMETs];
   float px[NMETs];
   float py[NMETs];
   float pz[NMETs];
   float phi[NMETs];
   float sumEt[NMETs];
};

struct PFJets{
   Int_t size;
   Int_t nco[NJets];
   float l1[NJets];
   float l1l2l3[NJets];
   float pt[NJets];
   float phi[NJets];
   float eta[NJets];
   float energy[NJets];
   float sigmapt[NJets];
   float sigmaphi[NJets];
   float neutralHadronFraction[NJets];
   float neutralEmFraction[NJets];
   float chargedHadronFraction[NJets];
   float chargedHadronMultiplicity[NJets];
   float chargedEmFraction[NJets];

   float puid_mva[NJets];
   int puid_idflag[NJets];
   bool puid_passloose[NJets];
   bool puid_passmedium[NJets];
   bool puid_passtight[NJets];
};

struct GenJets{
   Int_t size;
   float pt[NJets];
   float phi[NJets];
   float eta[NJets];
   float energy[NJets];
   float emEnergy[NJets];
   float hadEnergy[NJets];
   float invEnergy[NJets];
   float auxEnergy[NJets];
};

