//maximum number of particles and jets saved.
const Int_t NMuons=100;
const Int_t NElectrons=100;
const Int_t NParticles = 50000;
const Int_t NJets = 10000;  //for Wenu study it is more than enough
const Int_t NMETs = 10;
const Int_t NVertices = 100;
const Int_t NNeutrinos = 100;
const Int_t NBTags = 10000;

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

struct RecoElectron{
   int size;
   int charge[NElectrons];
   float pt[NElectrons];
   float p[NElectrons];
   float e[NElectrons];
   float phi[NElectrons];
   float eta[NElectrons];
   float px[NElectrons];
   float py[NElectrons];
   float pz[NElectrons];
   float supercluster_eta[NElectrons];
   bool IDveto[NElectrons];
   bool IDloose[NElectrons];
   bool IDmedium[NElectrons];
   bool IDtight[NElectrons];
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
   float dxx[NMETs];
   float dxy[NMETs];
   float dyy[NMETs];
   float sig[NMETs];
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

   float jetres_par0[NJets];
   float jetres_par1[NJets];
   float jetres_par2[NJets];
   float jetres_par3[NJets];
   float jetres_par4[NJets];
   float jetres_par5[NJets];
   float jetres_par6[NJets];
};

struct PFCandidates{
   Int_t size;
   Int_t type[NParticles];
   float pt[NParticles];
   float et[NParticles];
   float phi[NParticles];
   float eta[NParticles];
   Int_t jetIndex[NParticles];
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
   float nu_px[NJets];
   float nu_py[NJets];
   float nu_pt[NJets];
};

struct BTags{
   Int_t size;
   float pt[NBTags];
   float phi[NBTags];
   float eta[NBTags];
   float energy[NBTags];
   float discriminator[NBTags];
};

struct GenW{
   int id;
   float pt;
   float eta;
   float phi;
   float energy;
   int l_id;
   float l_pt;
   float l_eta;
   float l_phi;
   float l_energy;
   int nu_id;
   float nu_pt;
   float nu_eta;
   float nu_phi;
   float nu_energy;
};

struct GenNu{
   int size;
   int id[NNeutrinos];
   float pt[NNeutrinos];
   float eta[NNeutrinos];
   float phi[NNeutrinos];
   float energy[NNeutrinos];
   int status[NNeutrinos];
};

struct GenMu{
   int size;
   float pt[NMuons];
   float eta[NMuons];
   float phi[NMuons];
   float energy[NMuons];
   int status[NMuons];
};

struct GenInfo {
   Int_t pid;
   float pthat;
   float alphaQCD;
   float alphaQED;
   float scalePDF;
   Int_t id1;
   Int_t id2;
   float x1;
   float x2;
   float xPDF1;
   float xPDF2;
   float weight;
   float xsec;
   float eff;
};
