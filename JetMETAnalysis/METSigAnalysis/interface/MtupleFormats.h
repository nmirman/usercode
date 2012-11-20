//maximum number of particles and jets saved.
const Int_t NMuons=1000;
const Int_t NParticles = 50000;
const Int_t NJets = 10000;  //for Wenu study it is more than enough
const Int_t NMETs = 100;
const Int_t NVertices = 100;

struct METs{
    Int_t size;
    float et[NMETs];
    float phi[NMETs];
    float sig[NMETs];
    float sumEt[NMETs];
    float dxx[NMETs];
    float dxy[NMETs];
    float dyy[NMETs];
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
    float dpt[NJets];
    float dphi[NJets];
    float neutralHadronFraction[NJets];
    float neutralEmFraction[NJets];
    float chargedHadronFraction[NJets];
    float chargedHadronMultiplicity[NJets];
    float chargedEmFraction[NJets];
};

struct GenJets{
    Int_t size;
    float pt[NJets];
    float phi[NJets];
    float eta[NJets];
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

struct GenW{
    int id;
    int eid;
    int nid;
    double M;
    double we;
    double wpx;
    double wpy;
    double wpz;
    double wpt;
    double weta;
    double wphi;
    double wy;
    double m;
    double ee;
    double epx;
    double epy;
    double epz;
    double ept;
    double eeta;
    double ephi;
    double ne;
    double npx;
    double npy;
    double npz;
    double npt;
    double neta;
    double nphi;
};

struct Vertices{
    Int_t size;
    Bool_t isFake[NVertices];
    float ndof[NVertices];
    float x[NVertices];
    float y[NVertices];
    float z[NVertices];
    float Rho[NVertices];
};


struct RecoMuon{
    int size;
    int charge[NMuons];
    int isGlobal[NMuons];
    int isTracker[NMuons];
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
    float dr03TkSumPt[NMuons];
    float dr03EcalRecHitSumEt[NMuons];
    float dr03HcalTowerSumEt[NMuons];
    float drTO[20][NMuons];


    float gdr[NMuons];
    int   mid[NMuons];
    float mpx[NMuons];
    float mpy[NMuons];
    float mpz[NMuons];
    float men[NMuons];
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
