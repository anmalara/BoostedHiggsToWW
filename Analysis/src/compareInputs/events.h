//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jan 20 13:11:46 2019 by ROOT version 6.10/09
// from TTree events/events
// found on file: NN_HZ.root
//////////////////////////////////////////////////////////

#ifndef events_h
#define events_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class events {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxrun = 1;
   static constexpr Int_t kMaxevt = 1;
   static constexpr Int_t kMaxlumi = 1;
   static constexpr Int_t kMaxnVtx = 1;
   static constexpr Int_t kMaxnJets = 1;
   static constexpr Int_t kMaxnGenJets = 1;
   static constexpr Int_t kMaxnBJets = 1;
   static constexpr Int_t kMaxpvRho = 1;
   static constexpr Int_t kMaxpvz = 1;
   static constexpr Int_t kMaxpvchi2 = 1;
   static constexpr Int_t kMaxpvndof = 1;
   static constexpr Int_t kMaxrho = 1;
   static constexpr Int_t kMaxht = 1;
   static constexpr Int_t kMaxmet = 1;
   static constexpr Int_t kMaxmetGen = 1;
   static constexpr Int_t kMaxmetSig = 1;
   static constexpr Int_t kMaxmetGenSig = 1;
   static constexpr Int_t kMaxnpu = 1;
   static constexpr Int_t kMaxgenEvtWeight = 1;
   static constexpr Int_t kMaxlheOriginalXWGTUP = 1;

   // Declaration of leaf types
   Int_t           runNo;
   Int_t           evtNo;
   Int_t           lumi;
   Int_t           nvtx;
   Int_t           nJets;
   Int_t           nGenJets;
   Int_t           nBJets;
   Float_t         pvRho;
   Float_t         pvz;
   Float_t         pvchi2;
   Float_t         pvndof;
   Float_t         rho;
   Float_t         ht;
   Float_t         met;
   Float_t         metGen;
   Float_t         metSig;
   Float_t         metGenSig;
   vector<bool>    *jetIsBtag;
   vector<int>     *jetFlavor;
   vector<int>     *jetFlavorHadron;
   vector<bool>    *MatchedLeptons;
   vector<bool>    *MatchedHiggs;
   vector<bool>    *WMinusLep;
   vector<bool>    *WPlusLep;
   vector<float>   *jetPt;
   vector<float>   *jetCorr;
   vector<float>   *jetUnc;
   vector<float>   *jetBtag;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetMass;
   vector<float>   *jetEnergy;
   vector<float>   *jetChf;
   vector<float>   *jetChm;
   vector<float>   *jetNpr;
   vector<float>   *jetNhf;
   vector<float>   *jetPhf;
   vector<float>   *jetMuf;
   vector<float>   *jetElf;
   vector<float>   *jetMassSoftDrop;
   vector<float>   *jetEnergy;
   vector<float>   *jetChf;
   vector<float>   *jetNhf;
   vector<float>   *jetPhf;
   vector<float>   *jetMuf;
   vector<float>   *jetElf;
   vector<float>   *jetTau1;
   vector<float>   *jetTau2;
   vector<float>   *jetTau3;
   vector<float>   *jetTau4;
   vector<float>   *btagSub0;
   vector<float>   *btagSub1;
   vector<float>   *massSub0;
   vector<float>   *massSub1;
   vector<float>   *ptSub0;
   vector<float>   *ptSub1;
   vector<float>   *etaSub0;
   vector<float>   *etaSub1;
   vector<float>   *phiSub0;
   vector<float>   *phiSub1;
   vector<int>     *flavorSub0;
   vector<int>     *flavorSub1;
   vector<int>     *flavorHadronSub0;
   vector<int>     *flavorHadronSub1;
   vector<int>     *nSubJets;
   vector<int>     *nBSubJets;
   vector<float>   *CandPdgId;
   vector<float>   *CandPhi;
   vector<float>   *CandEta;
   vector<float>   *CandPt;
   vector<float>   *CandMass;
   vector<float>   *CandPx;
   vector<float>   *CandPy;
   vector<float>   *CandPz;
   vector<float>   *CandEnergy;
   vector<float>   *CandPuppiWeight;
   vector<float>   *CandPvAssociationQuality;
   vector<float>   *CandDXY;
   vector<float>   *CandDZ;
   vector<float>   *CandDZAssociatedPV;
   vector<float>   *CandSubJetPart;
   vector<int>     *ncandidates;
   vector<int>     *nGencandidates;
   vector<bool>    *triggerBit;
   vector<int>     *triggerPre;
   Int_t           npu;
   Float_t         genEvtWeight;
   Float_t         lheOriginalXWGTUP;
   vector<float>   *scaleWeights;
   vector<float>   *pdfWeights;
   vector<float>   *GenJetPt;
   vector<float>   *GenJetEta;
   vector<float>   *GenJetPhi;
   vector<float>   *GenJetEnergy;
   vector<float>   *GenJetMass;
   vector<bool>    *isBJetGen;
   vector<float>   *GenJetTau1;
   vector<float>   *GenJetTau2;
   vector<float>   *GenJetTau3;
   vector<float>   *GenJetTau4;
   vector<float>   *GenSoftDropMass;
   vector<float>   *GenCandPhi;
   vector<float>   *GenCandEta;
   vector<float>   *GenCandPt;
   vector<float>   *GenCandPx;
   vector<float>   *GenCandPy;
   vector<float>   *GenCandPz;
   vector<float>   *GenCandEnergy;
   vector<float>   *GenmassSub0;
   vector<float>   *GenmassSub1;
   vector<float>   *GenptSub0;
   vector<float>   *GenptSub1;
   vector<float>   *GenetaSub0;
   vector<float>   *GenetaSub1;
   vector<float>   *GenphiSub0;
   vector<float>   *GenphiSub1;

   // List of branches
   TBranch        *b_run_;   //!
   TBranch        *b_evt_;   //!
   TBranch        *b_lumi_;   //!
   TBranch        *b_nVtx_;   //!
   TBranch        *b_nJets_;   //!
   TBranch        *b_nGenJets_;   //!
   TBranch        *b_nBJets_;   //!
   TBranch        *b_pvRho_;   //!
   TBranch        *b_pvz_;   //!
   TBranch        *b_pvchi2_;   //!
   TBranch        *b_pvndof_;   //!
   TBranch        *b_rho_;   //!
   TBranch        *b_ht_;   //!
   TBranch        *b_met_;   //!
   TBranch        *b_metGen_;   //!
   TBranch        *b_metSig_;   //!
   TBranch        *b_metGenSig_;   //!
   TBranch        *b_jetIsBtag;   //!
   TBranch        *b_jetFlavor;   //!
   TBranch        *b_jetFlavorHadron;   //!
   TBranch        *b_MatchedLeptons;   //!
   TBranch        *b_MatchedHiggs;   //!
   TBranch        *b_WMinusLep;   //!
   TBranch        *b_WPlusLep;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetCorr;   //!
   TBranch        *b_jetUnc;   //!
   TBranch        *b_jetBtag;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetEnergy;   //!
   TBranch        *b_jetChf;   //!
   TBranch        *b_jetChm;   //!
   TBranch        *b_jetNpr;   //!
   TBranch        *b_jetNhf;   //!
   TBranch        *b_jetPhf;   //!
   TBranch        *b_jetMuf;   //!
   TBranch        *b_jetElf;   //!
   TBranch        *b_jetMassSoftDrop;   //!
   TBranch        *b_jetEnergy;   //!
   TBranch        *b_jetChf;   //!
   TBranch        *b_jetNhf;   //!
   TBranch        *b_jetPhf;   //!
   TBranch        *b_jetMuf;   //!
   TBranch        *b_jetElf;   //!
   TBranch        *b_jetTau1;   //!
   TBranch        *b_jetTau2;   //!
   TBranch        *b_jetTau3;   //!
   TBranch        *b_jetTau4;   //!
   TBranch        *b_btagSub0;   //!
   TBranch        *b_btagSub1;   //!
   TBranch        *b_massSub0;   //!
   TBranch        *b_massSub1;   //!
   TBranch        *b_ptSub0;   //!
   TBranch        *b_ptSub1;   //!
   TBranch        *b_etaSub0;   //!
   TBranch        *b_etaSub1;   //!
   TBranch        *b_phiSub0;   //!
   TBranch        *b_phiSub1;   //!
   TBranch        *b_flavorSub0;   //!
   TBranch        *b_flavorSub1;   //!
   TBranch        *b_flavorHadronSub0;   //!
   TBranch        *b_flavorHadronSub1;   //!
   TBranch        *b_nSubJets;   //!
   TBranch        *b_nBSubJets;   //!
   TBranch        *b_CandPdgId;   //!
   TBranch        *b_CandPhi;   //!
   TBranch        *b_CandEta;   //!
   TBranch        *b_CandPt;   //!
   TBranch        *b_CandMass;   //!
   TBranch        *b_CandPx;   //!
   TBranch        *b_CandPy;   //!
   TBranch        *b_CandPz;   //!
   TBranch        *b_CandEnergy;   //!
   TBranch        *b_CandPuppiWeight;   //!
   TBranch        *b_CandPvAssociationQuality;   //!
   TBranch        *b_CandDXY;   //!
   TBranch        *b_CandDZ;   //!
   TBranch        *b_CandDZAssociatedPV;   //!
   TBranch        *b_CandSubJetPart;   //!
   TBranch        *b_ncandidates;   //!
   TBranch        *b_nGencandidates;   //!
   TBranch        *b_triggerBit;   //!
   TBranch        *b_triggerPre;   //!
   TBranch        *b_npu_;   //!
   TBranch        *b_genEvtWeight_;   //!
   TBranch        *b_lheOriginalXWGTUP_;   //!
   TBranch        *b_scaleWeights;   //!
   TBranch        *b_pdfWeights;   //!
   TBranch        *b_GenJetPt;   //!
   TBranch        *b_GenJetEta;   //!
   TBranch        *b_GenJetPhi;   //!
   TBranch        *b_GenJetEnergy;   //!
   TBranch        *b_GenJetMass;   //!
   TBranch        *b_isBJetGen;   //!
   TBranch        *b_GenJetTau1;   //!
   TBranch        *b_GenJetTau2;   //!
   TBranch        *b_GenJetTau3;   //!
   TBranch        *b_GenJetTau4;   //!
   TBranch        *b_GenSoftDropMass;   //!
   TBranch        *b_GenCandPhi;   //!
   TBranch        *b_GenCandEta;   //!
   TBranch        *b_GenCandPt;   //!
   TBranch        *b_GenCandPx;   //!
   TBranch        *b_GenCandPy;   //!
   TBranch        *b_GenCandPz;   //!
   TBranch        *b_GenCandEnergy;   //!
   TBranch        *b_GenmassSub0;   //!
   TBranch        *b_GenmassSub1;   //!
   TBranch        *b_GenptSub0;   //!
   TBranch        *b_GenptSub1;   //!
   TBranch        *b_GenetaSub0;   //!
   TBranch        *b_GenetaSub1;   //!
   TBranch        *b_GenphiSub0;   //!
   TBranch        *b_GenphiSub1;   //!

   events(TTree *tree=0);
   virtual ~events();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef events_cxx
events::events(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("NN_HZ.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("NN_HZ.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("NN_HZ.root:/boostedAK8");
      dir->GetObject("events",tree);

   }
   Init(tree);
}

events::~events()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t events::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t events::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void events::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jetIsBtag = 0;
   jetFlavor = 0;
   jetFlavorHadron = 0;
   MatchedLeptons = 0;
   MatchedHiggs = 0;
   WMinusLep = 0;
   WPlusLep = 0;
   jetPt = 0;
   jetCorr = 0;
   jetUnc = 0;
   jetBtag = 0;
   jetEta = 0;
   jetPhi = 0;
   jetMass = 0;
   jetEnergy = 0;
   jetChf = 0;
   jetChm = 0;
   jetNpr = 0;
   jetNhf = 0;
   jetPhf = 0;
   jetMuf = 0;
   jetElf = 0;
   jetMassSoftDrop = 0;
   jetEnergy = 0;
   jetChf = 0;
   jetNhf = 0;
   jetPhf = 0;
   jetMuf = 0;
   jetElf = 0;
   jetTau1 = 0;
   jetTau2 = 0;
   jetTau3 = 0;
   jetTau4 = 0;
   btagSub0 = 0;
   btagSub1 = 0;
   massSub0 = 0;
   massSub1 = 0;
   ptSub0 = 0;
   ptSub1 = 0;
   etaSub0 = 0;
   etaSub1 = 0;
   phiSub0 = 0;
   phiSub1 = 0;
   flavorSub0 = 0;
   flavorSub1 = 0;
   flavorHadronSub0 = 0;
   flavorHadronSub1 = 0;
   nSubJets = 0;
   nBSubJets = 0;
   CandPdgId = 0;
   CandPhi = 0;
   CandEta = 0;
   CandPt = 0;
   CandMass = 0;
   CandPx = 0;
   CandPy = 0;
   CandPz = 0;
   CandEnergy = 0;
   CandPuppiWeight = 0;
   CandPvAssociationQuality = 0;
   CandDXY = 0;
   CandDZ = 0;
   CandDZAssociatedPV = 0;
   CandSubJetPart = 0;
   ncandidates = 0;
   nGencandidates = 0;
   triggerBit = 0;
   triggerPre = 0;
   scaleWeights = 0;
   pdfWeights = 0;
   GenJetPt = 0;
   GenJetEta = 0;
   GenJetPhi = 0;
   GenJetEnergy = 0;
   GenJetMass = 0;
   isBJetGen = 0;
   GenJetTau1 = 0;
   GenJetTau2 = 0;
   GenJetTau3 = 0;
   GenJetTau4 = 0;
   GenSoftDropMass = 0;
   GenCandPhi = 0;
   GenCandEta = 0;
   GenCandPt = 0;
   GenCandPx = 0;
   GenCandPy = 0;
   GenCandPz = 0;
   GenCandEnergy = 0;
   GenmassSub0 = 0;
   GenmassSub1 = 0;
   GenptSub0 = 0;
   GenptSub1 = 0;
   GenetaSub0 = 0;
   GenetaSub1 = 0;
   GenphiSub0 = 0;
   GenphiSub1 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNo", &runNo, &b_run_);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evt_);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi_);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nVtx_);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets_);
   fChain->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets_);
   fChain->SetBranchAddress("nBJets", &nBJets, &b_nBJets_);
   fChain->SetBranchAddress("pvRho", &pvRho, &b_pvRho_);
   fChain->SetBranchAddress("pvz", &pvz, &b_pvz_);
   fChain->SetBranchAddress("pvchi2", &pvchi2, &b_pvchi2_);
   fChain->SetBranchAddress("pvndof", &pvndof, &b_pvndof_);
   fChain->SetBranchAddress("rho", &rho, &b_rho_);
   fChain->SetBranchAddress("ht", &ht, &b_ht_);
   fChain->SetBranchAddress("met", &met, &b_met_);
   fChain->SetBranchAddress("metGen", &metGen, &b_metGen_);
   fChain->SetBranchAddress("metSig", &metSig, &b_metSig_);
   fChain->SetBranchAddress("metGenSig", &metGenSig, &b_metGenSig_);
   fChain->SetBranchAddress("jetIsBtag", &jetIsBtag, &b_jetIsBtag);
   fChain->SetBranchAddress("jetFlavor", &jetFlavor, &b_jetFlavor);
   fChain->SetBranchAddress("jetFlavorHadron", &jetFlavorHadron, &b_jetFlavorHadron);
   fChain->SetBranchAddress("MatchedLeptons", &MatchedLeptons, &b_MatchedLeptons);
   fChain->SetBranchAddress("MatchedHiggs", &MatchedHiggs, &b_MatchedHiggs);
   fChain->SetBranchAddress("WMinusLep", &WMinusLep, &b_WMinusLep);
   fChain->SetBranchAddress("WPlusLep", &WPlusLep, &b_WPlusLep);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetCorr", &jetCorr, &b_jetCorr);
   fChain->SetBranchAddress("jetUnc", &jetUnc, &b_jetUnc);
   fChain->SetBranchAddress("jetBtag", &jetBtag, &b_jetBtag);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetEnergy", &jetEnergy, &b_jetEnergy);
   fChain->SetBranchAddress("jetChf", &jetChf, &b_jetChf);
   fChain->SetBranchAddress("jetChm", &jetChm, &b_jetChm);
   fChain->SetBranchAddress("jetNpr", &jetNpr, &b_jetNpr);
   fChain->SetBranchAddress("jetNhf", &jetNhf, &b_jetNhf);
   fChain->SetBranchAddress("jetPhf", &jetPhf, &b_jetPhf);
   fChain->SetBranchAddress("jetMuf", &jetMuf, &b_jetMuf);
   fChain->SetBranchAddress("jetElf", &jetElf, &b_jetElf);
   fChain->SetBranchAddress("jetMassSoftDrop", &jetMassSoftDrop, &b_jetMassSoftDrop);
//    fChain->SetBranchAddress("jetEnergy", &jetEnergy, &b_jetEnergy);
//    fChain->SetBranchAddress("jetChf", &jetChf, &b_jetChf);
//    fChain->SetBranchAddress("jetNhf", &jetNhf, &b_jetNhf);
//    fChain->SetBranchAddress("jetPhf", &jetPhf, &b_jetPhf);
//    fChain->SetBranchAddress("jetMuf", &jetMuf, &b_jetMuf);
//    fChain->SetBranchAddress("jetElf", &jetElf, &b_jetElf);
   fChain->SetBranchAddress("jetTau1", &jetTau1, &b_jetTau1);
   fChain->SetBranchAddress("jetTau2", &jetTau2, &b_jetTau2);
   fChain->SetBranchAddress("jetTau3", &jetTau3, &b_jetTau3);
   fChain->SetBranchAddress("jetTau4", &jetTau4, &b_jetTau4);
   fChain->SetBranchAddress("btagSub0", &btagSub0, &b_btagSub0);
   fChain->SetBranchAddress("btagSub1", &btagSub1, &b_btagSub1);
   fChain->SetBranchAddress("massSub0", &massSub0, &b_massSub0);
   fChain->SetBranchAddress("massSub1", &massSub1, &b_massSub1);
   fChain->SetBranchAddress("ptSub0", &ptSub0, &b_ptSub0);
   fChain->SetBranchAddress("ptSub1", &ptSub1, &b_ptSub1);
   fChain->SetBranchAddress("etaSub0", &etaSub0, &b_etaSub0);
   fChain->SetBranchAddress("etaSub1", &etaSub1, &b_etaSub1);
   fChain->SetBranchAddress("phiSub0", &phiSub0, &b_phiSub0);
   fChain->SetBranchAddress("phiSub1", &phiSub1, &b_phiSub1);
   fChain->SetBranchAddress("flavorSub0", &flavorSub0, &b_flavorSub0);
   fChain->SetBranchAddress("flavorSub1", &flavorSub1, &b_flavorSub1);
   fChain->SetBranchAddress("flavorHadronSub0", &flavorHadronSub0, &b_flavorHadronSub0);
   fChain->SetBranchAddress("flavorHadronSub1", &flavorHadronSub1, &b_flavorHadronSub1);
   fChain->SetBranchAddress("nSubJets", &nSubJets, &b_nSubJets);
   fChain->SetBranchAddress("nBSubJets", &nBSubJets, &b_nBSubJets);
   fChain->SetBranchAddress("CandPdgId", &CandPdgId, &b_CandPdgId);
   fChain->SetBranchAddress("CandPhi", &CandPhi, &b_CandPhi);
   fChain->SetBranchAddress("CandEta", &CandEta, &b_CandEta);
   fChain->SetBranchAddress("CandPt", &CandPt, &b_CandPt);
   fChain->SetBranchAddress("CandMass", &CandMass, &b_CandMass);
   fChain->SetBranchAddress("CandPx", &CandPx, &b_CandPx);
   fChain->SetBranchAddress("CandPy", &CandPy, &b_CandPy);
   fChain->SetBranchAddress("CandPz", &CandPz, &b_CandPz);
   fChain->SetBranchAddress("CandEnergy", &CandEnergy, &b_CandEnergy);
   fChain->SetBranchAddress("CandPuppiWeight", &CandPuppiWeight, &b_CandPuppiWeight);
   fChain->SetBranchAddress("CandPvAssociationQuality", &CandPvAssociationQuality, &b_CandPvAssociationQuality);
   fChain->SetBranchAddress("CandDXY", &CandDXY, &b_CandDXY);
   fChain->SetBranchAddress("CandDZ", &CandDZ, &b_CandDZ);
   fChain->SetBranchAddress("CandDZAssociatedPV", &CandDZAssociatedPV, &b_CandDZAssociatedPV);
   fChain->SetBranchAddress("CandSubJetPart", &CandSubJetPart, &b_CandSubJetPart);
   fChain->SetBranchAddress("ncandidates", &ncandidates, &b_ncandidates);
   fChain->SetBranchAddress("nGencandidates", &nGencandidates, &b_nGencandidates);
   fChain->SetBranchAddress("triggerBit", &triggerBit, &b_triggerBit);
   fChain->SetBranchAddress("triggerPre", &triggerPre, &b_triggerPre);
   fChain->SetBranchAddress("npu", &npu, &b_npu_);
   fChain->SetBranchAddress("genEvtWeight", &genEvtWeight, &b_genEvtWeight_);
   fChain->SetBranchAddress("lheOriginalXWGTUP", &lheOriginalXWGTUP, &b_lheOriginalXWGTUP_);
   fChain->SetBranchAddress("scaleWeights", &scaleWeights, &b_scaleWeights);
   fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
   fChain->SetBranchAddress("GenJetPt", &GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetEta", &GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetPhi", &GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetEnergy", &GenJetEnergy, &b_GenJetEnergy);
   fChain->SetBranchAddress("GenJetMass", &GenJetMass, &b_GenJetMass);
   fChain->SetBranchAddress("isBJetGen", &isBJetGen, &b_isBJetGen);
   fChain->SetBranchAddress("GenJetTau1", &GenJetTau1, &b_GenJetTau1);
   fChain->SetBranchAddress("GenJetTau2", &GenJetTau2, &b_GenJetTau2);
   fChain->SetBranchAddress("GenJetTau3", &GenJetTau3, &b_GenJetTau3);
   fChain->SetBranchAddress("GenJetTau4", &GenJetTau4, &b_GenJetTau4);
   fChain->SetBranchAddress("GenSoftDropMass", &GenSoftDropMass, &b_GenSoftDropMass);
   fChain->SetBranchAddress("GenCandPhi", &GenCandPhi, &b_GenCandPhi);
   fChain->SetBranchAddress("GenCandEta", &GenCandEta, &b_GenCandEta);
   fChain->SetBranchAddress("GenCandPt", &GenCandPt, &b_GenCandPt);
   fChain->SetBranchAddress("GenCandPx", &GenCandPx, &b_GenCandPx);
   fChain->SetBranchAddress("GenCandPy", &GenCandPy, &b_GenCandPy);
   fChain->SetBranchAddress("GenCandPz", &GenCandPz, &b_GenCandPz);
   fChain->SetBranchAddress("GenCandEnergy", &GenCandEnergy, &b_GenCandEnergy);
   fChain->SetBranchAddress("GenmassSub0", &GenmassSub0, &b_GenmassSub0);
   fChain->SetBranchAddress("GenmassSub1", &GenmassSub1, &b_GenmassSub1);
   fChain->SetBranchAddress("GenptSub0", &GenptSub0, &b_GenptSub0);
   fChain->SetBranchAddress("GenptSub1", &GenptSub1, &b_GenptSub1);
   fChain->SetBranchAddress("GenetaSub0", &GenetaSub0, &b_GenetaSub0);
   fChain->SetBranchAddress("GenetaSub1", &GenetaSub1, &b_GenetaSub1);
   fChain->SetBranchAddress("GenphiSub0", &GenphiSub0, &b_GenphiSub0);
   fChain->SetBranchAddress("GenphiSub1", &GenphiSub1, &b_GenphiSub1);
   Notify();
}

Bool_t events::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void events::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t events::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef events_cxx
