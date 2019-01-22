#define AnalysisTree_cxx
#include "AnalysisTree_origin.h"
// #include "AnalysisTree_presel.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalysisTree::Loop() {

  TH1F *h_pt = new TH1F("jetPt", "; jetPt; A.U.", 100, 300., 500.);
  TH1F *h_eta = new TH1F("jetEta", "; jetEta; A.U.", 100, -5, 5.);
  TH1F *h_phi = new TH1F("jetPhi", "; jetPhi; A.U.", 100, -5, 5.);
  TH1F *h_mass = new TH1F("jetMass", "; jetMass; A.U.", 100, 0., 300.);
  TH1F *h_energy = new TH1F("jetEnergy", "; jetEnergy; A.U.", 100, 300., 1300.);
  TH1F *h_tau1 = new TH1F("jetTau1", "; jetTau1; A.U.", 100, 0., 1.);
  TH1F *h_tau2 = new TH1F("jetTau2", "; jetTau2; A.U.", 100, 0., 1.);
  TH1F *h_tau3 = new TH1F("jetTau3", "; jetTau3; A.U.", 100, 0., 1.);
  TH1F *h_cand = new TH1F("ncandidates", "; ncandidates; A.U.", 100, 0., 200.);
  TH1F *h_btag = new TH1F("jetBtag", "; jetBtag; A.U.", 100, 0., 1.);
  TH1F *h_tau21 = new TH1F("jetTau21", "; jetTau21; A.U.", 100, 0., 1.);
  TH1F *h_tau31 = new TH1F("jetTau31", "; jetTau31; A.U.", 100, 0., 1.);
  TH1F *h_tau32 = new TH1F("jetTau32", "; jetTau32; A.U.", 100, 0., 1.);


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << outname << " " << nentries << '\n';
// return true;
  // nentries = 900000;
  Long64_t nbytes = 0, nb = 0;
  int count = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%100000==0) std::cout << jentry << '\n';
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // std::cout << packedPatJetsAk8CHSJets_SoftDropCHS_ << '\n';
    // std::cout << packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[0] << " " << packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[1] << " " << packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[2] << '\n';
    // if (Cut(ientry) < 0) continue;

    for (size_t index = 0; index < packedPatJetsAk8CHSJets_SoftDropCHS_; index++) {
      if (packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[index]<300 || packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[index]>500) {
        continue;
      }
      count++;
      if (count%100==0)  std::cout << count << " " << index << " " << packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[index] << '\n';

      // std::cout << packedPatJetsAk8CHSJets_SoftDropCHS_m_numberOfDaughters[index]-(packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedMultiplicity[index]+packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralMultiplicity[index]) << '\n';
      TLorentzVector jet;
      jet.SetPtEtaPhiE(packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[index],packedPatJetsAk8CHSJets_SoftDropCHS_m_eta[index],packedPatJetsAk8CHSJets_SoftDropCHS_m_phi[index],packedPatJetsAk8CHSJets_SoftDropCHS_m_energy[index]);
      h_pt->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[index]);
      h_eta->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_eta[index]);
      h_phi->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_phi[index]);
      h_mass->Fill(jet.M());
      h_energy->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_energy[index]);
      h_tau1->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1[index]);
      h_tau2->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2[index]);
      h_tau3->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3[index]);
      h_cand->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_numberOfDaughters[index]);
      h_btag->Fill(packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertex[index]);
      h_tau21->Fill( (packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1[index]!=0) ? (packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2[index]/packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1[index]) :0);
      h_tau31->Fill( (packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1[index]!=0) ? (packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3[index]/packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1[index]) :0);
      h_tau32->Fill( (packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2[index]!=0) ? (packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3[index]/packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2[index]) :0);
    }
  }

  std::cout << "Writing" << '\n';
  TFile *MyFile = new TFile(outname,"RECREATE");
  h_pt->Write();
  h_eta->Write();
  h_phi->Write();
  h_mass->Write();
  h_energy->Write();
  h_tau1->Write();
  h_tau2->Write();
  h_tau3->Write();
  h_cand->Write();
  h_btag->Write();
  h_tau21->Write();
  h_tau31->Write();
  h_tau32->Write();
  MyFile->Close();
}


void MakeHistoNN(TString name, TString outname) {

  TFile *file = new TFile(name);
  TTree *tree   = (TTree*)file->Get("boostedAK8/events");

  TH1F *h_pt = new TH1F("jetPt", "; jetPt; A.U.", 100, 300., 500.);
  TH1F *h_eta = new TH1F("jetEta", "; jetEta; A.U.", 100, -5, 5.);
  TH1F *h_phi = new TH1F("jetPhi", "; jetPhi; A.U.", 100, -5, 5.);
  TH1F *h_mass = new TH1F("jetMass", "; jetMass; A.U.", 100, 0., 300.);
  TH1F *h_energy = new TH1F("jetEnergy", "; jetEnergy; A.U.", 100, 300., 1300.);
  TH1F *h_tau1 = new TH1F("jetTau1", "; jetTau1; A.U.", 100, 0., 1.);
  TH1F *h_tau2 = new TH1F("jetTau2", "; jetTau2; A.U.", 100, 0., 1.);
  TH1F *h_tau3 = new TH1F("jetTau3", "; jetTau3; A.U.", 100, 0., 1.);
  TH1F *h_cand = new TH1F("ncandidates", "; ncandidates; A.U.", 100, 0., 200.);
  TH1F *h_btag = new TH1F("jetBtag", "; jetBtag; A.U.", 100, 0., 1.);
  TH1F *h_tau21 = new TH1F("jetTau21", "; jetTau21; A.U.", 100, 0., 1.);
  TH1F *h_tau31 = new TH1F("jetTau31", "; jetTau31; A.U.", 100, 0., 1.);
  TH1F *h_tau32 = new TH1F("jetTau32", "; jetTau32; A.U.", 100, 0., 1.);

  vector<float>* jetPt = 0; tree->SetBranchAddress("jetPt", &jetPt); tree->SetBranchStatus("jetPt", 1);
  vector<float>* jetEta = 0; tree->SetBranchAddress("jetEta", &jetEta); tree->SetBranchStatus("jetEta", 1);
  vector<float>* jetPhi = 0; tree->SetBranchAddress("jetPhi", &jetPhi); tree->SetBranchStatus("jetPhi", 1);
  vector<float>* jetMass = 0; tree->SetBranchAddress("jetMass", &jetMass); tree->SetBranchStatus("jetMass", 1);
  vector<float>* jetEnergy = 0; tree->SetBranchAddress("jetEnergy", &jetEnergy); tree->SetBranchStatus("jetEnergy", 1);
  vector<float>* jetTau1 = 0; tree->SetBranchAddress("jetTau1", &jetTau1); tree->SetBranchStatus("jetTau1", 1);
  vector<float>* jetTau2 = 0; tree->SetBranchAddress("jetTau2", &jetTau2); tree->SetBranchStatus("jetTau2", 1);
  vector<float>* jetTau3 = 0; tree->SetBranchAddress("jetTau3", &jetTau3); tree->SetBranchStatus("jetTau3", 1);
  vector<float>* ncandidates = 0; tree->SetBranchAddress("ncandidates", &ncandidates); tree->SetBranchStatus("ncandidates", 1);
  vector<float>* jetBtag = 0; tree->SetBranchAddress("jetBtag", &jetBtag); tree->SetBranchStatus("jetBtag", 1);

  std::cout << tree->GetEntriesFast() << '\n';
  for (int i=0; i<tree->GetEntriesFast(); ++i) {
    tree->GetEvent(i);
    if (jetPt->size()==0 || jetPt->at(0)<300 || jetPt->at(0)>500) continue;

    h_pt->Fill(jetPt->at(0));
    h_eta->Fill(jetEta->at(0));
    h_phi->Fill(jetPhi->at(0));
    h_mass->Fill(jetMass->at(0));
    h_energy->Fill(jetEnergy->at(0));
    h_tau1->Fill(jetTau1->at(0));
    h_tau2->Fill(jetTau2->at(0));
    h_tau3->Fill(jetTau3->at(0));
    h_cand->Fill(ncandidates->at(0));
    h_btag->Fill(jetBtag->at(0));
    h_tau21->Fill( (jetTau1->at(0)!=0) ? (jetTau2->at(0)/jetTau1->at(0)) :0);
    h_tau31->Fill( (jetTau1->at(0)!=0) ? (jetTau3->at(0)/jetTau1->at(0)) :0);
    h_tau32->Fill( (jetTau2->at(0)!=0) ? (jetTau3->at(0)/jetTau2->at(0)) :0);
  }

  std::cout << "Writing" << '\n';
  TFile *MyFile = new TFile(outname,"RECREATE");
  h_pt->Write();
  h_eta->Write();
  h_phi->Write();
  h_mass->Write();
  h_energy->Write();
  h_tau1->Write();
  h_tau2->Write();
  h_tau3->Write();
  h_cand->Write();
  h_btag->Write();
  h_tau21->Write();
  h_tau31->Write();
  h_tau32->Write();
  MyFile->Close();
}



void MakeHisto() {
  AnalysisTree* t;
  TString path = "compareInputs/";
  // std::vector<TString> samples = {"presel_DY","presel_HZ","presel_TT","presel_ZZ"};
  // std::vector<TString> samples = {"origin_TT", "origin_HZ", "origin_DY", "origin_DY1", "origin_DY2", "origin_DY3", "origin_DY4", "origin_TTBarhad", "origin_TTBarlep", "origin_TTBarsemilep"};
  std::vector<TString> samples = {"origin_DY", "origin_DY1", "origin_DY2", "origin_DY3"};

  for (auto &sample: samples) {
    t = new AnalysisTree(path+sample+".root", path+"h_"+sample+".root");
    t->Loop();
  }

  // std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT"};
  // for (auto &sample: samples) {
  //   std::cout << sample << '\n';
  //   MakeHistoNN(path+sample+".root", path+"h_"+sample+".root");
  // }

}
