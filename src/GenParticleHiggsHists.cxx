#include "UHH2/BoostedHiggsToWW/include/GenParticleHiggsHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

GenParticleHiggsHists::GenParticleHiggsHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  //GenParticles
  number_muon = book<TH1F>("number_muon",     "number of muons",    7,-0.5,6.5);
  number_ele  = book<TH1F>("number_ele",      "number of electrons",7,-0.5,6.5);

  numbers     = book<TH1F>("numbers",     "numbers of leptons from Z",    8,-0.5,7.5);
  std::vector<string> labels = {"Z", "Z to mu", "Z to ele", "Z to else", "muons", "ele", "muon from Z", "ele from Z"};
  for (int i = 1; i <= 8; i++) {
    numbers->GetXaxis()->SetBinLabel(i,(labels[i-1]).c_str());
  }

  flavor_Z     = book<TH1F>("flavor_Z",     "Z decays particles",    50,-0.5, 49.5);
  flavor_H     = book<TH1F>("flavor_H",     "H decays particles",    50,-0.5, 49.5);
  flavor_W1     = book<TH1F>("flavor_W1",     "W1 decays particles",    50,-0.5, 49.5);
  flavor_W2     = book<TH1F>("flavor_W2",     "W2 decays particles",    50,-0.5, 49.5);
  flavor_muon     = book<TH1F>("flavor_muon",     "Muon mothers",    50,-0.5, 49.5);
  flavor_ele     = book<TH1F>("flavor_ele",     "Ele mothers",    50,-0.5, 49.5);

  // number_muon_Z = book<TH1F>("number_muon_Z",     "number of muons from Z",    7,-0.5,6.5);
  // number_ele_Z  = book<TH1F>("number_ele_Z",      "number of electrons from Z", 7,-0.5,6.5);
  //
  // number_Z_muon = book<TH1F>("number_Z_muon",     "number of Z to muons",    7,-0.5,6.5);
  // number_Z_ele  = book<TH1F>("number_Z_ele",      "number of Z to electrons", 7,-0.5,6.5);

  pt_muon     = book<TH1F>("pt_muon",         "p_{T} muon",         100,0,500);
  pt_muon_Z   = book<TH1F>("pt_muon_Z",       "p_{T} muon",         100,0,500);
  eta_muon    = book<TH1F>("eta_muon",        "#eta muon",          100,-3,3);
  phi_muon    = book<TH1F>("phi_muon",        "#phi muon",          100,-M_PI,M_PI);
  charge_muon = book<TH1F>("charge_muon",     "muon charge",        3,-1.5,1.5);

  pt_ele      = book<TH1F>("pt_ele",          "p_{T} electrons",    100,0,500);
  pt_ele_Z    = book<TH1F>("pt_ele_Z",          "p_{T} electrons",    100,0,500);
  eta_ele     = book<TH1F>("eta_ele",         "#eta electrons",     100,-3,3);
  phi_ele     = book<TH1F>("phi_ele",         "#phi electrons",     100,-M_PI,M_PI);
  charge_ele  = book<TH1F>("charge_ele",      "electrons charge",   3,-1.5,1.5);

  mass_higgs  = book<TH1F>("mass_Higgs",      "mass^{Higgs} [GeV]", 200, 120, 130);
  pt_higgs    = book<TH1F>("pt_Higgs",        "pT^{Higgs} [GeV]",   500, 0, 1000);
  eta_higgs   = book<TH1F>("eta_Higgs",       "#eta^{Higgs}",       50, -5, 5);
  phi_higgs   = book<TH1F>("phi_Higgs",       "#phi^{Higgs}",       100, -M_PI, M_PI);

  mass_w1     = book<TH1F>("mass_W_1",        "mass^{W^+} [GeV]",   120, 0, 120);
  pt_w1       = book<TH1F>("pt_W_1",          "pT^{W^+} [GeV]",     150, 0, 600);
  eta_w1      = book<TH1F>("eta_W_1",         "#eta^{W^+}",         50, -5, 5);
  phi_w1      = book<TH1F>("phi_W_1",         "#phi^{W^+}",         100, -M_PI, M_PI);

  mass_w2     = book<TH1F>("mass_W_2",        "mass^{W^-} [GeV]",   120, 0, 120);
  pt_w2       = book<TH1F>("pt_W_2",          "pT^{W^-} [GeV]",     150, 0, 600);
  eta_w2      = book<TH1F>("eta_W_2",         "#eta^{W^-}",         50, -5, 5);
  phi_w2      = book<TH1F>("phi_W_2",         "#phi^{W^-}",         100, -M_PI, M_PI);

  mass_ww     = book<TH1F>("mass_WW",         "mass^{WW} [GeV]",    200, 120, 130);
  pt_ww       = book<TH1F>("pt_WW",           "pT^{WW} [GeV]",      500, 0, 1000);
  eta_ww      = book<TH1F>("eta_WW",          "#eta^{WW}",          50, -5, 5);
  phi_ww      = book<TH1F>("phi_WW",          "#phi^{WW}",          100, -M_PI, M_PI);

  DeltaRW1H   = book<TH1F>("DeltaRW1Higgs",   "#Delta R(W1,H)",     200, 0, 4);
  DeltaRW2H   = book<TH1F>("DeltaRW2Higgs",   "#Delta R(W2,H)",     200, 0, 4);
  DeltaRWW    = book<TH1F>("DeltaR_WW",       "#Delta R(W,W)",      200, 0, 4);
  DeltaRW1jet = book<TH1F>("DeltaRW1jet",     "#Delta R(W1,jet)",   200, 0, 6);
  DeltaRW2jet = book<TH1F>("DeltaRW2jet",     "#Delta R(W1, jet)",  200, 0, 6);
  DeltaRHjet  = book<TH1F>("DeltaRHiggsJet",  "#Delta R(H, jet)",   200, 0, 6);

  DeltaRW1H_DeltaRW2H     = book<TH2F>("DeltaRW1H_DeltaRW2H",     "x=#Delta R(W1,H) y=#Delta R(W2,H)",      200, 0, 4, 200, 0, 4);
  DeltaRW1jet_DeltaRW2jet = book<TH2F>("DeltaRW1jet_DeltaRW2jet", "x=#Delta R(W1,jet) y=#Delta R(W2,jet)",  200, 0, 6, 200, 0, 6);
  DeltaRHjet_ptjet        = book<TH2F>("DeltaRHjet_ptjet",        "x=#Delta R(H,jet) y=pT^{jet}",           200, 0, 6, 500, 0, 1000);
  pTH_pTJet               = book<TH2F>("pTH_pTJet",               "x=pT^{H} y=pT^{jet}",                    500, 0, 1000, 500, 0, 1000);
  pt_muon_pt_muon         = book<TH2F>("pt_muon_pt_muon",         "x=pT^{#mu} y=pT^{#mu}",                  100, 0, 500, 100, 0, 500);
  pt_ele_pt_ele           = book<TH2F>("pt_ele_pt_ele",           "x=pT^{ele} y=pT^{ele}",                  100, 0, 500, 100, 0, 500);

}


void GenParticleHiggsHists::fill(const Event & event){

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  std::vector<TopJet>* jets = event.topjets;

  std::vector<GenParticle>* genParts = event.genparticles;

  GenParticle H, W1, W2;

  int n_muon = 0, n_ele = 0;

  for( auto & gp : *genParts){
    if(abs(gp.pdgId())==23){
      numbers->Fill("Z", weight);
      flavor_Z->Fill( abs(event.genparticles->at(gp.daughter1()).pdgId()), weight);
      flavor_Z->Fill( abs(event.genparticles->at(gp.daughter2()).pdgId()), weight);
      if (abs(event.genparticles->at(gp.daughter1()).pdgId())== 13 && abs(event.genparticles->at(gp.daughter2()).pdgId())== 13) {
        numbers->Fill("Z to mu", 2*weight);
        pt_muon_Z->Fill(event.genparticles->at(gp.daughter1()).pt(), weight);
        pt_muon_Z->Fill(event.genparticles->at(gp.daughter2()).pt(), weight);
        pt_muon_pt_muon->Fill(event.genparticles->at(gp.daughter1()).pt(), event.genparticles->at(gp.daughter2()).pt(), weight);
      }
      else if (abs(event.genparticles->at(gp.daughter1()).pdgId())== 11 && abs(event.genparticles->at(gp.daughter2()).pdgId())== 11){
        numbers->Fill("Z to ele", 2*weight);
        pt_ele_Z->Fill(event.genparticles->at(gp.daughter1()).pt(), weight);
        pt_ele_Z->Fill(event.genparticles->at(gp.daughter2()).pt(), weight);
        pt_ele_pt_ele->Fill(event.genparticles->at(gp.daughter1()).pt(), event.genparticles->at(gp.daughter2()).pt(), weight);
      }
      else numbers->Fill("Z to else", 2*weight);
    }

    if(abs(gp.pdgId())==13){
      pt_muon->Fill(gp.pt(),weight);
      eta_muon->Fill(gp.eta(),weight);
      phi_muon->Fill(gp.phi(),weight);
      charge_muon->Fill(gp.charge(), weight);
      n_muon++;
      numbers->Fill("muons", weight);
      if (gp.mother1() != 65535) flavor_muon->Fill( abs(event.genparticles->at(gp.mother1()).pdgId()), weight);
      if (gp.mother2() != 65535) flavor_muon->Fill( abs(event.genparticles->at(gp.mother2()).pdgId()), weight);
      // std::cout << "GEN: " << gp.pdgId() << " " << gp.status() << " " << gp.index() << " " << gp.mother1() << " " << gp.mother2() << " " << gp.daughter1() << " " << gp.daughter2() << '\n';
      if (abs(event.genparticles->at(gp.mother1()).pdgId())== 23 ) numbers->Fill("muon from Z", weight);
    }

    if(abs(gp.pdgId())==11){
      pt_ele->Fill(gp.pt(),weight);
      eta_ele->Fill(gp.eta(),weight);
      phi_ele->Fill(gp.phi(),weight);
      charge_ele->Fill(gp.charge(), weight);
      numbers->Fill("ele", weight);
      if (gp.mother1() != 65535) flavor_ele->Fill( abs(event.genparticles->at(gp.mother1()).pdgId()), weight);
      if (gp.mother2() != 65535) flavor_ele->Fill( abs(event.genparticles->at(gp.mother2()).pdgId()), weight);
      // std::cout << "GEN: " << gp.pdgId() << " " << gp.status() << " " << gp.index() << " " << gp.mother1() << " " << gp.mother2() << " " << gp.daughter1() << " " << gp.daughter2() << '\n';
      if (abs(event.genparticles->at(gp.mother1()).pdgId())== 23 ) numbers->Fill("ele from Z", weight);
      n_ele++;
    }

    if(gp.pdgId()==25){
      mass_higgs->Fill(gp.v4().M(),weight);
      pt_higgs->Fill(gp.pt(),weight);
      eta_higgs->Fill(gp.eta(),weight);
      phi_higgs->Fill(gp.phi(),weight);
      H = gp;
      flavor_H->Fill( abs(event.genparticles->at(gp.daughter1()).pdgId()), weight);
      flavor_H->Fill( abs(event.genparticles->at(gp.daughter2()).pdgId()), weight);
    }

    if(gp.pdgId()==24){
      mass_w1->Fill(gp.v4().M(),weight);
      pt_w1->Fill(gp.pt(),weight);
      eta_w1->Fill(gp.eta(),weight);
      phi_w1->Fill(gp.phi(),weight);
      flavor_W1->Fill( abs(event.genparticles->at(gp.daughter1()).pdgId()), weight);
      flavor_W1->Fill( abs(event.genparticles->at(gp.daughter2()).pdgId()), weight);
      W1 = gp;
    }

    if(gp.pdgId()== -24){
      mass_w2->Fill(gp.v4().M(),weight);
      pt_w2->Fill(gp.pt(),weight);
      eta_w2->Fill(gp.eta(),weight);
      phi_w2->Fill(gp.phi(),weight);
      flavor_W2->Fill( abs(event.genparticles->at(gp.daughter1()).pdgId()), weight);
      flavor_W2->Fill( abs(event.genparticles->at(gp.daughter2()).pdgId()), weight);
      W2 = gp;
    }
  }

  number_muon->Fill(n_muon, weight);
  number_ele->Fill(n_ele, weight);

  if(jets->size()>0 && H.pdgId()==25 && W1.pdgId()==24 && W2.pdgId()==-24){
    auto dibos = W1.v4() + W2.v4();
    mass_ww->Fill(dibos.M(),weight);
    pt_ww->Fill(dibos.pt(),weight);
    eta_ww->Fill(dibos.eta(),weight);
    phi_ww->Fill(dibos.phi(),weight);
    DeltaRWW->Fill(uhh2::deltaR(W1,W2),weight);
    DeltaRW1H->Fill(uhh2::deltaR(W1,H),weight);
    DeltaRW2H->Fill(uhh2::deltaR(W2,H),weight);
    DeltaRW1jet->Fill(uhh2::deltaR(jets->at(0), W1),weight);
    DeltaRW2jet->Fill(uhh2::deltaR(jets->at(0), W2),weight);
    DeltaRHjet->Fill(uhh2::deltaR(jets->at(0), H),weight);
    DeltaRW1H_DeltaRW2H->Fill(uhh2::deltaR(W1,H), uhh2::deltaR(W2,H), weight);
    DeltaRW1jet_DeltaRW2jet->Fill(uhh2::deltaR(jets->at(0), W1), uhh2::deltaR(jets->at(0), W2), weight);
    DeltaRHjet_ptjet->Fill(uhh2::deltaR(jets->at(0), H), jets->at(0).pt(), weight);
    pTH_pTJet->Fill(H.pt(), jets->at(0).pt(), weight);

  } else {
    std::cout << "\n @@@ WARNING -- GenParticleHiggsHists::fill -- unexpected number of GenParticles in the event --- check selection. \n";
    // std::cout << H.pdgId() << " " << W1.pdgId() << " " << W2.pdgId() << " " << jets->size() << '\n';
  }

}

GenParticleHiggsHists::~GenParticleHiggsHists(){}
