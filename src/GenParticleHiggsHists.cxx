#include "UHH2/BoostedHiggsToWW/include/GenParticleHiggsHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

GenParticleHiggsHists::GenParticleHiggsHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  //GenParticles
  mass_higgs  = book<TH1F>("mass_Higgs",      "mass^{Higgs} [GeV]", 200, 120, 130);
  pt_higgs    = book<TH1F>("pt_Higgs",        "pT^{Higgs} [GeV]",   500, 0, 1000);
  eta_higgs   = book<TH1F>("eta_Higgs",       "#eta^{Higgs}",       50, -5, 5);
  phi_higgs   = book<TH1F>("phi_Higgs",       "#phi^{Higgs}",       50, -M_PI, M_PI);

  mass_w1     = book<TH1F>("mass_W_1",        "mass^{W^+} [GeV]",   120, 0, 120);
  pt_w1       = book<TH1F>("pt_W_1",          "pT^{W^+} [GeV]",     150, 0, 600);
  eta_w1      = book<TH1F>("eta_W_1",         "#eta^{W^+}",         50, -5, 5);
  phi_w1      = book<TH1F>("phi_W_1",         "#phi^{W^+}",         50, -M_PI, M_PI);

  mass_w2     = book<TH1F>("mass_W_2",        "mass^{W^-} [GeV]",   120, 0, 120);
  pt_w2       = book<TH1F>("pt_W_2",          "pT^{W^-} [GeV]",     150, 0, 600);
  eta_w2      = book<TH1F>("eta_W_2",         "#eta^{W^-}",         50, -5, 5);
  phi_w2      = book<TH1F>("phi_W_2",         "#phi^{W^-}",         50, -M_PI, M_PI);

  mass_ww     = book<TH1F>("mass_WW",         "mass^{WW} [GeV]",    200, 120, 130);
  pt_ww       = book<TH1F>("pt_WW",           "pT^{WW} [GeV]",      500, 0, 1000);
  eta_ww      = book<TH1F>("eta_WW",          "#eta^{WW}",          50, -5, 5);
  phi_ww      = book<TH1F>("phi_WW",          "#phi^{WW}",          50, -M_PI, M_PI);

  DeltaRW1H   = book<TH1F>("DeltaRW1Higgs",   "#Delta R(W1,H)",     200, 0, 4);
  DeltaRW2H   = book<TH1F>("DeltaRW2Higgs",   "#Delta R(W2,H)",     200, 0, 4);
  DeltaRWW    = book<TH1F>("DeltaR_WW",       "#Delta R(W,W)",      200, 0, 4);
  DeltaRW1jet = book<TH1F>("DeltaRW1jet",     "#Delta R(W1,jet)",   200, 0, 6);
  DeltaRW2jet = book<TH1F>("DeltaRW2jet",     "#Delta R(W1, jet)",  200, 0, 6);
  DeltaRHjet  = book<TH1F>("DeltaRHiggsJet",  "#Delta R(H, jet)",   200, 0, 6);

  DeltaRW1H_DeltaRW2H     = book<TH2F>("DeltaRW1H_DeltaRW2H",     "x=#Delta R(W1,H) y=#Delta R(W2,H)",      200, 0, 4, 200, 0, 4);
  DeltaRW1jet_DeltaRW2jet = book<TH2F>("DeltaRW1jet_DeltaRW2jet", "x=#Delta R(W1,jet) y=#Delta R(W2,jet)",  200, 0, 6, 200, 0, 6);

}


void GenParticleHiggsHists::fill(const Event & event){

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  std::vector<TopJet>* jets = event.topjets;

  std::vector<GenParticle>* genParts = event.genparticles;

  GenParticle H, W1, W2;

  for( auto & gp : *genParts){
    if(gp.pdgId()==25){
      mass_higgs->Fill(gp.v4().M(),weight);
      pt_higgs->Fill(gp.pt(),weight);
      eta_higgs->Fill(gp.eta(),weight);
      phi_higgs->Fill(gp.phi(),weight);
      H = gp;
    }

    if(gp.pdgId()==24){
      mass_w1->Fill(gp.v4().M(),weight);
      pt_w1->Fill(gp.pt(),weight);
      eta_w1->Fill(gp.eta(),weight);
      phi_w1->Fill(gp.phi(),weight);
      W1 = gp;
    }

    if(gp.pdgId()== -24){
      mass_w2->Fill(gp.v4().M(),weight);
      pt_w2->Fill(gp.pt(),weight);
      eta_w2->Fill(gp.eta(),weight);
      phi_w2->Fill(gp.phi(),weight);
      W2 = gp;
    }
  }

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
  } else {
    std::cout << "\n @@@ WARNING -- GenParticleHiggsHists::fill -- unexpected number of GenParticles in the event --- check selection. \n";
    std::cout << H.pdgId() << " " << W1.pdgId() << " " << W2.pdgId() << " " << jets->size() << '\n';
  }

}

GenParticleHiggsHists::~GenParticleHiggsHists(){}
