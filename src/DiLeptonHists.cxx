#include "UHH2/BoostedHiggsToWW/include/DiLeptonHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

DiLeptonHists::DiLeptonHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  dilep_number  = book<TH1F>("number",        "number of dilepton",   20, 0, 20);
  // DI-ELECTRON
  diele_number  = book<TH1F>("diele_number",  "number of electrons",  20, 0, 20);
  diele_charge  = book<TH1F>("diele_charge",  "charge^{ee}",          5,-2,2);
  diele_m       = book<TH1F>("diele_m",       "M^{ee}",               100,0,100);
  diele_pt      = book<TH1F>("diele_pt",      "p_{T}^{ee}",           100,0,500);
  diele_eta     = book<TH1F>("diele_eta",     "#eta^{ee}",            100,-3,3);
  diele_phi     = book<TH1F>("diele_phi",     "#phi^{ee}",            100,-M_PI,M_PI);
  diele_deltaR  = book<TH1F>("diele__DR12",   "#Delta R(e,e)",        40, 0, 2.0);

  // DI-MUON
  dimuon_number  = book<TH1F>("dimuon_number",  "number of muons",    20, 0, 20);
  dimuon_charge  = book<TH1F>("dimuon_charge",  "charge^{ee}",        5,-2,2);
  dimuon_m       = book<TH1F>("dimuon_m",       "M^{ee}",             100,0,100);
  dimuon_pt      = book<TH1F>("dimuon_pt",      "p_{T}^{ee}",         100,0,500);
  dimuon_eta     = book<TH1F>("dimuon_eta",     "#eta^{ee}",          100,-3,3);
  dimuon_phi     = book<TH1F>("dimuon_phi",     "#phi^{ee}",          100,-M_PI,M_PI);
  dimuon_deltaR  = book<TH1F>("dimuon__DR12",   "#Delta R(e,e)",      40, 0, 2.0);
  pt_muon_2      = book<TH2F>("pt_muon_2", ";PT1;pT2", 100, 0, 500, 100, 0, 500.);
  pt_ele_2       = book<TH2F>("pt_ele_2",  ";PT1;pT2", 100, 0, 500, 100, 0, 500.);
}


void DiLeptonHists::fill(const Event & event){
  // fill the histograms.

  // Don't forget to always use the weight when filling.
  auto weight = event.weight;
  assert(event.electrons);
  assert(event.muons);

  if (event.muons->size()>1) {
    pt_muon_2->Fill( event.muons->at(0).pt(), event.muons->at(1).pt(), weight);
  }

  if (event.electrons->size()>1) {
    pt_ele_2->Fill( event.electrons->at(0).pt(), event.electrons->at(1).pt(), weight);
  }

  if (event.electrons->size()==2) {
    dilep_number->Fill(event.muons->size(), weight);
    diele_number->Fill(event.electrons->size(), weight);
    const auto& lep1 = event.electrons->at(0);
    const auto& lep2 = event.electrons->at(1);

    auto dilep = lep1.v4() + lep2.v4();

    diele_charge->Fill(lep1.charge()+lep2.charge(), weight);
    diele_deltaR->Fill(uhh2::deltaR(lep1, lep2)   , weight);

    diele_m->Fill(dilep.M()  , weight);
    diele_pt->Fill(dilep.Pt() , weight);
    diele_eta->Fill(dilep.Eta(), weight);
    diele_phi->Fill(dilep.Phi(), weight);

  } else if (event.muons->size()==2) {
    dilep_number->Fill(event.muons->size(), weight);
    dimuon_number->Fill(event.muons->size(), weight);
    const auto& lep1 = event.muons->at(0);
    const auto& lep2 = event.muons->at(1);

    auto dilep = lep1.v4() + lep2.v4();

    dimuon_charge->Fill(lep1.charge()+lep2.charge(), weight);
    dimuon_deltaR->Fill(uhh2::deltaR(lep1, lep2)   , weight);

    dimuon_m->Fill(dilep.M()  , weight);
    dimuon_pt->Fill(dilep.Pt() , weight);
    dimuon_eta->Fill(dilep.Eta(), weight);
    dimuon_phi->Fill(dilep.Phi(), weight);
  }

}

DiLeptonHists::~DiLeptonHists(){}
