#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

BoostedHiggsToWWHists::BoostedHiggsToWWHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  number_jet  = book<TH1F>("number",      "number of jets",         20, 0, 20);
  flavor_jet  = book<TH1F>("flavor_jets", "flavor of jets",         200, -100, 100);
  eta_jet     = book<TH1F>("eta_jet",     "#eta of jets",           40, -3.5, 3.5);
  pt_jet      = book<TH1F>("pt_jet",      "p_{T} of jets [GeV/c]",  100, 0, 1000);

  mass_jet1   = book<TH1F>("mass_jet1",   "#mass^{jet 1} [GeV]",    200, 0, 400);
  mass_jet2   = book<TH1F>("mass_jet2",   "#mass^{jet 2} [GeV]",    200, 0, 400);
  mass_jet3   = book<TH1F>("mass_jet3",   "#mass^{jet 3} [GeV]",    200, 0, 400);
  mass_jet4   = book<TH1F>("mass_jet4",   "#mass^{jet 4} [GeV]",    200, 0, 400);

  MT_jet1   = book<TH1F>("MT_jet1",       "M_T^{jet 1} [GeV]",      200, 0, 400);
  MT_jet2   = book<TH1F>("MT_jet2",       "M_T^{jet 2} [GeV]",      200, 0, 400);
  MT_jet3   = book<TH1F>("MT_jet3",       "M_T^{jet 3} [GeV]",      200, 0, 400);
  MT_jet4   = book<TH1F>("MT_jet4",       "M_T^{jet 4} [GeV]",      200, 0, 400);

  pt_jet1     = book<TH1F>("pt_jet1",     "p_{T}^{jet 1} [GeV/c]",  100, 0, 1000);
  pt_jet2     = book<TH1F>("pt_jet2",     "p_{T}^{jet 2} [GeV/c]",  100, 0, 1000);
  pt_jet3     = book<TH1F>("pt_jet3",     "p_{T}^{jet 3} [GeV/c]",  100, 0, 1000);
  pt_jet4     = book<TH1F>("pt_jet4",     "p_{T}^{jet 4} [GeV/c]",  100, 0, 1000);

  eta_jet1    = book<TH1F>("eta_jet1",    "#eta^{jet 1}",           40, -3.5, 3.5);
  eta_jet2    = book<TH1F>("eta_jet2",    "#eta^{jet 2}",           40, -3.5, 3.5);
  eta_jet3    = book<TH1F>("eta_jet3",    "#eta^{jet 3}",           40, -3.5, 3.5);
  eta_jet4    = book<TH1F>("eta_jet4",    "#eta^{jet 4}",           40, -3.5, 3.5);

  phi_jet1    = book<TH1F>("phi_jet1",    "#phi^{jet 1}",           40, -M_PI, M_PI);
  phi_jet2    = book<TH1F>("phi_jet2",    "#phi^{jet 2}",           40, -M_PI, M_PI);
  phi_jet3    = book<TH1F>("phi_jet3",    "#phi^{jet 3}",           40, -M_PI, M_PI);
  phi_jet4    = book<TH1F>("phi_jet4",    "#phi^{jet 4}",           40, -M_PI, M_PI);

  //jet substructure for the analysed jets
  tau1_jet1   = book<TH1F>("tau1_jet1",   "#tau1^{jet 1}",          20, 0.,1.);
  tau1_jet2   = book<TH1F>("tau1_jet2",   "#tau1^{jet 2}",          20, 0.,1.);
  tau1_jet3   = book<TH1F>("tau1_jet3",   "#tau1^{jet 3}",          20, 0.,1.);
  tau1_jet4   = book<TH1F>("tau1_jet4",   "#tau1^{jet 4}",          20, 0.,1.);

  tau2_jet1   = book<TH1F>("tau2_jet1",   "#tau2^{jet 1}",          20, 0.,1.);
  tau2_jet2   = book<TH1F>("tau2_jet2",   "#tau2^{jet 2}",          20, 0.,1.);
  tau2_jet3   = book<TH1F>("tau2_jet3",   "#tau2^{jet 3}",          20, 0.,1.);
  tau2_jet4   = book<TH1F>("tau2_jet4",   "#tau2^{jet 4}",          20, 0.,1.);

  tau3_jet1   = book<TH1F>("tau3_jet1",   "#tau3^{jet 1}",          20, 0.,1.);
  tau3_jet2   = book<TH1F>("tau3_jet2",   "#tau3^{jet 2}",          20, 0.,1.);
  tau3_jet3   = book<TH1F>("tau3_jet3",   "#tau3^{jet 3}",          20, 0.,1.);
  tau3_jet4   = book<TH1F>("tau3_jet4",   "#tau3^{jet 4}",          20, 0.,1.);

  tau21_jet1  = book<TH1F>("tau21_jet1",  "#tau21^{jet 1}",         20, 0.,1.);
  tau21_jet2  = book<TH1F>("tau21_jet2",  "#tau21^{jet 2}",         20, 0.,1.);
  tau21_jet3  = book<TH1F>("tau21_jet3",  "#tau21^{jet 3}",         20, 0.,1.);
  tau21_jet4  = book<TH1F>("tau21_jet4",  "#tau21^{jet 4}",         20, 0.,1.);

  tau31_jet1  = book<TH1F>("tau31_jet1",  "#tau31^{jet 1}",         20, 0.,1.);
  tau31_jet2  = book<TH1F>("tau31_jet2",  "#tau31^{jet 2}",         20, 0.,1.);
  tau31_jet3  = book<TH1F>("tau31_jet3",  "#tau31^{jet 3}",         20, 0.,1.);
  tau31_jet4  = book<TH1F>("tau31_jet4",  "#tau31^{jet 4}",         20, 0.,1.);

  tau32_jet1  = book<TH1F>("tau32_jet1",  "#tau32^{jet 1}",         20, 0.,1.);
  tau32_jet2  = book<TH1F>("tau32_jet2",  "#tau32^{jet 2}",         20, 0.,1.);
  tau32_jet3  = book<TH1F>("tau32_jet3",  "#tau32^{jet 3}",         20, 0.,1.);
  tau32_jet4  = book<TH1F>("tau32_jet4",  "#tau32^{jet 4}",         20, 0.,1.);

  SDmass_jet1 = book<TH1F>("SDmass_jet1", "#SDmass^{jet 1} [GeV]",  200, 0, 400);
  SDmass_jet2 = book<TH1F>("SDmass_jet2", "#SDmass^{jet 2} [GeV]",  200, 0, 400);
  SDmass_jet3 = book<TH1F>("SDmass_jet3", "#SDmass^{jet 3} [GeV]",  200, 0, 400);
  SDmass_jet4 = book<TH1F>("SDmass_jet4", "#SDmass^{jet 4} [GeV]",  200, 0, 400);

  // leptons
  number_mu   = book<TH1F>("N_mu",        "N^{#mu}",                10, 0, 10);
  pt_mu       = book<TH1F>("pt_mu",       "p_{T}^{#mu} [GeV/c]",    40, 0, 200);
  eta_mu      = book<TH1F>("eta_mu",      "#eta^{#mu}",             40, -3.5, 3.5);
  reliso_mu   = book<TH1F>("reliso_mu",   "#mu rel. Iso",           40, 0, 0.5);

  number_ele    = book<TH1F>("N_ele ",        "N^{ele}",                10, 0, 10);
  pt_ele        = book<TH1F>("pt_ele ",       "p_{T}^{ele} [GeV/c]",    40, 0, 200);
  eta_ele       = book<TH1F>("eta_ele ",      "#eta^{ele}",             40, -3.5, 3.5);
  reliso_ele    = book<TH1F>("reliso_ele ",   "ele rel. Iso",           40, 0, 0.5);

  // primary vertices
  number_pv   = book<TH1F>("N_pv",        "N^{PV}",                 50, 0, 50);
}


void BoostedHiggsToWWHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double w = event.weight;
  std::vector<TopJet>*  jets = event.topjets;
  int Njets = (jets)->size();

  number_jet->Fill(Njets, w);
  for(const auto & jet : *jets){
    pt_jet->Fill(jet.pt(), w);
    eta_jet->Fill(jet.eta(), w);
    flavor_jet->Fill(fabs(jet.pdgId()), w);
  }

  if(jets->size() > 0){
    const auto & jet = (*jets)[0];
    mass_jet1->Fill(jet.v4().M(), w);
    MT_jet1->Fill(jet.v4().Mt(), w);
    pt_jet1->Fill(jet.pt(), w);
    eta_jet1->Fill(jet.eta(), w);
    phi_jet1->Fill(jet.phi(), w);
    tau1_jet1->Fill(jet.tau1(), w);
    tau2_jet1->Fill(jet.tau2(), w);
    tau3_jet1->Fill(jet.tau3(), w);
    SDmass_jet1->Fill(jet.softdropmass(), w);
    if (jet.tau1()!=0) {
      tau21_jet1->Fill(jet.tau2()/jet.tau1(), w);
      tau31_jet1->Fill(jet.tau3()/jet.tau1(), w);
    }
    if (jet.tau2()!=0) {
      tau32_jet1->Fill(jet.tau3()/jet.tau2(), w);
    }
  }

  if(jets->size() > 1){
    const auto & jet = (*jets)[1];
    mass_jet2->Fill(jet.v4().M(), w);
    MT_jet1->Fill(jet.v4().Mt(), w);
    pt_jet2->Fill(jet.pt(), w);
    eta_jet2->Fill(jet.eta(), w);
    phi_jet2->Fill(jet.phi(), w);
    tau1_jet2->Fill(jet.tau1(), w);
    tau2_jet2->Fill(jet.tau2(), w);
    tau3_jet2->Fill(jet.tau3(), w);
    SDmass_jet2->Fill(jet.softdropmass(), w);
    if (jet.tau1()!=0) {
      tau21_jet2->Fill(jet.tau2()/jet.tau1(), w);
      tau31_jet2->Fill(jet.tau3()/jet.tau1(), w);
    }
    if (jet.tau2()!=0) {
      tau32_jet2->Fill(jet.tau3()/jet.tau2(), w);
    }
  }

  if(jets->size() > 2){
    const auto & jet = (*jets)[3];
    mass_jet3->Fill(jet.v4().M(), w);
    MT_jet3->Fill(jet.v4().Mt(), w);
    pt_jet3->Fill(jet.pt(), w);
    eta_jet3->Fill(jet.eta(), w);
    phi_jet3->Fill(jet.phi(), w);
    tau1_jet3->Fill(jet.tau1(), w);
    tau2_jet3->Fill(jet.tau2(), w);
    tau3_jet3->Fill(jet.tau3(), w);
    SDmass_jet3->Fill(jet.softdropmass(), w);
    if (jet.tau1()!=0) {
      tau21_jet3->Fill(jet.tau2()/jet.tau1(), w);
      tau31_jet3->Fill(jet.tau3()/jet.tau1(), w);
    }
    if (jet.tau2()!=0) {
      tau32_jet3->Fill(jet.tau3()/jet.tau2(), w);
    }
  }

  if(jets->size() > 3){
    const auto & jet = (*jets)[3];
    mass_jet4->Fill(jet.v4().M(), w);
    MT_jet4->Fill(jet.v4().Mt(), w);
    pt_jet4->Fill(jet.pt(), w);
    eta_jet4->Fill(jet.eta(), w);
    phi_jet4->Fill(jet.phi(), w);
    tau1_jet4->Fill(jet.tau1(), w);
    tau2_jet4->Fill(jet.tau2(), w);
    tau3_jet4->Fill(jet.tau3(), w);
    SDmass_jet4->Fill(jet.softdropmass(), w);
    if (jet.tau1()!=0) {
      tau21_jet4->Fill(jet.tau2()/jet.tau1(), w);
      tau31_jet4->Fill(jet.tau3()/jet.tau1(), w);
    }
    if (jet.tau2()!=0) {
      tau32_jet4->Fill(jet.tau3()/jet.tau2(), w);
    }
  }

  int Nmuons = event.muons->size();
  number_mu->Fill(Nmuons, w);
  for (const Muon & mu : *event.muons){
    pt_mu->Fill(mu.pt(), w);
    eta_mu->Fill(mu.eta(), w);
    reliso_mu->Fill(mu.relIso(), w);
  }

  int Nele = event.electrons->size();
  number_ele->Fill(Nele, w);
  for (const Electron & ele : *event.electrons){
    pt_ele->Fill(ele.pt(), w);
    eta_ele->Fill(ele.eta(), w);
    reliso_ele->Fill(ele.relIso(), w);
  }

  int Npvs = event.pvs->size();
  number_pv->Fill(Npvs, w);

}

BoostedHiggsToWWHists::~BoostedHiggsToWWHists(){}
