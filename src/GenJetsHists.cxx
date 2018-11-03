#include "UHH2/BoostedHiggsToWW/include/GenJetsHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

GenJetsHists::GenJetsHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  number_jet  = book<TH1F>("number",      "number of jets",         20, 0, 20);
  mass_jet    = book<TH1F>("mass_jet",    "#mass^{jet} [GeV]",      200, 0, 400);
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

  phi_jet1    = book<TH1F>("phi_jet1",    "#phi^{jet 1}",           100, -M_PI, M_PI);
  phi_jet2    = book<TH1F>("phi_jet2",    "#phi^{jet 2}",           100, -M_PI, M_PI);
  phi_jet3    = book<TH1F>("phi_jet3",    "#phi^{jet 3}",           100, -M_PI, M_PI);
  phi_jet4    = book<TH1F>("phi_jet4",    "#phi^{jet 4}",           100, -M_PI, M_PI);

  //jet substructure for the analysed jets
  tau1_jet1   = book<TH1F>("tau1_jet1",   "#tau1^{jet 1}",          20, -1.,1.);
  tau1_jet2   = book<TH1F>("tau1_jet2",   "#tau1^{jet 2}",          20, -1.,1.);
  tau1_jet3   = book<TH1F>("tau1_jet3",   "#tau1^{jet 3}",          20, -1.,1.);
  tau1_jet4   = book<TH1F>("tau1_jet4",   "#tau1^{jet 4}",          20, -1.,1.);

  tau2_jet1   = book<TH1F>("tau2_jet1",   "#tau2^{jet 1}",          20, -1.,1.);
  tau2_jet2   = book<TH1F>("tau2_jet2",   "#tau2^{jet 2}",          20, -1.,1.);
  tau2_jet3   = book<TH1F>("tau2_jet3",   "#tau2^{jet 3}",          20, -1.,1.);
  tau2_jet4   = book<TH1F>("tau2_jet4",   "#tau2^{jet 4}",          20, -1.,1.);

  tau3_jet1   = book<TH1F>("tau3_jet1",   "#tau3^{jet 1}",          20, -1.,1.);
  tau3_jet2   = book<TH1F>("tau3_jet2",   "#tau3^{jet 2}",          20, -1.,1.);
  tau3_jet3   = book<TH1F>("tau3_jet3",   "#tau3^{jet 3}",          20, -1.,1.);
  tau3_jet4   = book<TH1F>("tau3_jet4",   "#tau3^{jet 4}",          20, -1.,1.);

  tau21_jet1  = book<TH1F>("tau21_jet1",  "#tau21^{jet 1}",         20, -1.,1.);
  tau21_jet2  = book<TH1F>("tau21_jet2",  "#tau21^{jet 2}",         20, -1.,1.);
  tau21_jet3  = book<TH1F>("tau21_jet3",  "#tau21^{jet 3}",         20, -1.,1.);
  tau21_jet4  = book<TH1F>("tau21_jet4",  "#tau21^{jet 4}",         20, -1.,1.);

  tau31_jet1  = book<TH1F>("tau31_jet1",  "#tau31^{jet 1}",         20, -1.,1.);
  tau31_jet2  = book<TH1F>("tau31_jet2",  "#tau31^{jet 2}",         20, -1.,1.);
  tau31_jet3  = book<TH1F>("tau31_jet3",  "#tau31^{jet 3}",         20, -1.,1.);
  tau31_jet4  = book<TH1F>("tau31_jet4",  "#tau31^{jet 4}",         20, -1.,1.);

  tau32_jet1  = book<TH1F>("tau32_jet1",  "#tau32^{jet 1}",         20, -1.,1.);
  tau32_jet2  = book<TH1F>("tau32_jet2",  "#tau32^{jet 2}",         20, -1.,1.);
  tau32_jet3  = book<TH1F>("tau32_jet3",  "#tau32^{jet 3}",         20, -1.,1.);
  tau32_jet4  = book<TH1F>("tau32_jet4",  "#tau32^{jet 4}",         20, -1.,1.);


  DeltaRW1jet = book<TH1F>("DeltaRW1jet",     "#Delta R(W1,jet)",   200, 0, 6);
  DeltaRW2jet = book<TH1F>("DeltaRW2jet",     "#Delta R(W1, jet)",  200, 0, 6);
  DeltaRHjet  = book<TH1F>("DeltaRHiggsJet",  "#Delta R(H, jet)",   200, 0, 6);

  DeltaRW1jet_DeltaRW2jet = book<TH2F>("DeltaRW1jet_DeltaRW2jet", "x=#Delta R(W1,jet) y=#Delta R(W2,jet)",  200, 0, 6, 200, 0, 6);
  DeltaRHjet_ptjet        = book<TH2F>("DeltaRHjet_ptjet",        "x=#Delta R(H,jet) y=pT^{jet}",           200, 0, 6, 500, 0, 1000);
  DeltaRHjet_Mjet         = book<TH2F>("DeltaRHjet_Mjet",         "x=#Delta R(H,jet) y=M^{jet}",            200, 0, 6, 500, 0, 150);
  pTH_pTJet               = book<TH2F>("pTH_pTJet",               "x=pT^{H} y=pT^{jet}",                    500, 0, 1000, 500, 0, 1000);
  pTH_pTZ                 = book<TH2F>("pTH_pTZ",                 "x=pT^{H} y=pT^{Z}",                      500, 0, 1000, 500, 0, 1000);
  pTW1_pTW2               = book<TH2F>("pTW1_pTW2",               "x=pT^{W1} y=pT^{W2}",                    500, 0, 1000, 500, 0, 1000);



  DeltaRW1jet_1 = book<TH1F>("DeltaRW1jet_1",     "#Delta R(W1,jet_1)",   200, 0, 6);
  DeltaRW2jet_1 = book<TH1F>("DeltaRW2jet_1",     "#Delta R(W1, jet_1)",  200, 0, 6);
  DeltaRHjet_1  = book<TH1F>("DeltaRHiggsJet_1",  "#Delta R(H, jet_1)",   200, 0, 6);

  DeltaRW1jet_DeltaRW2jet_1 = book<TH2F>("DeltaRW1jet_DeltaRW2jet_1", "x=#Delta R(W1,jet_1) y=#Delta R(W2,jet_1)",  200, 0, 6, 200, 0, 6);
  DeltaRHjet_ptjet_1        = book<TH2F>("DeltaRHjet_ptjet_1",        "x=#Delta R(H,jet_1) y=pT^{jet_1}",           200, 0, 6, 500, 0, 1000);
  DeltaRHjet_Mjet_1         = book<TH2F>("DeltaRHjet_Mjet_1",         "x=#Delta R(H,jet_1) y=M^{jet_1}",            200, 0, 6, 500, 0, 150);
  pTH_pTJet_1               = book<TH2F>("pTH_pTJet_1",               "x=pT^{H} y=pT^{jet_1}",                    500, 0, 1000, 500, 0, 1000);

}


void GenJetsHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double w = event.weight;
  std::vector<GenTopJet>* jets = event.gentopjets;
  std::vector<GenParticle>* genParts = event.genparticles;
  GenParticle H, W1, W2, Z;

  for( auto & gp : *genParts){
    if(gp.pdgId()==25) H = gp;
    if(gp.pdgId()==24) W1 = gp;
    if(gp.pdgId()== -24) W2 = gp;
    if(gp.pdgId()==23) Z = gp;
  }

  if(jets->size() > 0 && H.pdgId()==25 && W1.pdgId()==24 && W2.pdgId()==-24){
    const auto & jet_1 = jets->at(0);
    DeltaRW1jet_1->Fill(uhh2::deltaR(jet_1, W1), w);
    DeltaRW2jet_1->Fill(uhh2::deltaR(jet_1, W2), w);
    DeltaRHjet_1->Fill(uhh2::deltaR(jet_1, H), w);
    DeltaRW1jet_DeltaRW2jet_1->Fill(uhh2::deltaR(jet_1, W1), uhh2::deltaR(jet_1, W2), w);
    DeltaRHjet_ptjet_1->Fill(uhh2::deltaR(jet_1, H), jet_1.pt(), w);
    DeltaRHjet_Mjet_1->Fill(uhh2::deltaR(jet_1, H), jet_1.v4().M(), w);
    pTH_pTJet_1->Fill(H.pt(), jet_1.pt(), w);
    pTH_pTZ->Fill(H.pt(), Z.pt(), w);
    pTW1_pTW2->Fill(W1.pt(), W2.pt(), w);

  } else {
    // std::cout << "\n @@@ WARNING -- GenParticleHiggsHists::fill -- unexpected number of GenParticles in the event --- check selection. \n";
    // std::cout << H.pdgId() << " " << W1.pdgId() << " " << W2.pdgId() << " " << jets->size() << '\n';
  }

  number_jet->Fill(jets->size(), w);
  for(const auto & jet : *jets){
    mass_jet->Fill(jet.v4().M(), w);
    pt_jet->Fill(jet.pt(), w);
    eta_jet->Fill(jet.eta(), w);

    if(H.pdgId()==25 && W1.pdgId()==24 && W2.pdgId()==-24){
      DeltaRW1jet->Fill(uhh2::deltaR(jet, W1), w);
      DeltaRW2jet->Fill(uhh2::deltaR(jet, W2), w);
      DeltaRHjet->Fill(uhh2::deltaR(jet, H), w);
      DeltaRW1jet_DeltaRW2jet->Fill(uhh2::deltaR(jet, W1), uhh2::deltaR(jet, W2), w);
      DeltaRHjet_ptjet->Fill(uhh2::deltaR(jet, H), jet.pt(), w);
      DeltaRHjet_Mjet->Fill(uhh2::deltaR(jet, H), jet.v4().M(), w);
      pTH_pTJet->Fill(H.pt(), jet.pt(), w);

    } else {
      std::cout << "\n @@@ WARNING -- GenParticleHiggsHists::fill -- unexpected number of GenParticles in the event --- check selection. \n";
      std::cout << H.pdgId() << " " << W1.pdgId() << " " << W2.pdgId() << " " << jets->size() << '\n';
    }

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
    if (jet.tau1()!=0) {
      tau21_jet4->Fill(jet.tau2()/jet.tau1(), w);
      tau31_jet4->Fill(jet.tau3()/jet.tau1(), w);
    }
    if (jet.tau2()!=0) {
      tau32_jet4->Fill(jet.tau3()/jet.tau2(), w);
    }
  }

}

GenJetsHists::~GenJetsHists(){}
