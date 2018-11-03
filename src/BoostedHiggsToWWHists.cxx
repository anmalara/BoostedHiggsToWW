#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

BoostedHiggsToWWHists::BoostedHiggsToWWHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  weights     = book<TH1F>("weights",     "weights",                2000, 0, 20);
  number_jet  = book<TH1F>("number",      "number of jets",         20, 0, 20);
  flavor_jet  = book<TH1F>("flavor_jets", "flavor of jets",         200, -100, 100);
  eta_jet     = book<TH1F>("eta_jet",     "#eta of jets",           40, -3.5, 3.5);
  pt_jet      = book<TH1F>("pt_jet",      "p_{T} of jets [GeV/c]",  100, 0, 1000);

  std::string name_temp, title_temp;
  #define SETTAUHISTOS(n,jet)\
  name_temp  = "tau"; name_temp += #n; name_temp += "_"; name_temp += #jet;\
  title_temp = "#tau"; title_temp = #n; title_temp = "^{"; title_temp += #jet; title_temp += "}";\
  tau##n##_##jet   = book<TH1F>(name_temp, title_temp, 20, 0.,1.);\

  #define SETBOOSTEDHISTOS(jet)\
  name_temp  = "mass_";  name_temp += #jet;\
  title_temp = "mass^{"; title_temp += #jet; title_temp += "} [GeV]";\
  mass_##jet   = book<TH1F>(name_temp, title_temp, 200, 0, 400);\
  name_temp  = "MT_";  name_temp += #jet;\
  title_temp = "M_T^{"; title_temp += #jet; title_temp += "} [GeV]";\
  MT_##jet   = book<TH1F>(name_temp, title_temp, 200, 0, 400);\
  name_temp  = "pt_";  name_temp += #jet;\
  title_temp = "p_T^{"; title_temp += #jet; title_temp += "} [GeV/c]";\
  pt_##jet   = book<TH1F>(name_temp, title_temp, 100, 0, 1000);\
  name_temp  = "eta_";  name_temp += #jet;\
  title_temp = "#eta^{"; title_temp += #jet; title_temp += "}";\
  eta_##jet   = book<TH1F>(name_temp, title_temp, 40, -3.5, 3.5);\
  name_temp  = "phi_";  name_temp += #jet;\
  title_temp = "#phi^{"; title_temp += #jet; title_temp += "}";\
  phi_##jet   = book<TH1F>(name_temp, title_temp, 40, -M_PI, M_PI);\
  name_temp  = "SDmass_";  name_temp += #jet;\
  title_temp = "#SDmass^{"; title_temp += #jet; title_temp += "} [GeV]";\
  SDmass_##jet   = book<TH1F>(name_temp, title_temp, 200, 0, 400);\
  SETTAUHISTOS(1,jet)\
  SETTAUHISTOS(2,jet)\
  SETTAUHISTOS(3,jet)\
  SETTAUHISTOS(21,jet)\
  SETTAUHISTOS(31,jet)\
  SETTAUHISTOS(32,jet)\


  SETBOOSTEDHISTOS(jet1)
  SETBOOSTEDHISTOS(jet2)
  SETBOOSTEDHISTOS(jet3)
  SETBOOSTEDHISTOS(jet4)

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
  weights->Fill(w,1);
  std::vector<TopJet>*  jets = event.topjets;
  int Njets = (jets)->size();

  number_jet->Fill(Njets, w);
  for(const auto & jet : *jets){
    pt_jet->Fill(jet.pt(), w);
    eta_jet->Fill(jet.eta(), w);
    flavor_jet->Fill(fabs(jet.pdgId()), w);
  }

  #define FILLBOOSTEDHISTOS(jet,n)\
  if(jets->size() > n){\
    const auto & Jet = (*jets)[n];\
    mass_##jet->Fill(Jet.v4().M(), w);\
    MT_##jet->Fill(Jet.v4().Mt(), w);\
    pt_##jet->Fill(Jet.pt(), w);\
    eta_##jet->Fill(Jet.eta(), w);\
    phi_##jet->Fill(Jet.phi(), w);\
    tau1_##jet->Fill(Jet.tau1(), w);\
    tau2_##jet->Fill(Jet.tau2(), w);\
    tau3_##jet->Fill(Jet.tau3(), w);\
    SDmass_##jet->Fill(Jet.softdropmass(), w);\
    if (Jet.tau1()!=0) {\
      tau21_##jet->Fill(Jet.tau2()/Jet.tau1(), w);\
      tau31_##jet->Fill(Jet.tau3()/Jet.tau1(), w);\
    }\
    if (Jet.tau2()!=0) tau32_##jet->Fill(Jet.tau3()/Jet.tau2(), w);\
  }\

  FILLBOOSTEDHISTOS(jet1,0)
  FILLBOOSTEDHISTOS(jet2,1)
  FILLBOOSTEDHISTOS(jet3,2)
  FILLBOOSTEDHISTOS(jet4,3)

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
