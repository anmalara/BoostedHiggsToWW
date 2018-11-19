#include "UHH2/core/include/Event.h"
#include "UHH2/BoostedHiggsToWW/include/DiLeptonHists.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

DiLeptonHists::DiLeptonHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  dilep_number  = book<TH1F>("number",        "number of dilepton",   20, 0, 20);

  std::string name_temp;
  #define SETDILEPTONHISTOS(lepton)                                                     \
  name_temp = "di"; name_temp += #lepton; name_temp += "_number";                       \
  di##lepton##_number   = book<TH1F>(name_temp, name_temp, 20, 0, 20);                  \
  name_temp = "di"; name_temp += #lepton; name_temp += "_charge";                       \
  di##lepton##_charge   = book<TH1F>(name_temp, name_temp, 5,-2,2);                     \
  name_temp = "di"; name_temp += #lepton; name_temp += "_m";                            \
  di##lepton##_m        = book<TH1F>(name_temp, name_temp, 100,0,100);                  \
  name_temp = "di"; name_temp += #lepton; name_temp += "_pt";                           \
  di##lepton##_pt       = book<TH1F>(name_temp, name_temp, 100,0,500);                  \
  name_temp = "di"; name_temp += #lepton; name_temp += "_eta";                          \
  di##lepton##_eta      = book<TH1F>(name_temp, name_temp, 100,-3,3);                   \
  name_temp = "di"; name_temp += #lepton; name_temp += "_phi";                          \
  di##lepton##_phi      = book<TH1F>(name_temp, name_temp, 100,-M_PI,M_PI);             \
  name_temp = "di"; name_temp += #lepton; name_temp += "_DR12";                         \
  di##lepton##_deltaR   = book<TH1F>(name_temp, name_temp, 40, 0, 2.0);                 \
  name_temp = "di"; name_temp += #lepton; name_temp += "_jet_Dphi";                     \
  di##lepton##_jet_Dphi = book<TH1F>(name_temp, name_temp, 100,-M_PI,M_PI);             \
  name_temp = "pt_";name_temp += #lepton; name_temp += "_2D";                           \
  pt_##lepton##_2D      = book<TH2F>(name_temp,  ";PT1;pT2", 100, 0, 500, 100, 0, 500.);\

  SETDILEPTONHISTOS(Electron)
  SETDILEPTONHISTOS(Muon)
}


void DiLeptonHists::fill(const Event & event){
  // fill the histograms.

  // Don't forget to always use the weight when filling.
  auto weight = event.weight;

  #define  FILLDILEPTONHISTOS(lepton,leptons)                       \
  dilep_number->Fill(leptons->size(), weight);                      \
  if (leptons->size()>=min_leptons) {                               \
    di##lepton##_number->Fill(leptons->size(), weight);             \
    const auto& lep1 = leptons->at(0);                              \
    const auto& lep2 = leptons->at(1);                              \
    auto dilep = lep1.v4() + lep2.v4();                             \
    di##lepton##_charge->Fill(lep1.charge()+lep2.charge(), weight); \
    di##lepton##_deltaR->Fill(uhh2::deltaR(lep1, lep2), weight);    \
    if ( event.JETVARIABLE->size()>0 ) {                            \
      const auto& jet = event.JETVARIABLE->at(0);                   \
      di##lepton##_jet_Dphi->Fill(deltaPhi(dilep, jet), weight);    \
    }                                                               \
    di##lepton##_m->Fill(dilep.M(), weight);                        \
    di##lepton##_pt->Fill(dilep.Pt(), weight);                      \
    di##lepton##_eta->Fill(dilep.Eta(), weight);                    \
    di##lepton##_phi->Fill(dilep.Phi(), weight);                    \
    pt_##lepton##_2D->Fill(lep1.pt(), lep2.pt(), weight);           \
  }                                                                 \

  FILLDILEPTONHISTOS(Electron,event.electrons)
  FILLDILEPTONHISTOS(Muon,event.muons)

}

DiLeptonHists::~DiLeptonHists(){}
