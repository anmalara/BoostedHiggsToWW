#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWSelections.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"

#include <stdexcept>
#include <set>

using namespace uhh2;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

#define SETNLeptonPtEtaIsoSelection(Lepton,leptons)                                                   \
N##Lepton##PtEtaIsoSelection:: N##Lepton##PtEtaIsoSelection(int n_lepton_, float pt_lepton_, float eta_lepton_, float min_isolation_): n_lepton(n_lepton_), pt_lepton(pt_lepton_), eta_lepton(eta_lepton_), min_isolation(min_isolation_) {}\
bool N##Lepton##PtEtaIsoSelection::passes(const Event & event){                                       \
  assert(event.leptons);                                                                              \
  int counter = 0;                                                                                    \
  for(const auto & lep : *event.leptons){                                                             \
    if (lep.pt() > pt_lepton && lep.eta() < eta_lepton && lep.relIso() < min_isolation) counter += 1; \
  }                                                                                                   \
  if (counter >= n_lepton) return true;                                                               \
  return false;                                                                                       \
}                                                                                                     \

SETNLeptonPtEtaIsoSelection(Muon,muons)
SETNLeptonPtEtaIsoSelection(Electron,electrons)

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

#define SETJetDiLeptonPhiAngularSelection(Lepton,leptons)\
JetDi##Lepton##PhiAngularSelection::JetDi##Lepton##PhiAngularSelection (float phi_min_, float phi_max_): phi_min(phi_min_), phi_max(phi_max_){}\
bool JetDi##Lepton##PhiAngularSelection::passes(const Event& event){\
  assert(event.JETVARIABLE && event.leptons);\
  if(event.JETVARIABLE->size() < 1 ) throw std::runtime_error("JetDiLeptonPhiAngularSelection::JetDiLeptonPhiAngularSelection -- unexpected number of jets");\
  if(event.leptons->size() < 2 ) throw std::runtime_error("JetDiLeptonPhiAngularSelection::JetDiLeptonPhiAngularSelection -- unexpected number of leptons");\
  for (unsigned int i = 0; i < event.leptons->size(); i++) {\
    const auto & lep1 = event.leptons->at(i);\
    for (unsigned int j = i+1; j < event.leptons->size(); j++) {\
      const auto & lep2 = event.leptons->at(j);\
      auto diLep = lep1.v4() + lep2.v4();\
      if( (fabs(diLep.M() - ZMASS) < ZWIDTH) && (lep1.charge() + lep2.charge() == 0) ) {\
        for(const auto & jet: *event.JETVARIABLE){\
          auto Dphi = deltaPhi(diLep, jet);\
          if( phi_min < Dphi  && Dphi< phi_max ) return true;\
        }\
      }\
    }\
  }\
  return false;\
}\

SETJetDiLeptonPhiAngularSelection(Muon,muons)
SETJetDiLeptonPhiAngularSelection(Electron,electrons)

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
