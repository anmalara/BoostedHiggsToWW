#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWSelections.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"
#include "UHH2/BoostedHiggsToWW/include/constants.hpp"

#include <stdexcept>
#include <set>

//using namespace uhh2examples;
using namespace uhh2;

#define SETLeptonPtIsoSelection(Lepton,leptons)\
Lepton##PtIsoSelection:: Lepton##PtIsoSelection(float pt_lepton_, float min_isolation_): pt_lepton(pt_lepton_), min_isolation(min_isolation_) {}\
bool Lepton##PtIsoSelection::passes(const Event & event){\
  assert(event.leptons);\
  std::vector<Lepton> result;\
  for(const auto & lep : *event.leptons){\
    if (lep.pt() > pt_lepton &&  lep.relIso() < min_isolation) result.push_back(lep);\
  }\
  std::swap(result, *event.leptons);\
  if (event.leptons->size() < 1) return false;\
  return true;\
}\

SETLeptonPtIsoSelection(Muon,muons)
SETLeptonPtIsoSelection(Electron,electrons)
////////////////////////////////////////////////////////

#define  SETDiLeptonSelection(Lepton,leptons)\
Di##Lepton##Selection::Di##Lepton##Selection(){}\
bool Di##Lepton##Selection::passes(const Event & event){\
  assert(event.leptons);\
  if(event.leptons->size() < 2) return false;\
  std::set <int> index;\
  for (unsigned int i = 0; i < event.leptons->size(); i++) {\
    const auto & lep1 = event.leptons->at(i);\
    for (unsigned int j = i+1; j < event.leptons->size(); j++) {\
      const auto & lep2 = event.leptons->at(j);\
      auto diLep = lep1.v4() + lep2.v4();\
      if ( (fabs(diLep.M() - ZMASS) < ZWIDTH) && (lep1.charge() + lep2.charge() == 0) ) {index.insert(i); index.insert(j);}\
    }\
  }\
  std::vector<Lepton> result;\
  for (auto i = index.begin(); i != index.end(); ++i) result.push_back( event.leptons->at(*i) );\
  std::swap(result, *event.leptons);\
  sort_by_pt<Lepton>(*event.leptons);\
  if (event.leptons->size() < 2) return false;\
  else return true;\
}\

SETDiLeptonSelection(Muon,muons)
SETDiLeptonSelection(Electron,electrons)

////////////////////////////////////////////////////////

#define SETJetLeptonPhiAngularSelection(Lepton,leptons)\
Jet##Lepton##PhiAngularSelection::Jet##Lepton##PhiAngularSelection (float phi_min_, float phi_max_): phi_min(phi_min_), phi_max(phi_max_){}\
bool Jet##Lepton##PhiAngularSelection::passes(const Event& event){\
  assert(event.JETVARIABLE && event.leptons);\
  if(event.JETVARIABLE->size() < 1 ) return false;\
  if(event.leptons->size() < 2 ) return false;\
  int pos1 = 0, pos2 = 0;\
  double min = ZWIDTH;\
  for (unsigned int i = 0; i < event.leptons->size(); i++) {\
    const auto & lep1 = event.leptons->at(i);\
    for (unsigned int j = i+1; j < event.leptons->size(); j++) {\
      const auto & lep2 = event.leptons->at(j);\
      auto diLep = lep1.v4() + lep2.v4();\
      for(const auto & jet: *event.JETVARIABLE){\
        auto Dphi = deltaPhi(diLep, jet);\
        if( (Dphi > phi_min && Dphi< phi_max) && (fabs(diLep.M() - ZMASS) < min) && (lep1.charge() + lep2.charge() == 0) ) {pos1 = i; pos2 = j; min = fabs(diLep.M() - ZMASS);};\
      }\
    }\
  }\
  std::vector<Lepton> result;\
  if (pos1!=0 || pos2!=0) { result.push_back(event.leptons->at(pos1)); result.push_back(event.leptons->at(pos2)); }\
  std::swap(result, *event.leptons);\
  if (event.leptons->size()<2) return false;\
  else return true;\
}\

SETJetLeptonPhiAngularSelection(Muon,muons)
SETJetLeptonPhiAngularSelection(Electron,electrons)

////////////////////////////////////////////////////////
