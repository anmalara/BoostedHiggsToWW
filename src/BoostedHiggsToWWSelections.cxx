#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWSelections.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"

#include <stdexcept>
#include <set>

//using namespace uhh2examples;
using namespace uhh2;

#define SETPtLeptonSelection(Lepton,lepton)\
Pt##Lepton##Selection::Pt##Lepton##Selection(float pt_lepton_, float min_isolation_): pt_lepton(pt_lepton_), min_isolation(min_isolation_) {}\
bool Pt##Lepton##Selection::passes(const Event & event){\
  assert(event.lepton);\
  std::vector<Lepton> result;\
  for(const auto & lep : *event.lepton){\
    if (lep.pt() > pt_lepton &&  lep.relIso() < min_isolation) result.push_back(lep);\
  }\
  if (result.size() < 1) return false;\
  std::swap(result, *event.lepton);\
  return true;\
}\

SETPtLeptonSelection(Muon,muons)
SETPtLeptonSelection(Electron,electrons)
////////////////////////////////////////////////////////
//
// PtElectronSelection::PtElectronSelection(float pt_electron_, float min_isolation_): pt_electron(pt_electron_), min_isolation(min_isolation_) {}
//
// bool PtElectronSelection::passes(const Event & event){
//   assert(event.electrons); // if this fails, it probably means electrons are not read in
//
//   std::vector<Electron> result;
//   for(const auto & ele : *event.electrons){
//     if (ele.pt() > pt_electron &&  ele.relIso() < min_isolation) result.push_back(ele);
//   }
//   if (result.size() < 1) return false;
//   std::swap(result, *event.electrons);
//   return true;
//
// }

////////////////////////////////////////////////////////
//
// PtMuonSelection::PtMuonSelection(float pt_muon_, float min_isolation_): pt_muon(pt_muon_), min_isolation(min_isolation_) {}
//
// bool PtMuonSelection::passes(const Event & event){
//   assert(event.muons); // if this fails, it probably means muons are not read in
//
//   std::vector<Muon> result;
//   for(const auto & muon : *event.muons){
//     if (muon.pt() > pt_muon &&  muon.relIso() < min_isolation) result.push_back(muon);
//   }
//   if (result.size() < 1) return false;
//   std::swap(result, *event.muons);
//   return true;
// }

////////////////////////////////////////////////////////

LeptonPtIsoSelection::LeptonPtIsoSelection(float pt_lep_, float min_isolation_, std::string lepton_class_): pt_lep(pt_lep_), min_isolation(min_isolation_), lepton_class(lepton_class_) {}

bool LeptonPtIsoSelection::passes(const Event & event){
  // std::vector< Particle >* leptons;
  // if (lepton_class=="muons") leptons = event.muons;
  // else if (lepton_class=="electrons") leptons = event.electrons;
  // else throw std::runtime_error("Invalid lepton class in LeptonPtIsoSelection.");
  // assert(leptons); // if this fails, it probably means muons are not read in
  //
  // std::vector<Particle> result;
  // for(const auto & lep : *leptons){
  //   if (lep.pt() > pt_lep &&  lep.relIso() < min_isolation) result.push_back(lep);
  // }
  // if (result.size() < 1) return false;
  // std::swap(result, *leptons);
  return true;
}

////////////////////////////////////////////////////////

DiElectronSelection::DiElectronSelection(){}

bool DiElectronSelection::passes(const Event & event){
  assert(event.electrons); // if this fails, it probably means electrons are not read in
  if(event.electrons->size() < 2) return false;
  float min = 10;
  std::set <int> index;
  for (unsigned int i = 0; i < event.electrons->size(); i++) {
    const auto & lep1 = event.electrons->at(i);
    for (unsigned int j = i+1; j < event.electrons->size(); j++) {
      const auto & lep2 = event.electrons->at(j);
      auto diLep = lep1.v4() + lep2.v4();
      if ( (fabs(diLep.M() - 91) < min) && (lep1.charge() + lep2.charge() == 0) ) {index.insert(i); index.insert(j);}
    }
  }
  std::vector<Electron> result;
  for (auto i = index.begin(); i != index.end(); ++i) result.push_back( event.electrons->at(*i) );
  std::swap(result, *event.electrons);
  sort_by_pt<Electron>(*event.electrons);
  if (event.electrons->size() < 2) return false;
  else return true;

}
////////////////////////////////////////////////////////

DiMuonSelection::DiMuonSelection() {}

bool DiMuonSelection::passes(const Event & event){
  assert(event.muons); // if this fails, it probably means muons are not read in
  if(event.muons->size() < 2) return false;
  float min = 10;
  std::set <int> index;
  for (unsigned int i = 0; i < event.muons->size(); i++) {
    const auto & lep1 = event.muons->at(i);
    for (unsigned int j = i+1; j < event.muons->size(); j++) {
      const auto & lep2 = event.muons->at(j);
      auto diLep = lep1.v4() + lep2.v4();
      if ( (fabs(diLep.M() - 91) < min) && (lep1.charge() + lep2.charge() == 0) ) {index.insert(i); index.insert(j);}
    }
  }
  std::vector<Muon> result;
  for (auto i = index.begin(); i != index.end(); ++i) result.push_back( event.muons->at(*i) );
  std::swap(result, *event.muons);
  sort_by_pt<Muon>(*event.muons);
  if (event.muons->size() < 2) return false;
  else return true;
}

////////////////////////////////////////////////////////

PhiAngularCut::PhiAngularCut (float phi_min_, float phi_max_): phi_min(phi_min_), phi_max(phi_max_){}
bool PhiAngularCut::passes(const Event& event){

  assert(event.topjets && event.muons);

  if(event.muons->size() < 2 ){
    throw std::runtime_error("\n @@@ WARNING -- PhiAngularCuts::passes -- unexpected number of muons in the event (!=2). returning 'false'\n");
  }

  if(!event.topjets->size()){
    throw std::runtime_error("\n @@@ WARNING -- PhiAngularCuts::passes -- unexpected number of jets in the event (==0). returning 'false'\n");
  }
  int pos1 = 0, pos2 = 0;
  double min = 10.;
  for (unsigned int i = 0; i < event.muons->size(); i++) {
    const auto & lep1 = event.muons->at(i);
    for (unsigned int j = i+1; j < event.muons->size(); j++) {
      const auto & lep2 = event.muons->at(j);
      auto diLep = lep1.v4() + lep2.v4();
      for(const auto & jet: *event.topjets){
        auto Dphi = deltaPhi(diLep, jet);
        if( (Dphi > phi_min && Dphi< phi_max) && (fabs(diLep.M() - 91) < min) && (lep1.charge() + lep2.charge() == 0) ) {pos1 = i; pos1 = j; min = fabs(diLep.M() - 91);};
      }
    }
  }
  std::vector<Muon> result;
  if (pos1!=0 && pos2!=0) { result.push_back(event.muons->at(pos1)); result.push_back(event.muons->at(pos2)); }
  std::swap(result, *event.muons);
  if (event.muons->size()<2) return false;
  else return true;
}

////////////////////////////////////////////////////////

PhiAngularSelectionElectron::PhiAngularSelectionElectron (float phi_min_, float phi_max_): phi_min(phi_min_), phi_max(phi_max_){}
bool PhiAngularSelectionElectron::passes(const Event& event){
  assert(event.topjets && event.electrons);
  if(event.electrons->size() < 2 ) return false;
  if(event.topjets->size() < 1 ) return false;

  int pos1 = 0, pos2 = 0;
  double min = 10.;
  for (unsigned int i = 0; i < event.electrons->size(); i++) {
    const auto & lep1 = event.electrons->at(i);
    for (unsigned int j = i+1; j < event.electrons->size(); j++) {
      const auto & lep2 = event.electrons->at(j);
      auto diLep = lep1.v4() + lep2.v4();
      for(const auto & jet: *event.topjets){
        auto Dphi = deltaPhi(diLep, jet);
        if( (Dphi > phi_min && Dphi< phi_max) && (fabs(diLep.M() - 91) < min) && (lep1.charge() + lep2.charge() == 0) ) {pos1 = i; pos1 = j; min = fabs(diLep.M() - 91);};
      }
    }
  }
  std::vector<Electron> result;
  if (pos1!=0 || pos2!=0) { result.push_back(event.electrons->at(pos1)); result.push_back(event.electrons->at(pos2)); }
  std::swap(result, *event.electrons);
  if (event.electrons->size()<2) return false;
  else return true;
}

////////////////////////////////////////////////////////

PhiAngularSelectionMuon::PhiAngularSelectionMuon (float phi_min_, float phi_max_): phi_min(phi_min_), phi_max(phi_max_){}
bool PhiAngularSelectionMuon::passes(const Event& event){

  assert(event.topjets && event.muons);
  if(event.muons->size() < 2 ) return false;
  if(event.topjets->size() < 1 ) return false;

  int pos1 = 0, pos2 = 0;
  double min = 10.;
  for (unsigned int i = 0; i < event.muons->size(); i++) {
    const auto & lep1 = event.muons->at(i);
    for (unsigned int j = i+1; j < event.muons->size(); j++) {
      const auto & lep2 = event.muons->at(j);
      auto diLep = lep1.v4() + lep2.v4();
      for(const auto & jet: *event.topjets){
        auto Dphi = deltaPhi(diLep, jet);
        if( (Dphi > phi_min && Dphi< phi_max) && (fabs(diLep.M() - 91) < min) && (lep1.charge() + lep2.charge() == 0) ) {pos1 = i; pos1 = j; min = fabs(diLep.M() - 91);};
      }
    }
  }
  std::vector<Muon> result;
  if (pos1!=0 || pos2!=0) { result.push_back(event.muons->at(pos1)); result.push_back(event.muons->at(pos2)); }
  std::swap(result, *event.muons);
  if (event.muons->size()<2) return false;
  else return true;
}

////////////////////////////////////////////////////////
