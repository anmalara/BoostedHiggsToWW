#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

//using namespace uhh2examples;
using namespace uhh2;

////////////////////////////////////////////////////////

BoostedJetSelection::BoostedJetSelection(float pt_min_): pt_min(pt_min_){}

bool BoostedJetSelection::passes(const Event & event){
  assert(event.topjets); // if this fails, it probably means topjets are not read in
  if(event.topjets->size() < 1) return false;
  const auto & jet1 = event.topjets->at(0);
  if(jet1.pt() < pt_min) return false;
  return true;
}

////////////////////////////////////////////////////////

DijetSelection::DijetSelection(float dphi_min_, float third_frac_max_): dphi_min(dphi_min_), third_frac_max(third_frac_max_){}

bool DijetSelection::passes(const Event & event){
  assert(event.topjets); // if this fails, it probably means topjets are not read in
  if(event.topjets->size() < 2) return false;
  const auto & jet0 = event.topjets->at(0);
  const auto & jet1 = event.topjets->at(1);
  auto dphi = deltaPhi(jet0, jet1);
  if(dphi < dphi_min) return false;
  if(event.topjets->size() == 2) return true;
  const auto & jet2 = event.topjets->at(2);
  auto third_jet_frac = jet2.pt() / (0.5 * (jet0.pt() + jet1.pt()));
  return third_jet_frac < third_frac_max;
}

////////////////////////////////////////////////////////

PtElecSelection::PtElecSelection(float pt_electron_, float min_isolation_): pt_electron(pt_electron_), min_isolation(min_isolation_) {}

bool PtElecSelection::passes(const Event & event){
  assert(event.electrons); // if this fails, it probably means electrons are not read in
  // std::sort_by_pt<Electron>(*event.electrons);
  // if (event.electrons->size() < 2) return false;
  // if ((event.electrons)->at(0)).pt() < pt_electron) return false;
  // if ((event.electrons)->at(1)).pt() < pt_electron) return false;

  std::vector<Electron> result;
  for(const auto & ele : *event.electrons){
    if (ele.pt() > pt_electron &&  ele.relIso() < min_isolation) result.push_back(ele);
  }
  if (result.size() < 1) return false;
  std::swap(result, *event.electrons);
  return true;

}

////////////////////////////////////////////////////////

PtMuonSelection::PtMuonSelection(float pt_muon_, float min_isolation_): pt_muon(pt_muon_), min_isolation(min_isolation_) {}

bool PtMuonSelection::passes(const Event & event){
  assert(event.muons); // if this fails, it probably means muons are not read in

  // sort_by_pt<Muon>(*event.muons);
  // if (event.muons->size() < 2) return false;
  // if ((event.muons)->at(0)).pt() < pt_muon) return false;
  // if ((event.muons)->at(1)).pt() < pt_muon) return false;

  std::vector<Muon> result;
  for(const auto & muon : *event.muons){
    if (muon.pt() > pt_muon &&  muon.relIso() < min_isolation) result.push_back(muon);
  }
  if (result.size() < 1) return false;
  std::swap(result, *event.muons);
  return true;
}

////////////////////////////////////////////////////////

DiElecSelection::DiElecSelection(){}

bool DiElecSelection::passes(const Event & event){
  assert(event.electrons); // if this fails, it probably means electrons are not read in
  if(event.electrons->size() < 2) return false;
  std::vector<Electron> result;

  for (unsigned int i = 0; i < event.electrons->size(); i++) {
    const auto & lep1 = event.electrons->at(i);
    bool isFound = false;
    for (unsigned int j = 0; j < event.electrons->size(); j++) {
      if (j != i) {
        const auto & lep2 = event.electrons->at(j);
        auto diLep = lep1.v4() + lep2.v4();
        if (diLep.M() > 81 && diLep.M() < 101 && (lep1.charge() + lep2.charge() == 0) ) isFound = true;
        if (isFound) continue;
      }
    }
    if (isFound) result.push_back(lep1);
  }
  std::swap(result, *event.electrons);
  if (event.electrons->size() > 2) {
    result.clear();
    float min = 10;
    int pos1 = 0, pos2 = 0;
    for (unsigned int i = 0; i < event.electrons->size(); i++) {
      const auto & lep1 = event.electrons->at(i);
      for (unsigned int j = 0; j < event.electrons->size(); j++) {
        if (j != i) {
          const auto & lep2 = event.electrons->at(j);
          auto diLep = lep1.v4() + lep2.v4();
          if ( fabs(diLep.M() - 91) < min ) {min = fabs(diLep.M() - 91); pos1 = i; pos2 = j;  }
        }
      }
    }
    if (pos1 != 0 || pos2 != 0){ result.push_back(event.electrons->at(pos1));  result.push_back(event.electrons->at(pos2));}
    std::swap(result, *event.electrons);
  }
  if (event.electrons->size() < 2) return false;
  if(event.electrons->size() == 2) return true;
  std::cout << "WARNING " << __LINE__  << '\n';
  return false;
}
////////////////////////////////////////////////////////

DiMuonSelection::DiMuonSelection() {}

bool DiMuonSelection::passes(const Event & event){
  assert(event.muons); // if this fails, it probably means muons are not read in
  if(event.muons->size() < 2) return false;
  std::vector<Muon> result;

  for (unsigned int i = 0; i < event.muons->size(); i++) {
    const auto & lep1 = event.muons->at(i);
    bool isFound = false;
    for (unsigned int j = 0; j < event.muons->size(); j++) {
      if (j != i) {
        const auto & lep2 = event.muons->at(j);
        auto diLep = lep1.v4() + lep2.v4();
        if (diLep.M() > 81 && diLep.M() < 101 && (lep1.charge() + lep2.charge() == 0) ) isFound = true;
        if (isFound) continue;
      }
    }
    if (isFound) result.push_back(lep1);
  }
  std::swap(result, *event.muons);
  if (event.muons->size() > 2) {
    result.clear();
    float min = 10;
    int pos1 = 0, pos2 = 0;
    for (unsigned int i = 0; i < event.muons->size(); i++) {
      const auto & lep1 = event.muons->at(i);
      for (unsigned int j = 0; j < event.muons->size(); j++) {
        if (j != i) {
          const auto & lep2 = event.muons->at(j);
          auto diLep = lep1.v4() + lep2.v4();
          if ( fabs(diLep.M() - 91) < min ) {min = fabs(diLep.M() - 91); pos1 = i; pos2 = j;  }
        }
      }
    }
    if (pos1 != 0 || pos2 != 0){ result.push_back(event.muons->at(pos1));  result.push_back(event.muons->at(pos2));}
    std::swap(result, *event.muons);
  }
  if (event.muons->size() < 2) return false;
  if(event.muons->size() == 2) return true;
  std::cout << "WARNING " << __LINE__  << '\n';
  return false;
}

////////////////////////////////////////////////////////

ZSelection::ZSelection(){}

bool ZSelection::passes(const Event & event){
  assert(event.muons || event.electrons); // if this fails, it probably means muons/electrons are not read in

  if (event.muons->size() != 2 && event.electrons->size() != 2 ){
    std::cout << event.muons->size() << " "  << event.electrons->size() << '\n';
    std::cout << "\n @@@ WARNING -- ZSelection::passes -- unexpected number of muons/electrons in the event (>2) --- check selection. returning 'false'\n";
    return false;
  }

  // if(!XOR(event.muons->size()==2, event.electrons->size()==2)){
  if(event.muons->size()==2 && event.electrons->size()==2){
    Electron ele1 = event.electrons->at(0);
    Electron ele2 = event.electrons->at(1);
    Muon muon1 = event.muons->at(0);
    Muon muon2 = event.muons->at(1);
    if ( fabs((event.muons->at(0).v4()+event.muons->at(1).v4()).M()-91) <= fabs((event.electrons->at(0).v4()+event.electrons->at(1).v4()).M()-91)){
      std::vector<Electron> result; std::swap(result, *event.electrons);
    } else{
      std::vector<Muon> result; std::swap(result, *event.muons);
    }
  }

  Particle lep1, lep2;

  if (event.electrons->size()==2){
    lep1= event.electrons->at(0);
    lep2= event.electrons->at(1);
  }

  if (event.muons->size()==2){
    lep1= event.muons->at(0);
    lep2= event.muons->at(1);
  }

  auto diLep = lep1.v4() + lep2.v4();
  if(diLep.M() < 81 || diLep.M() > 101 || (lep1.charge() + lep2.charge() != 0) ) return false;
  return true;
}

////////////////////////////////////////////////////////

PhiAngularCut::PhiAngularCut (float phi_min_, float phi_max_): phi_min(phi_min_), phi_max(phi_max_){}
bool PhiAngularCut::passes(const Event& event){

  assert(event.topjets && (event.muons || event.electrons));

  if(!XOR(event.muons->size()==2, event.electrons->size()==2)){
    std::cout << "\n @@@ WARNING -- PhiAngularCuts::passes -- unexpected number of muons/electrons in the event (!=2). returning 'false'\n";
    return false;
  }

  if(!event.topjets->size()){
    std::cout << "\n @@@ WARNING -- PhiAngularCuts::passes -- unexpected number of topjets in the event (==0). returning 'false'\n";
    return false;
  }

  Particle lep1, lep2;

  if (event.electrons->size()==2){
    lep1= event.electrons->at(0);
    lep2= event.electrons->at(1);
  }

  if (event.muons->size()==2){
    lep1= event.muons->at(0);
    lep2= event.muons->at(1);
  }

  auto diLep = lep1.v4() + lep2.v4();
  bool skip_jet = false;

  for(const auto & jet: *event.topjets){
    auto Dphi = deltaPhi(diLep, jet);
    skip_jet = skip_jet || (Dphi > phi_min && Dphi< phi_max);
    if(skip_jet) continue;
  }

  return skip_jet;
}

////////////////////////////////////////////////////////
