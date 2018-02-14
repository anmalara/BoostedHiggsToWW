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

PtElecSelection::PtElecSelection(float pt_electron_): pt_electron(pt_electron_){}

bool PtElecSelection::passes(const Event & event){
  assert(event.electrons); // if this fails, it probably means electrons are not read in

  bool skip_event = true;

  for(const auto & electron : *event.electrons){
    skip_event = skip_event && electron.pt() > pt_electron;
    if(!skip_event) continue;
  }

  return skip_event;
}

////////////////////////////////////////////////////////

PtMuonSelection::PtMuonSelection(float pt_muon_): pt_muon(pt_muon_){}

bool PtMuonSelection::passes(const Event & event){
  assert(event.muons); // if this fails, it probably means muons are not read in

  bool skip_event = true;

  for(const auto & muon: *event.muons){
    skip_event = skip_event && muon.pt() > pt_muon;
    if(!skip_event) continue;
  }

  return skip_event;
}

////////////////////////////////////////////////////////

DiElecSelection::DiElecSelection(float min_isolation_): min_isolation(min_isolation_){}

bool DiElecSelection::passes(const Event & event){
  assert(event.electrons); // if this fails, it probably means electron are not read in
  if(event.electrons->size() != 2) return false;
  const auto & ele1 = event.electrons->at(0);
  const auto & ele2 = event.electrons->at(1);
  auto isolation_1 = ele1.relIso();
  auto isolation_2 = ele2.relIso();
  return isolation_1 < min_isolation && isolation_2 < min_isolation;
}

////////////////////////////////////////////////////////

DiMuonSelection::DiMuonSelection(float min_isolation_): min_isolation(min_isolation_){}

bool DiMuonSelection::passes(const Event & event){
  assert(event.muons); // if this fails, it probably means muons are not read in
  if(event.muons->size() != 2) return false;
  const auto & muon1 = event.muons->at(0);
  const auto & muon2 = event.muons->at(1);
  auto isolation_1 = muon1.relIso();
  auto isolation_2 = muon2.relIso();
  return isolation_1 < min_isolation && isolation_2 < min_isolation;
}

////////////////////////////////////////////////////////

ZSelection::ZSelection(){}

bool ZSelection::passes(const Event & event){
  assert(event.muons || event.electrons); // if this fails, it probably means muons/electrons are not read in

  if(!XOR(event.muons->size()==2, event.electrons->size()==2)){
    std::cout << "\n @@@ WARNING -- ZSelection::passes -- unexpected number of muons/electrons in the event (!=2) --- check selection. returning 'false'\n";
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

  auto dilep = lep1.v4() + lep2.v4();
  if(dilep.M() < 81 || dilep.M() > 101) return false;
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
