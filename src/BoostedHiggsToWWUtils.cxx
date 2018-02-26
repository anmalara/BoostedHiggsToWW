#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"

using namespace uhh2;
using namespace std;

//================================================================================
bool XOR( bool a, bool b){
  return (!a&&b) || (!b && a);
}

//================================================================================

TopJetLeptonOverlapRemoval::TopJetLeptonOverlapRemoval(double deltaRmin_): deltaRmin(deltaRmin_){}

bool TopJetLeptonOverlapRemoval::process(Event & event){

  assert(event.topjets);
  std::vector<TopJet> result;

  if (event.muons) {
    for(const auto & jet : *event.topjets){
      bool skip_jet = true;
      for(const auto & lepton : *event.muons) skip_jet = skip_jet && deltaR(jet, lepton) > deltaRmin;
      if (skip_jet) result.push_back(jet);
    }
    std::swap(result, *event.topjets);
    result.clear();
  }

  if(event.electrons) {
    for(const auto & jet : *event.topjets){
      bool skip_jet = true;
      for(const auto & lepton : *event.electrons) skip_jet = skip_jet && deltaR(jet, lepton) > deltaRmin;
      if (skip_jet) result.push_back(jet);
    }
    std::swap(result, *event.topjets);
  }

  return true;
}

GenTopJetLeptonOverlapRemoval::GenTopJetLeptonOverlapRemoval(double deltaRmin_): deltaRmin(deltaRmin_){}
bool GenTopJetLeptonOverlapRemoval::process(Event & event){

  assert(event.gentopjets);
  std::vector<GenTopJet> result;

  if (event.genparticles) {
    for(const auto & jet : *event.gentopjets){
      bool skip_jet = true;
      for(const auto & lepton : *event.genparticles){ if ( XOR(abs(lepton.pdgId()) == 13, abs(lepton.pdgId()) == 11) ) skip_jet = skip_jet && deltaR(jet, lepton) > deltaRmin; }
      if (skip_jet) result.push_back(jet);
    }
    std::swap(result, *event.gentopjets);
  }

  return true;
}


//================================================================================

void sort_topjet_H(std::vector<TopJet> & jets, GenParticle H){
  std::sort(jets.begin(), jets.end(), [H](const TopJet &jet1, const TopJet &jet2 ){return uhh2::deltaR(H, jet1) < uhh2::deltaR(H, jet2);});
}

void sort_gentopjet_H(std::vector<GenTopJet> & jets, GenParticle H){
  std::sort(jets.begin(), jets.end(), [H](const GenTopJet &jet1, const GenTopJet &jet2 ){return uhh2::deltaR(H, jet1) < uhh2::deltaR(H, jet2);});
}

void sort_topjet(std::vector<TopJet> & jets, auto diLep){
  std::sort(jets.begin(), jets.end(), [diLep](const TopJet &jet1, const TopJet &jet2 ){return deltaR(diLep, jet1) > deltaR(diLep, jet2);});
}

void sort_topjet_by_dilepdist(uhh2::Event & event){

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
  sort_topjet(*event.topjets, diLep);
}
