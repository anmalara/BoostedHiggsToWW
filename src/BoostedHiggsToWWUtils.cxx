#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"

using namespace uhh2;
using namespace std;

//================================================================================
bool XOR( bool a, bool b){
  return (!a&&b) || (!b && a);
}

//================================================================================

JetLeptonOverlapRemoval::JetLeptonOverlapRemoval(double deltaRmin_): deltaRmin(deltaRmin_){}

bool JetLeptonOverlapRemoval::process(Event & event){

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

GenJetLeptonOverlapRemoval::GenJetLeptonOverlapRemoval(double deltaRmin_): deltaRmin(deltaRmin_){}
bool GenJetLeptonOverlapRemoval::process(Event & event){

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

void sort_jet_H(std::vector<Jet> & jets, GenParticle H){
  std::sort(jets.begin(), jets.end(), [H](const Jet &jet1, const Jet &jet2 ){return uhh2::deltaR(H, jet1) < uhh2::deltaR(H, jet2);});
}

void sort_genjet_H(std::vector<Particle> & jets, GenParticle H){
  std::sort(jets.begin(), jets.end(), [H](const Particle &jet1, const Particle &jet2 ){return uhh2::deltaR(H, jet1) < uhh2::deltaR(H, jet2);});
}

template<typename T>
void sort_myjet(std::vector<T> & jets, auto diLep){
  std::sort(jets.begin(), jets.end(), [diLep](const T &jet1, const T &jet2 ){return uhh2::deltaR(diLep, jet1) > uhh2::deltaR(diLep, jet2);});
}

void sort_jet_by_dilepdist(uhh2::Event & event){

  Electron ele1 = event.electrons->at(0);
  Electron ele2 = event.electrons->at(1);
  Muon muon1 = event.muons->at(0);
  Muon muon2 = event.muons->at(1);

  Particle lep1, lep2;

  if(event.muons->size()==2 && event.electrons->size()==2){
    if ( fabs((muon1.v4()+muon2.v4()).M()-91) <= fabs((ele1.v4()+ele2.v4()).M()-91)) { lep1= muon1; lep2= muon2;}
    else { lep1= ele1; lep2= ele2;}
  } else if (event.muons->size()==2)  { lep1= muon1; lep2= muon2; }
  else if (event.electrons->size()==2) { lep1= ele1; lep2= ele2; }
  // else return false;

  auto diLep = lep1.v4() + lep2.v4();
  sort_myjet<TopJet>(*event.topjets, diLep);
}
