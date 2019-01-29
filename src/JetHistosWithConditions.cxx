#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/BoostedHiggsToWW/include/HistsBase.hpp"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWJetHists.h"
#include "UHH2/BoostedHiggsToWW/include/JetHistosWithConditions.h"

using namespace std;
using namespace uhh2;


/*
█ ██████   ██████   ██████  ██   ██      █████  ███    ██ ██████      ███████ ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ████   ██ ██   ██     ██      ██ ██      ██
█ ██████  ██    ██ ██    ██ █████       ███████ ██ ██  ██ ██   ██     █████   ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ██  ██ ██ ██   ██     ██      ██ ██      ██
█ ██████   ██████   ██████  ██   ██     ██   ██ ██   ████ ██████      ██      ██ ███████ ███████
*/

bool JetHistosWithConditions::isLeptonic(int pdgId) { return (pdgId>= 11 && pdgId<=19) ? true : false; }

Decay JetHistosWithConditions::DobleDecay(int pdgId1, int pdgId2) {
  if (  isLeptonic(pdgId1) &&  isLeptonic(pdgId2) )  return leptonic;
  if (  isLeptonic(pdgId1) && !isLeptonic(pdgId2) )  return semileptonic;
  if ( !isLeptonic(pdgId1) &&  isLeptonic(pdgId2) )  return semileptonic;
  if ( !isLeptonic(pdgId1) && !isLeptonic(pdgId2) )  return hadronic;
  throw logic_error("Impossible case");
}

template <typename T>
bool JetHistosWithConditions::MatchGenPart(const Event & event, T jet, int pdgId, Decay decay) {
  for (auto& gp: *event.genparticles) {
    if (gp.pdgId()==pdgId && deltaR(jet,gp)<Dr) {
      if (decay == nodecay) return true;
      GenParticle W1 = *gp.daughter(event.genparticles, 1);
      GenParticle W2 = *gp.daughter(event.genparticles, 2);
      if ( fabs(W1.pdgId())==24 && fabs(W2.pdgId()) == 24 ){
        int pdgId1 = fabs(W1.daughter(event.genparticles, 1)->pdgId());
        int pdgId2 = fabs(W2.daughter(event.genparticles, 1)->pdgId());
        if (decay == DobleDecay(pdgId1, pdgId2) ) return true;
        return false;
      }
    }
  }
  return false;
}

void JetHistosWithConditions::fill(const Event & event){
  if (!isGen && !isTop) fill_internal(event, *event.jets);
  if (!isGen &&  isTop) fill_internal(event, *event.topjets);
  // if ( isGen && !isTop) fill_internal(*event.genjets, weight, isTop, isGen);
  // if ( isGen &&  isTop) fill_internal(*event.topgenjets, weight, isTop, isGen);

}

template <typename T>
void JetHistosWithConditions::fill_internal(const Event & event, vector<T> jets){
  auto weight = event.weight;
  fill_H1("number", jets.size(), weight);
  fill_H1("weights", weight, 1);
  for(unsigned int i = 0; i <jets.size(); i++){
    bool found = false;
    T& jet = jets[i];
    switch (condition) {
      case ZMatch:    found = MatchGenPart(event, jet, 23); break;
      case HMatch:    found = MatchGenPart(event, jet, 25); break;
      case WWfullLep: found = MatchGenPart(event, jet, 25,leptonic); break;
      case WWfullHad: found = MatchGenPart(event, jet, 25,hadronic); break;
      case WWsemiLep: found = MatchGenPart(event, jet, 25,semileptonic); break;
      default: found = false;
    }
    if (!found) continue;
    fill_jetHist(event, "_jet",jet);
    if(i<4) fill_jetHist(event, "_"+to_string(i+1),jet);
    else if (i<NumberOfPlottedJets) fill_jetHist(event, "_"+to_string(i+1),jet);
    if(i < 2){
      auto next_jet = closestParticle(jet, jets);
      auto drmin = next_jet ? deltaR(jet, *next_jet) : numeric_limits<float>::infinity();
      if(i==0) fill_H1("deltaRmin_1", drmin, weight);
      else fill_H1("deltaRmin_2", drmin, weight);
    }
  }
}
