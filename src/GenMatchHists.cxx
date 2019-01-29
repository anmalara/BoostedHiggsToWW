#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/BoostedHiggsToWW/include/GenMatchHists.h"
#include "UHH2/BoostedHiggsToWW/include/HistsBase.hpp"

#include "TH1F.h"
#include <type_traits>
#include <iostream>

using namespace std;
using namespace uhh2;


/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/


GenMatchHists::GenMatchHists(Context & ctx, const string & dname): HistsBase(ctx, dname) {

  vector<string> ParticleNames = {"Z", "Higgs", "Wplus", "Wminus"};
  for (auto& name : ParticleNames) {
    book_ParticleHist(name, name, 20, 1500);
    book_FamilyTree(name);
  }

}

void GenMatchHists::fill(const Event & event){
  string name;
  for( auto & gp : *event.genparticles){
    switch (gp.pdgId()) {
      case 23  : name = "Z"; break;
      case 25  : name = "Higgs"; break;
      case 24  : name = "Wplus"; break;
      case -24 : name = "Wminus"; break;
      default: continue;
    }
    fill_ParticleHist(event, name, gp);
    fill_FamilyTree(event, name, gp);
  }
}


/*
█ ██████   ██████   ██████  ██   ██      █████  ███    ██ ██████      ███████ ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ████   ██ ██   ██     ██      ██ ██      ██
█ ██████  ██    ██ ██    ██ █████       ███████ ██ ██  ██ ██   ██     █████   ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ██  ██ ██ ██   ██     ██      ██ ██      ██
█ ██████   ██████   ██████  ██   ██     ██   ██ ██   ████ ██████      ██      ██ ███████ ███████
*/

void GenMatchHists::book_ParticleHist(const string & histSuffix, const string & axisSuffix, double minPt, double maxPt){
  book_TH1F("N_"      +histSuffix, "Number of " +axisSuffix, 21, -0.5, 20.5);
  book_TH1F("pt_"     +histSuffix, "p_{T} "     +axisSuffix, 50,minPt,maxPt);
  book_TH1F("eta_"    +histSuffix, "#eta "      +axisSuffix, 100,-5,5);
  book_TH1F("phi_"    +histSuffix, "#phi "      +axisSuffix, 50,-M_PI,M_PI);
  book_TH1F("mass_"   +histSuffix, "M "         +axisSuffix, 100, 0, 200);
  book_TH1F("flavor_" +histSuffix, "flavor "    +axisSuffix, 101, -50.5, 50.5);
  return;
}

void GenMatchHists::fill_ParticleHist(const Event & event, const string & histSuffix, const GenParticle& gp) {
  auto weight = event.weight;
  fill_H1("N_"      +histSuffix, 1,           weight);
  fill_H1("pt_"     +histSuffix, gp.pt(),     weight);
  fill_H1("eta_"    +histSuffix, gp.eta(),    weight);
  fill_H1("phi_"    +histSuffix, gp.phi(),    weight);
  fill_H1("mass_"   +histSuffix, gp.v4().M(), weight);
  fill_H1("flavor_" +histSuffix, gp.pdgId(),  weight);
  return;
}


void GenMatchHists::book_FamilyTree(const string & ParticleName) {
  book_ParticleHist("mother1"   +ParticleName, "mother1"  +ParticleName, 20,1500);
  book_ParticleHist("mother2"   +ParticleName, "mother2"  +ParticleName, 20,1500);
  book_ParticleHist("daughter1" +ParticleName, "daughter1"+ParticleName, 20,1500);
  book_ParticleHist("daughter2" +ParticleName, "daughter2"+ParticleName, 20,1500);
  book_TH1F("deltaR_daughters"  +ParticleName, "#Delta R(daughters) "+ParticleName, 40, 0, M_PI);
  book_TH1F("deltaR_daughter1"  +ParticleName, "#Delta R(daughter1) "+ParticleName, 40, 0, M_PI);
  book_TH1F("deltaR_daughter2"  +ParticleName, "#Delta R(daughter2) "+ParticleName, 40, 0, M_PI);
}

void GenMatchHists::fill_FamilyTree(const Event & event, const string & ParticleName, const GenParticle& gp) {
  if (gp.mother(event.genparticles, 1))   fill_ParticleHist(event, "mother1"   +ParticleName, *gp.mother(event.genparticles, 1));
  if (gp.mother(event.genparticles, 2))   fill_ParticleHist(event, "mother2"   +ParticleName, *gp.mother(event.genparticles, 2));
  if (gp.daughter(event.genparticles, 1)) fill_ParticleHist(event, "daughter1" +ParticleName, *gp.daughter(event.genparticles, 1));
  if (gp.daughter(event.genparticles, 2)) fill_ParticleHist(event, "daughter2" +ParticleName, *gp.daughter(event.genparticles, 2));
  if (gp.daughter(event.genparticles, 1)) fill_H1("deltaR_daughter1" +ParticleName, deltaR(gp,*gp.daughter(event.genparticles, 1)), event.weight);
  if (gp.daughter(event.genparticles, 2)) fill_H1("deltaR_daughter2" +ParticleName, deltaR(gp,*gp.daughter(event.genparticles, 2)), event.weight);
  if (gp.daughter(event.genparticles, 1) && gp.daughter(event.genparticles, 2)) fill_H1("deltaR_daughters" +ParticleName, deltaR(*gp.daughter(event.genparticles, 1),*gp.daughter(event.genparticles, 2)), event.weight);
}
