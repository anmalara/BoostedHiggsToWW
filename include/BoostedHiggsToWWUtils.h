#pragma once
#include <iostream>
#include <string>
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"


//================================================================================
bool XOR( bool a, bool b);

//================================================================================

class JetLeptonOverlapRemoval: public uhh2::AnalysisModule {
public:
   explicit JetLeptonOverlapRemoval(double deltaRmin);
   virtual bool process(uhh2::Event & event) override;

private:
   double deltaRmin;
};

//================================================================================

class GenJetLeptonOverlapRemoval: public uhh2::AnalysisModule {
public:
   explicit GenJetLeptonOverlapRemoval(double deltaRmin);
   virtual bool process(uhh2::Event & event) override;

private:
   double deltaRmin;
};

//================================================================================

void sort_jet_H(std::vector<Jet> & jets, GenParticle H);
void sort_genjet_H(std::vector<Particle> & jets, GenParticle H);
template<typename T> void sort_myjet(std::vector<T> & jets, auto diLep);
void sort_jet_by_dilepdist(uhh2::Event & event);
