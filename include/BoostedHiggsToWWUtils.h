#pragma once
#include <iostream>
#include <string>
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"


//================================================================================
bool XOR( bool a, bool b);

//================================================================================

class TopJetLeptonOverlapRemoval: public uhh2::AnalysisModule {
public:
   explicit TopJetLeptonOverlapRemoval(double deltaRmin);
   virtual bool process(uhh2::Event & event) override;

private:
   double deltaRmin;
};

//================================================================================

void sort_topjet_H(std::vector<TopJet> & jets, GenParticle H);
void sort_topjet(std::vector<TopJet> & jets, auto diLep);
void sort_topjet_by_dilepdist(uhh2::Event & event);
