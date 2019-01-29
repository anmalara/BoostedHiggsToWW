#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/BoostedHiggsToWW/include/constants.hpp"
#include "UHH2/BoostedHiggsToWW/include/HistsBase.hpp"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWJetHists.h"
#include "TString.h"

#include <vector>
#include <string>


class JetHistosWithConditions: public BoostedHiggsToWWJetHists {
public:
  JetHistosWithConditions(uhh2::Context& ctx, const std::string& dname, const std::string& collection, Conditions condition_=ZMatch, float Dr_ = 1.6, const unsigned int NumberOfPlottedJets=4): BoostedHiggsToWWJetHists(ctx, dname, collection, NumberOfPlottedJets), condition(condition_), Dr(Dr_) {};
  virtual void fill(const uhh2::Event&) override;
protected:
  template <typename T>
  void fill_internal(const uhh2::Event&, std::vector<T>);
  bool isLeptonic(int);
  template <typename T>
  bool MatchGenPart(const uhh2::Event&, T , int , Decay decay=nodecay);
  Decay DobleDecay(int, int);
  Conditions condition;
  float Dr;
};
