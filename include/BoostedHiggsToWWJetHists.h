#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/BoostedHiggsToWW/include/constants.hpp"
#include "UHH2/BoostedHiggsToWW/include/HistsBase.hpp"
#include "TString.h"

#include <vector>
#include <string>


class BoostedHiggsToWWJetHists: public HistsBase {
public:
  BoostedHiggsToWWJetHists(uhh2::Context&, const std::string&, const std::string& ,const unsigned int NumberOfPlottedJets=4);
  virtual void fill(const uhh2::Event&) override;
protected:
  std::string collection;
  unsigned int NumberOfPlottedJets;
  std::vector<std::string> axis_suffix {"first jet","second jet","third jet","fourth jet"};
  bool isTop= false, isGen=false;
  void book_jetHist(const std::string&, const std::string&, double, double, bool, bool);
  template <typename T>
  void fill_jetHist(const uhh2::Event&,const std::string&, const T&);
  template <typename T>
  void fill_internal(const uhh2::Event&, std::vector<T>);
};
