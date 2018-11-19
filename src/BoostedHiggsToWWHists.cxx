#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWHists.h"
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


BoostedHiggsToWWHists::BoostedHiggsToWWHists(Context & ctx, const std::string & dname, const string & collection_, const unsigned int NumberOfPlottedJets_): HistsBase(ctx, dname), collection(collection_), NumberOfPlottedJets(NumberOfPlottedJets_){
  if (collection.find("top") != std::string::npos || collection.find("Top") != std::string::npos) isTop = true;
  if (collection.find("gen") != std::string::npos || collection.find("Gen") != std::string::npos) isGen = true;
  book_TH1F("number","number of jets",21, -.5, 20.5);
  book_TH1F("weights","weights",2000, 0, 20);
  book_jetHist("_jet", "jet",20,1500,isTop,isGen);
  vector<double> minPt {20,20,20,20};
  vector<double> maxPt {1500,1000,500,350};
  for(unsigned int i =0; i<NumberOfPlottedJets; i++){
    if(i<4) book_jetHist("_"+to_string(i+1),axis_suffix[i],minPt[i],maxPt[i],isTop,isGen);
    else book_jetHist("_"+to_string(i+1), to_string(i+1)+"-th jet",20,500,isTop,isGen);
  }
  book_TH1F("deltaRmin_1", "#Delta R_{min}(1st jet,nearest jet)", 40, 0, 8.0);
  book_TH1F("deltaRmin_2", "#Delta R_{min}(2nd jet,nearest jet)", 40, 0, 8.0);
}

void BoostedHiggsToWWHists::fill(const Event & event){
  if (!isGen && !isTop) fill_internal(event, *event.jets);
  if (!isGen &&  isTop) fill_internal(event, *event.topjets);
  // if ( isGen && !isTop) fill_internal(*event.genjets, weight, isTop, isGen);
  // if ( isGen &&  isTop) fill_internal(*event.topgenjets, weight, isTop, isGen);

}


/*
█ ██████   ██████   ██████  ██   ██      █████  ███    ██ ██████      ███████ ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ████   ██ ██   ██     ██      ██ ██      ██
█ ██████  ██    ██ ██    ██ █████       ███████ ██ ██  ██ ██   ██     █████   ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ██  ██ ██ ██   ██     ██      ██ ██      ██
█ ██████   ██████   ██████  ██   ██     ██   ██ ██   ████ ██████      ██      ██ ███████ ███████
*/


template <typename T>
void BoostedHiggsToWWHists::fill_internal(const Event & event, vector<T> jets){
  auto weight = event.weight;
  fill_H1("number", jets.size(), weight);
  fill_H1("weights", weight, 1);
  for(unsigned int i = 0; i <jets.size(); i++){
    T& jet = jets[i];
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

void BoostedHiggsToWWHists::book_jetHist(const string & histSuffix, const string & axisSuffix, double minPt, double maxPt, bool isTop, bool isGen){
  book_TH1F("mass"+histSuffix,"M "+axisSuffix+" [GeV/c^{2}]",100,0,300);
  book_TH1F("mT"+histSuffix,"m_{T} "+axisSuffix+" [GeV/c^{2}]",100,0,300);
  book_TH1F("pt"+histSuffix,"p_{T} "+axisSuffix+" [GeV]",50,minPt,maxPt);
  book_TH1F("eta"+histSuffix,"#eta "+axisSuffix,100,-5,5);
  book_TH1F("phi"+histSuffix,"#phi "+axisSuffix,50,-M_PI,M_PI);
  book_TH1F("csv"+histSuffix,"csv-disriminator "+axisSuffix,50,0,1);
  book_TH1F("flavor"+histSuffix,"flavor "+axisSuffix,200,-100,100);
  book_TH1F("deltaR_lepton"+histSuffix,"#Delta R(lepton,"+axisSuffix+")",40, 0, 8.0);
  if (isTop) {
    book_TH1F("SDmass"+histSuffix,"SDmass^{"+axisSuffix+"} [GeV/c^{2}]",100,0,300);
    book_TH1F("tau1"+histSuffix,"#tau1^{"+axisSuffix+"}",20,0,1);
    book_TH1F("tau2"+histSuffix,"#tau2^{"+axisSuffix+"}",20,0,1);
    book_TH1F("tau3"+histSuffix,"#tau3^{"+axisSuffix+"}",20,0,1);
    book_TH1F("tau21"+histSuffix,"#tau21^{"+axisSuffix+"}",20,0,1);
    book_TH1F("tau31"+histSuffix,"#tau31^{"+axisSuffix+"}",20,0,1);
    book_TH1F("tau32"+histSuffix,"#tau32^{"+axisSuffix+"}",20,0,1);
    book_TH1F("nsubjet"+histSuffix,"nsubjet^{"+axisSuffix+"}",41, -.5, 40.5);
  }
  return;
}

template <>
void BoostedHiggsToWWHists::fill_jetHist<Jet>(const Event & event, const string& histSuffix, const Jet& jet){
  auto weight = event.weight;
  fill_H1("mass"+histSuffix, jet.v4().M(), weight);
  fill_H1("mT"+histSuffix, jet.v4().Mt(), weight);
  fill_H1("pt"+histSuffix, jet.pt(), weight);
  fill_H1("eta"+histSuffix, jet.eta(), weight);
  fill_H1("phi"+histSuffix, jet.phi(), weight);
  fill_H1("csv"+histSuffix, jet.btag_combinedSecondaryVertex(), weight);
  fill_H1("flavor"+histSuffix, jet.pdgId(), weight);
  for (const auto & muon : *event.muons ) fill_H1("deltaR_lepton"+histSuffix, deltaR(jet, muon), weight);
  for (const auto & electron : *event.electrons ) fill_H1("deltaR_lepton"+histSuffix, deltaR(jet, electron), weight);
  return;
}

template <>
void BoostedHiggsToWWHists::fill_jetHist<TopJet>(const Event & event, const string& histSuffix, const TopJet& jet){
  auto weight = event.weight;
  fill_jetHist<Jet>(event, histSuffix,jet);
  fill_H1("SDmass"+histSuffix, jet.softdropmass(), weight);
  fill_H1("tau1"+histSuffix, jet.tau1(), weight);
  fill_H1("tau2"+histSuffix, jet.tau2(), weight);
  fill_H1("tau3"+histSuffix, jet.tau3(), weight);
  if (jet.tau1()!=0) {
    fill_H1("tau21"+histSuffix, jet.tau2()/jet.tau1(), weight);
    fill_H1("tau31"+histSuffix, jet.tau3()/jet.tau1(), weight);
  }
  if (jet.tau2()!=0) {
    fill_H1("tau32"+histSuffix, jet.tau3()/jet.tau2(), weight);
  }
  fill_H1("nsubjet"+histSuffix, jet.subjets().size(), weight);
  return;
}
