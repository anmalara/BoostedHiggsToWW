  #include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/BoostedHiggsToWW/include/ModuleBase.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWSelections.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWJetHists.h"
#include "UHH2/BoostedHiggsToWW/include/DiLeptonHists.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"
#include "UHH2/BoostedHiggsToWW/include/constants.hpp"
#include "UHH2/BoostedHiggsToWW/include/Macros.hpp"

using namespace std;

/*
█ ██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
█ ██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
█ ██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
█ ██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
█ ██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class GenLevelMatchModule: public ModuleBASE {

public:

  explicit GenLevelMatchModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  // Define variables
  bool is_mc;
  bool muonchannel, electronchannel;
  Event::Handle< float > h_weights_GLP;

  std::shared_ptr<Selection> TopJetHiggsMatch;

};


void GenLevelMatchModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag;
    mytag = "nTopJet_" + tag;   book_HFolder(mytag, new BoostedHiggsToWWJetHists(ctx,mytag,"topjets"));
    mytag = "nJet_" + tag;      book_HFolder(mytag, new BoostedHiggsToWWJetHists(ctx,mytag,"jets"));
    mytag = "ele_" + tag;       book_HFolder(mytag, new ElectronHists(ctx,mytag));
    mytag = "muon_" + tag;      book_HFolder(mytag, new MuonHists(ctx,mytag));
    mytag = "event_" + tag;     book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = "diLepton_" + tag;  book_HFolder(mytag, new DiLeptonHists(ctx,mytag));
  }
}

void GenLevelMatchModule::fill_histograms(uhh2::Event& event, string tag){
  string mytag;
  mytag = "nTopJet_" + tag;   HFolder(mytag)->fill(event);
  mytag = "nJet_" + tag;      HFolder(mytag)->fill(event);
  mytag = "ele_" + tag;       HFolder(mytag)->fill(event);
  mytag = "muon_" + tag;      HFolder(mytag)->fill(event);
  mytag = "event_" + tag;     HFolder(mytag)->fill(event);
  mytag = "diLepton_" + tag;  HFolder(mytag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

GenLevelMatchModule::GenLevelMatchModule(uhh2::Context& ctx){

  // Set up histograms:

  std::vector<std::string> histogram_tags({ "init", "after"});
  book_histograms(ctx, histogram_tags);

  // Set up variables
  is_mc = ctx.get("dataset_type") == "MC";

  muonchannel = string2bool(ctx.get("muonchannel", "true"));
  electronchannel = string2bool(ctx.get("electronchannel", "false"));

  if (muonchannel == electronchannel) throw std::runtime_error("In GenLevelMatchModule.cxx: Choose exactly one lepton channel.");

  TopJetHiggsMatch.reset(new TopJetHiggsMatching());

  h_weights_GLP = ctx.declare_event_input< float >("weights_GLP");

}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool GenLevelMatchModule::process(uhh2::Event& event) {

  event.weight = event.get(h_weights_GLP);

  fill_histograms(event, "init");

  if(TopJetHiggsMatch->passes(event)) fill_histograms(event, "after");

  for (auto& jet: *event.topjets) {
    // jet.set_tag(NN_IsHiggs, 0.);
    // jet.set_tag(NN_IsQCD, 0.);
    // jet.set_tag(NN_IsTop, 0.);
  }


  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the GenLevelMatchModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(GenLevelMatchModule)
