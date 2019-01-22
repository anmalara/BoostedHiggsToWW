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
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWHists.h"
#include "UHH2/BoostedHiggsToWW/include/DiLeptonHists.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"
#include "UHH2/BoostedHiggsToWW/include/constants.hpp"
#include "UHH2/BoostedHiggsToWW/include/Macros.hpp"
#include "UHH2/BoostedHiggsToWW/include/NeuralNetworkBase.hpp"


using namespace std;

/*
█ ██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
█ ██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
█ ██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
█ ██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
█ ██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class NeuralNetworkModule: public ModuleBASE {

public:

  explicit NeuralNetworkModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  // Define variables
  bool is_mc;
  bool muonchannel, electronchannel;
  Event::Handle< float > h_lumi_weights, h_my_weights, h_pu_weights;
  Event::Handle< float > h_lumi_weights_in, h_my_weights_in, h_pu_weights_in;
  Event::Handle< float > h_isHiggs, h_isQCD, h_isTop;

  std::shared_ptr<NeuralNetworkBase> NeuralNetwork;

};


void NeuralNetworkModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag;
    mytag = "nTopJet_" + tag;   book_HFolder(mytag, new BoostedHiggsToWWHists(ctx,mytag,"topjets"));
    mytag = "nJet_" + tag;      book_HFolder(mytag, new BoostedHiggsToWWHists(ctx,mytag,"jets"));
    mytag = "ele_" + tag;       book_HFolder(mytag, new ElectronHists(ctx,mytag));
    mytag = "muon_" + tag;      book_HFolder(mytag, new MuonHists(ctx,mytag));
    mytag = "event_" + tag;     book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = "diLepton_" + tag;  book_HFolder(mytag, new DiLeptonHists(ctx,mytag));
  }
}

void NeuralNetworkModule::fill_histograms(uhh2::Event& event, string tag){
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

NeuralNetworkModule::NeuralNetworkModule(uhh2::Context& ctx){

  // Set up histograms:

  std::vector<std::string> histogram_tags({ "beforeNN", "afterNN"});
  book_histograms(ctx, histogram_tags);

  // Set up variables
  is_mc = ctx.get("dataset_type") == "MC";

  muonchannel = string2bool(ctx.get("muonchannel", "true"));
  electronchannel = string2bool(ctx.get("electronchannel", "false"));

  if (muonchannel == electronchannel) throw std::runtime_error("In NeuralNetworkModule.cxx: Choose exactly one lepton channel.");

  // if (muonchannel) JetDiLeptonPhiAngularSel.reset(new JetDiMuonPhiAngularSelection(min_delta_phi));
  // if (electronchannel) JetDiLeptonPhiAngularSel.reset(new JetDiElectronPhiAngularSelection(min_delta_phi));

  // unsorted_jets = ctx.declare_event_output< std::vector< TopJet > >("unsorted_jets");

  h_lumi_weights_in = ctx.declare_event_input< float >("lumi_weights");
  h_my_weights_in = ctx.declare_event_input< float >("my_weights");
  h_pu_weights_in = ctx.declare_event_input< float >("weight_pu");

  h_lumi_weights = ctx.declare_event_output< float >("lumi_weights_out");
  h_my_weights = ctx.declare_event_output< float >("my_weights_out");
  h_pu_weights = ctx.declare_event_output< float >("weight_pu_out");


  h_isHiggs = ctx.declare_event_output< float >("isHiggs");
  h_isQCD = ctx.declare_event_output< float >("isQCD");
  h_isTop = ctx.declare_event_output< float >("isTop");

  NeuralNetwork.reset(new NeuralNetworkBase("/nfs/dust/cms/user/amalara/WorkingArea/File/NeuralNetwork/model_300epochs_1500k/"));


}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool NeuralNetworkModule::process(uhh2::Event& event) {

  float lumi_weights_in = event.get(h_lumi_weights_in);
  float my_weights_in = event.get(h_my_weights_in);
  float pu_weights_in = event.get(h_pu_weights_in);

  event.weight = my_weights_in;

  fill_histograms(event, "beforeNN");

  //
  // Matrix2D inputs (1, {0.17794486,   -1.02782187,    0.25901741,    0.49207153, 0.12045357,    0.96785651,    0.66286739,    0.46307791, 2.85150328,  104.        ,    0.23433149,    0.17218957, 0.23612824});
  //
  // std::cout << "Input" << '\n';
  // for (size_t i = 0; i < inputs[0].size(); i++) {
  //   std::cout << inputs[0][i] << '\t';
  // }
  // std::cout << '\n';

  std::cout << "New Event " <<  (*event.topjets).size() << '\n';
  for(const auto & jet : *event.topjets){
    std::cout << "jet pt " << jet.pt() << '\n';
    NeuralNetwork->Apply(jet);
    Matrix2D outputs = NeuralNetwork->getOutputs();
    event.set(h_isHiggs, outputs[0][0]);
    event.set(h_isQCD, outputs[0][1]);
    event.set(h_isTop, outputs[0][2]);
  }
  //
  // std::cout << "Output" << '\n';
  // for (size_t i = 0; i < inputs[0].size(); i++) {
  //   std::cout << inputs[0][i] << '\t';
  // }
  // std::cout << '\n \n';
  //
  // for(const auto & jet : *event.topjets){
  //   Matrix2D inputs (1, {jet.pt(), jet.eta(), jet.phi(), jet.v4().M(),jet.energy(), jet.tau1(), jet.tau2(), jet.tau3(), jet.chargedMultiplicity()+jet.neutralMultiplicity(), jet.btag_combinedSecondaryVertex(), jet.tau2()/jet.tau1(), jet.tau3()/jet.tau1(), jet.tau3()/jet.tau2()});
  //   std::cout << "Input" << '\n';
  //   for (size_t i = 0; i < inputs[0].size(); i++) {
  //     std::cout << inputs[0][i] << '\t';
  //   }
  //   std::cout << '\n';
  //
  //   NeuralNetwork->Apply(inputs);
  //
  //   std::cout << "Output" << '\n';
  //   for (size_t i = 0; i < inputs[0].size(); i++) {
  //     std::cout << inputs[0][i] << '\t';
  //   }
  //   std::cout << '\n \n';
  // }

  fill_histograms(event, "afterNN");

  event.set(h_lumi_weights, lumi_weights_in);
  event.set(h_my_weights, my_weights_in);
  event.set(h_pu_weights, pu_weights_in);


  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the NeuralNetworkModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(NeuralNetworkModule)
