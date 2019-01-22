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

#include <UHH2/BoostedHiggsToWW/include/ModuleBase.h>
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWSelections.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWHists.h"
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

class GenericCleaningModule: public ModuleBASE {

public:

  explicit GenericCleaningModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  // Define variables
  bool is_mc, lumisel, mclumiweight, mcpileupreweight, jec, topjec, jersmear, topjersmear, eleid, muid, tauid, jetpfidcleaner, topjetpfidcleaner, metfilters, jetlepcleaner, topjetlepcleaner, jetid, topjetid, do_metcorrection;
  bool muonchannel, electronchannel;
  string SysType_PU, JEC_Version, JEC_LABEL;


  // Define common modules
  std::unique_ptr<uhh2::Selection> lumi_selection;
  std::vector<std::unique_ptr<AnalysisModule>> modules;
  std::unique_ptr<uhh2::AndSelection> metfilters_selection;
  // std::unique_ptr<uhh2::Selection> trigger_sel;

  #define  DEFINEJEC(type)                                        \
  std::vector<std::string> JEC_corr_##type##_AK4PFchs;            \
  std::vector<std::string> JEC_corr_##type##_AK8PFchs;            \
  std::vector<std::string> JEC_corr_##type##_AK4PFPuppi;          \
  std::vector<std::string> JEC_corr_##type##_AK8PFPuppi;          \
  std::unique_ptr<JetCorrector> jet_corrector_##type;             \
  std::unique_ptr<TopJetCorrector> topjet_corrector_##type;       \
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> JLC_##type;    \
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> TopJLC_##type; \

  DEFINEJEC(MC)
  DEFINEJEC(B)
  DEFINEJEC(C)
  DEFINEJEC(D)
  DEFINEJEC(E)
  DEFINEJEC(F)

  std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;
  std::unique_ptr<GenericJetResolutionSmearer> topjet_resolution_smearer;
  std::unique_ptr<JetCleaner> jet_cleaner;
  std::unique_ptr<TopJetCleaner> topjet_cleaner;

  // Define selections

  std::shared_ptr<Selection> NBoostedJetSel, NoLeptonSel;
  std::shared_ptr<VetoSelection> VetoLeptonSel;

  Event::Handle< float > h_weights_lumi, h_weights_GLP, h_weight_pu;

};


void GenericCleaningModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag;
    mytag = "nTopJet_" + tag;  book_HFolder(mytag, new BoostedHiggsToWWHists(ctx,mytag,"topjets"));
    mytag = "nJet_" + tag;     book_HFolder(mytag, new BoostedHiggsToWWHists(ctx,mytag,"jets"));
    mytag = "ele_" + tag;      book_HFolder(mytag, new ElectronHists(ctx,mytag));
    mytag = "muon_" + tag;     book_HFolder(mytag, new MuonHists(ctx,mytag));
    mytag = "event_" + tag;    book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = "diLepton_" + tag; book_HFolder(mytag, new DiLeptonHists(ctx,mytag));
  }
}

void GenericCleaningModule::fill_histograms(uhh2::Event& event, string tag){
  string mytag;
  mytag = "nTopJet_" + tag;  HFolder(mytag)->fill(event);
  mytag = "nJet_" + tag;     HFolder(mytag)->fill(event);
  mytag = "ele_" + tag;      HFolder(mytag)->fill(event);
  mytag = "muon_" + tag;     HFolder(mytag)->fill(event);
  mytag = "event_" + tag;    HFolder(mytag)->fill(event);
  mytag = "diLepton_" + tag; HFolder(mytag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

GenericCleaningModule::GenericCleaningModule(uhh2::Context& ctx){

  // Set up histograms:

  std::vector<std::string> histogram_tags({ "nocuts", "cleaned", "Veto", "NBoostedJet"});
  book_histograms(ctx, histogram_tags);

  // Set up variables
  is_mc = ctx.get("dataset_type") == "MC";
  JEC_Version = ctx.get("JEC_Version", "Fall17_17Nov2017_V6");

  lumisel = string2bool(ctx.get("lumisel", "true"));
  mclumiweight = string2bool(ctx.get("mclumiweight", "true"));
  mcpileupreweight = string2bool(ctx.get("mcpileupreweight", "true"));
  SysType_PU = ctx.get("SysType_PU");
  jec = string2bool(ctx.get("jec", "true"));
  topjec = string2bool(ctx.get("topjec", "true"));
  jersmear = string2bool(ctx.get("jersmear", "true"));
  topjersmear = string2bool(ctx.get("topjersmear", "true"));
  do_metcorrection = string2bool(ctx.get("do_metcorrection", "true"));
  eleid = string2bool(ctx.get("eleid", "true"));
  muid = string2bool(ctx.get("muid", "true"));
  tauid = string2bool(ctx.get("tauid", "false"));
  jetid = string2bool(ctx.get("jetid", "true"));
  topjetid = string2bool(ctx.get("topjetid", "true"));
  jetpfidcleaner = true; // WP fixed TODO
  topjetpfidcleaner = true; // WP fixed TODO
  metfilters = string2bool(ctx.get("metfilters", "true"));
  jetlepcleaner = string2bool(ctx.get("jetlepcleaner", "true"));
  topjetlepcleaner = string2bool(ctx.get("topjetlepcleaner", "true"));

  muonchannel = string2bool(ctx.get("muonchannel", "true"));
  electronchannel = string2bool(ctx.get("electronchannel", "false"));

  if (muonchannel == electronchannel) throw std::runtime_error("In GenericCleaningModule.cxx: Choose exactly one lepton channel.");

  const MuonId muId(AndId<Muon> (MuonID(Muon::CutBasedIdTight), PtEtaCut(min_pt_lepton, max_eta_lepton), MuonIso(iso_lepton)));
  const ElectronId eleId(AndId<Electron>(ElectronID_Fall17_loose, PtEtaSCCut(min_pt_lepton, max_eta_lepton)));
  const JetId jetId(AndId<Jet> (JetPFID(JetPFID::WP_TIGHT), PtEtaCut(30, max_eta_lepton)));
  const TopJetId jetIdBoosted(AndId<TopJet> (JetPFID(JetPFID::WP_TIGHT), PtEtaCut(300, max_eta_lepton)));

  MAKE_JEC(Fall17_17Nov2017_V6, AK4PFchs)
  MAKE_JEC(Fall17_17Nov2017_V6, AK8PFchs)
  // MAKE_JEC(Fall17_17Nov2017_V6, AK4PFPuppi)
  // MAKE_JEC(Fall17_17Nov2017_V6, AK8PFPuppi)


  PrimaryVertexId pvid = StandardPrimaryVertexId(); // TODO
  modules.emplace_back(new PrimaryVertexCleaner(pvid));
  if(is_mc) {
    ctx.declare_event_input<std::vector<Particle> >(ctx.get("TopJetCollectionGEN"), "topjetsGEN");
    if(mclumiweight)  modules.emplace_back(new MCLumiWeight(ctx));
    if(mcpileupreweight) modules.emplace_back(new MCPileupReweight(ctx,SysType_PU));
    if(jec) jet_corrector_MC.reset(new JetCorrector(ctx, JEC_corr_MC_AK4PFchs));
    if(topjec) topjet_corrector_MC.reset(new TopJetCorrector(ctx, JEC_corr_MC_AK8PFchs));
    if(jersmear) jet_resolution_smearer.reset(new JetResolutionSmearer(ctx, JERSmearing::SF_13TeV_Summer16_25nsV1));
    if(topjersmear) topjet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "topjetsGEN", JERSmearing::SF_13TeV_Summer16_25nsV1, "Fall17_25nsV1_MC_PtResolution_AK8PFchs.txt"));
  } else {
    if(lumisel) lumi_selection.reset(new LumiSelection(ctx));
    if(jec){ SET_JET_CORRECTION(JetCorrector, AK8PFchs, jet) }
    if(topjec){ SET_JET_CORRECTION(TopJetCorrector, AK8PFchs, topjet) }
  }

  if(eleid) modules.emplace_back(new ElectronCleaner(eleId));
  if(muid)  modules.emplace_back(new MuonCleaner(muId));
  // if(tauid) modules.emplace_back(new TauCleaner(tauId));
  if(jetpfidcleaner) modules.emplace_back(new JetCleaner(ctx, JetPFID(JetPFID::WP_TIGHT))); // TODO
  if(topjetpfidcleaner) modules.emplace_back(new TopJetCleaner(ctx, JetPFID(JetPFID::WP_TIGHT), "topjets"));

  // modules.emplace_back(new HTCalculator(ctx,HT_jetid));

  if(!is_mc && metfilters){ SET_MET_FILTERS() }

  if(jetlepcleaner){ SET_JETLEPTON_CLEANER(JLC, AK4PFchs, jets) }
  if(topjetlepcleaner){ SET_JETLEPTON_CLEANER(TopJLC, AK8PFchs, topjets) }

  if(jetid) jet_cleaner.reset(new JetCleaner(ctx, jetId));
  if(topjetid) topjet_cleaner.reset(new TopJetCleaner(ctx, jetIdBoosted, "topjets"));


  // Set up selections

  NBoostedJetSel.reset(new NTopJetSelection(1));
  if (muonchannel) NoLeptonSel.reset(new NElectronSelection(1));
  if (electronchannel) NoLeptonSel.reset(new NMuonSelection(1));

  VetoLeptonSel.reset(new VetoSelection(NoLeptonSel));

  h_weights_lumi = ctx.declare_event_output< float >("weights_lumi");
  if (is_mc) h_weight_pu = ctx.get_handle< float >("weight_pu");
  else h_weight_pu = ctx.declare_event_output< float >("weight_pu");
  h_weights_GLP = ctx.declare_event_output< float >("weights_GLP");

}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool GenericCleaningModule::process(uhh2::Event& event) {

  bool apply_B = false, apply_C = false, apply_D = false, apply_E = false, apply_F = false;

  if(event.isRealData){
    if(event.run <= s_runnr_A) throw std::runtime_error("In GenericCleaningModule.cxx: Run number for RunA used.");
    else if(event.run <= s_runnr_B) apply_B = true;
    else if(event.run <= s_runnr_C) apply_C = true;
    else if(event.run <= s_runnr_D) apply_D = true;
    else if(event.run <= s_runnr_E) apply_E = true;
    else if(event.run <= s_runnr_F) apply_F = true;
    else throw std::runtime_error("In GenericCleaningModule.cxx: Run number not covered by if-statements in process-routine.");
    if(apply_B+apply_C+apply_D+apply_E+apply_F != 1) throw std::runtime_error("In GenericCleaningModule.cxx: Sum of apply_* when applying JECs is not == 1. Fix this.");
  }

  fill_histograms(event, "nocuts");

  // COMMON MODULES

  auto weight_gen = event.weight;
  if(event.isRealData && lumisel) if(!lumi_selection->passes(event)) return false;

  for(auto & m : modules){
    m->process(event);
  }

  double weight_pu = 1;
  if (is_mc) weight_pu = event.get(h_weight_pu);
  else event.set(h_weight_pu, weight_pu);
  event.set(h_weights_lumi, event.weight/(weight_gen*weight_pu));

  if(event.isRealData && metfilters) if(!metfilters_selection->passes(event)) return false;

  if(jetlepcleaner){    APPLY_PROCESS(JLC,process) }
  if(topjetlepcleaner){ APPLY_PROCESS(TopJLC,process) }
  if(jec){              APPLY_PROCESS(jet_corrector,process) }
  if(topjec){           APPLY_PROCESS(topjet_corrector,process) }
  if(jersmear && is_mc) jet_resolution_smearer->process(event);
  if(topjersmear && is_mc) topjet_resolution_smearer->process(event);
  if(do_metcorrection){ APPLY_PROCESS(jet_corrector,correct_met) }
  if(jetid) jet_cleaner->process(event);
  if(topjetid) topjet_cleaner->process(event);

  fill_histograms(event, "cleaned");

  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);
  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);

  if(! VetoLeptonSel->passes(event)) return false;
  fill_histograms(event, "Veto");

  if(! NBoostedJetSel->passes(event)) return false;
  fill_histograms(event, "NBoostedJet");

  event.set(h_weights_GLP, event.weight);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the GenericCleaningModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(GenericCleaningModule)
