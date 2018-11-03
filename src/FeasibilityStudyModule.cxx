#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWSelections.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWHists.h"
#include "UHH2/BoostedHiggsToWW/include/DiLeptonHists.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"

using namespace std;
using namespace uhh2;

class FeasibilityStudyModule: public AnalysisModule {
public:

  explicit FeasibilityStudyModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // Define variables
  string SysType_PU;
  bool isMC;

  // Define histograms
  #define  DEFINEHISTOS(extention)\
  std::unique_ptr<Hists> h_nJet_##extention, h_Electron_##extention, h_Muon_##extention, h_DiLepton_##extention;\

  DEFINEHISTOS(nocuts)
  DEFINEHISTOS(cleaned)
  DEFINEHISTOS(trigger)

  // CommonModules
  std::unique_ptr<CommonModules> common;
  // std::unique_ptr<JetLeptonOverlapRemoval> overlap_removal;

  // Define selections
  std::unique_ptr<uhh2::Selection> trigger_sel;

  std::shared_ptr<Selection> NBoostedJetSel;

  DEFINEHISTOS(NBoostedJetSel)

  #define DEFINELEPTONSELECTION(Lepton)\
  std::unique_ptr<uhh2::AndSelection> Lepton##Sel;\
  std::shared_ptr<Selection> N##Lepton##Sel, Pt##Lepton##Sel, Di##Lepton##Sel, PhiAngularSel##Lepton;\
  DEFINEHISTOS(N##Lepton)\
  DEFINEHISTOS(Pt##Lepton)\
  DEFINEHISTOS(Di##Lepton)\
  DEFINEHISTOS(Jet##Lepton)\
  DEFINEHISTOS(PhiAngularSel##Lepton)\


  DEFINELEPTONSELECTION(Muon)
  DEFINELEPTONSELECTION(Electron)
  DEFINEHISTOS(DiLepton)
  DEFINEHISTOS(eventSel)


  // DEFINEHISTOS(PhiAngularSel)

  // Event::Handle< std::vector< TopJet > > unsorted_jets;

};


FeasibilityStudyModule::FeasibilityStudyModule(Context & ctx){

  // Set up histograms:
  std::string name_temp;
  #define SETHISTOS(extension)\
  name_temp = "nJet_"; name_temp += #extension;\
  h_nJet_##extension.reset(new BoostedHiggsToWWHists(ctx, name_temp.c_str()));\
  name_temp = "ele_"; name_temp += #extension;\
  h_Electron_##extension.reset(new ElectronHists(ctx, name_temp));\
  name_temp = "muon_"; name_temp += #extension;\
  h_Muon_##extension.reset(new MuonHists(ctx, name_temp));\
  name_temp = "diLepton_"; name_temp += #extension;\
  h_DiLepton_##extension.reset(new DiLeptonHists(ctx, name_temp));\

  SETHISTOS(nocuts)

  // Set up variables
  SysType_PU = ctx.get("SysType_PU");
  isMC = (ctx.get("dataset_type") == "MC");
  const std::string& trigger = ctx.get("trigger", "NULL");

  const JetId jetId(AndId<Jet> (JetPFID(JetPFID::WP_TIGHT), PtEtaCut(30, 2.4)));
  const TopJetId JetIdBoosted(AndId<TopJet> (JetPFID(JetPFID::WP_TIGHT), PtEtaCut(300,2.4)));
  const MuonId muonId(AndId<Muon> (MuonID(Muon::CutBasedIdTight), PtEtaCut(1., 4.)));
  const ElectronId eleId(AndId<Electron>(ElectronID_MVA_Fall17_loose_iso, PtEtaSCCut(1., 4.)));

  // CommonModules
  common.reset(new CommonModules());
  common->disable_lumisel();
  // common->disable_mclumiweight();
  // common->disable_jec();
  common->disable_jersmear();
  if (isMC) common->disable_metfilters();
  // common->disable_mcpileupreweight();
  common->switch_metcorrection();
  common->disable_jetpfidfilter();
  common->set_jet_id(jetId);
  common->set_muon_id(muonId);
  common->set_electron_id(eleId);
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->init(ctx,SysType_PU);

  SETHISTOS(cleaned)

  // overlap_removal.reset(new JetLeptonOverlapRemoval(0.8));

  // Set up selections
  if(!isMC && trigger != "NULL") trigger_sel.reset(new TriggerSelection(trigger));
  else trigger_sel.reset(new AndSelection(ctx));
  SETHISTOS(trigger)

  NBoostedJetSel.reset(new NTopJetSelection(1,-1,JetIdBoosted));
  SETHISTOS(NBoostedJetSel)

  #define SETLEPTONSELECTION(Lepton)\
  N##Lepton##Sel.reset(new N##Lepton##Selection(2));\
  SETHISTOS(N##Lepton)\
  Pt##Lepton##Sel.reset(new Lepton##PtIsoSelection(30, 0.2));\
  SETHISTOS(Pt##Lepton)\
  Di##Lepton##Sel.reset(new Di##Lepton##Selection());\
  SETHISTOS(Di##Lepton)\
  name_temp = #Lepton; name_temp += "_filters";\
  Lepton##Sel.reset(new uhh2::AndSelection(ctx, name_temp));\
  name_temp = "2 "; name_temp += #Lepton;\
  Lepton##Sel->add(name_temp, N##Lepton##Sel);\
  name_temp = "pt_"; name_temp += #Lepton;\
  Lepton##Sel->add(name_temp, Pt##Lepton##Sel);\
  name_temp = "di"; name_temp += #Lepton;\
  Lepton##Sel->add(name_temp, Di##Lepton##Sel);\
  Lepton##Sel->add("nboostedjet >=2 ", NBoostedJetSel);\
  SETHISTOS(Jet##Lepton)\
  PhiAngularSel##Lepton.reset(new Jet##Lepton##PhiAngularSelection(2.7));\
  SETHISTOS(PhiAngularSel##Lepton)\
  Lepton##Sel->add("#Delta#phi(jet,dilep) >2.7", PhiAngularSel##Lepton);\

  SETLEPTONSELECTION(Muon)
  SETLEPTONSELECTION(Electron)

  SETHISTOS(DiLepton)
  SETHISTOS(eventSel)

  // unsorted_jets = ctx.declare_event_output< std::vector< TopJet > >("unsorted_jets");

}


bool FeasibilityStudyModule::process(Event & event) {

  #define FILLHISTOS(extension)\
  h_nJet_##extension->fill(event);\
  h_Electron_##extension->fill(event);\
  h_Muon_##extension->fill(event);\
  h_DiLepton_##extension->fill(event);\

  FILLHISTOS(nocuts)
  // 1. run all modules other modules.
  common->process(event);
  FILLHISTOS(cleaned)
  // overlap_removal->process(event);

  sort_by_pt<TopJet>(*event.topjets);
  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);

  const bool pass_trigger = trigger_sel->passes(event);
  if(!pass_trigger) return false;
  FILLHISTOS(trigger)

  // boosted Jet Selection
  bool JetPass = true;
  bool NBoostedJetPass = NBoostedJetSel->passes(event);
  JetPass = JetPass && NBoostedJetPass;
  if(JetPass) FILLHISTOS(NBoostedJetSel)
  sort_by_pt<TopJet>(*event.topjets);

  // Lepton Selection

  #define PASSLEPTONSELECTION(Lepton)\
  bool Lepton##Pass = true;\
  bool N##Lepton##Pass = N##Lepton##Sel->passes(event);\
  Lepton##Pass = Lepton##Pass && N##Lepton##Pass;\
  if(Lepton##Pass) {FILLHISTOS(N##Lepton)}\
  bool Pt##Lepton##Pass = Pt##Lepton##Sel->passes(event);\
  Lepton##Pass = Lepton##Pass && Pt##Lepton##Pass;\
  if(Lepton##Pass) {FILLHISTOS(Pt##Lepton)}\
  bool Di##Lepton##Pass = Di##Lepton##Sel->passes(event);\
  Lepton##Pass = Lepton##Pass && Di##Lepton##Pass;\
  if(Lepton##Pass) {FILLHISTOS(Di##Lepton)}\
  Lepton##Pass = Lepton##Pass && JetPass;\
  if(Lepton##Pass) {FILLHISTOS(Jet##Lepton)}\
  bool PhiAngularSel##Lepton##Pass = PhiAngularSel##Lepton->passes(event);\
  Lepton##Pass = Lepton##Pass && PhiAngularSel##Lepton##Pass;\
  if(Lepton##Pass) {FILLHISTOS(PhiAngularSel##Lepton)}\
  bool Total##Lepton##Pass = Lepton##Sel->passes(event);\
  if (Lepton##Pass != Total##Lepton##Pass) throw std::runtime_error("Recheck selection!");\

  PASSLEPTONSELECTION(Muon)
  PASSLEPTONSELECTION(Electron)

  bool eventPass = MuonPass || ElectronPass;
  if (eventPass) FILLHISTOS(DiLepton)

  // std::vector< Jet > my_unsorted_jets = *event.topjets;
  // event.set(unsorted_jets, std::move(my_unsorted_jets));
  // sort_topjet_by_dilepdist(event);

  eventPass = XOR(MuonPass,ElectronPass);
  if (!eventPass) return false;
  FILLHISTOS(eventSel)
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the FeasibilityStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(FeasibilityStudyModule)
