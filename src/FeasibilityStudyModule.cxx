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

namespace uhh2examples {

  /** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
  *
  * This is the central class which calls other AnalysisModules, Hists or Selection classes.
  * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
  */
  class FeasibilityStudyModule: public AnalysisModule {
  public:

    explicit FeasibilityStudyModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<JetCleaner> 	jetcleaner;
    std::unique_ptr<MuonCleaner>	muoncleaner;
    std::unique_ptr<ElectronCleaner>	elecleaner;
    std::unique_ptr<TopJetLeptonOverlapRemoval> overlap_removal;

    // Selections
    std::unique_ptr<Selection> nJet_sel, boostedJet_sel;
    std::unique_ptr<Selection> ptMuon_sel, ptEle_sel, nMuon_sel, nEle_sel, diMuon_sel, diEle_sel, Z_sel;
    std::unique_ptr<Selection> phiAngular_sel;

    //std::unique_ptr<uhh2::Selection> trigger_sel;

    // Histos
    std::unique_ptr<Hists> h_nJet_nocuts, h_ele_nocuts, h_muon_nocuts, h_diLepton_nocuts;
    std::unique_ptr<Hists> h_nJet_cleaned, h_ele_cleaned, h_muon_cleaned, h_diLepton_cleaned;
    std::unique_ptr<Hists> h_nJet, h_boostedJet_sel;
    std::unique_ptr<Hists> h_ptmuon_sel, h_ptele_sel, h_nmuon_sel, h_nele_sel, h_dimuon_sel, h_diele_sel, h_dilep_sel;
    std::unique_ptr<Hists> h_dimuon_Z_sel, h_diele_Z_sel, h_dilep_Z_sel;
    std::unique_ptr<Hists> h_jet_angular_sel, h_dilep_angular_sel;

    Event::Handle< std::vector< TopJet > > unsorted_jets;
  };


  FeasibilityStudyModule::FeasibilityStudyModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.

    // If needed, access the configuration of the module here, e.g.:
    // string testvalue = ctx.get("TestKey", "<not set>");
    // cout << "TestKey in the configuration was: " << testvalue << endl;

    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    // for(auto & kv : ctx.get_all()){
    //   cout << " " << kv.first << " = " << kv.second << endl;
    // }

    // 1. setup other modules. CommonModules and the JetCleaner:
    common.reset(new CommonModules());
    // TODO: configure common here, e.g. by
    // calling common->set_*_id or common->disable_*
    common->disable_lumisel();
    common->disable_mclumiweight();
    // common->disable_jec();
    common->disable_mcpileupreweight();
    // const bool isMC = (ctx.get("dataset_type") == "MC");
    if (ctx.get("dataset_type") == "MC") {
      common->disable_metfilters();
    }
    common->init(ctx);

    jetcleaner.reset(new JetCleaner(ctx, 30.0, 2.4));
    muoncleaner.reset(new MuonCleaner (AndId<Muon>(MuonIDMedium(), PtEtaCut(1., 4))));
    elecleaner.reset(new ElectronCleaner (AndId<Electron>(ElectronID_Spring16_medium, PtEtaCut(1., 4))));
    overlap_removal.reset(new TopJetLeptonOverlapRemoval(0.8));
    // note that the JetCleaner is only kept for the sake of example;
    // instead of constructing a jetcleaner explicitly,
    // the cleaning can also be achieved with less code via CommonModules with:
    // common->set_jet_id(PtEtaCut(30.0, 2.4));
    // before the 'common->init(ctx)' line.

    // 2. set up selections
    nJet_sel.reset(new NJetSelection(1)); // see common/include/NSelections.h
    boostedJet_sel.reset(new BoostedJetSelection(300));
    ptMuon_sel.reset(new PtMuonSelection(30));
    ptEle_sel.reset(new PtElecSelection(30));
    nMuon_sel.reset(new NMuonSelection(2,2));
    nEle_sel.reset(new NElectronSelection(2,2));
    diMuon_sel.reset(new DiMuonSelection());
    diEle_sel.reset(new DiElecSelection());
    Z_sel.reset(new ZSelection());
    phiAngular_sel.reset(new PhiAngularCut(2.7));

    //const std::string& trigger = ctx.get("trigger", "NULL");
    //if(trigger != "NULL") trigger_sel.reset(new TriggerSelection(trigger));

    // 3. Set up Hists classes:
    h_nJet_nocuts.reset(new BoostedHiggsToWWHists(ctx, "nJet_noCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "ele_noCuts"));
    h_muon_nocuts.reset(new MuonHists(ctx, "muon_noCuts"));
    h_diLepton_nocuts.reset(new DiLeptonHists(ctx, "diLepton_noCuts"));

    h_nJet_cleaned.reset(new BoostedHiggsToWWHists(ctx, "nJet_cleaned"));
    h_ele_cleaned.reset(new ElectronHists(ctx, "ele_cleaned"));
    h_muon_cleaned.reset(new MuonHists(ctx, "muon_cleaned"));
    h_diLepton_cleaned.reset(new DiLeptonHists(ctx, "diLepton_cleaned"));

    h_nJet.reset(new BoostedHiggsToWWHists(ctx, "nJet_sel"));
    h_boostedJet_sel.reset(new BoostedHiggsToWWHists(ctx, "boostedJet_sel"));

    h_ptmuon_sel.reset(new MuonHists(ctx, "ptMuon_sel"));
    h_ptele_sel.reset(new ElectronHists(ctx, "ptEle_sel"));
    h_nmuon_sel.reset(new MuonHists(ctx, "nmuon_sel"));
    h_nele_sel.reset(new ElectronHists(ctx, "nele_sel"));
    h_dimuon_sel.reset(new MuonHists(ctx, "diMuon_sel"));
    h_diele_sel.reset(new ElectronHists(ctx, "diEle_sel"));
    h_dilep_sel.reset(new DiLeptonHists(ctx, "diLepton_sel"));

    h_dimuon_Z_sel.reset(new MuonHists(ctx, "diMuon_Z_sel"));
    h_diele_Z_sel.reset(new ElectronHists(ctx, "diEle_Z_sel"));
    h_dilep_Z_sel.reset(new DiLeptonHists(ctx, "diLepton_Z_sel"));

    h_jet_angular_sel.reset(new BoostedHiggsToWWHists(ctx, "boostedJet_angular_sel"));
    h_dilep_angular_sel.reset(new DiLeptonHists(ctx, "diLepton_angular_sel"));

    unsorted_jets = ctx.declare_event_output< std::vector< TopJet > >("unsorted_jets");

  }


  bool FeasibilityStudyModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.

    cout << "Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    //std::cout << "EVENT.WEIGHT "<< event.weight << std::endl;

    h_nJet_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_muon_nocuts->fill(event);
    h_diLepton_nocuts->fill(event);

    // 1. run all modules other modules.
    common->process(event);
    jetcleaner->process(event);
    muoncleaner->process(event);
    elecleaner->process(event);
    overlap_removal->process(event);

    sort_by_pt(*event.topjets);
    sort_by_pt(*event.muons);
    sort_by_pt(*event.electrons);

    //const bool pass_trigger = trigger_sel->passes(event);
    //if(!pass_trigger) return false;

    // 2. test selections and fill histograms
    h_nJet_cleaned->fill(event);
    h_ele_cleaned->fill(event);
    h_muon_cleaned->fill(event);
    h_diLepton_cleaned->fill(event);

    // boosted Jet Selection
    bool nJet_selection = nJet_sel->passes(event);
    bool Jet_selection = nJet_selection;
    if(Jet_selection) h_nJet->fill(event);
    bool boostedJet_selection = boostedJet_sel->passes(event);
    Jet_selection = Jet_selection && boostedJet_selection;
    if(Jet_selection) h_boostedJet_sel->fill(event);


    // Muons Selection
    bool ptMuon_selection = ptMuon_sel->passes(event);
    bool Muon_selection = ptMuon_selection;
    if(Muon_selection) h_ptmuon_sel->fill(event);
    bool nMuon_selection = nMuon_sel->passes(event);
    Muon_selection = Muon_selection && nMuon_selection;
    if(Muon_selection) h_nmuon_sel->fill(event);
    bool diMuon_selection = diMuon_sel->passes(event);
    Muon_selection = Muon_selection && diMuon_selection;
    if(Muon_selection) h_dimuon_sel->fill(event);


    // Elecs Selection
    bool ptEle_selection = ptEle_sel->passes(event);
    bool Ele_selection = ptEle_selection;
    if(Ele_selection) h_ptele_sel->fill(event);
    bool nEle_selection = nEle_sel->passes(event);
    Ele_selection = Ele_selection && nEle_selection;
    if(Ele_selection) h_nele_sel->fill(event);
    bool diEle_selection = diEle_sel->passes(event);
    Ele_selection = Ele_selection && diEle_selection;
    if(Ele_selection) h_diele_sel->fill(event);
    bool event_selection = XOR(Muon_selection, Ele_selection);
    if(event_selection) h_dilep_sel->fill(event);

    if (!event_selection || (event.muons->size()==2 && event.electrons->size() == 2) ) return false;
    // Z Selection
    bool Z_selection = Z_sel->passes(event);
    Muon_selection = Muon_selection && Z_selection;
    Ele_selection = Ele_selection && Z_selection;
    if(Muon_selection) h_dimuon_Z_sel->fill(event);
    if(Ele_selection) h_diele_Z_sel->fill(event);
    event_selection = XOR(Muon_selection, Ele_selection);
    if(event_selection) h_dilep_Z_sel->fill(event);

    // Phi Angular diLep-Jet Selection
    event_selection = event_selection && Jet_selection;
    if (!event_selection) return false;
    bool phiAngular_selection = phiAngular_sel->passes(event);
    event_selection = event_selection && phiAngular_selection;

    std::vector< TopJet > my_unsorted_jets = *event.topjets;
    event.set(unsorted_jets, std::move(my_unsorted_jets));

    if(phiAngular_selection){
      sort_topjet_by_dilepdist(event);
      h_jet_angular_sel->fill(event);
      h_dilep_angular_sel->fill(event);
    }



    // 3. decide whether or not to keep the current event in the output:
    return event_selection;
  }

  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the FeasibilityStudyModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(FeasibilityStudyModule)

}
