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
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/BoostedHiggsToWW/include/GenParticleHiggsHists.h"
#include "UHH2/BoostedHiggsToWW/include/GenJetsHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  /** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
  *
  * This is the central class which calls other AnalysisModules, Hists or Selection classes.
  * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
  */
  class MCTruthMatchingModule: public AnalysisModule {
  public:

    explicit MCTruthMatchingModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<JetCleaner> 	jetcleaner;
    std::unique_ptr<MuonCleaner>	muoncleaner;
    std::unique_ptr<ElectronCleaner>	elecleaner;
    std::unique_ptr<TopJetLeptonOverlapRemoval> overlap_removal;
    std::unique_ptr<GenTopJetLeptonOverlapRemoval> overlap_removal_genjet;

    // Selections
    std::unique_ptr<Selection> nJet_sel, boostedJet_sel;
    std::unique_ptr<Selection> ptMuon_sel, ptEle_sel, nMuon_sel, nEle_sel, diMuon_sel, diEle_sel, Z_sel;
    std::unique_ptr<Selection> phiAngular_sel;

    std::unique_ptr<uhh2::Selection> trigger_sel;

    std::unique_ptr<uhh2::AndSelection> muon_sel, ele_sel;

    // Histos
    std::unique_ptr<Hists> h_genpart_noCuts, h_genpart, h_genpart1, h_genpart2, h_genpart3, h_genpart4, h_genpart5, h_genpart6;
    std::unique_ptr<Hists> h_genjets_noCuts, h_genjets, h_genjets_ordered;
    std::unique_ptr<Hists> h_nJet_nocuts, h_ele_nocuts, h_muon_nocuts, h_diLepton_nocuts;
    std::unique_ptr<Hists> h_nJet_trigger, h_ele_trigger, h_muon_trigger, h_diLepton_trigger;
    std::unique_ptr<Hists> h_nJet_cleaned, h_ele_cleaned, h_muon_cleaned, h_diLepton_cleaned;
    std::unique_ptr<Hists> h_nJet, h_boostedJet_sel;
    std::unique_ptr<Hists> h_ptmuon_sel, h_ptele_sel, h_nmuon_sel, h_nele_sel, h_dimuon_sel, h_diele_sel, h_dilep_sel;
    std::unique_ptr<Hists> h_dimuon_Z_sel, h_diele_Z_sel, h_dilep_Z_sel;
    std::unique_ptr<Hists> h_jet_angular_sel, h_dilep_angular_sel;
    std::unique_ptr<Hists> h_jet_angular_match, h_dilep_angular_match;
    std::unique_ptr<Hists> h_jet_angular_matching_H;

    Event::Handle< std::vector< TopJet > > unsorted_jets, matched_jets;
  };


  MCTruthMatchingModule::MCTruthMatchingModule(Context & ctx){
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
    common->switch_metcorrection();
    const bool isMC = (ctx.get("dataset_type") == "MC");
    if (isMC) {
      common->disable_metfilters();
    }
    common->init(ctx);

    jetcleaner.reset(new JetCleaner(ctx, 30.0, 2.4));
    muoncleaner.reset(new MuonCleaner (AndId<Muon>(MuonID(Muon::CutBasedIdMedium), PtEtaCut(1., 4))));
    elecleaner.reset(new ElectronCleaner (AndId<Electron>(ElectronID_Spring16_medium, PtEtaCut(1., 4))));
    overlap_removal.reset(new TopJetLeptonOverlapRemoval(0.8));
    overlap_removal_genjet.reset(new GenTopJetLeptonOverlapRemoval(0.8));

    // note that the JetCleaner is only kept for the sake of example;
    // instead of constructing a jetcleaner explicitly,
    // the cleaning can also be achieved with less code via CommonModules with:
    // common->set_jet_id(PtEtaCut(30.0, 2.4));
    // before the 'common->init(ctx)' line.

    // 2. set up selections
    nMuon_sel.reset(new NMuonSelection(2));
    nEle_sel.reset(new NElectronSelection(2));
    ptMuon_sel.reset(new PtMuonSelection(30, 0.2));
    ptEle_sel.reset(new PtElecSelection(30, 0.2));
    diMuon_sel.reset(new DiMuonSelection());
    diEle_sel.reset(new DiElecSelection());
    Z_sel.reset(new ZSelection());
    nJet_sel.reset(new NJetSelection(1)); // see common/include/NSelections.h
    boostedJet_sel.reset(new BoostedJetSelection(300));
    phiAngular_sel.reset(new PhiAngularCut(2.7));

    const std::string& trigger = ctx.get("trigger", "NULL");
    if(trigger != "NULL") trigger_sel.reset(new TriggerSelection(trigger));
    else trigger_sel.reset(new AndSelection(ctx));

    muon_sel.reset(new uhh2::AndSelection(ctx, "muon_filters"));
    muon_sel->add<NMuonSelection>("2 muons", 2);
    muon_sel->add<PtMuonSelection>("pt_muon", 30);
    muon_sel->add<DiMuonSelection>("diMuon");
    muon_sel->add<ZSelection>("Z_selection");
    muon_sel->add<NJetSelection>("1 Jet", 1);
    muon_sel->add<BoostedJetSelection>("boosted Jet", 300);
    muon_sel->add<PhiAngularCut>("PhiAngularCut", 2.7);

    ele_sel.reset(new uhh2::AndSelection(ctx, "ele_filters"));
    ele_sel->add<NElectronSelection>("2 electrons", 2);
    ele_sel->add<PtElecSelection>("pt_electron", 30);
    ele_sel->add<DiElecSelection>("diEle");
    ele_sel->add<ZSelection>("Z_selection");
    ele_sel->add<NJetSelection>("1 Jet", 1);
    ele_sel->add<BoostedJetSelection>("boosted Jet", 300);
    ele_sel->add<PhiAngularCut>("PhiAngularCut", 2.7);

    // 3. Set up Hists classes:
    h_genpart_noCuts.reset(new GenParticleHiggsHists(ctx, "gen_particles_noCuts"));
    h_genpart.reset(new GenParticleHiggsHists(ctx, "gen_particles"));
    h_genpart1.reset(new GenParticleHiggsHists(ctx, "gen_particles1"));
    h_genpart2.reset(new GenParticleHiggsHists(ctx, "gen_particles2"));
    h_genpart3.reset(new GenParticleHiggsHists(ctx, "gen_particles3"));
    h_genpart4.reset(new GenParticleHiggsHists(ctx, "gen_particles4"));
    h_genpart5.reset(new GenParticleHiggsHists(ctx, "gen_particles5"));
    h_genpart6.reset(new GenParticleHiggsHists(ctx, "gen_particles6"));
    h_genjets_noCuts.reset(new GenJetsHists(ctx, "gen_jets_noCuts"));
    h_genjets.reset(new GenJetsHists(ctx, "gen_jets"));
    h_genjets_ordered.reset(new GenJetsHists(ctx, "gen_jets_ordered"));

    h_nJet_nocuts.reset(new BoostedHiggsToWWHists(ctx, "nJet_noCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "ele_noCuts"));
    h_muon_nocuts.reset(new MuonHists(ctx, "muon_noCuts"));
    h_diLepton_nocuts.reset(new DiLeptonHists(ctx, "diLepton_noCuts"));

    h_nJet_trigger.reset(new BoostedHiggsToWWHists(ctx, "nJet_trigger"));
    h_ele_trigger.reset(new ElectronHists(ctx, "ele_trigger"));
    h_muon_trigger.reset(new MuonHists(ctx, "muon_trigger"));
    h_diLepton_trigger.reset(new DiLeptonHists(ctx, "diLepton_trigger"));

    h_nJet_cleaned.reset(new BoostedHiggsToWWHists(ctx, "nJet_cleaned"));
    h_ele_cleaned.reset(new ElectronHists(ctx, "ele_cleaned"));
    h_muon_cleaned.reset(new MuonHists(ctx, "muon_cleaned"));
    h_diLepton_cleaned.reset(new DiLeptonHists(ctx, "diLepton_cleaned"));

    h_nmuon_sel.reset(new MuonHists(ctx, "nmuon_sel"));
    h_nele_sel.reset(new ElectronHists(ctx, "nele_sel"));
    h_ptmuon_sel.reset(new MuonHists(ctx, "ptMuon_sel"));
    h_ptele_sel.reset(new ElectronHists(ctx, "ptEle_sel"));
    h_dimuon_sel.reset(new MuonHists(ctx, "diMuon_sel"));
    h_diele_sel.reset(new ElectronHists(ctx, "diEle_sel"));
    h_dilep_sel.reset(new DiLeptonHists(ctx, "diLepton_sel"));

    h_dimuon_Z_sel.reset(new MuonHists(ctx, "diMuon_Z_sel"));
    h_diele_Z_sel.reset(new ElectronHists(ctx, "diEle_Z_sel"));
    h_dilep_Z_sel.reset(new DiLeptonHists(ctx, "diLepton_Z_sel"));

    h_nJet.reset(new BoostedHiggsToWWHists(ctx, "nJet_sel"));
    h_boostedJet_sel.reset(new BoostedHiggsToWWHists(ctx, "boostedJet_sel"));

    h_jet_angular_sel.reset(new BoostedHiggsToWWHists(ctx, "boostedJet_angular_sel"));
    h_dilep_angular_sel.reset(new DiLeptonHists(ctx, "diLepton_angular_sel"));

    h_jet_angular_match.reset(new BoostedHiggsToWWHists(ctx, "boostedJet_angular_match"));
    h_dilep_angular_match.reset(new DiLeptonHists(ctx, "diLepton_angular_match"));

    h_jet_angular_matching_H.reset(new BoostedHiggsToWWHists(ctx, "boostedJet_angular_matching_H"));

    unsorted_jets = ctx.declare_event_output< std::vector< TopJet > >("unsorted_jets");
    matched_jets = ctx.declare_event_output< std::vector< TopJet > >("matched_jets");

  }


  bool MCTruthMatchingModule::process(Event & event) {
    // cout << "Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;


    // for(const auto & muon : *event.muons) std::cout << muon.get_bool(Muon::global) << '\n';
    // std::cout << "PF" << '\n';
    // for(const auto & muon : *event.genparticles){if (fbs(muon.pdgId()) == 13) std::cout << muon.get_bool(Muon::global) << '\n'; }

    muon_sel->passes(event);
    ele_sel->passes(event);
    h_nJet_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_muon_nocuts->fill(event);
    h_diLepton_nocuts->fill(event);
    if (event.gentopjets->size() >0 && event.topjets->size() >0) {
      h_genpart_noCuts->fill(event);
      h_genjets_noCuts->fill(event);
    }

    GenParticle H, muon1, muon2;
    for(const auto & gp : *event.genparticles){
      if (gp.pdgId() == 25)  H = gp;
      if (gp.pdgId() == 13)  muon1 = gp;
      if (gp.pdgId() == -13)  muon2 = gp;
    }

    // 1. run all modules other modules.
    common->process(event);
    jetcleaner->process(event);
    muoncleaner->process(event);
    elecleaner->process(event);
    overlap_removal->process(event);
    overlap_removal_genjet->process(event);

    sort_by_pt<TopJet>(*event.topjets);
    sort_by_pt<Muon>(*event.muons);
    sort_by_pt<Electron>(*event.electrons);
    sort_by_pt<GenTopJet>(*event.gentopjets);

    const bool pass_trigger = trigger_sel->passes(event);
    if(!pass_trigger) return false;

    h_nJet_trigger->fill(event);
    h_ele_trigger->fill(event);
    h_muon_trigger->fill(event);
    h_diLepton_trigger->fill(event);


    // 2. test selections and fill histograms
    h_nJet_cleaned->fill(event);
    h_ele_cleaned->fill(event);
    h_muon_cleaned->fill(event);
    h_diLepton_cleaned->fill(event);
    if (event.gentopjets->size() >0) h_genjets->fill(event);

    sort_gentopjet_H(*event.gentopjets, H);
    if (event.gentopjets->size() >0) h_genjets_ordered->fill(event);

    // Muons Selection
    bool nMuon_selection = nMuon_sel->passes(event);
    bool Muon_selection = nMuon_selection;
    if(Muon_selection) h_nmuon_sel->fill(event);
    bool ptMuon_selection = ptMuon_sel->passes(event);
    Muon_selection = Muon_selection && ptMuon_selection;
    if(Muon_selection) h_ptmuon_sel->fill(event);
    bool diMuon_selection = diMuon_sel->passes(event);
    Muon_selection = Muon_selection && diMuon_selection;
    if(Muon_selection) h_dimuon_sel->fill(event);


    // Elecs Selection
    bool nEle_selection = nEle_sel->passes(event);
    bool Ele_selection = nEle_selection;
    if(Ele_selection) h_ptele_sel->fill(event);
    bool ptEle_selection = ptEle_sel->passes(event);
    Ele_selection = Ele_selection && ptEle_selection;
    if(Ele_selection) h_nele_sel->fill(event);
    bool diEle_selection = diEle_sel->passes(event);
    Ele_selection = Ele_selection && diEle_selection;
    if(Ele_selection) h_diele_sel->fill(event);
    bool event_selection = XOR(Muon_selection, Ele_selection);
    if(event_selection) h_dilep_sel->fill(event);

    if (!event_selection) return false;
    // Z Selection
    bool Z_selection = Z_sel->passes(event);
    Muon_selection = Muon_selection && Z_selection;
    Ele_selection = Ele_selection && Z_selection;
    if(Muon_selection) h_dimuon_Z_sel->fill(event);
    if(Ele_selection) h_diele_Z_sel->fill(event);
    event_selection = XOR(Muon_selection, Ele_selection);
    if(event_selection) h_dilep_Z_sel->fill(event);

    // boosted Jet Selection
    bool nJet_selection = nJet_sel->passes(event);
    bool Jet_selection = nJet_selection;
    if(Jet_selection) h_nJet->fill(event);
    bool boostedJet_selection = boostedJet_sel->passes(event);
    Jet_selection = Jet_selection && boostedJet_selection;
    if(Jet_selection) h_boostedJet_sel->fill(event);

    // Phi Angular diLep-Jet Selection
    event_selection = event_selection && Jet_selection;
    if (!event_selection) return false;
    bool phiAngular_selection = phiAngular_sel->passes(event);
    event_selection = event_selection && phiAngular_selection;

    if (!event_selection) return false;

    std::vector< TopJet > my_unsorted_jets = *event.topjets;
    event.set(unsorted_jets, std::move(my_unsorted_jets));
    sort_topjet_by_dilepdist(event);
    h_jet_angular_sel->fill(event);
    h_dilep_angular_sel->fill(event);
    h_genpart->fill(event);

    // GenParticle H;
    // for(const auto & gp : *event.genparticles){
    //   // std::cout << "GEN: " << gp.pdgId() << " " << gp.status() << " " << gp.index() << " " << gp.mother1() << " " << gp.mother2() << " " << gp.daughter1() << " " << gp.daughter2() << '\n';
    //   // abs(mother().daughter(0).pdgId()) != 24
    //   if (gp.pdgId() == 25)  H = gp;
    // }

    bool match_selection = uhh2::deltaR((event.topjets)->at(0), H) < 0.8;
    if (match_selection) {
      h_jet_angular_match->fill(event);
      h_dilep_angular_match->fill(event);
    }

    if (event.gentopjets->size() >0) {
      bool match_hadron_selection = uhh2::deltaR((event.topjets)->at(0), (event.gentopjets)->at(0)) < 0.8;
      if (match_hadron_selection) h_jet_angular_matching_H->fill(event);
    }
    std::vector< TopJet > my_matched_jets = *event.topjets;
    event.set(matched_jets, std::move(my_matched_jets));
    sort_topjet_H(*event.topjets, H);


    // 3. decide whether or not to keep the current event in the output:
    return event_selection;
  }

  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the MCTruthMatchingModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(MCTruthMatchingModule)

}
