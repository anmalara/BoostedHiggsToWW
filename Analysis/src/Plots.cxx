#include <vector>

#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/Analysis/src/Plotter.cxx"

void Plots() {
  Plotter* plotter = new Plotter("/nfs/dust/cms/user/amalara/sframe_all/GenericCleaning_MC/muonchannel/");
  plotter->SetNameFiles({"uhh2.AnalysisModuleRunner.MC.MC_HZ_HiggsToWWZToLL"});
  plotter->SetNameFolders({"nTopJet", "nTopJet_ZMatch", "nTopJet_WWfullHad", "nTopJet_WWsemiLep", "nTopJet_WWfullLep"});
  // plotter->SetHistoTags({"IDs", "metfilters", "jetlep", "JEC", "JER", "MET", "jetIDnoboost", "jetID", "cleaned", "Veto", "NBoostedJet", "NLeptonSel", "JetDiLeptonPhiAngular"});

  plotter->SetHistoTags({"IDs", "jetlep", "JEC", "MET", "jetID", "cleaned", "Veto", "NBoostedJet"});

  plotter->LoadMap();
  plotter->PrintMap();

  plotter->PlotHistosEvolution(1, 0);


  // plotter->SetHistoTags({"IDs", "jetlep", "cleaned", "JetDiLeptonPhiAngular"});
  plotter->SetNameFolders({"nTopJet_ZMatch", "nTopJet_WWfullHad", "nTopJet_WWsemiLep", "nTopJet_WWfullLep"});
  plotter->PlotHistosComposition("uhh2.AnalysisModuleRunner.MC.MC_HZ_HiggsToWWZToLL", "nTopJet", 1);


  // plotter->PlotMaps(1, 0);

  // plotter->PlotSingleHistoComposition("uhh2.AnalysisModuleRunner.MC.MC_HZ_HiggsToWWZToLL", "IDs", "nTopJet", "mass_jet", 1);

  // plotter->SetNameFolders({"nTopJet_WWfullLep", "nTopJet_WWsemiLep", "nTopJet_WWfullHad" });
  // plotter->PlotSingleHistoComposition("uhh2.AnalysisModuleRunner.MC.MC_HZ_HiggsToWWZToLL", "IDs", "nTopJet_HMatch", "mass_jet", 1);


  // plotter->PlotHistosComposition("uhh2.AnalysisModuleRunner.MC.MC_HZ_HiggsToWWZToLL", "nTopJet_HMatch", 1);
  // plotter->SetNameFolders({"nTopJet_ZMatch", "nTopJet_WWfullHad", "nTopJet_WWsemiLep", "nTopJet_WWfullLep"});
  // plotter->PlotHistosComposition("uhh2.AnalysisModuleRunner.MC.MC_HZ_HiggsToWWZToLL", "nTopJet", 1);

  // plotter->PlotAllHistos(1, 0);

  // plotter->PlotMaps("nTopJet_ZMatch_", "mass_jet", 1, 0);
  // plotter->PlotMaps("nTopJet_HMatch_", "mass_jet", 1, 0);
  // plotter->PlotMaps("nTopJet_WWfullHad_", "mass_jet", 1, 0);
  // plotter->PlotMaps("nTopJet_WWfullLep_", "mass_jet", 1, 0);
  // plotter->PlotMaps("nTopJet_WWsemiLep_", "mass_jet", 1, 0);


}
//
//
//
//
// for (size_t i = 1; i < 95; i++) {
//   TString fname = "/pnfs/desy.de/cms/tier2/store/user/pgunnell/RunII_94X_v1/HiggsZProduction/crab_HZ_HiggsToWWZToLL_Pythia8_TuneCP5/180214_092805/0000/Ntuple_";
//   fname += i; fname += ".root";
//   TFile *file_ = TFile::Open(fname);
//   TTree* tt = (TTree*)file_->Get("AnalysisTree");
//   std::cout << "FILE " << fname << '\n';
//   tt->Scan("run:event:luminosityBlock:updatedPatJetsPatJetsAK8PFPUPPI.m_pt:updatedPatJetsPatJetsAK8PFPUPPI.m_eta", "(event>145000)*(event<145100)");
//   file_->Close();
// }
//
//
//
//
//
//
// TString fname = "/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/Preselection/muonchannel/uhh2.AnalysisModuleRunner.MC_HZ.root";
// TFile *file_ = TFile::Open(fname);
// TTree* tt = (TTree*)file_->Get("AnalysisTree");
// std::cout << "FILE " << fname << '\n';
// tt->Scan("run:event:luminosityBlock:packedPatJetsAk8CHSJets_SoftDropCHS.m_pt:packedPatJetsAk8CHSJets_SoftDropCHS.m_eta:packedPatJetsAk8CHSJets_SoftDropCHS.m_phi", "(event>145000)*(event<145100)");
// file_->Close();
//
//
// TString fname1 = "/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/Preselection/electronchannel/uhh2.AnalysisModuleRunner.MC_HZ.root";
// TFile *file1 = TFile::Open(fname1);
// TTree* tt1 = (TTree*)file1->Get("AnalysisTree");
// std::cout << "FILE " << fname1 << '\n';
// tt1->Scan("run:event:luminosityBlock:packedPatJetsAk8CHSJets_SoftDropCHS.m_pt:packedPatJetsAk8CHSJets_SoftDropCHS.m_eta:packedPatJetsAk8CHSJets_SoftDropCHS.m_phi", "(event>145000)*(event<145100)");
// file1->Close();
//
//
//
// for (size_t i = 1; i < 100; i++) {
//   TString fname = "/pnfs/desy.de/cms/tier2/store/user/pgunnell/HiggsZProduction/MC-CandidatesP8-Higgs-GenInfos-v5-ext1/180518_072651/0000/flatTreeFileHiggs-nano_";
//   fname += i; fname += ".root";
//   TFile *file_ = TFile::Open(fname);
//   TTree* tt = (TTree*)file_->Get("boostedAK8/events");
//   std::cout << "FILE " << fname << '\n';
//   tt->Scan("runNo.run:evtNo.evt:lumi.lumi:jetPt:jetEta", "(evtNo.evt>145000)*(evtNo.evt<145100)");
//   file_->Close();
// }
//
//
//
//
//
// for (size_t i = 1; i < 1000; i++) {
//   TString fname = "/pnfs/desy.de/cms/tier2/store/user/pgunnell/HiggsZProduction/MC-CandidatesP8-Higgs-GenInfos-v4/180413_072857/0000/flatTreeFileHiggs-nano_";
//   fname += i; fname += ".root";
//   TFile *file_ = TFile::Open(fname);
//   TTree* tt = (TTree*)file_->Get("boostedAK8/events");
//   std::cout << "FILE " << fname << '\n';
//   tt->Scan("runNo.run:evtNo.evt:lumi.lumi:jetPt:jetEta", "(evtNo.evt>145000)*(evtNo.evt<145100)");
//   file_->Close();
// }





//
//
//
//
// t->Scan("runNo.run:evtNo.evt:lumi.lumi:jetPt:jetEta:jetPhi", "evtNo.evt>17000")
//
//
//
// TH1F* h0 = (TH1F*)_file0->Get("pt_ave_wNext_trg80")
// TH1F* h1 = (TH1F*)_file1->Get("pt_ave_wNext_trg140")
// TH1F* h2 = (TH1F*)_file2->Get("pt_ave_wNext_trg200")
// TH1F* h3 = (TH1F*)_file3->Get("pt_ave_wNext_trg260")
// TH1F* h4 = (TH1F*)_file4->Get("pt_ave_wNext_trg320")
// TH1F* h5 = (TH1F*)_file5->Get("pt_ave_wNext_trg400")
// TH1F* h6 = (TH1F*)_file6->Get("pt_ave_wNext_trg450")
// TH1F* h7 = (TH1F*)_file7->Get("pt_ave_wNext_trg40")
// TH1F* h8 = (TH1F*)_file8->Get("pt_ave_wNext_trg60")
//
//
// Double_t SmoothFit(Double_t *v, Double_t *par){
//   Double_t fitval  = 0.;
//   if(par[2] != 0.){
//     fitval = 0.5 * par[2] * (1. + TMath::Erf((v[0]-par[0]) / (TMath::Power(2, 0.5) * par[1] ) ) );
//   }
//
//   return fitval;
// }
//
// TF1* f0 = new TF1("myfunc0",SmoothFit,100,600,3)
// TF1* f1 = new TF1("myfunc1",SmoothFit,100,600,3)
// TF1* f2 = new TF1("myfunc2",SmoothFit,100,600,3)
// TF1* f3 = new TF1("myfunc3",SmoothFit,250,600,3)
// TF1* f4 = new TF1("myfunc4",SmoothFit,250,600,3)
// TF1* f5 = new TF1("myfunc5",SmoothFit,100,600,3)
// TF1* f6 = new TF1("myfunc6",SmoothFit,10,1200,3)
// TF1* f7 = new TF1("myfunc7",SmoothFit,100,600,3)
// TF1* f8 = new TF1("myfunc8",SmoothFit,100,600,3)
// //
// // new TCanvas
// // h0->Draw()
// // h1->Draw()
// // h2->Draw()
// // h3->Draw()
// // h4->Draw()
// // h5->Draw()
// // h6->Draw()
// // h7->Draw()
// // h8->Draw()
//
// new TCanvas
// h0->Fit(f0,"R")
//
// new TCanvas
// h1->Fit(f1,"R")
//
// new TCanvas
// h2->Fit(f2,"R")
//
// new TCanvas
// h3->Fit(f3,"R")
//
// new TCanvas
// h4->Fit(f4,"R")
//
// new TCanvas
// h5->Fit(f5,"R")
//
// new TCanvas
// h6->Fit(f6,"R")
//
// new TCanvas
// h7->Fit(f7,"R")
//
// new TCanvas
// h8->Fit(f8,"R")
