#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <TPaveText.h>
#include <THStack.h>
#include <TLatex.h>
#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/tdrstyle_all.C"

TCanvas* overlaphistos(std::map<TString,TH1*> histos, TH1* histo_data, std::vector<TString> names,TString name, std::map<TString,double> options);
TH1F* CutFlow(TString file_name, std::vector<TString> name_dir, TString name = "histo_cut");
void MergeHistos(std::map<TString,TH1*> & histos, std::vector<TString> names, TString name);

void LoadVariables(TString bkg, std::vector<TString> name_dir, std::map<TString,TH1*> & histos_cuts, std::map<TString,TH1*> & histos_SD1) {
  TString inputdir = "/nfs/dust/cms/user/amalara/sframe_all/FeasibilityStudy_MC/";
  TString file_name = inputdir+"uhh2.AnalysisModuleRunner.MC.MC_"; file_name += bkg; file_name += ".root";
  TH1F* histo_cut = CutFlow(file_name, name_dir, bkg);
  histo_cut->SetName(bkg);
  histos_cuts[bkg] = histo_cut;

  TFile *file = new TFile(file_name);
  TH1F* histo_tot = (TH1F*)file  ->Get("nJet_cleaned/numbers");
  TH1F* histo = (TH1F*)file  ->Get("nJet_JetElectron/SDmass_jet1");
  histos_SD1[bkg] = histo;
}

void DisplayCutFlow(TString outdir = "/nfs/dust/cms/user/amalara/sframe_all/FeasibilityStudy_MC/", bool quick = true) {

  #define LOADVARIABLE(bkg,var,color)\
  TString file_name_##bkg = "/nfs/dust/cms/user/amalara/sframe_all/FeasibilityStudy_MC/uhh2.AnalysisModuleRunner.MC.MC_"; file_name_##bkg += #bkg; file_name_##bkg += ".root";\
  TFile *file_##bkg = new TFile(file_name_##bkg);\
  TH1F* histo_##var##_##bkg##_tot = (TH1F*)file_##bkg  ->Get("nJet_cleaned/numbers");\
  TH1F* histo_##var##_##bkg = (TH1F*)file_##bkg  ->Get("nJet_JetElectron/SDmass_jet1");\
  tdrDraw(histo_##var##_##bkg, "hist p", kFullCircle, color, kSolid, color, 0);\
  text = #bkg; text += ": "; text += histo_##var##_##bkg->GetEntries();  text += " of "; text += histo_##var##_##bkg##_tot->GetEntries(); text += " eff: "; text += round(100.*histo_##var##_##bkg->GetEntries()/histo_##var##_##bkg##_tot->GetEntries()); text += "\%"; \
  leg->AddEntry(histo_##var##_##bkg, text, "p");\
  histos_SD1[#bkg] = histo_##var##_##bkg;\

  // #define LOADBKG(bkg,name_dir)\
  // TString file_name_##bkg = "/nfs/dust/cms/user/amalara/sframe_all/FeasibilityStudy_MC/uhh2.AnalysisModuleRunner.MC.MC_"; file_name_##bkg += #bkg; file_name_##bkg += ".root";\
  // TH1F* histo_cut_##bkg = CutFlow(file_name_##bkg, name_dir, #bkg);\
  // histo_cut_##bkg->SetName(#bkg);\
  // histos_cuts[#bkg] = histo_cut_##bkg;\

  std::map<TString,TH1*> histos_cuts, histos_SD1;
  std::vector<TString> name_dir = {"nocuts", "cleaned", "NBoostedJetSel", "NMuon", "NElectron", "PtMuon", "PtElectron", "DiMuon", "DiElectron", "JetMuon", "JetElectron", "PhiAngularSelMuon", "PhiAngularSelElectron", "eventSel","DiLepton"};

  LoadVariables("HZ_HiggsToWWZToLL", name_dir, histos_cuts,histos_SD1);
  LoadVariables("DY1JetsToLL", name_dir, histos_cuts,histos_SD1);
  LoadVariables("DY2JetsToLL", name_dir, histos_cuts,histos_SD1);
  LoadVariables("DY3JetsToLL", name_dir, histos_cuts,histos_SD1);
  LoadVariables("DY4JetsToLL", name_dir, histos_cuts,histos_SD1);
  LoadVariables("TTToHadronic", name_dir, histos_cuts,histos_SD1);
  LoadVariables("TTbarSemiLeptonic", name_dir, histos_cuts,histos_SD1);
  LoadVariables("TTTo2L2Nu", name_dir, histos_cuts,histos_SD1);
  LoadVariables("WZ", name_dir, histos_cuts,histos_SD1);
  LoadVariables("ZZ", name_dir, histos_cuts,histos_SD1);


  std::vector<TString> names = {"HZ_HiggsToWWZToLL"};
  MergeHistos(histos_cuts, names, "HWW");
  MergeHistos(histos_SD1, names, "HWW");
  names = {"TTToHadronic","TTbarSemiLeptonic","TTTo2L2Nu"};
  MergeHistos(histos_cuts, names, "TT");
  MergeHistos(histos_SD1, names, "TT");
  names.clear();
  names = {"DY1JetsToLL","DY2JetsToLL","DY3JetsToLL","DY4JetsToLL"};
  MergeHistos(histos_cuts, names, "DY");
  MergeHistos(histos_SD1, names, "DY");
  names.clear();
  // names = {"HZ_HiggsToWWZToLL","DY1JetsToLL","TTToHadronic","TTbarSemiLeptonic","TTTo2L2Nu","WZ","ZZ"};
  names = {"HWW","DY","TT", "WZ", "ZZ" };


  std::map<TString,double> options = { {"rangemin",1e-7}, {"rangemax",1e5}, {"drawoption", 1 }, {"rebin",1} } ;
  TCanvas* c_cutflow = overlaphistos(histos_cuts, NULL, names, "CutFlow", options);
  c_cutflow->SetLogy();
  // c_cutflow->SaveAs((outdir+"cutflow.png").c_str());
  // c_cutflow->SaveAs((outdir+"cutflow.pdf").c_str());


  options["rangemin"] = 1e-10;
  options["rangemax"] = 1.;
  options["drawoption"] = 0;
  options["rebin"] = 4;
  // TCanvas* c = tdrCanvas("SDmass", 0, 200, 1e-10, 1E2, "SDmass_jet1" , "Events");
  TCanvas* c = overlaphistos(histos_SD1, NULL, names, "SDmass_jet1", options);
  c->SetLogy(1);

  // TString text;
  // TLegend* leg = tdrLeg(0.3,0.65,0.9,0.9); leg->SetTextFont(42);  leg->SetTextSize(0.035);

  // LOADVARIABLE(HZ_HiggsToWWZToLL,SD1,kRed+1)
  // LOADVARIABLE(DY1JetsToLL,SD1,kGreen-2)
  // LOADVARIABLE(DY2JetsToLL,SD1,kGreen-2)
  // LOADVARIABLE(DY3JetsToLL,SD1,kGreen-2)
  // LOADVARIABLE(DY4JetsToLL,SD1,kGreen-2)
  // LOADVARIABLE(TTToHadronic,SD1,kMagenta)
  // LOADVARIABLE(TTbarSemiLeptonic,SD1,kViolet-3)
  // LOADVARIABLE(TTTo2L2Nu,SD1,kViolet)
  // LOADVARIABLE(WZ,SD1,kOrange)
  // LOADVARIABLE(ZZ,SD1,kBlue-4)


  return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

TCanvas* overlaphistos(std::map<TString,TH1*> histos, TH1* histo_data, std::vector<TString> names, TString name, std::map<TString,double> options) {
  gStyle->SetOptStat(false);
  std::string drawoption;
  if ((int)options["drawoption"]== 0) drawoption = "p";
  else if ((int)options["drawoption"]== 1) drawoption = "hist";
  TCanvas* c_histo = new TCanvas(name,name, 800, 600);
  TLegend *leg = tdrLeg(0.40,0.70,0.9,0.9, 0.035);
  THStack* hs = new THStack(name, name);

  std::map<TString,int> colors= { {"HWW",kRed+1}, {"DY",kGreen-2}, {"TT",kBlue}, {"WZ",kOrange}, {"ZZ",kViolet-3}, {"DY1JetsToLL",kSpring-9}, {"DY1JetsToLL",kGreen-2},{"DY1JetsToLL",kGreen+1},{"DY1JetsToLL",kGreen+3},{"TTToHadronic",kAzure-4},{"TTbarSemiLeptonic",kCyan},{"TTTo2L2Nu",kAzure+6}};

  for (auto entry : histos) {
    if(std::find(names.begin(), names.end(), entry.first) == names.end()) continue;
    c_histo->cd();
    entry.second->Rebin((int)options["rebin"]);
    entry.second->SetLineWidth(4);
    entry.second->GetYaxis()->SetRangeUser(options["rangemin"], options["rangemax"]);
    tdrDraw(entry.second, drawoption, colors[entry.first],colors[entry.first], 1, colors[entry.first], 0, colors[entry.first]);
    hs->Add(entry.second);
    leg->AddEntry(entry.second, entry.second->GetName(),"lep");
  }

  if (histo_data) {
    tdrDraw(histo_data, drawoption+" nostack", kBlack,kBlack, 1, kBlack, 1001, kBlack);
    histo_data->SetLineWidth(2);
    hs->Add(histo_data);
    leg->AddEntry(histo_data, histo_data->GetName(),"lep");
  }

  hs->Draw( "hist nostack");
  leg->SetLineColorAlpha(1, 0.7);
  leg->Draw("same");
  return c_histo;
}



TH1F* CutFlow(TString file_name, std::vector<TString> name_dir, TString name = "histo_cut") {
  std::vector<TH1F*> histos;
  TFile *file = new TFile(file_name);
  for (int i = 0; i < name_dir.size(); i++) {
    TString folder;
    if (name_dir[i].Contains("Jet") )folder = "nJet_";
    else if (name_dir[i].Contains("Muon") )folder = "muon_";
    else if (name_dir[i].Contains("Electron") ) folder = "ele_";
    else folder = "nJet_";
    histos.push_back((TH1F*)file->Get(folder+name_dir[i]+"/number"));
  }

  TH1F* histo_cut = new TH1F(name,name, name_dir.size(), 0, name_dir.size());
  for (int i = 0; i < name_dir.size(); i++) {
    histo_cut->SetBinContent(histo_cut->FindBin(i), histos[i]->GetSumOfWeights());
    histo_cut->GetXaxis()->SetBinLabel(i+1,name_dir[i]);
  }

  return histo_cut;
}

void MergeHistos(std::map<TString,TH1*> & histos, std::vector<TString> names, TString name) {
  bool first = true;
  TH1F* histo;
  for (auto entry : histos) {
    if(std::find(names.begin(), names.end(), entry.first) == names.end()) continue;
    if (first) {histo = (TH1F*)entry.second->Clone(name); first = false;}
    else histo->Add(entry.second);
  }
  histos[name] = histo;
}
