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

//bool do_merge = true;
bool do_merge = false;

void DisplayPlots(std::vector<TString> name_vars, std::vector<TString> name_samples, std::vector<TString> name_subdirs, TString inputdir, TString outdir);

TCanvas* overlaphistos(std::map<TString,TH1*> histos, TH1* histo_data, std::vector<TString> names,TString name, std::map<TString,double> options);
TH1F* CutFlow(TString file_name, std::vector<TString> name_subdirs, TString name);
void MergeHistos(std::map<TString,TH1*>& histos_cuts, std::map<std::pair<TString,TString>,TH1*>& histos_vars, std::vector<TString> name_vars, std::vector<TString> names, TString name);
void LoadVariables(TString name_var, TString bkg, std::vector<TString> name_subdirs, std::map<TString,TH1*>& histos_cuts, std::map<std::pair<TString,TString>,TH1*>& histos_vars, TString inputdir);

void PlotHistos(std::map<TString,TH1*> histos_cuts, std::map<std::pair<TString,TString>,TH1*> histos_vars, std::vector<TString> names, std::vector<TString> name_vars, TString outdir);

void DisplayCutFlow(){
  std::vector<TString> channels = {"muonchannel", "electronchannel"};
  for(auto channel : channels) {
    std::cout << channel << '\n';

    TString outdir = "/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/plot/"+channel+"/";
    std::vector<TString> name_vars = {"nTopJet_JetDiLeptonPhiAngular/SDmass_1", "diLepton_JetDiLeptonPhiAngular/diMuon_DR12", "diLepton_JetDiLeptonPhiAngular/diElectron_DR12", "event_JetDiLeptonPhiAngular/MET"};
    std::vector<TString> name_subdirs = {"cleaned", "Veto", "NBoostedJetLepton", "NLeptonSel", "JetDiLeptonPhiAngular"};
    if (do_merge) {
      TString inputdir = "/nfs/dust/cms/user/amalara/sframe_all/FeasibilityStudy_MC/"+channel+"/";
      std::vector<TString> name_samples = {"HZ_HiggsToWWZToLL", "DY1JetsToLL", "DY2JetsToLL", "DY3JetsToLL", "DY4JetsToLL", "TTToHadronic", "TTbarSemiLeptonic", "TTTo2L2Nu", "WZ", "ZZ"};
      DisplayPlots( name_vars, name_samples, name_subdirs, inputdir, outdir);
    } else {
      TString inputdir = "/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/Preselection/"+channel+"/";
      std::vector<TString> name_samples = {"HZ", "DYJets", "TTbar", "WZ", "ZZ"};
      DisplayPlots( name_vars, name_samples, name_subdirs, inputdir, outdir);
    }
  }
}

void DisplayPlots(std::vector<TString> name_vars, std::vector<TString> name_samples, std::vector<TString> name_subdirs, TString inputdir, TString outdir) {

  std::map<TString,TH1*> histos_cuts;
  std::map<std::pair<TString,TString>,TH1*> histos_vars;

  for(auto name_var : name_vars) {
    for(auto sample : name_samples) {
      LoadVariables(name_var, sample, name_subdirs, histos_cuts,histos_vars, inputdir);
    }
  }

  // for (std::map<std::pair<TString,TString>,TH1*>::iterator it=histos_vars.begin(); it!=histos_vars.end(); ++it)
  // std::cout << it->first.first << " " << it->first.second << " => " << it->second << '\n';

  std::vector<TString> names;

  names = {"HZ_HiggsToWWZToLL"};
  if (!do_merge) names = {"HZ"};
  MergeHistos(histos_cuts, histos_vars, name_vars, names, "HZ");
  names = {"WZ"};
  MergeHistos(histos_cuts, histos_vars, name_vars, names, "WZ");
  names = {"ZZ"};
  MergeHistos(histos_cuts, histos_vars, name_vars, names, "ZZ");
  names = {"TTToHadronic","TTbarSemiLeptonic","TTTo2L2Nu"};
  if (!do_merge) names = {"TTbar"};
  MergeHistos(histos_cuts, histos_vars, name_vars, names, "TTbar");
  names.clear();
  names = {"DY1JetsToLL","DY2JetsToLL","DY3JetsToLL","DY4JetsToLL"};
  if (!do_merge) names = {"DYJets"};
  MergeHistos(histos_cuts, histos_vars, name_vars, names, "DYJets");
  names.clear();

  // for (auto entry : histos_cuts) {
  //   std::cout << entry.first << " " << entry.second << " " << entry.second->GetName() << '\n';
  // }
  //
  // for (std::map<std::pair<TString,TString>,TH1*>::iterator it=histos_vars.begin(); it!=histos_vars.end(); ++it)
  //  std::cout << it->first.first << " " << it->first.second << " => " << it->second << '\n';


  // names = {"HZ_HiggsToWWZToLL","DY1JetsToLL","TTToHadronic","TTbarSemiLeptonic","TTTo2L2Nu","WZ","ZZ"};
  names = {"HZ","DYJets","TTbar", "WZ", "ZZ" };
  // names = {"HZ","DYJets","TTbar" };

  PlotHistos(histos_cuts, histos_vars, names, name_vars, outdir);


  return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

TCanvas* overlaphistos(std::map<TString,TH1*> histos, TH1* histo_data, std::vector<TString> names, TString name, std::map<TString,double> options) {
  gStyle->SetOptStat(false);
  std::cout << "option " << options["drawoption"] << '\n';
  std::string drawoption, name_signal;
  if ((int)options["drawoption"]== 0) drawoption = "p";
  else if ((int)options["drawoption"]== 1) drawoption = "hist";

  if ((int)options["name"]== 1) name_signal = "HZ";
  else name_signal = "HZ_HiggsToWWZToLL";

  TCanvas* c_histo = new TCanvas(name,name, 800, 600);
  TLegend *leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);
  THStack* hs = new THStack(name, name);

  std::map<TString,int> colors= { {"HZ",kRed+1}, {"DYJets",kGreen-2}, {"TTbar",kBlue}, {"WZ",kOrange}, {"ZZ",kViolet-3}, {"DY1JetsToLL",kSpring-9}, {"DY1JetsToLL",kGreen-2},{"DY1JetsToLL",kGreen+1},{"DY1JetsToLL",kGreen+3},{"TTToHadronic",kAzure-4},{"TTbarSemiLeptonic",kCyan},{"TTTo2L2Nu",kAzure+6}};

  for (auto entry : histos) {
    if(std::find(names.begin(), names.end(), entry.first) == names.end()) continue;
    c_histo->cd();
    entry.second->Rebin((int)options["rebin"]);
    entry.second->SetLineWidth(4);
    entry.second->GetYaxis()->SetRangeUser(options["rangemin"], options["rangemax"]);
    tdrDraw(entry.second, drawoption, colors[entry.first],colors[entry.first], 1, colors[entry.first], 0, colors[entry.first]);
    hs->Add(entry.second);
    if ((int)options["drawoption"]== 1) {
      TString name = Form("%s : 10e-%d ", entry.second->GetName(), (int)round(TMath::Log10(entry.second->GetBinContent(1)/entry.second->GetBinContent(entry.second->GetNbinsX()))) );
      if (entry.second->GetName() == name_signal) name = Form("%s : %.f %% ", entry.second->GetName(), entry.second->GetBinContent(entry.second->GetNbinsX())*100./entry.second->GetBinContent(1));
      leg->AddEntry(entry.second, name ,"lep");
    } else leg->AddEntry(entry.second, entry.second->GetName(),"lep");
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



TH1F* CutFlow(TString file_name, std::vector<TString> name_subdirs, TString name) {
  std::vector<TH1F*> histos;
  TFile *file = new TFile(file_name);
  for (int i = 0; i < name_subdirs.size(); i++) {
    TString folder;
    if (name_subdirs[i].Contains("Jet") )folder = "nJet_";
    else if (name_subdirs[i].Contains("Muon") )folder = "muon_";
    else if (name_subdirs[i].Contains("Electron") ) folder = "ele_";
    else folder = "nJet_";
    histos.push_back((TH1F*)file->Get(folder+name_subdirs[i]+"/number"));
  }

  TH1F* histo_cut = new TH1F(name,name, name_subdirs.size(), 0, name_subdirs.size());
  for (int i = 0; i < name_subdirs.size(); i++) {
    histo_cut->SetBinContent(histo_cut->FindBin(i), histos[i]->GetSumOfWeights());
    histo_cut->GetXaxis()->SetBinLabel(i+1,name_subdirs[i]);
  }

  return histo_cut;
}

void MergeHistos(std::map<TString,TH1*>& histos_cuts, std::map<std::pair<TString,TString>,TH1*>& histos_vars, std::vector<TString> name_vars, std::vector<TString> bkgs, TString name) {
  bool first = true;
  TH1F* histo_cut;
  std::map<TString, TH1F*> histo_var;
  for (auto bkg : bkgs) {
    if (first) {
      first = false;
      histo_cut = (TH1F*)histos_cuts[bkg]->Clone(name);
      for (auto name_var : name_vars) histo_var[name_var] = (TH1F*) histos_vars[make_pair(name_var,bkg)]->Clone(name);
    }
    else {
      histo_cut->Add(histos_cuts[bkg]);
      for (auto name_var : name_vars) histo_var[name_var]->Add(histos_vars[make_pair(name_var,bkg)]);
    }
  }
  histos_cuts[name] = histo_cut;
  for (auto name_var : name_vars) histos_vars[make_pair(name_var,name)] = histo_var[name_var];

}

void LoadVariables(TString name_var, TString bkg, std::vector<TString> name_subdirs, std::map<TString,TH1*> & histos_cuts, std::map<std::pair<TString,TString>,TH1*> & histos_vars, TString inputdir) {
  TString prename = "uhh2.AnalysisModuleRunner.MC.";
  if (do_merge) prename += "MC_";
  TString file_name = inputdir+prename+bkg; file_name += ".root";
  TH1F* histo_cut = CutFlow(file_name, name_subdirs, bkg);
  histo_cut->SetName(bkg);
  histos_cuts[bkg] = histo_cut;

  TFile *file = new TFile(file_name);
  // std::cout << file_name << " " << file << '\n';
  TH1F* histo = (TH1F*)file  ->Get(name_var);
  // histos_vars[bkg] = histo;
  histos_vars[make_pair(name_var,bkg)] = histo;
}


void PlotHistos(std::map<TString,TH1*> histos_cuts, std::map<std::pair<TString,TString>,TH1*> histos_vars, std::vector<TString> names, std::vector<TString> name_vars, TString outdir) {
  std::map<TString,double> options = { {"rangemin",1e-10}, {"rangemax",1e5}, {"drawoption", 1 }, {"rebin",1}, {"name",1} } ;
  TCanvas* c_cutflow;
  c_cutflow = overlaphistos(histos_cuts, NULL, names, "CutFlow", options);
  c_cutflow->SetLogy();
  c_cutflow->SaveAs(outdir+"cutflow.png");
  c_cutflow->SaveAs(outdir+"cutflow.pdf");

  options["rangemin"] = 1e-10;
  options["rangemax"] = 100.;
  options["drawoption"] = 0;
  options["rebin"] = 4;

  for(auto name_var : name_vars) {
    TString name = name_var(name_var.First("/")+1,name_var.Length());
    std::map<TString,TH1*> histos_vars_plot;
    for(auto sample : names) histos_vars_plot[sample] = histos_vars[make_pair(name_var, sample)];
    TCanvas* c = overlaphistos(histos_vars_plot, NULL, names, name, options);
    c->SetLogy(1);
    c->SaveAs(outdir+name+".png");
    c->SaveAs(outdir+name+".pdf");
  }
}
