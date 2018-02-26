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
#include "tdrstyle.cxx"
// #include "../src/BoostedHiggsToWWUtils.cxx"

TCanvas* overlap2histos(TH1* histo1, TH1* histo2, bool full=false, bool norm=false);
TCanvas* overlaphistos(std::vector<TH1*> histos, TH1* histo_data, bool data = false, bool norm = false, bool full = false);
TH1F* CutFlow(string file_name, std::vector<string> name_dir, string histo_name = "histo_cut", bool norm = false, float cross_sec = 1);
bool OrderHistos(TH1* h1, TH1* h2);

void DisplayCutFlow(string inputdir = "../file/outputfile/", string outdir = "../file/outputfile/", bool quick = true) {

  // string file_name_signal = inputdir+"uhh2.AnalysisModuleRunner.MC.signal.root";
  string file_name_signal = "../file/outputfile/MCTruthMatching/uhh2.AnalysisModuleRunner.MC.signal_new.root";
  //string file_name_signal = "../file/outputfile/feasibilitystudy_signalnew/uhh2.AnalysisModuleRunner.MC.signal_new.root";
  string file_name_bkg_DY = inputdir+"uhh2.AnalysisModuleRunner.MC.bkg_DYJetsToLL_Pt-250To400.root";
  string file_name_bkg_TT = inputdir+"uhh2.AnalysisModuleRunner.MC.bkg_TT.root";
  string file_name_bkg_ZZ = inputdir+"uhh2.AnalysisModuleRunner.MC.bkg_ZZ.root";
  string file_name_bkg_ZW = inputdir+"uhh2.AnalysisModuleRunner.MC.bkg_ZW.root";

  TFile *file_signal = new TFile((file_name_signal).c_str());
  TFile *file_bkg_DY = new TFile((file_name_bkg_DY).c_str());
  TFile *file_bkg_TT = new TFile((file_name_bkg_TT).c_str());
  TFile *file_bkg_ZZ = new TFile((file_name_bkg_ZZ).c_str());
  TFile *file_bkg_ZW = new TFile((file_name_bkg_ZW).c_str());

  std::vector<string> name_dir;
  std::vector<string> name_histo;


  TList *dirlist = file_bkg_DY->GetListOfKeys();
  TIterator *iter = dirlist->MakeIterator();
  TObject *key = iter->Next();

  while (key) {
    string name_avoid = key->GetName();
    if (name_avoid == "SFrame") {
      break;
    }
    if (key->IsFolder()) {
      TDirectory *dir = file_bkg_DY->GetDirectory(key->GetName());
      string name = dir->GetName();
      name_dir.push_back(name);
      gSystem->Exec(("mkdir -p "+outdir+name).c_str());
      TList *dirlist_1 = dir->GetListOfKeys();
      TIterator *iter_1 = dirlist_1->MakeIterator();
      TObject *key_1 = iter_1->Next();
      while (key_1) {
        TObject* h= dir->Get( key_1->GetName() );
        if (h->InheritsFrom( "TH1" )){
          string name = key->GetName();
          string name_1 = key_1->GetName();
          name_histo.push_back(name+"/"+name_1);
        }
        key_1 = iter_1->Next();
      }
    }
    key = iter->Next();
  }

  if(quick) {
    name_histo = {"boostedJet_angular_sel/SDmass_jet1"};
    name_dir = {"nJet_noCuts", "boostedJet_sel", "ptMuon_sel", "ptEle_sel", "diMuon_sel", "diEle_sel", "diLepton_Z_sel", "boostedJet_angular_sel" };
  }

  bool cut_norm = true;
  TH1F* histo_cut_signal = CutFlow(file_name_signal, name_dir, "signal", cut_norm, 0.00022); //0.752
  TH1F* histo_cut_bkg_DY = CutFlow(file_name_bkg_DY, name_dir, "DY", cut_norm, 2.991);
  TH1F* histo_cut_bkg_TT = CutFlow(file_name_bkg_TT, name_dir, "TT", cut_norm, 832);
  TH1F* histo_cut_bkg_ZZ = CutFlow(file_name_bkg_ZZ, name_dir, "ZZ", cut_norm, 10.32);
  TH1F* histo_cut_bkg_ZW = CutFlow(file_name_bkg_ZW, name_dir, "ZW", cut_norm, 22.82);

  std::vector<TH1*> histo_cuts;
  histo_cuts.push_back(histo_cut_signal);
  histo_cuts.push_back(histo_cut_bkg_TT);
  histo_cuts.push_back(histo_cut_bkg_DY);
  histo_cuts.push_back(histo_cut_bkg_ZZ);
  histo_cuts.push_back(histo_cut_bkg_ZW);


  TCanvas* c_cutflow = overlaphistos(histo_cuts, histo_cut_signal);
  c_cutflow->SaveAs((outdir+"cutflow.png").c_str());



  std::vector<TCanvas*> canvases;
  TH1F *histo1, *histo2, *histotemp;
  TCanvas *c;
  /*
  for (int i = 0; i < name_histo.size(); i++) {
  histotemp=(TH1F*)file_signal->Get((name_histo[i]).c_str());
  string name = histotemp->GetName();
  histo1= (TH1F*)histotemp->Clone((name+" Signal").c_str() );
  histotemp=(TH1F*)file_bkg_DY->Get((name_histo[i]).c_str());
  histo2= (TH1F*)histotemp->Clone((name+" Bkg").c_str() );
  c = overlap2histos(histo1, histo2, false, true);
  name = outdir+name_histo[i]+".png";
  c->SaveAs((name).c_str());
}
*/
for (int i = 0; i < name_histo.size(); i++) {
  TH1F* histo_signal = (TH1F*)file_signal->Get((name_histo[i]).c_str());
  TH1F* histo_bkg_DY = (TH1F*)file_bkg_DY->Get((name_histo[i]).c_str());
  TH1F* histo_bkg_TT = (TH1F*)file_bkg_TT->Get((name_histo[i]).c_str());
  TH1F* histo_bkg_ZZ = (TH1F*)file_bkg_ZZ->Get((name_histo[i]).c_str());
  TH1F* histo_bkg_ZW = (TH1F*)file_bkg_ZW->Get((name_histo[i]).c_str());
  histo_signal->SetName(("Signal "+to_string(histo_signal->Integral())).c_str());
  histo_signal->GetSumw2();
  histo_signal->Scale(1./histo_signal->Integral());
  histo_bkg_DY->SetName(("DY "+to_string(histo_bkg_DY->Integral())).c_str());
  histo_bkg_DY->GetSumw2();
  histo_bkg_DY->Scale(1./histo_bkg_DY->Integral());
  histo_bkg_TT->SetName(("TT "+to_string(histo_bkg_TT->Integral())).c_str());
  histo_bkg_TT->GetSumw2();
  histo_bkg_TT->Scale(1./histo_bkg_TT->Integral());
  histo_bkg_ZZ->SetName(("ZZ "+to_string(histo_bkg_ZZ->Integral())).c_str());
  histo_bkg_ZZ->GetSumw2();
  histo_bkg_ZZ->Scale(1./histo_bkg_ZZ->Integral());
  histo_bkg_ZW->SetName(("ZW "+to_string(histo_bkg_ZW->Integral())).c_str());
  histo_bkg_ZW->GetSumw2();
  histo_bkg_ZW->Scale(1./histo_bkg_ZW->Integral());
  std::vector<TH1*> histos;
  histos.push_back(histo_signal);
  histos.push_back(histo_bkg_TT);
  histos.push_back(histo_bkg_DY);
  histos.push_back(histo_bkg_ZZ);
  histos.push_back(histo_bkg_ZW);
  TCanvas* c = overlaphistos(histos, histo_signal);
  c->SetLogy();
  c->SaveAs((outdir+name_histo[i]+".png").c_str());
}



}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

TCanvas* overlaphistos(std::vector<TH1*> histos, TH1* histo_data, bool data = false, bool norm = false, bool full = false){
  setTDRStyle();
  gStyle->SetOptStat(false);

  TCanvas* c_histo = new TCanvas((histos[0]->GetName()),"histo overlap", 800, 600);
  std::vector<int> colors = {kAzure-8, kRed+1, kOrange -4, kCyan, kSpring-6, kRed, kOrange, kGreen+1, kMagenta-10};

  if (norm) {
    for (int i = 0; i < histos.size(); i++) {
      histos[i]->GetSumw2();
      histos[i]->Scale(1./histos[i]->Integral());
    }
    if (data){
      histo_data->GetSumw2();
      histo_data->Scale(1./histo_data->Integral());
    }
  }

  for (int i = 0; i < histos.size(); i++) {
    histos[i]->SetFillColor(colors[i]);
    histos[i]->SetFillStyle(3004);
    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetLineStyle(1);
    histos[i]->SetLineWidth(2);
  }

  histos[0]->SetFillColorAlpha(kAzure-8, 0.85);
  histos[0]->SetFillStyle(1001);
  histos[0]->SetLineColor(kBlue);
  histos[0]->SetLineStyle(1);
  histos[0]->SetLineWidth(2);

  if (data) {
    histo_data->SetFillColor(kBlack);
    histo_data->SetFillStyle(0);
    histo_data->SetLineColor(kBlack);
    histo_data->SetLineStyle(1);
    histo_data->SetLineWidth(2);
  }

  if (full) {
    for (int i = 0; i < histos.size(); i++) histos[i]->SetFillStyle(1001);
  }

  THStack* hs = new THStack("hs","Cut Flow ");

  for (int i = 0; i < histos.size(); i++) hs->Add(histos[i]);
  if (data) hs->Add(histo_data);
  hs->Draw("hist nostack");
  //hs->Draw("nostack");


  TLegend *leg;
  if (!norm) {
    leg = new TLegend(0.60,0.70,0.9,0.9,NULL,"brNDC");
    for (int i = 0; i < histos.size(); i++) leg->AddEntry(histos[i], histos[i]->GetName(),"lep");
    if (data) leg->AddEntry(histo_data, histo_data->GetName(),"lep");

  }
  else{
    leg = new TLegend(0.60,0.75,0.9,0.9,NULL,"brNDC");
    for (int i = 0; i < histos.size(); i++){
      string name = histos[i]->GetName();
      name = name + " "+  to_string(histos[i]->GetEntries());
      leg->AddEntry(histos[i], name.c_str(),"lep");
    }
    if (data){
      string name = histo_data->GetName();
      name = name + " "+  to_string(histo_data->GetEntries());
      leg->AddEntry(histo_data, name.c_str(),"lep");
    }
  }
  leg->SetTextSize(0.035);
  leg->SetLineColorAlpha(1, 0.7);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillStyle(0);
  leg->Draw("same");
  return c_histo;
}



TCanvas* overlap2histos(TH1* histo1, TH1* histo2, bool full = false, bool norm = false){

  setTDRStyle();
  gStyle->SetOptStat(false);

  TCanvas* c_histo = new TCanvas("histo overlap","histo overlap", 800, 600);

  histo1->SetFillColorAlpha(kAzure-8, 0.85);
  histo1->SetFillStyle(1001);
  histo1->SetLineColor(kBlue);
  histo1->SetLineStyle(1);
  histo1->SetLineWidth(2);

  histo2->SetFillColorAlpha(kRed, 1);
  histo2->SetFillStyle(3004);
  histo2->SetLineColor(kRed);
  histo2->SetLineStyle(1);
  histo2->SetLineWidth(2);

  if (norm) {
    histo1->GetSumw2();
    histo2->GetSumw2();
    histo1->Scale(1./histo1->Integral());
    histo2->Scale(1./histo2->Integral());
  }
  if (full) {
    histo2->SetFillStyle(1001);
  }

  if (histo1->GetMaximum() > histo2->GetMaximum()) {
    histo1->Draw();
    histo2->Draw("same");
  }
  else{
    histo2->Draw();
    histo1->Draw("same");
  }

  TLegend *leg;
  if (!norm) {
    leg = new TLegend(0.70,0.75,0.9,0.9,NULL,"brNDC");
    leg->AddEntry(histo1, histo1->GetName(),"lep");
    leg->AddEntry(histo2, histo2->GetName(), "lep");
  }
  else{
    leg = new TLegend(0.60,0.75,0.9,0.9,NULL,"brNDC");
    string name = histo1->GetName();
    name = name + " "+  to_string(histo1->GetEntries());
    leg->AddEntry(histo1, name.c_str(), "lep");
    name = histo2->GetName();
    name = name + " "+  to_string(histo2->GetEntries());
    leg->AddEntry(histo2, name.c_str(), "lep");
  }
  leg->SetTextSize(0.035);
  leg->SetLineColorAlpha(1, 0.7);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillStyle(0);
  leg->Draw("same");
  return c_histo;
}


TH1F* CutFlow(string file_name, std::vector<string> name_dir, string histo_name = "histo_cut", bool norm = false, float cross_sec = 1) {
  std::vector<TH1F*> histos;
  TFile *file = new TFile((file_name).c_str());
  for (int i = 0; i < name_dir.size(); i++) {
    histos.push_back((TH1F*)file->Get((name_dir[i]+"/number").c_str()));
  }

  TH1F* histo_cut = new TH1F(histo_name.c_str(),histo_name.c_str(), name_dir.size(), 0, name_dir.size());
  for (int i = 0; i < name_dir.size(); i++) {
    histo_cut->SetBinContent(histo_cut->FindBin(i), histos[i]->GetSumOfWeights());
    histo_cut->GetXaxis()->SetBinLabel(i+1,(name_dir[i]).c_str());
  }

  int Norm = histo_cut->GetBinContent(1);
  histo_cut->SetName((histo_name+" "+to_string(Norm)).c_str());
  if (norm){
    histo_cut->GetSumw2();
    histo_cut->Scale(cross_sec/histo_cut->GetBinContent(1));
    std::cout << "norm " << Norm << std::endl;
  }

  return histo_cut;
}


bool OrderHistos(TH1* h1, TH1* h2){
  return h1->GetMaximum() > h2->GetMaximum();
}
