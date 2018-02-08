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


TH1F* CutFlow(string file_name, std::vector<string> name_dir) {
  // TH1F* CutFlow(string file_name, string histo_ref_name_, std::vector<string> name_dir) {

  // string histo_ref_name = histo_ref_name_+"/numbers";
  std::vector<TH1F*> histos;
  TFile *file = new TFile((file_name).c_str());
  // TH1F* histo_ref = (TH1F*)file->Get((histo_ref_name).c_str());
  for (int i = 0; i < name_dir.size(); i++) {
    histos.push_back((TH1F*)file->Get((name_dir[i]+"/number").c_str()));
  }

  TH1F* histo_cut = new TH1F("histo_cut","histo_cut", name_dir.size(), 0, name_dir.size());

  for (int i = 0; i < name_dir.size(); i++) {
    histo_cut->SetBinContent(histo_cut->FindBin(i), histos[i]->GetEntries());
    histo_cut->GetXaxis()->SetBinLabel(i+1,(name_dir[i]).c_str());
    // histo_cut->SetBinContent(histo_cut->FindBin(i), histos[i]->GetEntries());
    // h->GetXaxis()->SetBinLabel(i,people[i-1]);
  }
  return histo_cut;
}


TCanvas* overlapHistos(std::string file_name1, std::string histo_name1,
  std::string file_name2, std::string histo_name2,
  std::string file_name3, std::string histo_name3,
  std::string file_name4, std::string histo_name4){

    setTDRStyle();
    TFile *file_1=new TFile((file_name1).c_str());
    TH1F *histo_1=(TH1F*)file_1->Get((histo_name1).c_str());

    TFile *file_2=new TFile((file_name2).c_str());
    TH1F *histo_2=(TH1F*)file_2->Get((histo_name2).c_str());

    TFile *file_3=new TFile((file_name3).c_str());
    TH1F *histo_3=(TH1F*)file_3->Get((histo_name3).c_str());

    TFile *file_4=new TFile((file_name4).c_str());
    TH1F *histo_4=(TH1F*)file_4->Get((histo_name4).c_str());

    TCanvas* c_histo = new TCanvas(histo_name1.c_str(), histo_name1.c_str(), 800, 600);

    histo_2->SetLineColor(kRed);
    histo_4->SetLineColor(kRed);
    histo_3->SetLineStyle(kDashed);
    histo_4->SetLineStyle(kDashed);
    histo_1->Draw();
    histo_2->Draw("same");
    histo_3->Draw("same");
    histo_4->Draw("same");

    return c_histo;
  }


  void DisplayCuts() {

    std::vector<string> name_dir = { "nJet_noCuts", "nJet_sel", "nJet_sel_1", "ele_sel", "ele_sel_2","muon_sel", "muon_sel_2", "diLepton_sel"} ;

    // string file_name_signal = "/Users/andrea/Desktop/untitledfolder/uhh2.AnalysisModuleRunner.MC.Signal.root";
    // string file_name_bkg = "/Users/andrea/Desktop/untitledfolder/uhh2.AnalysisModuleRunner.MC.Backgroung.root";
    string file_name_signal = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.Signal.root";
    string file_name_bkg = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.Backgroung.root";

    string histo_ref_name = "nJet_noCuts/number_jets";

    TH1F* histo_cut_signal = CutFlow(file_name_signal, name_dir);
    TH1F* histo_cut_bkg = CutFlow(file_name_bkg, name_dir);
    // TH1F* histo_cut_signal = CutFlow(file_name_signal, histo_ref_name, name_dir);
    // TH1F* histo_cut_bkg = CutFlow(file_name_bkg, histo_ref_name, name_dir);

    TCanvas* c_histo = new TCanvas("c_histo", "c_histo", 800, 600);

    histo_cut_bkg->SetLineColor(kRed);
    histo_cut_bkg->Draw();
    histo_cut_signal->Draw("same");

    std::vector<string> name_histo_1 = { "nJet_noCuts/pt_jet", "nJet_noCuts/pt_jet", "ele_noCuts/pt", "ele_noCuts/pt", "muon_noCuts/pt", "muon_noCuts/pt", "diLepton_noCuts/diele_pt"} ;
    std::vector<string> name_histo_2 = { "nJet_sel/pt_jet", "nJet_sel_1/pt_jet", "ele_sel/pt", "ele_sel_2/pt", "muon_sel/pt", "muon_sel_2/pt", "diLepton_sel/diele_pt"} ;

    std::vector<TCanvas*> canvases;
    for (int i = 0; i < name_histo_1.size(); i++) {
      // canvases.push_back(overlapHistos(file_name_signal, name_histo_1[i], file_name_bkg, name_histo_1[i]));
      // canvases.push_back(overlapHistos(file_name_signal, name_histo_1[i], file_name_signal, name_histo_2[i], file_name_bkg, name_histo_1[i], file_name_bkg, name_histo_2[i]));
    }
  }
