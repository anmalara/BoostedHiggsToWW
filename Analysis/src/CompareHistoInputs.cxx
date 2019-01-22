#include <vector>
#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/tdrstyle_all.C"

#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/Analysis/macros/AnalysisTree.h"

void PlotHistos(std::vector<TString> samples, TString& dir) {
  dir = "compareInputs/plots/"+dir;
  gSystem->Exec("mkdir -p "+dir);
  std::vector<TString> variables = {"jetPt", "jetEta", "jetPhi", "jetMass", "jetEnergy", "jetTau1", "jetTau2", "jetTau3", "ncandidates", "jetBtag", "jetTau21", "jetTau31", "jetTau32"};
  std::map<TString,int> colors= { {"NN_HZ",kRed+1}, {"NN_QCD",kGreen-2}, {"NN_TT",kBlue},
  {"presel_HZ",kRed}, {"presel_DY",kGreen}, {"presel_TT",kViolet-3}, {"origin_HZ",kRed-1},
  {"origin_DY",kSpring-9}, {"origin_DY1",kSpring-7}, {"origin_DY2",kSpring+10},{"origin_DY3",kGreen+1},{"origin_DY4",kGreen+3},
  {"origin_TT",kAzure-4},{"origin_TTBarhad",kAzure-4},{"origin_TTBarlep",kCyan},{"origin_TTBarsemilep",kAzure+6}};

  std::map<TString, std::map<TString, TH1F*> > map_histos;

  for (auto &sample : samples) {
    TH1::AddDirectory(kFALSE);
    TFile *file = new TFile("compareInputs/h_"+sample+".root");
    for (auto & var : variables) map_histos[sample][var] = (TH1F*) ((TH1F*) file->Get(var))->Clone();
    file->Close();
    std::cout << sample << " " << map_histos[sample]["jetPt"]->GetEntries() << '\n';
  }

  for (auto & var : variables) {
    double xmin = 10.;
    double ymin = 10.;
    double xmax = 0;
    double ymax = 0;
    for (auto &sample : samples) {
      map_histos[sample][var]->Scale(1./map_histos[sample][var]->Integral());
      xmin = std::min(xmin, map_histos[sample][var]->GetXaxis()->GetXmin());
      ymin = std::min(ymin, 0.00000001);
      xmax = std::max(xmax, map_histos[sample][var]->GetXaxis()->GetXmax());
      ymax = std::max(ymax, map_histos[sample][var]->GetMaximum());
    }

    if (var== "jetEnergy" || var == "jetPt") xmin = 300;
    if (var== "jetTau2" || var == "jetTau3") xmax = 0.5;

    TCanvas* c_histo = tdrCanvas(var, xmin, xmax, ymin, 1.5*ymax, var, "A.U.");
    TLegend *leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);

    for (auto &sample : samples) {
      tdrDraw(map_histos[sample][var], "same", kFullCircle, colors[sample], 1, colors[sample], 0, colors[sample]);
      leg->AddEntry(map_histos[sample][var], sample ,"lep");
    }
    leg->SetLineColorAlpha(1, 0.7);
    leg->Draw("same");

    c_histo->SaveAs(dir+"/"+var+".pdf", "pdf");
    c_histo->SaveAs(dir+"/"+var+".png", "png");
  }

}



void CompareHistoInputs() {

  std::vector<TString> samples_presel = {"presel_DY","presel_HZ","presel_TT"};
  std::vector<TString> samples_origin = {"origin_TT", "origin_HZ", "origin_DY", "origin_DY1", "origin_DY2", "origin_DY3", "origin_DY4", "origin_TTBarhad", "origin_TTBarlep", "origin_TTBarsemilep"};
  std::vector<TString> samples_NN = {"NN_HZ","NN_QCD","NN_TT"};

  { TString dir = "all";    std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "presel_HZ", "presel_DY", "presel_TT", "origin_HZ", "origin_DY", "origin_TT", "origin_DY1", "origin_DY2", "origin_DY3", "origin_DY4", "origin_TTBarhad", "origin_TTBarlep", "origin_TTBarsemilep"}; PlotHistos(samples,dir); }
  { TString dir = "origin"; std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "origin_HZ", "origin_DY", "origin_TT"}; PlotHistos(samples,dir); }
  { TString dir = "presel"; std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "presel_HZ", "presel_DY", "presel_TT"}; PlotHistos(samples,dir); }
  { TString dir = "HZ";     std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "presel_HZ", "origin_HZ"}; PlotHistos(samples,dir); }
  { TString dir = "DY";     std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "presel_DY", "origin_DY"}; PlotHistos(samples,dir); }
  { TString dir = "TT";     std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "presel_TT", "origin_TT"}; PlotHistos(samples,dir); }
  { TString dir = "DY_1";   std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "presel_DY", "origin_DY", "origin_DY1", "origin_DY2", "origin_DY3", "origin_DY4", "origin_TTBarhad", "origin_TTBarlep", "origin_TTBarsemilep"}; PlotHistos(samples,dir); }
  { TString dir = "TT_1";   std::vector<TString> samples = {"NN_HZ","NN_QCD","NN_TT", "presel_TT", "origin_TT", "origin_TTBarhad", "origin_TTBarlep", "origin_TTBarsemilep"}; PlotHistos(samples,dir); }


}


void PlotNNOut() {

  TString dir = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/Analysis/src/compareInputs/plots/";
  TTree* tH = (TTree*)_file1->Get("AnalysisTree");
  TTree* tD = (TTree*)_file0->Get("AnalysisTree");
  TTree* tT = (TTree*)_file3->Get("AnalysisTree");
  tH->Draw("isHiggs>>hHH");
  tH->Draw("isQCD>>hHQ");
  tH->Draw("isTop>>hHT");
  tD->Draw("isHiggs>>hDH");
  tD->Draw("isQCD>>hDQ");
  tD->Draw("isTop>>hDT");
  tT->Draw("isHiggs>>hTH");
  tT->Draw("isQCD>>hTQ");
  tT->Draw("isTop>>hTT");

  std::map<TString,int> colors= { {"HZ",kRed+1}, {"QCD",kGreen-2}, {"TT",kBlue}};
  TCanvas* c_histo;
  TLegend *leg;

  c_histo = tdrCanvas("NNout", -0.01, 1.01, 0.0001, 0.5, "isHiggs", "A.U."); c_histo->SetLogy();
  leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);
  hHH->Scale(1./hHH->Integral()); tdrDraw(hHH, "same", kFullCircle, colors["HZ"], 1, colors["HZ"], 0, colors["HZ"]); leg->AddEntry(hHH, "HZ" ,"lep");
  hDH->Scale(1./hDH->Integral()); tdrDraw(hDH, "same", kFullCircle, colors["QCD"], 1, colors["QCD"], 0, colors["QCD"]); leg->AddEntry(hDH, "DY" ,"lep");
  hTH->Scale(1./hTH->Integral()); tdrDraw(hTH, "same", kFullCircle, colors["TT"], 1, colors["TT"], 0, colors["TT"]); leg->AddEntry(hTH, "TT" ,"lep");
  leg->Draw();
  c_histo->SaveAs(dir+"isHiggs.pdf", "pdf");
  c_histo->SaveAs(dir+"isHiggs.png", "png");

  c_histo = tdrCanvas("NNout", -0.01, 1.01, 0.0001, 0.5, "isQCD", "A.U."); c_histo->SetLogy();
  leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);
  hHQ->Scale(1./hHQ->Integral()); tdrDraw(hHQ, "same", kFullCircle, colors["HZ"], 1, colors["HZ"], 0, colors["HZ"]); leg->AddEntry(hHQ, "HZ" ,"lep");
  hDQ->Scale(1./hDQ->Integral()); tdrDraw(hDQ, "same", kFullCircle, colors["QCD"], 1, colors["QCD"], 0, colors["QCD"]); leg->AddEntry(hDQ, "DY" ,"lep");
  hTQ->Scale(1./hTQ->Integral()); tdrDraw(hTQ, "same", kFullCircle, colors["TT"], 1, colors["TT"], 0, colors["TT"]); leg->AddEntry(hTQ, "TT" ,"lep");
  leg->Draw();
  c_histo->SaveAs(dir+"isQCD.pdf", "pdf");
  c_histo->SaveAs(dir+"isQCD.png", "png");

  c_histo = tdrCanvas("NNout", -0.01, 1.01, 0.0001, 0.5, "isTop", "A.U."); c_histo->SetLogy();
  leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);
  hHT->Scale(1./hHT->Integral()); tdrDraw(hHT, "same", kFullCircle, colors["HZ"], 1, colors["HZ"], 0, colors["HZ"]); leg->AddEntry(hHT, "HZ" ,"lep");
  hDT->Scale(1./hDT->Integral()); tdrDraw(hDT, "same", kFullCircle, colors["QCD"], 1, colors["QCD"], 0, colors["QCD"]); leg->AddEntry(hDT, "DY" ,"lep");
  hTT->Scale(1./hTT->Integral()); tdrDraw(hTT, "same", kFullCircle, colors["TT"], 1, colors["TT"], 0, colors["TT"]); leg->AddEntry(hTT, "DTT" ,"lep");
  leg->Draw();
  c_histo->SaveAs(dir+"isTop.pdf", "pdf");
  c_histo->SaveAs(dir+"isTop.png", "png");

}
