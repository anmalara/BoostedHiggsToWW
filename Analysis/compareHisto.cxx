#include <iostream>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TF1.h>



Double_t pol_line(Double_t *x, Double_t *par)
{
  Double_t x_0,x_1,x_2,x_3;

  x_0 = par[0];
  x_1 = par[1];
  x_2 = par[2];
  x_3 = par[3];

  return x_0 + x_1*x[0]+ x_2*x[0]*x[0]+ x_3*x[0]*x[0]*x[0];
}

// std::vector<TH1F> checkHistos(std::string file_name1, std::string histo_name1, std::string file_name2, std::string histo_name2, bool fit=false){
//
//   TFile *file_1=new TFile((file_name1).c_str());
//   TH1F *histo_1=(TH1F*)file_1->Get((histo_name1).c_str());
//
//   TFile *file_2=new TFile((file_name2).c_str());
//   TH1F *histo_2=(TH1F*)file_2->Get((histo_name2).c_str());
//
//   TH1F *histo_Ratio= (TH1F*)histo_1->Clone();
//   TH1F *histo_Diff= (TH1F*)histo_1->Clone();
//
//   histo_Ratio->Divide(histo_2);
//   histo_Diff->Add(histo_2,-1);
//
//   if (fit) {
//     TF1 *fit_pol = new TF1("fit_pol",pol_line,0,200,1);
//     fit_pol->SetParameter(0, 1.);
//     histo_Ratio->Fit(fit_pol,"Quiet","", 0, 200);
//
//     std::cout << std::endl;
//     std::cout << "##################################################"<< std::endl;
//     std::cout << histo_name1 << std::endl;
//     std::cout << "Ratio plot " << std::endl;
//     std::cout << "Fit param = " << fit_pol->GetParameter(0) << " +- " << fit_pol->GetParError(0) << std::endl;
//     std::cout << "ChiSquare = " << fit_pol->GetChisquare()/fit_pol->GetNDF()<< std::endl;
//     std::cout << std::endl;
//
//     fit_pol->SetParameter(0, 0.);
//     histo_Diff->Fit(fit_pol,"Quiet","", 0, 200);
//
//     std::cout << "Diff plot " << '\n';
//     std::cout << "Fit param = " << fit_pol->GetParameter(0) << " +- " << fit_pol->GetParError(0) << std::endl;
//     std::cout << "ChiSquare = " << fit_pol->GetChisquare()<< std::endl;
//     std::cout << "##################################################"<< std::endl;
//     std::cout << std::endl;
//   }
//
//   std::vector<TH1F> *histos;
//   histos.push_back(histo_Ratio);
//   histos.push_back(histo_Diff);
//   return histos;
// }
//
//
//
// TCanvas* overlapHistos(std::string file_name1, std::string histo_name1, std::string file_name2, std::string histo_name2){
//
//   TFile *file_1=new TFile((file_name1).c_str());
//   TH1F *histo_1=(TH1F*)file_1->Get((histo_name1).c_str());
//
//   TFile *file_2=new TFile((file_name2).c_str());
//   TH1F *histo_2=(TH1F*)file_2->Get((histo_name2).c_str());
//
//   TCanvas* c_histo = new TCanvas("c_histo", "c_histo", 800, 600);
//
//   histo_1->Draw();
//   histo_2->SetLineColor(kRed);
//   histo_2->Draw();
//
//   return c_histo;
// }
//
//
//
//
// void pullHisto(std::string dir_name1, std::string dir_name2, std::string file_name1, std::string file_name2, std::string folder_name1, std::string folder_name2){
//
//   gROOT->SetStyle("Plain");
//   gStyle->SetPadGridX(0);
//   gStyle->SetPadGridY(0);
//   gStyle->SetOptStat(0000);
//
//   gSystem->Exec(("mkdir -p ./output/"+file_name1+"_"+file_name2).c_str());
//   gSystem->Exec(("mkdir -p ./output/"+file_name1+"_"+file_name2+"/"+folder_name1+"_"+folder_name2 ).c_str());
//
//   TFile *file_1=new TFile((dir_name1+"/"+file_name1).c_str());
//   TDirectory *dir_1;
//
//   TList *dirlist = file_1->GetListOfKeys();
//   TIterator *iter = dirlist->MakeIterator();
//   TObject *key = iter->Next();
//
//   bool found = false;
//   while (key && !found) {
//     if (key->IsFolder() && key->GetName() == folder_name1) {
//       dir_1 = file_1->GetDirectory(key->GetName());
//       found = true;
//     }
//     key = iter->Next();
//   }
//
//   TFile *file_2=new TFile((dir_name2+"/"+file_name2).c_str());
//   TDirectory *dir_2;
//
//   dirlist = file_2->GetListOfKeys();
//   iter = dirlist->MakeIterator();
//   key = iter->Next();
//
//   found = false;
//   while (key && !found) {
//     if (key->IsFolder() && key->GetName() == folder_name2) {
//       dir_2 = file_2->GetDirectory(key->GetName());
//       found = true;
//     }
//     key = iter->Next();
//   }
//
//   TList *dirlist_1 = dir_1->GetListOfKeys();
//   TIterator *iter_1 = dirlist_1->MakeIterator();
//   TObject *key_1 = iter_1->Next();
//
//   TList *dirlist_2 = dir_2->GetListOfKeys();
//   TIterator *iter_2 = dirlist_2->MakeIterator();
//   TObject *key_2 = iter_2->Next();
//
//   TFile *outfile=new TFile(("./output/"+file_name1+"_"+file_name2+"/"+folder_name1+"_"+folder_name2+"/ratio_plot.root").c_str(), "recreate");
//
//   while (key_1 && key_2) {
//     TObject* h_1 = dir_1->Get( key_1->GetName() );
//     TObject* h_2 = dir_2->Get( key_2->GetName() );
//
//     if ( h_1->InheritsFrom( "TH1F" ) && h_2->InheritsFrom( "TH1F" ) ){
//
//       TH1F *h1=(TH1F*)file_1->Get((folder_name1+"/"+h_1->GetName()).c_str());
//       TH1F *h2=(TH1F*)file_2->Get((folder_name2+"/"+h_2->GetName()).c_str());
//
//       TH1F *h3 = (TH1F*)h2->Clone("h3");
//       h3->Divide(h1);
//       h3->SetLineColor(kBlack);
//
//       h1->SetLineColor(kBlack);
//       h2->SetLineColor(kRed);
//       h1->SetLineWidth(2);
//       h2->SetLineWidth(2);
//       h1->GetYaxis()->SetTitleSize(20);
//       h1->GetYaxis()->SetTitleFont(43);
//       h1->GetYaxis()->SetTitleOffset(1.55);
//
//       // Ratio plot (h3) settings
//       h3->SetTitle(""); // Remove the ratio title
//       h3->GetXaxis()->SetTitle(h1->GetTitle());
//       h3->GetXaxis()->SetTitleSize(20);
//       h3->GetXaxis()->SetTitleFont(43);
//       h3->GetXaxis()->SetTitleOffset(4.);
//       h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
//       h3->GetXaxis()->SetLabelSize(15);
//
//       h3->GetYaxis()->SetTitle("ratio h2/h1");
//       h3->GetYaxis()->SetTitleSize(20);
//       h3->GetYaxis()->SetTitleFont(43);
//       h3->GetYaxis()->SetTitleOffset(1.55);
//       h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
//       h3->GetYaxis()->SetLabelSize(15);
//       h3->SetStats(0);      // No statistics on lower plot
//
//       double xPad = 0.3;
//       TCanvas *c=new TCanvas("c", "c", 800, 800);
//       c->SetFillStyle(4000);
//       c->SetFrameFillColor(0);
//
//       TPad *p_1=new TPad("p_1", "p_1", 0, xPad, 1, 1);
//       p_1->SetFillStyle(4000);
//       p_1->SetFrameFillColor(0);
//       p_1->SetBottomMargin(0.01);
//
//       TPad* p_2 = new TPad("p_2", "p_2",0,0,1,xPad);
//       p_2->SetTopMargin(0.01);
//       p_2->SetBottomMargin((1.-xPad)/xPad*0.09);
//       p_2->SetFillColor(0);
//       p_2->SetBorderMode(0);
//       p_2->SetBorderSize(2);
//       p_2->SetFrameBorderMode(0);
//       p_2->SetFrameBorderMode(0);
//
//       p_1->Draw();
//       p_2->Draw();
//       p_1->cd();
//       h1->Draw();
//       h2->Draw("same");
//       p_2->cd();
//       h3->Draw();       // Draw the ratio plot
//       // break;
//
//       c->Write();
//       std::string name_pdf = h_1->GetName();
//       c->SaveAs(("./output/"+file_name1+"_"+file_name2+"/"+folder_name1+"_"+folder_name2+"/"+name_pdf+".pdf").c_str());
//       c->SaveAs(("./output/"+file_name1+"_"+file_name2+"/"+folder_name1+"_"+folder_name2+"/"+name_pdf+".png").c_str());
//
//       delete c;
//     }
//     key_1 = iter_1->Next();
//     key_2 = iter_2->Next();
//   }
//
//   outfile->Close();
//
// }
//
//
//
// void ScrollFile(std::string dir_name1, std::string dir_name2, std::string file_name1, std::string file_name2, std::string folder_name1, std::string folder_name2, std::string histo_name1, std::string histo_name2, double rebin){
//
//   gROOT->SetStyle("Plain");
//   gStyle->SetPadGridX(0);
//   gStyle->SetPadGridY(0);
//   gStyle->SetOptStat(0000);
//
//   TFile *file_1 = TFile::Open("../OutputFile/uhh2.AnalysisModuleRunner.Data.Output_19_01.root");
//   TFile *file_2 = TFile::Open("../OutputFile/uhh2.AnalysisModuleRunner.Data.19_01.root");
//   TList *dirlist_1 = file_1->GetListOfKeys();
//   TIterator *iter_1 = dirlist_1->MakeIterator();
//   TObject *key_1 = iter_1->Next();
//   std::vector<string> dirs;
//   while (key_1) {
//     if (key_1->IsFolder()) {
//       TDirectory *dir = file_1->GetDirectory(key_1->GetName());
//       if (dir) {
//         TList *dirlist_2 = dir->GetListOfKeys();
//         TIterator *iter_2 = dirlist_2->MakeIterator();
//         TObject *key_2 = iter_2->Next();
//
//         while (key_2) {
//           TObject* histo = dir->Get( key_2->GetName() );
//           if (histo->InheritsFrom( "TH1F" )) {
//             std::string dir_name = dir->GetName();
//             std::string histo_name = key_2->GetName();
//             // std::cout << dir_name+"/"+histo_name << '\n';
//             TH1F *histo_2=(TH1F*)file_2->Get((dir_name+"/"+histo_name).c_str());
//             histo->Draw();
//             histo_2->SetLineColor(kRed);
//             histo_2->Scale(1./10);
//             histo_2->Draw("same");
//             break;
//           }
//           key_2 = iter_2->Next();
//         }
//       }
//     }
//     key_1 = iter_1->Next();
//   }
// }
//
//




void MCTruthMatching(){
  std::string file_name_det = "../file/outputfile/feasibilitystudy_signalnew/uhh2.AnalysisModuleRunner.MC.signal_new.root";
  std::string file_name_had = "../file/outputfile/MCTruthMatching/uhh2.AnalysisModuleRunner.MC.signal_new.root";
  std::string histo_name1 = "boostedJet_angular_matching_H/SDmass_jet1";
  std::string histo_name2 = "gen_jets_ordered/mass_jet1";
  std::string histo_name3 = "boostedJet_angular_match/SDmass_jet1";
  std::string histo_name4 = "boostedJet_angular_sel/SDmass_jet1";

  // TFile *file_det = new TFile((file_name_det).c_str());
  TFile *file_had = new TFile((file_name_had).c_str());
  TH1F *histo_1 = (TH1F*)file_had->Get((histo_name1).c_str());
  TH1F *histo_2 = (TH1F*)file_had->Get((histo_name2).c_str());
  TH1F *histo_3 = (TH1F*)file_had->Get((histo_name3).c_str());
  TH1F *histo_4 = (TH1F*)file_had->Get((histo_name4).c_str());

  TH1F *histo_acc= (TH1F*)histo_1->Clone();
  TH1F *histo_bkg= (TH1F*)histo_4->Clone();

  histo_acc->Divide(histo_2);
  histo_bkg->Add(histo_3,-1);
  histo_bkg->Divide(histo_4);
  histo_bkg->SetLineColor(kRed);

  bool fit = false;
  if (fit) {
    TF1 *fit_pol = new TF1("fit_pol",pol_line,0,200,1);
    fit_pol->SetParameter(0, 1.);
    histo_acc->Fit(fit_pol,"Quiet","", 0, 200);

    std::cout << std::endl;
    std::cout << "##################################################"<< std::endl;
    std::cout << histo_name1 << std::endl;
    std::cout << "Ratio plot " << std::endl;
    std::cout << "Fit param = " << fit_pol->GetParameter(0) << " +- " << fit_pol->GetParError(0) << std::endl;
    std::cout << "ChiSquare = " << fit_pol->GetChisquare()/fit_pol->GetNDF()<< std::endl;
    std::cout << std::endl;

    fit_pol->SetParameter(0, 0.);
    histo_bkg->Fit(fit_pol,"Quiet","", 0, 200);

    std::cout << "Diff plot " << '\n';
    std::cout << "Fit param = " << fit_pol->GetParameter(0) << " +- " << fit_pol->GetParError(0) << std::endl;
    std::cout << "ChiSquare = " << fit_pol->GetChisquare()<< std::endl;
    std::cout << "##################################################"<< std::endl;
    std::cout << std::endl;
  }

  new TCanvas;
  histo_acc->Draw();
  histo_bkg->Draw("same");
}
