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
#include "../src/BoostedHiggsToWWUtils.cxx"

TCanvas* overlap2histos(TH1* histo1, TH1* histo2, bool full=false, bool norm=false);
TCanvas* overlap3histos(TH1* histo1, TH1* histo2, TH1* histo3, bool full=false);
TCanvas* overlap4histos(TH1* histo1, TH1* histo2, TH1* histo3, TH1* histo4, bool full=false);
TCanvas* overlapHistos(std::string file_name1, std::string histo_name1, std::string file_name2, std::string histo_name2, std::string file_name3, std::string histo_name3, std::string file_name4, std::string histo_name4);
TCanvas* overlaphistos(std::vector<TH1*> histos, TH1* histo_data, bool data = false, bool norm = false, bool full = false);
TH1F* CutFlow(string file_name, std::vector<string> name_dir, string histo_name = "histo_cut", bool norm = false);

void HistoFunction(bool quick = true, string outdir = "./") {


    // string file_name_signal = "/Users/andrea/Desktop/untitledfolder/uhh2.AnalysisModuleRunner.MC.Signal.root";
    // string file_name_bkg = "/Users/andrea/Desktop/untitledfolder/uhh2.AnalysisModuleRunner.MC.Backgroung.root";
    //string file_name_signal = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.Signal.root";
    //string file_name_bkg = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.Backgroung.root";
    //string file_name_bkg_DY = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.Backgroung.root";
    string file_name_signal = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.signal.root";
    string file_name_bkg_DY = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.bkg_DYJetsToLL_Pt-250To400.root";
    string file_name_bkg_ZZ = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.bkg_ZZ.root";
    string file_name_bkg_ZW = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.bkg_ZW.root";
    string file_name_bkg_TT = "../file/outputfile/uhh2.AnalysisModuleRunner.MC.bkg_TT.root";

    TFile *file_signal=new TFile((file_name_signal).c_str());
    TFile *file_bkg_DY=new TFile((file_name_bkg_TT).c_str());
    TFile *file_bkg_ZZ=new TFile((file_name_bkg_TT).c_str());
    TFile *file_bkg_ZW=new TFile((file_name_bkg_TT).c_str());
    TFile *file_bkg_TT=new TFile((file_name_bkg_TT).c_str());


    std::vector<string> name_dir;
    std::vector<string> name_histo;


        TList *dirlist = file_signal->GetListOfKeys();
        TIterator *iter = dirlist->MakeIterator();
        TObject *key = iter->Next();

        while (key) {
          string name_avoid = key->GetName();
          if (name_avoid == "SFrame") {
            break;
          }
          if (key->IsFolder()) {
            TDirectory *dir = file_signal->GetDirectory(key->GetName());
            string name = dir->GetName();
            name_dir.push_back(name);
            gSystem->Exec(("mkdir -p ../file/outputfile/"+name).c_str());
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


    TH1F* histo_cut_signal = CutFlow(file_name_signal, name_dir, "signal", true );
    //TH1F* histo_cut_bkg = CutFlow(file_name_bkg, name_dir, "DY", true);
    TH1F* histo_cut_bkg_DY = CutFlow(file_name_bkg_DY, name_dir, "DY", true);
    TH1F* histo_cut_bkg_ZZ = CutFlow(file_name_bkg_ZZ, name_dir, "ZZ", true);
    TH1F* histo_cut_bkg_ZW = CutFlow(file_name_bkg_ZW, name_dir, "ZW", true);
    TH1F* histo_cut_bkg_TT = CutFlow(file_name_bkg_TT, name_dir, "TT", true);

    //histo_cut_signal->Scale(1./10000000);

    std::vector<TH1*> histo_cuts;
    histo_cuts.push_back(histo_cut_signal);
    //histo_cuts.push_back(histo_cut_bkg);
    histo_cuts.push_back(histo_cut_bkg_DY);
    histo_cuts.push_back(histo_cut_bkg_ZZ);
    histo_cuts.push_back(histo_cut_bkg_ZW);
    histo_cuts.push_back(histo_cut_bkg_TT);

    string name = "../file/outputfile/"+outdir+"cutflow.png";
    //TCanvas* c_cutflow = overlap2histos(histo_cut_signal, histo_cut_bkg);
    //TCanvas* c_cutflow = overlaphistos(histo_cuts, histo_cut_signal);
    TCanvas* c_cutflow = overlaphistos(histo_cuts, histo_cut_signal);
    c_cutflow->SaveAs((name).c_str());



    std::vector<TCanvas*> canvases;
    TH1F *histo1, *histo2, *histotemp;
    TCanvas *c;
    if(quick) name_histo = {};
    for (int i = 0; i < name_histo.size(); i++) {
      histotemp=(TH1F*)file_signal->Get((name_histo[i]).c_str());
      name = histotemp->GetName();
      histo1= (TH1F*)histotemp->Clone((name+" Signal").c_str() );
      histotemp=(TH1F*)file_signal->Get((name_histo[i]).c_str());
      histo2= (TH1F*)histotemp->Clone((name+" Bkg").c_str() );
      c = overlap2histos(histo1, histo2, false, true);
      name = "../file/outputfile/"+outdir+name_histo[i]+".png";
      c->SaveAs((name).c_str());
    }
  }


void A_M(int max1,int max2,int max3,int max4){
  if (max1 > max2) {
    if (max1 > max3) {
      if (max1 > max4) {
        std::cout << max1 << '\n';
      } else {
        std::cout << max4 << '\n';
      }
    } else if (max3 > max4) {
      std::cout << max3 << '\n';
    } else {
      std::cout << max4 << '\n';
    }
  } else if (max2 > max3) {
    if (max2 > max4) {
      std::cout << max2 << '\n';
    } else {
      std::cout << max4 << '\n';
    }
  } else {
    if (max3 > max4) {
      std::cout << max3 << '\n';
    } else {
      std::cout << max4 << '\n';
    }
  }
}




TCanvas* overlaphistos(std::vector<TH1*> histos, TH1* histo_data, bool data = false, bool norm = false, bool full = false){
    setTDRStyle();
    gStyle->SetOptStat(false);

    TCanvas* c_histo = new TCanvas("histo overlap","histo overlap", 800, 600);
    std::vector<int> colors = {kAzure-8, kRed, kOrange, kGreen+1, kViolet +2};

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

    THStack* hs = new THStack("hs","test stacked histograms");

    for (int i = 0; i < histos.size(); i++) hs->Add(histos[i]);
    if (data) hs->Add(histo_data);
    hs->Draw("hist nostack");


    TLegend *leg;
    if (!norm) {
      leg = new TLegend(0.70,0.75,0.9,0.9,NULL,"brNDC");
      for (int i = 0; i < histos.size(); i++) leg->AddEntry(histos[i], histos[i]->GetName(),"lep");
      if (data) leg->AddEntry(histo_data, histo_data->GetName(),"lep");

    }
    else{
      leg = new TLegend(0.60,0.75,0.9,0.9,NULL,"brNDC");
      for (int i = 0; i < histos.size(); i++){
	string name = histos[i]->GetName();
	name = name + " "+  itoa(histos[i]->GetEntries());
	leg->AddEntry(histos[i], name.c_str(),"lep");
      }
      if (data){
        string name = histo_data->GetName();
        name = name + " "+  itoa(histo_data->GetEntries());
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
      name = name + " "+  itoa(histo1->GetEntries());
      leg->AddEntry(histo1, name.c_str(), "lep");
      name = histo2->GetName();
      name = name + " "+  itoa(histo2->GetEntries());
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

TCanvas* overlap3histos(TH1* histo1, TH1* histo2, TH1* histo3, bool full=false){
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

        histo3->SetFillColorAlpha(kOrange, 1);
        histo3->SetFillStyle(3004);
        histo3->SetLineColor(kOrange);
        histo3->SetLineStyle(1);
        histo3->SetLineWidth(2);

        if (full) {
          histo2->SetFillStyle(1001);
          histo3->SetFillStyle(1001);
        }

        double max1 = histo1->GetMaximum();
        double max2 = histo2->GetMaximum();
        double max3 = histo3->GetMaximum();

        if (max1 > max2) {
          if (max1 > max3) {
            histo1->Draw("same");
            histo2->Draw("same");
            histo3->Draw("same");
          } else {
            histo3->Draw("same");
            histo1->Draw("same");
            histo2->Draw("same");
          }
        } else {
          if (max2 > max3) {
            histo2->Draw("same");
            histo1->Draw("same");
            histo3->Draw("same");
          }
          else{
            histo3->Draw("same");
            histo1->Draw("same");
            histo2->Draw("same");
          }
        }


        TLegend *leg = new TLegend(0.70,0.75,0.9,0.9,NULL,"brNDC");
        leg->SetTextSize(0.035);
        leg->SetLineColorAlpha(1, 0.7);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillStyle(0);
        leg->AddEntry(histo1, histo1->GetName(),"lep");
        leg->AddEntry(histo2, histo2->GetName(), "lep");
        leg->AddEntry(histo3, histo3->GetName(), "lep");
        leg->Draw("same");
        return c_histo;
      }

    TCanvas* overlap4histos(TH1* histo1, TH1* histo2, TH1* histo3, TH1* histo4, bool full=false){

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

        histo3->SetFillColorAlpha(kOrange, 1);
        histo3->SetFillStyle(3004);
        histo3->SetLineColor(kOrange);
        histo3->SetLineStyle(1);
        histo3->SetLineWidth(2);

        histo4->SetFillColorAlpha(kGreen+1, 1);
        histo4->SetFillStyle(3003);
        histo4->SetLineColor(kGreen+1);
        histo4->SetLineStyle(1);
        histo4->SetLineWidth(2);

        if (full) {
          histo2->SetFillStyle(1001);
          histo3->SetFillStyle(1001);
          histo4->SetFillStyle(1001);
        }

        double max1 = histo1->GetMaximum();
        double max2 = histo2->GetMaximum();
        double max3 = histo3->GetMaximum();
        double max4 = histo3->GetMaximum();

        if (max1 > max2) {
          if (max1 > max3) {
            if (max1 > max4) { histo1->Draw("same"); histo2->Draw("same"); histo3->Draw("same"); histo4->Draw("same");
            } else { histo4->Draw("same"); histo1->Draw("same"); histo2->Draw("same"); histo3->Draw("same");}
          } else if (max3 > max4) { histo3->Draw("same"); histo1->Draw("same"); histo2->Draw("same"); histo4->Draw("same");
          } else { histo4->Draw("same"); histo1->Draw("same"); histo2->Draw("same"); histo3->Draw("same");}
        } else if (max2 > max3) {
          if (max2 > max4) { histo2->Draw("same"); histo1->Draw("same"); histo3->Draw("same"); histo4->Draw("same");
        } else {histo4->Draw("same"); histo1->Draw("same"); histo2->Draw("same"); histo3->Draw("same");}
        } else {
          if (max3 > max4) { histo3->Draw("same"); histo1->Draw("same"); histo2->Draw("same"); histo4->Draw("same");
          } else { histo4->Draw("same"); histo1->Draw("same"); histo2->Draw("same"); histo3->Draw("same");}
        }


        TLegend *leg = new TLegend(0.70,0.75,0.9,0.9,NULL,"brNDC");
        leg->SetTextSize(0.035);
        leg->SetLineColorAlpha(1, 0.7);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillStyle(0);
        leg->AddEntry(histo1, histo1->GetName(),"lep");
        leg->AddEntry(histo2, histo2->GetName(), "lep");
        leg->AddEntry(histo3, histo3->GetName(), "lep");
        leg->AddEntry(histo4, histo4->GetName(), "lep");
        leg->Draw("same");
        return c_histo;
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

          TCanvas* c_histo = new TCanvas("histo overlap", "histo overlap", 800, 600);

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

        TH1F* CutFlow(string file_name, std::vector<string> name_dir, string histo_name = "histo_cut", bool norm = false) {
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
	  
	  if (norm){
                int Norm = histo_cut->GetBinContent(1);
                histo_cut->SetName((histo_name+" "+itoa(Norm)).c_str());
		histo_cut->GetSumw2();
                histo_cut->Scale(1./histo_cut->GetBinContent(1));
                std::cout << "norm " << Norm << std::endl;
            }

          return histo_cut;
        }
