#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/Analysis/include/Plotter.hpp"

// ClassImp(Plotter);




/*
██       ██████   █████  ██████       █████  ███    ██ ██████      ██████  ██████  ██ ███    ██ ████████
██      ██    ██ ██   ██ ██   ██     ██   ██ ████   ██ ██   ██     ██   ██ ██   ██ ██ ████   ██    ██
██      ██    ██ ███████ ██   ██     ███████ ██ ██  ██ ██   ██     ██████  ██████  ██ ██ ██  ██    ██
██      ██    ██ ██   ██ ██   ██     ██   ██ ██  ██ ██ ██   ██     ██      ██   ██ ██ ██  ██ ██    ██
███████  ██████  ██   ██ ██████      ██   ██ ██   ████ ██████      ██      ██   ██ ██ ██   ████    ██
*/



void Plotter::LoadMap(){
  TH1::AddDirectory(kFALSE);
  MapHistos.clear();
  for (auto &FileName: NameFiles) {
    TFile *file = new TFile(InputDir+FileName+".root");
    for (auto &Tag: HistoTags) {
      for (auto &FoldName: RetrieveKeys(NameHistTypeHistMap)) {
        for (auto &VarName: TypeHistVarNameMap[NameHistTypeHistMap[FoldName]]) {
          if (MapHistos[FileName][Tag][FoldName][VarName] == 0) std::cout << FileName << "\t" << Tag << "\t" << FoldName << "\t" << VarName << "\t" << MapHistos[FileName][Tag][FoldName][VarName] << '\n';
          MapHistos[FileName][Tag][FoldName][VarName] = (TH1F*) ((TH1F*)file->Get(FoldName+"_"+Tag+"/"+VarName))->Clone();
        }
      }
    }
    file->Close();
  }
  std::cout << "LOAD END" << '\n';
}

void Plotter::PrintMap(bool debug){
  for (auto const& FileName : MapHistos)
  for (auto const& Tag : FileName.second)
  for (auto const& FoldName : Tag.second)
  for (auto const& VarName : FoldName.second)
  if (debug || !VarName.second)
  std::cout << FileName.first << "\t" << Tag.first << "\t" << FoldName.first << "\t" << VarName.first <<  "\t"<< VarName.second << '\n';
}





/*
██████  ██       ██████  ████████ ████████ ██ ███    ██  ██████      ██       ██████   ██████  ██████  ███████
██   ██ ██      ██    ██    ██       ██    ██ ████   ██ ██           ██      ██    ██ ██    ██ ██   ██ ██
██████  ██      ██    ██    ██       ██    ██ ██ ██  ██ ██   ███     ██      ██    ██ ██    ██ ██████  ███████
██      ██      ██    ██    ██       ██    ██ ██  ██ ██ ██    ██     ██      ██    ██ ██    ██ ██           ██
██      ███████  ██████     ██       ██    ██ ██   ████  ██████      ███████  ██████   ██████  ██      ███████
*/








void Plotter::PlotHistosComposition(TString FileName, TString FoldName_, bool isLog) {
  for (auto &Tag: HistoTags) {
    for (auto &VarName: TypeHistVarNameMap[NameHistTypeHistMap[FoldName_]]) {
      PlotSingleHistoComposition(FileName,Tag,FoldName_,VarName,isLog);
    }
  }
}

void Plotter::SortByFolderName(std::vector<TString> & v, TString FileName, TString Tag, TString VarName){
  for (size_t i = 0; i < v.size(); i++) for (size_t j = 0; j < v.size(); j++) {
    if (MapHistos[FileName][Tag][v[i]][VarName]->Integral() < MapHistos[FileName][Tag][v[j]][VarName]->Integral())
    iter_swap(v.begin() + i, v.begin() + j);
  }
}

void Plotter::PlotSingleHistoComposition(TString FileName, TString Tag, TString FoldName_, TString VarName, bool isLog) {
  SetPlotLimits(NameFolders, FileName, Tag, VarName, isLog);

  TCanvas* c_histo = tdrCanvas(Tag+FoldName_+VarName, xmin, xmax, ymin, factor*ymax, VarName, "A.U.");
  c_histo->SetLogy(isLog);
  TLegend *leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);

  THStack* hs = new THStack(FoldName_+VarName, FoldName_+VarName);
  hs->GetYaxis()->SetRangeUser(ymin, factor*ymax);

  SortByFolderName(NameFolders, FileName, Tag, VarName);
  for (auto &FoldName : NameFolders) {
    hs->Add(MapHistos[FileName][Tag][FoldName][VarName]);
    SetLegend(leg, MapHistos[FileName][Tag][FoldName][VarName], FoldName, colors[FoldName], 1001);
  }

  hs->Draw( "hist");
  SetLegend(leg, MapHistos[FileName][Tag][FoldName_][VarName], FoldName_, colors[FoldName_], 1001);

  leg->SetLineColorAlpha(1, 0.7);
  leg->Draw("same");

  if (isSave) SaveCanvas(c_histo, OutDir+FileName+"/Composition/"+Tag+"/"+FoldName_+"/", VarName);

}






void Plotter::PlotMaps(bool isLog, bool isNorm) {
  for (auto &FileName: NameFiles) {
    for (auto &FoldName: NameFolders) {
      for (auto &VarName: NameVars) {
        PlotMaps(FileName, FoldName, VarName, isLog, isNorm);
      }
    }
  }
}


void Plotter::PlotAllHistos(bool isLog, bool isNorm) {
  for (auto &FileName: NameFiles) {
    for (auto &Tag: HistoTags) {
      for (auto &FoldName: RetrieveKeys(NameHistTypeHistMap)) {
        for (auto &VarName: TypeHistVarNameMap[NameHistTypeHistMap[FoldName]]) {
          PlotSingleHisto(FileName, Tag, FoldName, VarName, isLog, isNorm);
        }
      }
    }
  }
}

void Plotter::PlotHistosEvolution(bool isLog, bool isNorm) {
  for (auto &FileName: NameFiles) {
    for (auto &FoldName: RetrieveKeys(NameHistTypeHistMap)) {
      for (auto &VarName: TypeHistVarNameMap[NameHistTypeHistMap[FoldName]]) {
        PlotMaps(FileName, FoldName, VarName, isLog, isNorm);
      }
    }
  }
}



void Plotter::PlotMaps(TString FileName, TString FoldName, TString VarName, bool isLog, bool isNorm) {
  xmin = 10.;
  ymin = 1e-03;
  xmax = 0;
  ymax = 0;
  factor = isLog ? 1e03 : 2;

  for (auto &Tag : HistoTags) {
    if (isNorm)   MapHistos[FileName][Tag][FoldName][VarName]->Scale(1./  MapHistos[FileName][Tag][FoldName][VarName]->Integral());
    xmin = std::min(xmin,   MapHistos[FileName][Tag][FoldName][VarName]->GetXaxis()->GetXmin());
    ymin = std::min(ymin,   MapHistos[FileName][Tag][FoldName][VarName]->GetMinimum());
    xmax = std::max(xmax,   MapHistos[FileName][Tag][FoldName][VarName]->GetXaxis()->GetXmax());
    ymax = std::max(ymax,   MapHistos[FileName][Tag][FoldName][VarName]->GetMaximum());
  }

  if (ymax==0 || ymax<ymin) ymax = 1;
  if (ymin<=1e-10 ) ymin = 1e-10;


  TCanvas* c_histo = tdrCanvas(FoldName+VarName, xmin, xmax, ymin, factor*ymax, VarName, "A.U.");
  c_histo->SetLogy(isLog);
  TLegend *leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);

  for (auto &Tag : HistoTags) {
    SetLegend(leg, MapHistos[FileName][Tag][FoldName][VarName], Tag, colors[Tag], 0);
  }

  leg->SetLineColorAlpha(1, 0.7);
  leg->Draw("same");

  if (isSave) SaveCanvas(c_histo, OutDir+FileName+"/Evolution/"+FoldName+"/", VarName);

}


void Plotter::PlotSingleHisto(TString FileName, TString Tag, TString FoldName, TString VarName, bool isLog, bool isNorm) {
  SetPlotLimits({FoldName}, FileName, Tag, VarName, isLog);
  if (isNorm)   MapHistos[FileName][Tag][FoldName][VarName]->Scale(1./  MapHistos[FileName][Tag][FoldName][VarName]->Integral());


  TCanvas* c_histo = tdrCanvas(Tag+FoldName+VarName, xmin, xmax, ymin, factor*ymax, VarName, "A.U.");
  c_histo->SetLogy(isLog);
  TLegend *leg = tdrLeg(0.50,0.70,0.9,0.9, 0.035);

  for (auto &Tag : HistoTags) {
    tdrDraw(  MapHistos[FileName][Tag][FoldName][VarName], "same", kFullCircle, colors[Tag], 1, colors[Tag], 0, colors[Tag]);
    leg->AddEntry(  MapHistos[FileName][Tag][FoldName][VarName], Tag ,"lep");
  }
  leg->SetLineColorAlpha(1, 0.7);
  leg->Draw("same");

  if (isSave) SaveCanvas(c_histo, OutDir+FileName+"/All/"+Tag+"/"+FoldName+"/", VarName);

}



/*
██████  ██       ██████  ████████ ████████ ██ ███    ██  ██████      ██    ██  █████  ██████  ██  █████  ██████  ██      ███████ ███████
██   ██ ██      ██    ██    ██       ██    ██ ████   ██ ██           ██    ██ ██   ██ ██   ██ ██ ██   ██ ██   ██ ██      ██      ██
██████  ██      ██    ██    ██       ██    ██ ██ ██  ██ ██   ███     ██    ██ ███████ ██████  ██ ███████ ██████  ██      █████   ███████
██      ██      ██    ██    ██       ██    ██ ██  ██ ██ ██    ██      ██  ██  ██   ██ ██   ██ ██ ██   ██ ██   ██ ██      ██           ██
██      ███████  ██████     ██       ██    ██ ██   ████  ██████        ████   ██   ██ ██   ██ ██ ██   ██ ██████  ███████ ███████ ███████
*/




void Plotter::SaveCanvas(TCanvas* c, TString outdir, TString name) {
  gSystem->Exec("mkdir -p "+outdir);
  c->SaveAs(outdir+name+".pdf", "pdf");
}

void Plotter::SetPlotLimits(std::vector<TString> vecloop, TString FileName, TString Tag, TString VarName, bool isLog) {
  xmin = 10.;
  ymin = 1e-03;
  xmax = 0;
  ymax = 0;
  factor = isLog ? 1e03 : 2;

  for (auto &FoldName : vecloop) {
    xmin = std::min(xmin,   MapHistos[FileName][Tag][FoldName][VarName]->GetXaxis()->GetXmin());
    ymin = std::min(ymin,   MapHistos[FileName][Tag][FoldName][VarName]->GetMinimum());
    xmax = std::max(xmax,   MapHistos[FileName][Tag][FoldName][VarName]->GetXaxis()->GetXmax());
    ymax = std::max(ymax,   MapHistos[FileName][Tag][FoldName][VarName]->GetMaximum());
  }

  if (ymax==0 || ymax<ymin) ymax = 1;
  if (ymin<=1e-10 ) ymin = 1e-10;

}

void Plotter::SetLegend(TLegend* leg, TH1* h_, TString name, int color, int fill) {
  tdrDraw(h_, "same", kFullCircle, color, 1, color, fill, color);
  leg->AddEntry(h_, name ,"lep");
}




//
// void Plotter::LoadMap(TString FoldName, TString VarName){
//   for (auto &FileName: NameFiles) {
//     TH1::AddDirectory(kFALSE);
//     TFile *file = new TFile(InputDir+FileName+".root");
//     for (auto &Tag: HistoTags)  MapHistos[FileName][Tag][FoldName][VarName] = (TH1F*) ((TH1F*)file->Get(FoldName+Tag+"/"+VarName))->Clone();
//     file->Close();
//   }
// }
