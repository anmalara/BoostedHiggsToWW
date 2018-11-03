#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/tdrstyle_all.C"

void plot(string inputdir = "/nfs/dust/cms/user/amalara/sframe_all/FeasibilityStudy_MC/", string outdir = "./") {

  TCanvas* c = tdrCanvas("SDmass", 0, 200, 1e-10, 1E2, "SDmass_jet1" , "Events");
  c->SetLogy(1);

  TString text;
  TLegend* leg = tdrLeg(0.3,0.65,0.9,0.9); leg->SetTextFont(42);  leg->SetTextSize(0.035);

  double f;

  #define LOADBKG(bkg,color,xsec,bkg_)\
  string file_name_##bkg = inputdir+"uhh2.AnalysisModuleRunner.MC.MC_"; file_name_##bkg += #bkg_; file_name_##bkg += ".root";\
  TFile *file_##bkg = new TFile((file_name_##bkg).c_str());\
  TH1F* histo_SD2_##bkg##_tot = (TH1F*)file_##bkg  ->Get("nJet_nocuts/number");\
  TH1F* histo_SD1_##bkg##_tot = (TH1F*)file_##bkg  ->Get("nJet_nocuts/SDmass_jet1");\
  TH1F* histo_SD2_##bkg = (TH1F*)file_##bkg  ->Get("nJet_JetElectron/number");\
  TH1F* histo_SD1_##bkg = (TH1F*)file_##bkg  ->Get("nJet_JetElectron/SDmass_jet1");\
  std::cout << "new " <<histo_SD2_##bkg##_tot->Integral() << " " << histo_SD2_##bkg##_tot->GetEntries() << " " << histo_SD2_##bkg##_tot->GetSumOfWeights() << " " << '\n';\
  std::cout << histo_SD1_##bkg##_tot->Integral() << " " << histo_SD1_##bkg##_tot->GetEntries() << " " << histo_SD1_##bkg##_tot->GetSumOfWeights() << " " << '\n';\
  std::cout << histo_SD2_##bkg->Integral() << " " << histo_SD2_##bkg->GetEntries() << " " << histo_SD2_##bkg->GetSumOfWeights() << " " << '\n';\
  std::cout << histo_SD1_##bkg->Integral() << " " << histo_SD1_##bkg->GetEntries() << " " << histo_SD1_##bkg->GetSumOfWeights() << " " << '\n';\
  tdrDraw(histo_SD1_##bkg, "hist p", kFullCircle, color, kSolid, color, 0);\
  text = #bkg; text += ": "; text += histo_SD1_##bkg->GetEntries();  text += " of "; text += histo_SD1_##bkg##_tot->GetEntries(); text += " eff: "; text += round(100.*histo_SD1_##bkg->GetEntries()/histo_SD1_##bkg##_tot->GetEntries()); text += "\%"; \
  leg->AddEntry(histo_SD1_##bkg, text, "p");\


  #define LOADBKG1(bkg,color,xsec,bkg_)\
  string file_name_##bkg = inputdir+"uhh2.AnalysisModuleRunner.MC.MC_"; file_name_##bkg += #bkg_; file_name_##bkg += ".root";\
  TFile *file_##bkg = new TFile((file_name_##bkg).c_str());\
  TH1F* histo_SD2_##bkg##_tot = (TH1F*)file_##bkg  ->Get("nJet_nocuts/weights");\
  TH1F* histo_SD1_##bkg##_tot = (TH1F*)file_##bkg  ->Get("nJet_nocuts/SDmass_jet1");\
  TH1F* histo_SD2_##bkg = (TH1F*)file_##bkg  ->Get("nJet_JetElectron/number");\
  TH1F* histo_SD1_##bkg = (TH1F*)file_##bkg  ->Get("nJet_JetElectron/SDmass_jet1");\
  f = xsec*histo_SD1_##bkg##_tot->GetSumOfWeights()/histo_SD2_##bkg##_tot->GetSumOfWeights();\
  std::cout << f << " " << xsec << " " << histo_SD1_##bkg##_tot->GetSumOfWeights() << " " << histo_SD1_##bkg->GetSumOfWeights() << '\n';\
  histo_SD1_##bkg->Scale(f);\
  tdrDraw(histo_SD1_##bkg, "hist p", kFullCircle, color, kSolid, color, 0);\
  text = #bkg; text += ": "; text += histo_SD1_##bkg->GetEntries();  text += " of "; text += histo_SD1_##bkg##_tot->GetEntries(); text += " eff: "; text += round(100.*histo_SD1_##bkg->GetEntries()/histo_SD1_##bkg##_tot->GetEntries()); text += "\%"; \
  leg->AddEntry(histo_SD1_##bkg, text, "p");\


  #define LOADBKG2(bkg,color,xsec,bkg_)\
  string file_name_##bkg = inputdir+"uhh2.AnalysisModuleRunner.MC.MC_"; file_name_##bkg += #bkg_; file_name_##bkg += ".root";\
  TFile *file_##bkg = new TFile((file_name_##bkg).c_str());\
  TH1F* histo_SD1_##bkg##_tot = (TH1F*)file_##bkg  ->Get("nJet_nocuts/number");\
  TH1F* histo_SD1_##bkg = (TH1F*)file_##bkg  ->Get("nJet_JetElectron/SDmass_jet1");\
  histo_SD1_##bkg->Scale(1./(xsec*xsec*histo_SD1_##bkg##_tot->GetSumOfWeights()));\
  tdrDraw(histo_SD1_##bkg, "hist p", kFullCircle, color, kSolid, color, 0);\
  text = #bkg; text += ": "; text += histo_SD1_##bkg->GetEntries();  text += " of "; text += histo_SD1_##bkg##_tot->GetEntries(); text += " eff: "; text += round(100.*histo_SD1_##bkg->GetEntries()/histo_SD1_##bkg##_tot->GetEntries()); text += "\%"; \
  leg->AddEntry(histo_SD1_##bkg, text, "p");\

  // LOADBKG(HZJ_HToWW_ZTo2L,kRed+1,0.00022*1000)
  // LOADBKG(HZ_HiggsToWWZToLL,kRed+1,1,HZ_HiggsToWWZToLL)
  LOADBKG(DY1JetsToLL,kGreen-2,1,DY1JetsToLL)
  // LOADBKG(DY2JetsToLL,kGreen-2,1,DY2JetsToLL)
  // LOADBKG(DY3JetsToLL,kGreen-2,1,DY3JetsToLL)
  // LOADBKG(DY4JetsToLL,kGreen-2,1,DY4JetsToLL)
  // LOADBKG(TTToHadronic,kMagenta,1,TTToHadronic)
  // LOADBKG(TTbarSemiLeptonic,kViolet-3,1,TTbarSemiLeptonic)
  // LOADBKG(TTTo2L2Nu,kViolet,1,TTTo2L2Nu)
  // LOADBKG(WZ,kOrange,1,WZ)
  // LOADBKG(ZZ,kBlue-4,1,ZZ)

  // LOADBKG1(HZ_HiggsToWWZToLL1,kRed,0.00022*1000,HZ_HiggsToWWZToLL)
  // LOADBKG1(HZ_HiggsToWWZToLL1,kRed,1,HZ_HiggsToWWZToLL)

  LOADBKG1(DY1JetsToLL1,kRed,1./37342.76571,DY1JetsToLL)
  // LOADBKG1(TTToHadronic1,kRed,1./100683950.1,TTToHadronic)
  // LOADBKG1(TTbarSemiLeptonic1,kRed,1./91335604.46,TTbarSemiLeptonic)
  // LOADBKG1(TTTo2L2Nu1,kRed,1./56758143.39,TTTo2L2Nu)
  // LOADBKG1(WZ1,kRed,1./21005.792,WZ)
  // LOADBKG1(ZZ1,kRed,1./60107.365,ZZ)

  //
  // LOADBKG2(HZJ_HToWW_ZTo2L2,kRed,0.00022,HZJ_HToWW_ZTo2L)

  LOADBKG2(DY1JetsToLL2,kBlue,1./37342.76571,DY1JetsToLL)
  // LOADBKG2(TTToHadronic2,kBlue,377.96,TTToHadronic)
  // LOADBKG2(TTbarSemiLeptonic2,kBlue,365.34,TTbarSemiLeptonic)
  // LOADBKG2(TTTo2L2Nu2,kBlue,88.29,TTTo2L2Nu)
  // LOADBKG2(WZ2,kBlue,1./21005.792,WZ)
  // LOADBKG2(ZZ2,kBlue,1./60107.365,ZZ)




  TH1F* histo_SD1_check = (TH1F*)histo_SD1_DY1JetsToLL->Clone();
  histo_SD1_check->Divide(histo_SD1_DY1JetsToLL1);
  new TCanvas;
  tdrDraw(histo_SD1_check, "hist p", kFullCircle, kBlue, kSolid, kBlue, 0);
  //
  // leg->Draw();

  return true;


}
