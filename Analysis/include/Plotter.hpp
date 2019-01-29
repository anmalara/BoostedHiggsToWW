#include <vector>
#include <map>

#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/tdrstyle_all.C"
#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/Analysis/include/Tools.h"


class Plotter {
public:

  /** \brief Constructor and Settings
  *
  * Create a Plotter Class that handle the in-output paths,
  * as well as file/histos/folder names.
  * The idea of this Class is to load all histos saved in root files,
  * creating a map to re-access arterwords, allowing whatsoever loop you might want to.
  * Name of plots expected for each Class of histo folders is hard coded below.
  * There is a commented script to print iteratively the names for a give input example.
  * Color code is also hard coded in this header.
  *
  */

  Plotter(TString path, TString outpath="../plots/temp/"): InputDir(path), OutDir(outpath) { MakeMapTypeHist(); gSystem->Exec("mkdir -p "+OutDir); }
  void MakeMapTypeHist();
  void SetNameFiles(std::vector<TString> v)   { NameFiles.clear(); std::swap(v, NameFiles);}
  void SetNameVars(std::vector<TString> v)    { NameVars.clear(); std::swap(v, NameVars);}
  void SetNameFolders(std::vector<TString> v) { NameFolders.clear(); std::swap(v, NameFolders);}
  void SetHistoTags(std::vector<TString> v)   { HistoTags.clear(); std::swap(v, HistoTags);}
  void SetisSave(bool set) {isSave=set;};

  /** \brief Load and Print stuff
  * Those functions has to be called setting HistoTags and NameFiles.
  * This is to limit time and memory load process.
  */

  void LoadMap();
  void PrintMap(bool debug=0);

  /**
  * Generic Plotting Rules
  */

  void SaveCanvas(TCanvas* c, TString outdir, TString name);
  void SetPlotLimits(std::vector<TString> vecloop, TString FileName, TString Tag, TString VarName, bool isLog);
  void SetLegend(TLegend* leg, TH1* h_, TString name, int color, int fill);

  /** \brief Plot Macros
  * Create and add your own Plotting loop functions
  */

  void PlotHistosComposition(TString FileName, TString FoldName_, bool isLog);
  void PlotSingleHistoComposition(TString FileName, TString Tag, TString FoldName_, TString VarName, bool isLog);
  void SortByFolderName(std::vector<TString>& vect, TString FileName, TString Tag, TString VarName);

  void PlotHistosEvolution(bool isLog = 0, bool isNorm=0);
  void PlotMaps(TString FileName, TString FoldName, TString VarName, bool isLog = 0, bool isNorm=1);

  void PlotMaps(bool isLog = 0, bool isNorm=0);
  void PlotAllHistos(bool isLog = 0, bool isNorm=0);
  void PlotSingleHisto(TString FileName, TString Tag, TString FoldName, TString VarName, bool isLog=0, bool isNorm=1);
  // void LoadMap(TString foldName, TString var);
protected:
  // General Variables
  TString InputDir, OutDir;
  std::vector<TString> NameFiles;
  std::vector<TString> NameFolders;
  std::vector<TString> NameVars;
  std::vector<TString> HistoTags;

  // Internal Maps
  std::map<TString, int> colors;
  std::map<TString, TString > NameHistTypeHistMap;
  std::map<TString, std::vector<TString > > TypeHistVarNameMap;
  std::map<TString, std::map<TString, std::map<TString, std::map<TString, TH1* > > > > MapHistos;

  // Plot variables
  bool isSave = true;
  double xmin;
  double ymin;
  double xmax;
  double ymax;
  double factor;

  // ClassDef(Plotter, 0);

};



/** \brief Internal Map Definitions
*
* This function is Called at Constructor Time
* Simply create maps for Histo names and colors.
*
*/

void Plotter::MakeMapTypeHist() {

  NameHistTypeHistMap["gen"] = "GenMatchHists";
  NameHistTypeHistMap["nTopJet_ZMatch"] = "JetHistosWithConditions";
  NameHistTypeHistMap["nTopJet_HMatch"] = "JetHistosWithConditions";
  NameHistTypeHistMap["nTopJet_WWfullHad"] = "JetHistosWithConditions";
  NameHistTypeHistMap["nTopJet_WWfullLep"] = "JetHistosWithConditions";
  NameHistTypeHistMap["nTopJet_WWsemiLep"] = "JetHistosWithConditions";
  NameHistTypeHistMap["nTopJet"] = "BoostedHiggsToWWTopJetHists";
  NameHistTypeHistMap["nJet"] = "BoostedHiggsToWWJetHists";
  NameHistTypeHistMap["ele"] = "ElectronHists";
  NameHistTypeHistMap["muon"] = "MuonHists";
  NameHistTypeHistMap["event"] = "EventHists";
  NameHistTypeHistMap["diLepton"] = "DiLeptonHists";

  TypeHistVarNameMap["GenMatchHists"] = {"N_Z", "pt_Z", "eta_Z", "phi_Z", "mass_Z", "flavor_Z",
  "N_mother1Z", "pt_mother1Z", "eta_mother1Z", "phi_mother1Z", "mass_mother1Z", "flavor_mother1Z",
  "N_mother2Z", "pt_mother2Z", "eta_mother2Z", "phi_mother2Z", "mass_mother2Z", "flavor_mother2Z",
  "N_daughter1Z", "pt_daughter1Z", "eta_daughter1Z", "phi_daughter1Z", "mass_daughter1Z", "flavor_daughter1Z",
  "N_daughter2Z", "pt_daughter2Z", "eta_daughter2Z", "phi_daughter2Z", "mass_daughter2Z", "flavor_daughter2Z",
  "deltaR_daughtersZ", "deltaR_daughter1Z", "deltaR_daughter2Z",
  "N_Higgs", "pt_Higgs", "eta_Higgs", "phi_Higgs", "mass_Higgs", "flavor_Higgs",
  "N_mother1Higgs", "pt_mother1Higgs", "eta_mother1Higgs", "phi_mother1Higgs", "mass_mother1Higgs", "flavor_mother1Higgs",
  "N_mother2Higgs", "pt_mother2Higgs", "eta_mother2Higgs", "phi_mother2Higgs", "mass_mother2Higgs", "flavor_mother2Higgs",
  "N_daughter1Higgs", "pt_daughter1Higgs", "eta_daughter1Higgs", "phi_daughter1Higgs", "mass_daughter1Higgs", "flavor_daughter1Higgs",
  "N_daughter2Higgs", "pt_daughter2Higgs", "eta_daughter2Higgs", "phi_daughter2Higgs", "mass_daughter2Higgs", "flavor_daughter2Higgs",
  "deltaR_daughtersHiggs", "deltaR_daughter1Higgs", "deltaR_daughter2Higgs",
  "N_Wplus", "pt_Wplus", "eta_Wplus", "phi_Wplus", "mass_Wplus", "flavor_Wplus",
  "N_mother1Wplus", "pt_mother1Wplus", "eta_mother1Wplus", "phi_mother1Wplus", "mass_mother1Wplus", "flavor_mother1Wplus",
  "N_mother2Wplus", "pt_mother2Wplus", "eta_mother2Wplus", "phi_mother2Wplus", "mass_mother2Wplus", "flavor_mother2Wplus",
  "N_daughter1Wplus", "pt_daughter1Wplus", "eta_daughter1Wplus", "phi_daughter1Wplus", "mass_daughter1Wplus", "flavor_daughter1Wplus",
  "N_daughter2Wplus", "pt_daughter2Wplus", "eta_daughter2Wplus", "phi_daughter2Wplus", "mass_daughter2Wplus", "flavor_daughter2Wplus",
  "deltaR_daughtersWplus", "deltaR_daughter1Wplus", "deltaR_daughter2Wplus",
  "N_Wminus", "pt_Wminus", "eta_Wminus", "phi_Wminus", "mass_Wminus", "flavor_Wminus",
  "N_mother1Wminus", "pt_mother1Wminus", "eta_mother1Wminus", "phi_mother1Wminus", "mass_mother1Wminus", "flavor_mother1Wminus",
  "N_mother2Wminus", "pt_mother2Wminus", "eta_mother2Wminus", "phi_mother2Wminus", "mass_mother2Wminus", "flavor_mother2Wminus",
  "N_daughter1Wminus", "pt_daughter1Wminus", "eta_daughter1Wminus", "phi_daughter1Wminus", "mass_daughter1Wminus", "flavor_daughter1Wminus",
  "N_daughter2Wminus", "pt_daughter2Wminus", "eta_daughter2Wminus", "phi_daughter2Wminus", "mass_daughter2Wminus", "flavor_daughter2Wminus",
  "deltaR_daughtersWminus", "deltaR_daughter1Wminus", "deltaR_daughter2Wminus"};

  TypeHistVarNameMap["BoostedHiggsToWWJetHists"] = {"number", "weights", "mass_jet", "mT_jet", "pt_jet", "eta_jet", "phi_jet", "csv_jet", "flavor_jet", "deltaR_lepton_jet","deltaRmin_1", "deltaRmin_2",
  "mass_1", "mT_1", "pt_1", "eta_1", "phi_1", "csv_1", "flavor_1", "deltaR_lepton_1",
  "mass_2", "mT_2", "pt_2", "eta_2", "phi_2", "csv_2", "flavor_2", "deltaR_lepton_2",
  "mass_3", "mT_3", "pt_3", "eta_3", "phi_3", "csv_3", "flavor_3", "deltaR_lepton_3",
  "mass_4", "mT_4", "pt_4", "eta_4", "phi_4", "csv_4", "flavor_4", "deltaR_lepton_4"};

  TypeHistVarNameMap["BoostedHiggsToWWTopJetHists"] = {"number", "weights", "mass_jet", "mT_jet", "pt_jet", "eta_jet", "phi_jet", "csv_jet", "flavor_jet",
  "SDmass_jet", "tau1_jet", "tau2_jet", "tau3_jet", "tau21_jet", "tau31_jet", "tau32_jet", "nsubjet_jet", "deltaR_lepton_jet", "deltaRmin_1", "deltaRmin_2",
  "mass_1", "mT_1", "pt_1", "eta_1", "phi_1", "csv_1", "flavor_1", "deltaR_lepton_1", "SDmass_1", "tau1_1", "tau2_1", "tau3_1", "tau21_1", "tau31_1", "tau32_1", "nsubjet_1",
  "mass_2", "mT_2", "pt_2", "eta_2", "phi_2", "csv_2", "flavor_2", "deltaR_lepton_2", "SDmass_2", "tau1_2", "tau2_2", "tau3_2", "tau21_2", "tau31_2", "tau32_2", "nsubjet_2",
  "mass_3", "mT_3", "pt_3", "eta_3", "phi_3", "csv_3", "flavor_3", "deltaR_lepton_3", "SDmass_3", "tau1_3", "tau2_3", "tau3_3", "tau21_3", "tau31_3", "tau32_3", "nsubjet_3",
  "mass_4", "mT_4", "pt_4", "eta_4", "phi_4", "csv_4", "flavor_4", "deltaR_lepton_4", "SDmass_4", "tau1_4", "tau2_4", "tau3_4", "tau21_4", "tau31_4", "tau32_4", "nsubjet_4"};

  TypeHistVarNameMap["ElectronHists"] = {"number", "pt", "eta", "phi", "isolation", "charge", "ptrel", "deltaRmin", "deltaRmin_ptrel", "eff_sub", "eff_tot", "pt_response",
  "pt_1", "eta_1", "phi_1", "isolation_1", "charge_1", "ptrel_1", "deltaRmin_1", "deltaRmin_ptrel_1",
  "pt_2", "eta_2", "phi_2", "isolation_2", "ptrel_2", "charge_2", "deltaRmin_2", "deltaRmin_ptrel_2"};

  TypeHistVarNameMap["MuonHists"] = TypeHistVarNameMap["ElectronHists"];

  TypeHistVarNameMap["EventHists"] = {"N_PrimVertices", "N_TrueInteractions", "Weights", "MET", "HT", "HTLep", "ST", "WeightsLogBins"};

  TypeHistVarNameMap["DiLeptonHists"] = {"number", "diMuon_number", "diMuon_charge", "diMuon_m", "diMuon_pt", "diMuon_eta", "diMuon_phi", "diMuon_DR12", "diMuon_jet_Dphi", "pt_Muon_2D",
  "diElectron_number", "diElectron_charge", "diElectron_m", "diElectron_pt", "diElectron_eta", "diElectron_phi", "diElectron_DR12", "diElectron_jet_Dphi", "pt_Electron_2D"};

  TypeHistVarNameMap["JetHistosWithConditions"] = TypeHistVarNameMap["BoostedHiggsToWWTopJetHists"];

  colors["nocuts"] = kRed+1;
  colors["IDs"] = kGreen-2;
  colors["metfilters"] = kBlue;
  colors["jetlep"] = kRed;
  colors["JEC"] = kGreen;
  colors["JER"] = kViolet-3;
  colors["MET"] = kRed-1;
  colors["jetIDnoboost"] = kOrange+1;
  colors["jetID"] = kSpring-9;
  colors["cleaned"] = kSpring-7;
  colors["Veto"] = kSpring+10;
  colors["NBoostedJet"] = kGreen+1;
  colors["NLeptonSel"] = kOrange-1;
  colors["JetDiLeptonPhiAngular"] = kOrange;

  colors["origin_DY4"] = kGreen+3;
  colors["origin_TT"] = kAzure-4;
  colors["origin_TTBarhad"] = kAzure-4;
  colors["origin_TTBarlep"] = kCyan;
  colors["origin_TTBarsemilep"] = kAzure+6;
  colors["nTopJet"] = kBlack;
  colors["nTopJet_ZMatch"] = kRed+1;
  colors["nTopJet_HMatch"] = kBlue;
  colors["nTopJet_WWsemiLep"] = kGreen+3;
  colors["nTopJet_WWfullLep"] = kOrange;
  colors["nTopJet_WWfullHad"] = kViolet-3;

  //
  // TFile *file = new TFile("/nfs/dust/cms/user/amalara/sframe_all/GenericCleaning_MC/muonchannel/uhh2.AnalysisModuleRunner.MC.MC_HZ_HiggsToWWZToLL.root");
  // TIterator *iter = (file->GetListOfKeys())->MakeIterator();
  // TObject *key = iter->Next();
  // while (key) {
  //   if ( "SFrame" == (string)key->GetName() || "uhh2_meta" == (string)key->GetName() || "AnalysisTree" == (string)key->GetName() ) break;
  //   if (key->IsFolder()) {
  //     std::vector<TString> Name; TString tok; int from = 0;
  //     while (((TString)key->GetName()).Tokenize(tok, from, "_")) Name.push_back(tok);
  //     if (Name.size()==2) Name.insert(Name.begin()+1, "nominal");
  //     else if (Name.size()!=3) throw std::invalid_argument( "Case not implemented" ) ;
  //     for (size_t i = 0; i < Name.size(); i++) std::cout << Name[i] << '\t' << '\t'; std::cout << '\n';
  //     TIterator *iter_1 = ((file->GetDirectory(key->GetName()))->GetListOfKeys())->MakeIterator();
  //     TObject *key_1 = iter_1->Next();
  //     while (key_1) { std::cout << '"' << key_1->GetName() << '"' << ", "; key_1 = iter_1->Next(); } std::cout << '\n' << '\n';
  //   }
  //   key = iter->Next();
  // }

}
