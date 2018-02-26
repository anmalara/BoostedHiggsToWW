/* ************************************
*  *              Library               *
*   **************************************/

#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>  // needed to writw down on file
#include <string>

/* ************************************
*  *          General Variables         *
*   **************************************/

using namespace std;
void xmlWriterRecursive(string name_file, string path, int first, int last, bool add);
std::string itoa(int i);

/* ************************************
*  *           Main Program             *
*   **************************************/

void xmlfileWriter() {

  string name_file, path;
  int first, last;
  bool append;

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  for (int i = 0; i <= 8; i++) {
    name_file = "../config/xmlfile/SingleMu_C_2_1_18.xml";
    path = "/pnfs/desy.de/cms/tier2/store/user/raggleto/SingleMuon/crab_SingleMu_C_2_1_18/180202_202316/000"+itoa(i)+"/Ntuple_";

    if (i==0){ append = false; first = 1; last = 1000;}
    else {append = true; first = last; last = last+1000; }
    std::cout << first << " " << last << " " << path << '\n';
    xmlWriterRecursive(name_file, path, first, last, append);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////



  name_file = "../config/xmlfile/Run2017C-PromptReco-v1.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/areimers/NTuplesRun2/2017/DATA/SingleMuon/crab_Run2017C-PromptReco-v1/170809_094309/0000/Ntuple_";
  first = 1;
  last = 1000;
  append = false;
  // xmlWriterRecursive(name_file, path, first, last, append);

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  name_file = "../config/xmlfile/HZJ_HToWW_ZTo2L_M125_13TeV_powheg_pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/pgunnell/RunII_80X_v3/HZJ_HToWW_ZTo2L_M125_13TeV_powheg_pythia8/crab_HZJ_HToWW_ZTo2L_M125_13TeV_powheg_pythia8/180126_073414/0000/Ntuple_";
  first = 1;
  last = 1000;
  append = false;
  //xmlWriterRecursive(name_file, path, first, last, append);

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  name_file = "../config/xmlfile/HZ_HiggsToWWZToLL_Pythia8_TuneCP5.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/pgunnell/RunII_94X_v1/HiggsZProduction/crab_HZ_HiggsToWWZToLL_Pythia8_TuneCP5/180214_092805/0000/Ntuple_";
  first = 1;
  last = 1000;
  append = false;
  // xmlWriterRecursive(name_file, path, first, last, append);

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  name_file = "../config/xmlfile/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/akaravdi/RunII_80X_v3_Dep2016Campaign/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/170130_164539/0000/Ntuple_";
  first = 1;
  last = 1000;
  append = false;
  //xmlWriterRecursive(name_file, path, first, last, append);

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  name_file = "../config/xmlfile/WZ_TuneCUETP8M1_13TeV-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/isandova/NTuples_Moriond17/MC/WZ_TuneCUETP8M1_13TeV-pythia8/crab_MC_WZ/170115_174107/0000/Ntuple_";
  first = 1;
  last = 103;
  append = false;
  //xmlWriterRecursive(name_file, path, first, last, append);

  ////////////////////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////////////////////////

  name_file = "../config/xmlfile/ZZ_TuneCUETP8M1_13TeV-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/isandova/NTuples_Moriond17/MC/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_MC_ZZ/170115_174410/0000/Ntuple_";
  first = 1;
  last = 100;
  append = false;
  //xmlWriterRecursive(name_file, path, first, last, append);

  name_file = "../config/xmlfile/ZZ_TuneCUETP8M1_13TeV-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/isandova/NTuples_Moriond17/MC/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_MC_ZZ_v2/170126_214507/0000/Ntuple_";
  first = 1;
  last = 300;
  append = true;
  //xmlWriterRecursive(name_file, path, first, last, append);

  name_file = "../config/xmlfile/ZZ_TuneCUETP8M1_13TeV-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/isandova/NTuples_Moriond17/MC/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_MC_ZZ_v3/170129_041858/0000/Ntuple_";
  first = 1;
  last = 300;
  append = true;
  //xmlWriterRecursive(name_file, path, first, last, append);


  ////////////////////////////////////////////////////////////////////////////////////////////////////

  name_file = "../config/xmlfile/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/drberry/SFrameNtuples/RunII_80X_v3_Background/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170830_055134/0000/Ntuple_";
  first = 1;
  last = 53;
  append = false;
  //xmlWriterRecursive(name_file, path, first, last, append);

  name_file = "../config/xmlfile/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/drberry/SFrameNtuples/RunII_80X_v3_Background/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext2-v1/170830_055436/0000/Ntuple_";
  first = 1;
  last = 60;
  append = true;
  //xmlWriterRecursive(name_file, path, first, last, append);

  name_file = "../config/xmlfile/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/drberry/SFrameNtuples/RunII_80X_v3_Background/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1/170830_055305/0000/Ntuple_";
  first = 1;
  last = 60;
  append = true;
  //xmlWriterRecursive(name_file, path, first, last, append);

  name_file = "../config/xmlfile/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/drberry/SFrameNtuples/RunII_80X_v3_Background/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext5-v1/170830_055614/0000/Ntuple_";
  first = 1;
  last = 1000;
  append = true;
  //xmlWriterRecursive(name_file, path, first, last, append);

  name_file = "../config/xmlfile/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.xml";
  path = "/pnfs/desy.de/cms/tier2/store/user/drberry/SFrameNtuples/RunII_80X_v3_Background/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext5-v1/170830_055614/0001/Ntuple_";
  first = 1000;
  last = 2000;
  append = true;
  //xmlWriterRecursive(name_file, path, first, last, append);
  ////////////////////////////////////////////////////////////////////////////////////////////////////

}




/**********************************************************/

void xmlWriterRecursive(string name_file, string path, int first, int last, bool add=false) {
  char x[500];
  /******************* Write on txt *************************/

  std::ofstream myfile;
  if(!add) myfile.open(name_file);
  else myfile.open(name_file, ios::app);

  for (int i = first; i <= last; i++) {
    sprintf(x, "%s%d.root", path.c_str(), i);
    ifstream ifile(x);
    sprintf(x, "<In FileName=\"%s%d.root\" Lumi=\"0.0\"/>", path.c_str(), i);
    if(ifile) myfile << x << endl;
  }

  myfile.close();
}



std::string itoa(int i)
{
  char res[10];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}
