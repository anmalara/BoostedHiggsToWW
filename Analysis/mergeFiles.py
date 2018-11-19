import sys
from glob import glob
sys.path.append("/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/")
from parallelise import *


def MergeFiles(list_processes, srcfiles, inputDir, outputDir):
    for key in srcfiles:
        command = ["hadd", "-f", outputDir+"/"+key]
        for el in glob(inputDir+"/"+srcfiles[key]):
            command.append(el)
        list_processes.append(command)


inputDir  = "/nfs/dust/cms/user/amalara/sframe_all/"
outputDir = "/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/Preselection/"
channels = ["electronchannel","muonchannel"]
Data_MC = {"MC": "FeasibilityStudy_MC", "DATA": "FeasibilityStudy_DATA",}

prename = "uhh2.AnalysisModuleRunner"

srcfiles_MC   = {prename+".MC.HZ.root"       : prename+".MC.MC_HZ_HiggsToWWZToLL.root",
                 prename+".MC.DYJets.root"   : prename+".MC.MC_DY*.root",
                 prename+".MC.TTbar.root"    : prename+".MC.MC_TT*.root",
                 prename+".MC.WZ.root"       : prename+".MC.MC_WZ*.root",
                 prename+".MC.ZZ.root"       : prename+".MC.MC_ZZ*.root",};

srcfiles_DATA = {prename+".DATA.DATA.root" : prename+".DATA.Data_*.root",}

srcfiles = {"MC": srcfiles_MC, "DATA": srcfiles_DATA}


list_processes = []
channel = "electronchannel"
for channel in channels:
    for key in Data_MC:
        MergeFiles(list_processes, srcfiles[key], inputDir+Data_MC[key]+"/"+channel, outputDir+"/"+channel)

parallelise(list_processes, 10)
