import os
import sys

sys.path.append("/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/")
sys.path.append("../python")
from parallelise import *

from ModuleRunnerBase import ModuleRunnerBase
from CreateConfigFiles import *
from functions import *

#This file contains all the classes and functions to clean up the steer.py code


class ModuleRunner(ModuleRunnerBase):
    def __init__(self):
        self.Base = ModuleRunnerBase()
        self.Module = ""
        self.ConfigFile = ""
        subprocess.call("COND", shell=True, executable="/bin/zsh")
        self.defineModules()
        self.CompileModules()

    def defineModules(self):
        self.ModuleSamples  = {"GenericCleaning"    : self.Base.Samples_original,
                               "GenLevelMatch"      : self.Base.Samples_original,
                               "FeasibilityStudy"   : self.Base.Samples_original,
                               "NeuralNetwork"      : self.Base.Samples}
        self.XMLFiles       = {}
        self.ModuleFiles    = {}
        self.ModuleStorage  = {}
        for mod in self.ModuleSamples:
            self.XMLFiles[mod]      = mod + "Config.xml"
            self.ModuleFiles[mod]   = mod + "Module.cxx"
            self.ModuleStorage[mod] = self.Base.Path_STORAGE+"Analysis/"+mod+"/"

    def CompileModules(self):
        cwd = os.getcwd()
        os.chdir(self.Base.Path_ANALYSIS);
        process = subprocess.Popen("make -j 20", shell=True)
        process.wait()
        os.chdir(cwd)

    def SetModule(self,module):
        self.Module     = module
        self.ModuleFile = self.ModuleFiles[module]
        self.ConfigFile = self.XMLFiles[module]
        self.Samples    = self.ModuleSamples[module]
        print "Module Name: \t", self.Module, "\nRunning \t", self.ConfigFile, "\nUsing \t \t", self.ModuleFile, "\nSamples:\t", self.Samples, "\n"

    def DeleteWorkdirs(self):
        delete_workdir(self.Base.SubmitDir+self.Module+"/")
        delete_workdir(self.Base.Path_SFRAME+self.Module+"*/")

    def CreateConfigFiles(self):
        a = os.system("rm -fr " +self.Base.SubmitDir+self.Module+"/")
        CreateConfigFiles(self.Samples, self.Base.Channels, self.Base.Path_ANALYSIS, self.Base.SubmitDir+self.Module+"/", self.ConfigFile)

    def CondorControl(self,option=""):
        condor_control(self.Base.SubmitDir+self.Module+"/", option)

    def RunLocal(self,option):
        local_run(self.Base.SubmitDir+self.Module+"/", option)

    def MergePreselection(self):
        list_processes = []
        for channel in self.Base.channels :
            for sample_sel in self.Base.samples_sel:
                command = ["hadd", "-f", self.Base.Path_ANALYSISPreSel+channel+"channel/"+self.Base.PREFIX+sample_sel+".root"]
                mode = "MC" if "MC" in sample_sel else "DATA"
                for sample_presel in self.Base.samples_dict[sample_sel]:
                    command.append(self.Base.Path_SFRAME+self.Module+"_"+mode+"/"+channel+"channel/"+self.Base.PREFIX+mode+"."+sample_presel+".root")
                    if len(self.Base.samples_dict[sample_sel])==1:
                        command = ["cp", command[3], command[2]]
                list_processes.append(command)
        parallelise(list_processes, 16)
