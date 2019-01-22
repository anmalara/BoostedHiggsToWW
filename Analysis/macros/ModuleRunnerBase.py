class ModuleRunnerBase:
    def __init__(self):
        self.defineDirectories()
        self.defineSamples()
        self.PrefixrootFile     = "uhh2.AnalysisModuleRunner."

    def defineDirectories(self):
        self.Path_NFS           = "/nfs/dust/cms/user/amalara/"
        self.Path_SFRAME        = self.Path_NFS+"sframe_all/"
        self.Path_ANALYSIS      = self.Path_NFS+"WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/"
        self.Path_STORAGE       = self.Path_NFS+"WorkingArea/File/"
        self.ConfigDir          = self.Path_ANALYSIS+"config/"
        self.SubmitDirName      = "SubmittedJobs/"
        self.SubmitDir          = self.ConfigDir+self.SubmitDirName

    def defineSamples(self):
        self.Channels = ["muon", "electron"]
        self.Samples_dict = {"MC_DYJets" : ["MC_DY1JetsToLL","MC_DY2JetsToLL","MC_DY3JetsToLL","MC_DY4JetsToLL"],
                             "MC_TTbar"  : ["MC_TTbarSemiLeptonic", "MC_TTTo2L2Nu", "MC_TTToHadronic"],
                             "MC_WZ"     : ["MC_WZ"],
                             "MC_ZZ"     : ["MC_ZZ"],
                             "MC_HZ"     : ["MC_HZ_HiggsToWWZToLL"],
                             "DATA"      : ["Data_2017C"]}
        self.Samples_original = []
        self.Samples = []
        for key in self.Samples_dict:
            self.Samples_original.extend(self.Samples_dict[key])
            self.Samples.append(key)
