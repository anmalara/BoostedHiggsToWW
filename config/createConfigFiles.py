import sys
sys.path.append("/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/")
from fileManipulation import *

@timeit
def createConfigFiles(samples=["MC_TT", "Data_2017C"], channels= ["muon"], original_dir = "./submittedJobs/", original_file = "FeasibilityStudy.xml"):
    outdir = original_file[:len(original_file)-4]
    for channel in channels:
        path = original_dir+channel+"channel/"
        if not os.path.exists(path):
            os.makedirs(path)
        for process in samples:
            filename = original_file[:len(original_file)-4]+"_"+process+".xml"
            cmd = "cp %s %s" % (original_file, path+filename)
            a = os.system(cmd)
            cmd = "cp %s %s" % ("JobConfig.dtd", path)
            a = os.system(cmd)
            comments = []
            for el in list(set(samples) - set([process])):
                if "MC" in el:
                    comments.append(["<InputData", "Type", "MC", '"'+el+'"'])
                elif "DATA".lower() in el.lower():
                    comments.append(["<InputData", "Type", "DATA", '"'+el+'"'])
            comment_lines(path, filename, comments, remove=True)
            if "MC" in process:
                sample = "MC"
            elif "DATA".lower() in process.lower():
                sample = "DATA"
            changes = []
            changes.append(["<ConfigSGE", "Workdir", "workdir_"+outdir, "workdir_"+outdir+"_"+process])
            changes.append(["<!ENTITY", "OUTDIR", outdir , outdir+"_"+sample+"/"+channel+"channel"])
            if "muon".lower() in channel.lower():
                changes.append(["<!ENTITY", "electronchannel", '"true"', '"false"'])
            if "electron".lower() in channel.lower():
                changes.append(["<!ENTITY", "muonchannel", '"true"', '"false"'])
            if process in ["MC_WZ", "MC_ZZ", "MC_HZJ_HToWW_ZTo2L"]:
                changes.append(["<Item", "TopJetCollection", '"packedPatJetsAk8CHSJets_SoftDropCHS"', '"slimmedJetsAK8_SoftDrop"'])
            if "DATA".lower() in process.lower():
                changes.append(["<Item", "readTrigger", '"false"', '"true"'])
            change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
