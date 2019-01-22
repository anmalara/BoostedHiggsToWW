import sys
sys.path.append("/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/")
from fileManipulation import *

@timeit
def CreateConfigFiles(samples, channels, original_dir, SubmitDir,  original_file):
    outdir = original_file.replace("Config.xml", "")
    for channel in channels:
        path = SubmitDir+channel+"channel/"
        if not os.path.exists(path):
            os.makedirs(path)
        for process in samples:
            filename = outdir+"_"+process+".xml"
            a = os.system("cp "+original_dir+"config/"+original_file+" "+path+filename)
            a = os.system("cp "+original_dir+"JobConfig.dtd "+path)
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
            changes.append(["<!ENTITY", "CHANNEL", '"channel"' , '"'+channel+"channel"+'"'])
            if "muon".lower() in channel.lower():
                changes.append(["<!ENTITY", "electronchannel", '"true"', '"false"'])
            if "electron".lower() in channel.lower():
                changes.append(["<!ENTITY", "muonchannel", '"true"', '"false"'])
            if "DATA".lower() in process.lower():
                changes.append(["<Item", "readTrigger", '"false"', '"true"'])
            change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
