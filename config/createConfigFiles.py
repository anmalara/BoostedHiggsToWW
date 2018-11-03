import sys
sys.path.append("/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/PersonalCode/")
from fileManipulation import *

@timeit
def createConfigFiles(processes=["MC_TT", "Data_2017C"], original_dir = "./submittedJobs/", original_file = "FeasibilityStudy.xml"):
  outdir = original_file[:len(original_file)-4]
  path = original_dir
  if not os.path.exists(path):
    os.makedirs(path)
  for process in processes:
    filename = original_file[:len(original_file)-4]+"_"+process+".xml"
    cmd = "cp %s %s" % (original_file, path+filename)
    a = os.system(cmd)
    cmd = "cp %s %s" % ("JobConfig.dtd", path)
    a = os.system(cmd)
    controls = []
    for el in list(set(processes) - set([process])):
      if "MC" in el:
        controls.append(["<InputData", "Type", "MC", '"'+el+'"'])
      elif "DATA".lower() in el.lower():
        controls.append(["<InputData", "Type", "DATA", '"'+el+'"'])
    comment_lines(path, filename, controls, remove=True)
    if "MC" in process:
      sample = "MC"
    elif "DATA".lower() in process.lower():
      sample = "DATA"
    controls = []
    controls.append(["<ConfigSGE", "Workdir", "workdir_"+outdir, "workdir_"+outdir+"_"+process])
    controls.append(["<!ENTITY", "OUTDIR", outdir , outdir+"_"+sample])
    if process in ["MC_WZ", "MC_ZZ", "MC_HZJ_HToWW_ZTo2L"]:
      controls.append(["<Item", "TopJetCollection", '"packedPatJetsAk8CHSJets_SoftDropCHS"', '"slimmedJetsAK8_SoftDrop"'])
    if "DATA".lower() in process.lower():
      controls.append(["<Item", "readTrigger", '"false"', '"true"'])
    change_lines(path, filename, [el[0:2] for el in controls ], [el[2:3] for el in controls ], [el[3:4] for el in controls ])
