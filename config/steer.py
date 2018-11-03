from createConfigFiles import *
from glob import glob

def cont_event(original_dir ="/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/config/SubmittedJobs/"):
  count = 0
  for sample in glob(original_dir+"*"):
    if not ".xml" in sample:
      continue
    count += 1
  return count

@timeit
def condor_control(original_dir ="/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/config/SubmittedJobs/", option=""):
    internal_option= {"list": "-l", "submit": "-s", "resubmit": "-r", "add": "-f", "merge": "-f"}
    count = 0
    all_events = cont_event(original_dir)
    for el in glob(original_dir+"*"):
        if not ".xml" in el:
            continue
        path   = el[:el.rfind("/")+1 ]
        sample = el[ el.rfind("/")+1:]
        count += 1
        os.chdir(path)
        if option in internal_option:
            command = ['sframe_batch.py', internal_option[option], sample]
        else:
            command = ['sframe_batch.py', sample]
        process = subprocess.Popen(command)
        process.wait()
        print "Already completed "+str(count)+" out of "+str(all_events)+" jobs --> "+str(round(float(count)/float(all_events)*100,2))+"%."
        if "submit" in option:
            time.sleep(10)
        os.chdir(original_dir)


@timeit
def delete_workdir(original_dir ="/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/BoostedHiggsToWW/config/SubmittedJobs/"):
  for el in glob(original_dir+"*"):
      if os.path.isdir(el) and "workdir" in el:
          cmd = "rm -fr %s" % (el)
          a = os.system(cmd)



##################################################
#                                                #
#                   MAIN Program                 #
#                                                #
##################################################

def main_program(option="", samples=[], original_dir="./SubmittedJobs/", original_file="FeasibilityStudy.xml"):
    if option == "new":
      cmd = "rm -fr %s" % (original_dir)
      a = os.system(cmd)
      createConfigFiles(samples, original_dir, original_file)
    elif option == "remove" or option == "delete":
      delete_workdir(original_dir)
    else:
      condor_control(original_dir, option)


QCD_samples= []
QCD_samples.append("MC_DY1JetsToLL")
QCD_samples.append("MC_DY2JetsToLL")
QCD_samples.append("MC_DY3JetsToLL")
QCD_samples.append("MC_DY4JetsToLL")
QCD_samples.append("MC_TTbarSemiLeptonic")
QCD_samples.append("MC_TTTo2L2Nu")
QCD_samples.append("MC_TTToHadronic")
QCD_samples.append("MC_WZ")
QCD_samples.append("MC_ZZ")
QCD_samples.append("MC_HZ_HiggsToWWZToLL")

Data_samples= []
Data_samples.append("Data_2017C")

samples = QCD_samples+Data_samples

original_file = "FeasibilityStudy.xml"
original_dir_ = os.getcwd()

original_dir = original_dir_
original_dir += "/SubmittedJobs/"

try:
  option = sys.argv[1]
except:
  option = ""


main_program(option, samples, original_dir, original_file)
