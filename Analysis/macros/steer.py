#!/usr/bin/env python
from ModuleRunner import *

"""
This macro and guide describes and steers the complete progress of the BoostedHiggsToWW analysis.
Modify these items, everything else will work. The files given here must already be adapted to new settings in case there are any.
"""


Modules = ModuleRunner()

#################################################
#                                               #
#                GenericCleaning                #
#                                               #
#################################################

# Modules.SetModule("GenericCleaning")
# Modules.DeleteWorkdirs()
# Modules.CreateConfigFiles()
# Modules.CondorControl("Submit")
# Modules.CondorControl("Merge")
# Modules.CondorControl("List")

# Modules.CondorControl("ForceMerge")
# Modules.MergePreselection()

#################################################
#                                               #
#                GenLevelMatch                  #
#                                               #
#################################################

Modules.SetModule("GenLevelMatch")
# Modules.DeleteWorkdirs()
# Modules.CreateConfigFiles()
# Modules.CondorControl()
# Modules.CondorControl("List")
# Modules.CondorControl("Submit")

# Modules.CondorControl("ForceMerge")
# Modules.MergePreselection()



#################################################
#                                               #
#                   PRESELECTION                #
#                                               #
#################################################


# Modules.CompileModules()
# Modules.SetModule("FeasibilityStudy")
# Modules.RunLocal("Check")
# Modules.RunLocal("Local")

# Modules.DeleteWorkdirs()
# Modules.CreateConfigFiles()
# Modules.CondorControl("Submit")
# Modules.CondorControl("List")

# Modules.CondorControl("ForceMerge")
# Modules.MergePreselection()



#################################################
#                                               #
#                   NEURAL NETWORK              #
#                                               #
#################################################
# #
# Modules.CompileModules()
# Modules.SetModule("NeuralNetwork")
# Modules.DeleteWorkdirs()
# Modules.CreateConfigFiles()
# Modules.RunLocal("NoSplit")
# Modules.CondorControl("Submit")


#
#
# t->Draw("isHiggs>>hist","isHiggs>0.1");
# ((TH1F*)gDirectory->Get("hist"))->GetEntries()/t->GetEntries()*100
