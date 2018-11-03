#pragma once

#include <cmath>
#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"

namespace uhh2 {

  /* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
  * below 20% of the average of the leading two jets, where the minimum deltaphi and
  * maximum third jet pt fraction can be changed in the constructor.
  * The jets are assumed to be sorted in pt.
  */

  /////
  #define DEFINELeptonPtIsoSelection(Lepton)\
  class Lepton##PtIsoSelection: public Selection {\
  public:\
    Lepton##PtIsoSelection(float pt_lepton = 5, float min_isolation = 0.2);\
    virtual bool passes(const Event & event) override;\
  private:\
    float pt_lepton, min_isolation;\
  };\

  DEFINELeptonPtIsoSelection(Muon)
  DEFINELeptonPtIsoSelection(Electron)

  /////

  #define  DEFINEDiLeptonSelection(Lepton)\
  class Di##Lepton##Selection: public Selection {\
  public:\
    Di##Lepton##Selection();\
    virtual bool passes(const Event & event) override;\
  };\

  DEFINEDiLeptonSelection(Muon)
  DEFINEDiLeptonSelection(Electron)

  /////

  #define DEFINEJetLeptonPhiAngularSelection(Lepton)\
  class Jet##Lepton##PhiAngularSelection: public Selection {\
  public:\
    Jet##Lepton##PhiAngularSelection(float phi_min = 0, float phi_max = M_PI);\
    virtual bool passes(const Event & event) override;\
  private:\
    float phi_min, phi_max;\
  };\

  DEFINEJetLeptonPhiAngularSelection(Muon)
  DEFINEJetLeptonPhiAngularSelection(Electron)

}
