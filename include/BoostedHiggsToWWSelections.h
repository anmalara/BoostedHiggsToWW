#pragma once

#include <cmath>
#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"
#include "UHH2/BoostedHiggsToWW/include/constants.hpp"

namespace uhh2 {

  /* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
  * below 20% of the average of the leading two jets, where the minimum deltaphi and
  * maximum third jet pt fraction can be changed in the constructor.
  * The jets are assumed to be sorted in pt.
  */

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  #define DEFINENLeptonPtEtaIsoSelection(Lepton)\
  class N##Lepton##PtEtaIsoSelection: public Selection {\
    public:\
    N##Lepton##PtEtaIsoSelection(int n_lepton = 2, float pt_lepton = 5, float eta_lepton = 4., float min_isolation = 0.2);\
    virtual bool passes(const Event & event) override;\
    private:\
    float n_lepton, pt_lepton, eta_lepton, min_isolation;\
  };\

  DEFINENLeptonPtEtaIsoSelection(Muon)
  DEFINENLeptonPtEtaIsoSelection(Electron)

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  #define DEFINEJetDiLeptonPhiAngularSelection(Lepton)\
  class JetDi##Lepton##PhiAngularSelection: public Selection {\
    public:\
    JetDi##Lepton##PhiAngularSelection(float phi_min = 0, float phi_max = M_PI);\
    virtual bool passes(const Event & event) override;\
    private:\
    float phi_min, phi_max;\
  };\

  DEFINEJetDiLeptonPhiAngularSelection(Muon)
  DEFINEJetDiLeptonPhiAngularSelection(Electron)

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  class TopJetHiggsMatching: public Selection {
  public:
    TopJetHiggsMatching(float deltaR_min = 0.8 );
    virtual bool passes(const Event & event) override;
  private:
    float deltaR_min;
  };

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

}
