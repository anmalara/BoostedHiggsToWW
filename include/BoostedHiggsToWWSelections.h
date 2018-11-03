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
  #define DEFINEPtLeptonSelection(Lepton)\
  class Pt##Lepton##Selection: public Selection {\
  public:\
    Pt##Lepton##Selection(float pt_lepton = 5, float min_isolation = 0.2);\
    virtual bool passes(const Event & event) override;\
  private:\
    float pt_lepton, min_isolation;\
  };\

  DEFINEPtLeptonSelection(Muon)
  DEFINEPtLeptonSelection(Electron)
  /////
  //
  // class PtElectronSelection: public Selection {
  // public:
  //   PtElectronSelection(float pt_electron = 5, float min_isolation = 0.2);
  //   virtual bool passes(const Event & event) override;
  // private:
  //   float pt_electron, min_isolation;
  // };

  /////
  //
  // class PtMuonSelection: public Selection {
  // public:
  //   PtMuonSelection(float pt_muon = 5, float min_isolation = 0.2);
  //   virtual bool passes(const Event & event) override;
  // private:
  //   float pt_muon, min_isolation;
  // };

  /////

  class LeptonPtIsoSelection: public Selection {
  public:
    LeptonPtIsoSelection(float pt_lep = 5, float min_isolation = 0.2, std::string lepton_class = "muon" );
    virtual bool passes(const Event & event) override;
  private:
    float pt_lep, min_isolation;
    std::string lepton_class;
  };

  /////

  class DiElectronSelection: public Selection {
  public:
    DiElectronSelection();
    virtual bool passes(const Event & event) override;
  };

  /////

  class DiMuonSelection: public Selection {
  public:
    DiMuonSelection();
    virtual bool passes(const Event & event) override;
  };

  /////

  class PhiAngularCut: public Selection {
  public:
    PhiAngularCut(float phi_min = 0, float phi_max = M_PI);
    virtual bool passes(const Event & event) override;
  private:
    float phi_min, phi_max;
  };


  /////

  class PhiAngularSelectionElectron: public Selection {
  public:
    PhiAngularSelectionElectron(float phi_min = 0, float phi_max = M_PI);
    virtual bool passes(const Event & event) override;
  private:
    float phi_min, phi_max;
  };

  /////

  class PhiAngularSelectionMuon: public Selection {
  public:
    PhiAngularSelectionMuon(float phi_min = 0, float phi_max = M_PI);
    virtual bool passes(const Event & event) override;
  private:
    float phi_min, phi_max;
  };


}
