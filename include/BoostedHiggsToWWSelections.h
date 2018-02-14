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

  class BoostedJetSelection: public Selection {
  public:
    BoostedJetSelection(float pt_min = 300);
    virtual bool passes(const Event & event) override;
  private:
    float pt_min;
  };

  /////

  class DijetSelection: public Selection {
  public:
    DijetSelection(float dphi_min = 2.7f, float third_frac_max = 0.2f);
    virtual bool passes(const Event & event) override;
  private:
    float dphi_min, third_frac_max;
  };

  /////

  class PtElecSelection: public Selection {
  public:
    PtElecSelection(float pt_electron = 5);
    virtual bool passes(const Event & event) override;
  private:
    float pt_electron;
  };

  /////

  class PtMuonSelection: public Selection {
  public:
    PtMuonSelection(float pt_muon = 5);
    virtual bool passes(const Event & event) override;
  private:
    float pt_muon;
  };

  /////

  class DiElecSelection: public Selection {
  public:
    DiElecSelection(float min_isolation = 0.2);
    virtual bool passes(const Event & event) override;
  private:
    float min_isolation;
  };

  /////

  class DiMuonSelection: public Selection {
  public:
    DiMuonSelection(float min_isolation = 0.2);
    virtual bool passes(const Event & event) override;
  private:
    float min_isolation;
  };

  /////

  class ZSelection: public Selection {
  public:
    ZSelection();
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


}
