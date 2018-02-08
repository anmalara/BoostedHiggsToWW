#pragma once

#include "UHH2/core/include/Hists.h"

namespace uhh2examples {

  /**  \brief Example class for booking and filling histograms
  *
  * NOTE: This class uses the 'hist' method to retrieve histograms.
  * This requires a string lookup and is therefore slow if you have
  * many histograms. Therefore, it is recommended to use histogram
  * pointers as member data instead, like in 'common/include/ElectronHists.h'.
  */
  class BoostedHiggsToWWHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    BoostedHiggsToWWHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~BoostedHiggsToWWHists();

  protected:

    // declare all histograms as members. Note that one could also use get_hist
    // as in the example's ExampleHists instead of saving the histograms here. However,
    // that would entail quite a runtime overhead and it is much faster to declare the histograms
    // here directly.

    //jets
    TH1F *number_jet, *flavor_jet, *eta_jet, *pt_jet;
    TH1F *mass_jet1, *mass_jet2, *mass_jet3, *mass_jet4;
    TH1F *MT_jet1, *MT_jet2, *MT_jet3, *MT_jet4;
    TH1F *pt_jet1, *pt_jet2, *pt_jet3, *pt_jet4;
    TH1F *eta_jet1, *eta_jet2, *eta_jet3, *eta_jet4;
    TH1F *phi_jet1, *phi_jet2, *phi_jet3, *phi_jet4;
    //jet substructure
    TH1F *tau1_jet1, *tau1_jet2, *tau1_jet3, *tau1_jet4;
    TH1F *tau2_jet1, *tau2_jet2, *tau2_jet3, *tau2_jet4;
    TH1F *tau3_jet1, *tau3_jet2, *tau3_jet3, *tau3_jet4;
    TH1F *tau21_jet1, *tau21_jet2, *tau21_jet3, *tau21_jet4;
    TH1F *tau31_jet1, *tau31_jet2, *tau31_jet3, *tau31_jet4;
    TH1F *tau32_jet1, *tau32_jet2, *tau32_jet3, *tau32_jet4;
    TH1F *SDmass_jet1, *SDmass_jet2, *SDmass_jet3, *SDmass_jet4;
    // leptons
    TH1F *number_mu, *pt_mu, *eta_mu, *reliso_mu;
    TH1F *number_ele, *pt_ele, *eta_ele, *reliso_ele;
    // primary vertices
    TH1F *number_pv;
  };

}
