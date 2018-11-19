#pragma once

#include "UHH2/core/include/Hists.h"
#include "TString.h"
#include "UHH2/BoostedHiggsToWW/include/constants.hpp"

/**  \brief Example class for booking and filling histograms
*
* NOTE: This class uses the 'hist' method to retrieve histograms.
* This requires a string lookup and is therefore slow if you have
* many histograms. Therefore, it is recommended to use histogram
* pointers as member data instead, like in 'common/include/ElectronHists.h'.
*/
class GenericHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  GenericHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~GenericHists();

protected:

  // declare all histograms as members. Note that one could also use get_hist
  // as in the example's ExampleHists instead of saving the histograms here. However,
  // that would entail quite a runtime overhead and it is much faster to declare the histograms
  // here directly.

  //jets
  TH1F *weights, *number_jet, *flavor_jet, *eta_jet, *pt_jet;

  #define DEFINEBOOSTEDHISTOS(jet)\
  TH1F *mass_##jet, *MT_##jet, *pt_##jet, *eta_##jet, *phi_##jet, *SDmass_##jet;\
  TH1F *tau1_##jet, *tau2_##jet, *tau3_##jet, *tau21_##jet, *tau31_##jet, *tau32_##jet;\

  DEFINEBOOSTEDHISTOS(jet1)
  DEFINEBOOSTEDHISTOS(jet2)
  DEFINEBOOSTEDHISTOS(jet3)
  DEFINEBOOSTEDHISTOS(jet4)

  // leptons
  TH1F *number_mu, *pt_mu, *eta_mu, *reliso_mu;
  TH1F *number_ele, *pt_ele, *eta_ele, *reliso_ele;
  // primary vertices
  TH1F *number_pv;
};
