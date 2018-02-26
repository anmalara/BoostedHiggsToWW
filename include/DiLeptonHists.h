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
class DiLeptonHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    DiLeptonHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~DiLeptonHists();

  private:
    TH1F *dilep_number;
    TH1F *diele_number, *diele_charge, *diele_m, *diele_pt, *diele_eta, *diele_phi, *diele_deltaR;
    TH1F *dimuon_number, *dimuon_charge, *dimuon_m, *dimuon_pt, *dimuon_eta, *dimuon_phi, *dimuon_deltaR;
    TH2F *pt_muon_2, *pt_ele_2;
};

}
