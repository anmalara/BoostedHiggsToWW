#pragma once

#include "UHH2/core/include/Hists.h"

/**  \brief Example class for booking and filling histograms
*
* NOTE: This class uses the 'hist' method to retrieve histograms.
* This requires a string lookup and is therefore slow if you have
* many histograms. Therefore, it is recommended to use histogram
* pointers as member data instead, like in 'common/include/ElectronHists.h'.
*/
class GenParticleHiggsHists: public uhh2::Hists {
public:
	// use the same constructor arguments as Hists for forwarding:
	GenParticleHiggsHists(uhh2::Context & ctx, const std::string & dirname);

	virtual void fill(const uhh2::Event & ev) override;
	virtual ~GenParticleHiggsHists();

private:
	TH1F *numbers, *flavor_Z, *flavor_H, *flavor_W1, *flavor_W2, *flavor_muon, *flavor_ele;
	TH1F *number_muon, *pt_muon, *pt_muon_Z, *eta_muon, *phi_muon, *charge_muon;
	TH1F *number_ele, *pt_ele, *pt_ele_Z, *eta_ele, *phi_ele, *charge_ele;
	TH1F *mass_higgs, *pt_higgs, *eta_higgs, *phi_higgs;
	TH1F *mass_w1, *pt_w1, *eta_w1, *phi_w1;
	TH1F *mass_w2, *pt_w2, *eta_w2, *phi_w2;
	TH1F *mass_ww, *pt_ww, *eta_ww, *phi_ww;
	TH1F *DeltaRW1H, *DeltaRW2H, *DeltaRWW, *DeltaRW1jet, *DeltaRW2jet, *DeltaRHjet;
	TH2F *DeltaRW1H_DeltaRW2H, *DeltaRW1jet_DeltaRW2jet, *DeltaRHjet_ptjet, *pTH_pTJet, *pt_muon_pt_muon, *pt_ele_pt_ele;

};
