#pragma once

#include "UHH2/core/include/Hists.h"

/**  \brief Example class for booking and filling histograms
*
* NOTE: This class uses the 'hist' method to retrieve histograms.
* This requires a string lookup and is therefore slow if you have
* many histograms. Therefore, it is recommended to use histogram
* pointers as member data instead, like in 'common/include/ElectronH
* pointers as member data instead, like in 'common/include/MuonHists.h'.
*/
class DiLeptonHists: public uhh2::Hists {
public:
	// use the same constructor arguments as Hists for forwarding:
	DiLeptonHists(uhh2::Context & ctx, const std::string & dirname);

	virtual void fill(const uhh2::Event & ev) override;
	virtual ~DiLeptonHists();

private:
	TH1F *dilep_number;

	#define DEFINEDILEPTONHISTOS(lepton)\
	TH1F *di##lepton##_number, *di##lepton##_charge, *di##lepton##_m, *di##lepton##_pt, *di##lepton##_eta, *di##lepton##_phi, *di##lepton##_deltaR;\
	TH2F *pt_##lepton##_2D;\

	DEFINEDILEPTONHISTOS(Electron)
	DEFINEDILEPTONHISTOS(Muon)

};
