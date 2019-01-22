//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 18 15:57:23 2019 by ROOT version 6.10/09
// from TTree AnalysisTree/Format: User, data type: MC
// found on file: uhh2.AnalysisModuleRunner.MC_HZ.root
//////////////////////////////////////////////////////////

#ifndef AnalysisTree_h
#define AnalysisTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/FlavorParticle.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/GenInfo.h"
#include "vector"
#include "UHH2/core/include/PrimaryVertex.h"
#include "vector"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/JetBTagInfo.h"
#include "UHH2/core/include/Tags.h"
#include "UHH2/core/include/TopJet.h"
#include "vector"
#include "UHH2/core/include/RecParticle.h"
#include "UHH2/core/include/Electron.h"
#include "vector"
#include "vector"
#include "UHH2/core/include/MET.h"
#include "vector"
#include "UHH2/core/include/Muon.h"
#include "vector"
#include "UHH2/core/include/Tau.h"

class AnalysisTree {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  TString outname;

  // Fixed size dimensions of array or collections stored in the TTree if any.
  static constexpr Int_t kMaxGenParticles = 14;
  static constexpr Int_t kMaxofflineSlimmedPrimaryVertices = 89;
  static constexpr Int_t kMaxpackedPatJetsAk8CHSJets_SoftDropCHS = 3;
  static constexpr Int_t kMaxslimmedElectronsUSER = 3;
  static constexpr Int_t kMaxslimmedGenJets = 16;
  static constexpr Int_t kMaxslimmedJets = 9;
  static constexpr Int_t kMaxslimmedMuonsUSER = 1;
  static constexpr Int_t kMaxslimmedTaus = 9;

  // Declaration of leaf types
  Int_t           GenParticles_;
  Short_t         GenParticles_m_charge[kMaxGenParticles];   //[GenParticles_]
  Float_t         GenParticles_m_pt[kMaxGenParticles];   //[GenParticles_]
  Float_t         GenParticles_m_eta[kMaxGenParticles];   //[GenParticles_]
  Float_t         GenParticles_m_phi[kMaxGenParticles];   //[GenParticles_]
  Float_t         GenParticles_m_energy[kMaxGenParticles];   //[GenParticles_]
  Int_t           GenParticles_m_pdgId[kMaxGenParticles];   //[GenParticles_]
  Short_t         GenParticles_m_status[kMaxGenParticles];   //[GenParticles_]
  UShort_t        GenParticles_m_index[kMaxGenParticles];   //[GenParticles_]
  UShort_t        GenParticles_m_mother1[kMaxGenParticles];   //[GenParticles_]
  UShort_t        GenParticles_m_mother2[kMaxGenParticles];   //[GenParticles_]
  UShort_t        GenParticles_m_daughter1[kMaxGenParticles];   //[GenParticles_]
  UShort_t        GenParticles_m_daughter2[kMaxGenParticles];   //[GenParticles_]
  Short_t         GenParticles_m_spin[kMaxGenParticles];   //[GenParticles_]
  Float_t         beamspot_x0;
  Float_t         beamspot_y0;
  Float_t         beamspot_z0;
  Int_t           event;
  //GenInfo         *genInfo;
  vector<float>   m_binningValues;
  vector<float>   m_weights;
  vector<float>   m_systweights;
  Float_t         m_originalXWGTUP;
  Float_t         m_alphaQCD;
  Float_t         m_alphaQED;
  Float_t         m_qScale;
  Int_t           m_pdf_id1;
  Int_t           m_pdf_id2;
  Float_t         m_pdf_x1;
  Float_t         m_pdf_x2;
  Float_t         m_pdf_xPDF1;
  Float_t         m_pdf_xPDF2;
  Float_t         m_pdf_scalePDF;
  Int_t           m_pileup_NumInteractions_intime;
  Int_t           m_pileup_NumInteractions_ootbefore;
  Int_t           m_pileup_NumInteractions_ootafter;
  Float_t         m_pileup_TrueNumInteractions;
  Float_t         m_PU_pT_hat_max;
  Bool_t          isRealData;
  Float_t         lumi_weights;
  Int_t           luminosityBlock;
  Float_t         my_weights;
  Int_t           offlineSlimmedPrimaryVertices_;
  Float_t         offlineSlimmedPrimaryVertices_m_x[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_y[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_z[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  UInt_t          offlineSlimmedPrimaryVertices_m_nTracks[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_chi2[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_ndof[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_;
  Short_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_charge[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_pt[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_eta[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_phi[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_energy[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_pdgId[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_jetArea[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_numberOfDaughters[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralEmEnergyFraction[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronEnergyFraction[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedEmEnergyFraction[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedHadronEnergyFraction[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_muonEnergyFraction[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_photonEnergyFraction[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_muonMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_electronMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_photonMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_puppiMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralPuppiMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronPuppiMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_photonPuppiMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_HFHadronPuppiMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_HFEMPuppiMultiplicity[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertex[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertexMVA[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probb[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probbb[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexAK8[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexCA15[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_factor_raw[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_JER_factor_raw[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_L1factor_raw[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_genjet_index[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_hadronFlavor[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackMomentum[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEta[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEtaRel[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDeltaR[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dVal[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSig[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dVal[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dSig[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDecayLenVal[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackChi2[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNTotalHits[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNPixelHits[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRel[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPPar[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRatio[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPParRatio[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistVal[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistSig[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistVal[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistSig[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackWeight[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dVal[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dSig[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dVal[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dSig[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexJetDeltaR[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_JetNSecondaryVertices[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  vector<int>     packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNTracks[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<TLorentzVector> packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_SecondaryVertex[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexChi2[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNdof[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<float>   packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNormalizedChi2[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  Int_t           packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexCategoryJTC[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexMassJTC[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexEnergyRatioJTC[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  map<int,float>  packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_tags_tagdata[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<long>    packedPatJetsAk8CHSJets_SoftDropCHS_m_lepton_keys[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  map<int,float>  packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  vector<Jet>     packedPatJetsAk8CHSJets_SoftDropCHS_m_subjets[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_qjets_volatility[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1_groomed[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2_groomed[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3_groomed[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4_groomed[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta1[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta2[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta1[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta2[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_mvahiggsdiscr[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_prunedmass[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  Float_t         packedPatJetsAk8CHSJets_SoftDropCHS_m_softdropmass[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];   //[packedPatJetsAk8CHSJets_SoftDropCHS_]
  // map<int,float>  packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
  Float_t         rho;
  Int_t           run;
  Int_t           slimmedElectronsUSER_;
  Short_t         slimmedElectronsUSER_m_charge[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_pt[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_eta[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_phi[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_energy[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_ptError[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_etaError[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_phiError[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_supercluster_eta[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_supercluster_phi[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_dB[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_neutralHadronIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_chargedHadronIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_photonIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_trackIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_puChargedHadronIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Int_t           slimmedElectronsUSER_m_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_gsfTrack_px[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_gsfTrack_py[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_gsfTrack_pz[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_gsfTrack_vx[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_gsfTrack_vy[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_gsfTrack_vz[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Bool_t          slimmedElectronsUSER_m_passconversionveto[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_dEtaIn[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_dPhiIn[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_sigmaIEtaIEta[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_HoverE[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_fbrem[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_EoverPIn[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_EcalEnergy[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_hcalOverEcal[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_ecalPFClusterIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_hcalPFClusterIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_dr03TkSumPt[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_mvaIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_mvaNoIso[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_AEff[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_pfMINIIso_CH[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_pfMINIIso_NH[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_pfMINIIso_Ph[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_pfMINIIso_PU[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_pfMINIIso_NH_pfwgt[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_pfMINIIso_Ph_pfwgt[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  vector<source_candidate> slimmedElectronsUSER_m_source_candidates[kMaxslimmedElectronsUSER];
  map<int,float>  slimmedElectronsUSER_tags_tagdata[kMaxslimmedElectronsUSER];
  Int_t           slimmedElectronsUSER_m_Nclusters[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Int_t           slimmedElectronsUSER_m_Class[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Bool_t          slimmedElectronsUSER_m_isEcalDriven[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_dxy[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_dEtaInSeed[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_full5x5_e1x5[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_full5x5_e2x5Max[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Float_t         slimmedElectronsUSER_m_full5x5_e5x5[kMaxslimmedElectronsUSER];   //[slimmedElectronsUSER_]
  Int_t           slimmedGenJets_;
  Short_t         slimmedGenJets_m_charge[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_pt[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_eta[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_phi[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_energy[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Int_t           slimmedJets_;
  Short_t         slimmedJets_m_charge[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_pt[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_eta[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_phi[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_energy[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_pdgId[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_jetArea[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_numberOfDaughters[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_neutralEmEnergyFraction[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_neutralHadronEnergyFraction[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_chargedEmEnergyFraction[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_chargedHadronEnergyFraction[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_muonEnergyFraction[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_photonEnergyFraction[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_chargedMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_neutralMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_muonMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_electronMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_photonMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_puppiMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_neutralPuppiMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_neutralHadronPuppiMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_photonPuppiMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_HFHadronPuppiMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_HFEMPuppiMultiplicity[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btag_combinedSecondaryVertex[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btag_combinedSecondaryVertexMVA[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btag_DeepCSV_probb[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btag_DeepCSV_probbb[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btag_BoostedDoubleSecondaryVertexAK8[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btag_BoostedDoubleSecondaryVertexCA15[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_JEC_factor_raw[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_JER_factor_raw[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_JEC_L1factor_raw[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_genjet_index[kMaxslimmedJets];   //[slimmedJets_]
  Int_t           slimmedJets_m_hadronFlavor[kMaxslimmedJets];   //[slimmedJets_]
  vector<float>   slimmedJets_m_btaginfo_m_TrackMomentum[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackEta[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackEtaRel[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackDeltaR[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackSip3dVal[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackSip3dSig[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackSip2dVal[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackSip2dSig[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackDecayLenVal[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackChi2[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackNTotalHits[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackNPixelHits[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackPtRel[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackPPar[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackPtRatio[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackPParRatio[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackJetDistVal[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackJetDistSig[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackGhostTrackDistVal[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackGhostTrackDistSig[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_TrackGhostTrackWeight[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_FlightDistance2dVal[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_FlightDistance2dSig[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_FlightDistance3dVal[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_FlightDistance3dSig[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_VertexJetDeltaR[kMaxslimmedJets];
  Int_t           slimmedJets_m_btaginfo_m_JetNSecondaryVertices[kMaxslimmedJets];   //[slimmedJets_]
  vector<int>     slimmedJets_m_btaginfo_m_VertexNTracks[kMaxslimmedJets];
  vector<TLorentzVector> slimmedJets_m_btaginfo_m_SecondaryVertex[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_VertexChi2[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_VertexNdof[kMaxslimmedJets];
  vector<float>   slimmedJets_m_btaginfo_m_VertexNormalizedChi2[kMaxslimmedJets];
  Int_t           slimmedJets_m_btaginfo_m_VertexCategoryJTC[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btaginfo_m_VertexMassJTC[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btaginfo_m_VertexEnergyRatioJTC[kMaxslimmedJets];   //[slimmedJets_]
  Float_t         slimmedJets_m_btaginfo_m_TrackSip3dSigAboveCharmJTC[kMaxslimmedJets];   //[slimmedJets_]
  map<int,float>  slimmedJets_m_btaginfo_tags_tagdata[kMaxslimmedJets];
  vector<long>    slimmedJets_m_lepton_keys[kMaxslimmedJets];
  map<int,float>  slimmedJets_tags_tagdata[kMaxslimmedJets];
  //MET             *slimmedMETs;
  Float_t         m_pt;
  Float_t         m_phi;
  Float_t         m_mEtSig;
  Float_t         m_shiftedPx_JetEnUp;
  Float_t         m_shiftedPx_JetEnDown;
  Float_t         m_shiftedPx_JetResUp;
  Float_t         m_shiftedPx_JetResDown;
  Float_t         m_shiftedPx_UnclusteredEnUp;
  Float_t         m_shiftedPx_UnclusteredEnDown;
  Float_t         m_shiftedPx_ElectronEnUp;
  Float_t         m_shiftedPx_ElectronEnDown;
  Float_t         m_shiftedPx_TauEnUp;
  Float_t         m_shiftedPx_TauEnDown;
  Float_t         m_shiftedPx_MuonEnDown;
  Float_t         m_shiftedPx_MuonEnUp;
  Float_t         m_shiftedPy_JetEnUp;
  Float_t         m_shiftedPy_JetEnDown;
  Float_t         m_shiftedPy_JetResUp;
  Float_t         m_shiftedPy_JetResDown;
  Float_t         m_shiftedPy_UnclusteredEnUp;
  Float_t         m_shiftedPy_UnclusteredEnDown;
  Float_t         m_shiftedPy_ElectronEnUp;
  Float_t         m_shiftedPy_ElectronEnDown;
  Float_t         m_shiftedPy_TauEnUp;
  Float_t         m_shiftedPy_TauEnDown;
  Float_t         m_shiftedPy_MuonEnDown;
  Float_t         m_shiftedPy_MuonEnUp;
  Float_t         m_rawCHS_px;
  Float_t         m_rawCHS_py;
  Float_t         m_uncorr_pt;
  Float_t         m_uncorr_phi;
  Int_t           slimmedMuonsUSER_;
  Short_t         slimmedMuonsUSER_m_charge[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_pt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_eta[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_phi[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_energy[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  ULong_t         slimmedMuonsUSER_sel_bits[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_dxy[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_dxy_error[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_dz[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_dz_error[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_globalTrack_normalizedChi2[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_globalTrack_numberOfValidMuonHits[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_numberOfMatchedStations[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_innerTrack_trackerLayersWithMeasurement[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_innerTrack_numberOfValidPixelHits[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_innerTrack_validFraction[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_combinedQuality_chi2LocalPosition[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_combinedQuality_trkKink[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_segmentCompatibility[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_sumChargedHadronPt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_sumNeutralHadronEt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_sumPhotonEt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_sumPUPt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_pfMINIIso_CH[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_pfMINIIso_NH[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_pfMINIIso_Ph[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_pfMINIIso_PU[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_pfMINIIso_NH_pfwgt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_pfMINIIso_Ph_pfwgt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_simType[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_simFlavor[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_simPdgId[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_simMotherPdgId[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_simHeaviestMotherFlavor[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_tunePTrackPt[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_tunePTrackEta[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Float_t         slimmedMuonsUSER_m_tunePTrackPhi[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  Int_t           slimmedMuonsUSER_m_tunePTrackType[kMaxslimmedMuonsUSER];   //[slimmedMuonsUSER_]
  vector<source_candidate> slimmedMuonsUSER_m_source_candidates[kMaxslimmedMuonsUSER];
  map<int,float>  slimmedMuonsUSER_tags_tagdata[kMaxslimmedMuonsUSER];
  Int_t           slimmedTaus_;
  Short_t         slimmedTaus_m_charge[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_pt[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_eta[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_phi[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_energy[kMaxslimmedTaus];   //[slimmedTaus_]
  ULong_t         slimmedTaus_id_bits[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_byCombinedIsolationDeltaBetaCorrRaw3Hits[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_byIsolationMVA3newDMwoLTraw[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_byIsolationMVArun2v1DBnewDMwLTraw[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_chargedIsoPtSum[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_neutralIsoPtSum[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_puCorrPtSum[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_neutralIsoPtSumWeight[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_footprintCorrection[kMaxslimmedTaus];   //[slimmedTaus_]
  Float_t         slimmedTaus_m_photonPtSumOutsideSignalCone[kMaxslimmedTaus];   //[slimmedTaus_]
  Int_t           slimmedTaus_m_decayMode[kMaxslimmedTaus];   //[slimmedTaus_]
  map<int,float>  slimmedTaus_tags_tagdata[kMaxslimmedTaus];
  Float_t         weight_pu;
  Float_t         weight_pu_down;
  Float_t         weight_pu_up;

  // List of branches
  TBranch        *b_GenParticles_;   //!
  TBranch        *b_GenParticles_m_charge;   //!
  TBranch        *b_GenParticles_m_pt;   //!
  TBranch        *b_GenParticles_m_eta;   //!
  TBranch        *b_GenParticles_m_phi;   //!
  TBranch        *b_GenParticles_m_energy;   //!
  TBranch        *b_GenParticles_m_pdgId;   //!
  TBranch        *b_GenParticles_m_status;   //!
  TBranch        *b_GenParticles_m_index;   //!
  TBranch        *b_GenParticles_m_mother1;   //!
  TBranch        *b_GenParticles_m_mother2;   //!
  TBranch        *b_GenParticles_m_daughter1;   //!
  TBranch        *b_GenParticles_m_daughter2;   //!
  TBranch        *b_GenParticles_m_spin;   //!
  TBranch        *b_beamspot_x0;   //!
  TBranch        *b_beamspot_y0;   //!
  TBranch        *b_beamspot_z0;   //!
  TBranch        *b_event;   //!
  TBranch        *b_genInfo_m_binningValues;   //!
  TBranch        *b_genInfo_m_weights;   //!
  TBranch        *b_genInfo_m_systweights;   //!
  TBranch        *b_genInfo_m_originalXWGTUP;   //!
  TBranch        *b_genInfo_m_alphaQCD;   //!
  TBranch        *b_genInfo_m_alphaQED;   //!
  TBranch        *b_genInfo_m_qScale;   //!
  TBranch        *b_genInfo_m_pdf_id1;   //!
  TBranch        *b_genInfo_m_pdf_id2;   //!
  TBranch        *b_genInfo_m_pdf_x1;   //!
  TBranch        *b_genInfo_m_pdf_x2;   //!
  TBranch        *b_genInfo_m_pdf_xPDF1;   //!
  TBranch        *b_genInfo_m_pdf_xPDF2;   //!
  TBranch        *b_genInfo_m_pdf_scalePDF;   //!
  TBranch        *b_genInfo_m_pileup_NumInteractions_intime;   //!
  TBranch        *b_genInfo_m_pileup_NumInteractions_ootbefore;   //!
  TBranch        *b_genInfo_m_pileup_NumInteractions_ootafter;   //!
  TBranch        *b_genInfo_m_pileup_TrueNumInteractions;   //!
  TBranch        *b_genInfo_m_PU_pT_hat_max;   //!
  TBranch        *b_isRealData;   //!
  TBranch        *b_lumi_weights;   //!
  TBranch        *b_luminosityBlock;   //!
  TBranch        *b_my_weights;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_x;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_y;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_z;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_nTracks;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_chi2;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_ndof;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_charge;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_pt;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_eta;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_phi;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_energy;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_pdgId;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_jetArea;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_numberOfDaughters;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralEmEnergyFraction;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronEnergyFraction;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedEmEnergyFraction;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedHadronEnergyFraction;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_muonEnergyFraction;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_photonEnergyFraction;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_muonMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_electronMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_photonMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_puppiMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralPuppiMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronPuppiMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_photonPuppiMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_HFHadronPuppiMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_HFEMPuppiMultiplicity;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertex;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertexMVA;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probb;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probbb;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexAK8;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexCA15;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_factor_raw;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_JER_factor_raw;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_L1factor_raw;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_genjet_index;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_hadronFlavor;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackMomentum;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEta;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEtaRel;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDeltaR;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dVal;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSig;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dVal;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dSig;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDecayLenVal;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackChi2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNTotalHits;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNPixelHits;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRel;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPPar;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRatio;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPParRatio;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistVal;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistSig;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistVal;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistSig;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackWeight;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dVal;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dSig;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dVal;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dSig;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexJetDeltaR;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_JetNSecondaryVertices;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNTracks;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_SecondaryVertex;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_lepton_keys;   //!
  // TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_subjets;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_qjets_volatility;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1_groomed;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2_groomed;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3_groomed;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4_groomed;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta1;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta1;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_mvahiggsdiscr;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_prunedmass;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_softdropmass;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_run;   //!
  TBranch        *b_slimmedElectronsUSER_;   //!
  TBranch        *b_slimmedElectronsUSER_m_charge;   //!
  TBranch        *b_slimmedElectronsUSER_m_pt;   //!
  TBranch        *b_slimmedElectronsUSER_m_eta;   //!
  TBranch        *b_slimmedElectronsUSER_m_phi;   //!
  TBranch        *b_slimmedElectronsUSER_m_energy;   //!
  TBranch        *b_slimmedElectronsUSER_m_ptError;   //!
  TBranch        *b_slimmedElectronsUSER_m_etaError;   //!
  TBranch        *b_slimmedElectronsUSER_m_phiError;   //!
  TBranch        *b_slimmedElectronsUSER_m_supercluster_eta;   //!
  TBranch        *b_slimmedElectronsUSER_m_supercluster_phi;   //!
  TBranch        *b_slimmedElectronsUSER_m_dB;   //!
  TBranch        *b_slimmedElectronsUSER_m_neutralHadronIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_chargedHadronIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_photonIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_trackIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_puChargedHadronIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_gsfTrack_trackerExpectedHitsInner_numberOfLostHits;   //!
  TBranch        *b_slimmedElectronsUSER_m_gsfTrack_px;   //!
  TBranch        *b_slimmedElectronsUSER_m_gsfTrack_py;   //!
  TBranch        *b_slimmedElectronsUSER_m_gsfTrack_pz;   //!
  TBranch        *b_slimmedElectronsUSER_m_gsfTrack_vx;   //!
  TBranch        *b_slimmedElectronsUSER_m_gsfTrack_vy;   //!
  TBranch        *b_slimmedElectronsUSER_m_gsfTrack_vz;   //!
  TBranch        *b_slimmedElectronsUSER_m_passconversionveto;   //!
  TBranch        *b_slimmedElectronsUSER_m_dEtaIn;   //!
  TBranch        *b_slimmedElectronsUSER_m_dPhiIn;   //!
  TBranch        *b_slimmedElectronsUSER_m_sigmaIEtaIEta;   //!
  TBranch        *b_slimmedElectronsUSER_m_HoverE;   //!
  TBranch        *b_slimmedElectronsUSER_m_fbrem;   //!
  TBranch        *b_slimmedElectronsUSER_m_EoverPIn;   //!
  TBranch        *b_slimmedElectronsUSER_m_EcalEnergy;   //!
  TBranch        *b_slimmedElectronsUSER_m_hcalOverEcal;   //!
  TBranch        *b_slimmedElectronsUSER_m_ecalPFClusterIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_hcalPFClusterIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_dr03TkSumPt;   //!
  TBranch        *b_slimmedElectronsUSER_m_mvaIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_mvaNoIso;   //!
  TBranch        *b_slimmedElectronsUSER_m_AEff;   //!
  TBranch        *b_slimmedElectronsUSER_m_pfMINIIso_CH;   //!
  TBranch        *b_slimmedElectronsUSER_m_pfMINIIso_NH;   //!
  TBranch        *b_slimmedElectronsUSER_m_pfMINIIso_Ph;   //!
  TBranch        *b_slimmedElectronsUSER_m_pfMINIIso_PU;   //!
  TBranch        *b_slimmedElectronsUSER_m_pfMINIIso_NH_pfwgt;   //!
  TBranch        *b_slimmedElectronsUSER_m_pfMINIIso_Ph_pfwgt;   //!
  TBranch        *b_slimmedElectronsUSER_m_source_candidates;   //!
  TBranch        *b_slimmedElectronsUSER_tags_tagdata;   //!
  TBranch        *b_slimmedElectronsUSER_m_Nclusters;   //!
  TBranch        *b_slimmedElectronsUSER_m_Class;   //!
  TBranch        *b_slimmedElectronsUSER_m_isEcalDriven;   //!
  TBranch        *b_slimmedElectronsUSER_m_dxy;   //!
  TBranch        *b_slimmedElectronsUSER_m_dEtaInSeed;   //!
  TBranch        *b_slimmedElectronsUSER_m_full5x5_e1x5;   //!
  TBranch        *b_slimmedElectronsUSER_m_full5x5_e2x5Max;   //!
  TBranch        *b_slimmedElectronsUSER_m_full5x5_e5x5;   //!
  TBranch        *b_slimmedGenJets_;   //!
  TBranch        *b_slimmedGenJets_m_charge;   //!
  TBranch        *b_slimmedGenJets_m_pt;   //!
  TBranch        *b_slimmedGenJets_m_eta;   //!
  TBranch        *b_slimmedGenJets_m_phi;   //!
  TBranch        *b_slimmedGenJets_m_energy;   //!
  TBranch        *b_slimmedJets_;   //!
  TBranch        *b_slimmedJets_m_charge;   //!
  TBranch        *b_slimmedJets_m_pt;   //!
  TBranch        *b_slimmedJets_m_eta;   //!
  TBranch        *b_slimmedJets_m_phi;   //!
  TBranch        *b_slimmedJets_m_energy;   //!
  TBranch        *b_slimmedJets_m_pdgId;   //!
  TBranch        *b_slimmedJets_m_jetArea;   //!
  TBranch        *b_slimmedJets_m_numberOfDaughters;   //!
  TBranch        *b_slimmedJets_m_neutralEmEnergyFraction;   //!
  TBranch        *b_slimmedJets_m_neutralHadronEnergyFraction;   //!
  TBranch        *b_slimmedJets_m_chargedEmEnergyFraction;   //!
  TBranch        *b_slimmedJets_m_chargedHadronEnergyFraction;   //!
  TBranch        *b_slimmedJets_m_muonEnergyFraction;   //!
  TBranch        *b_slimmedJets_m_photonEnergyFraction;   //!
  TBranch        *b_slimmedJets_m_chargedMultiplicity;   //!
  TBranch        *b_slimmedJets_m_neutralMultiplicity;   //!
  TBranch        *b_slimmedJets_m_muonMultiplicity;   //!
  TBranch        *b_slimmedJets_m_electronMultiplicity;   //!
  TBranch        *b_slimmedJets_m_photonMultiplicity;   //!
  TBranch        *b_slimmedJets_m_puppiMultiplicity;   //!
  TBranch        *b_slimmedJets_m_neutralPuppiMultiplicity;   //!
  TBranch        *b_slimmedJets_m_neutralHadronPuppiMultiplicity;   //!
  TBranch        *b_slimmedJets_m_photonPuppiMultiplicity;   //!
  TBranch        *b_slimmedJets_m_HFHadronPuppiMultiplicity;   //!
  TBranch        *b_slimmedJets_m_HFEMPuppiMultiplicity;   //!
  TBranch        *b_slimmedJets_m_btag_combinedSecondaryVertex;   //!
  TBranch        *b_slimmedJets_m_btag_combinedSecondaryVertexMVA;   //!
  TBranch        *b_slimmedJets_m_btag_DeepCSV_probb;   //!
  TBranch        *b_slimmedJets_m_btag_DeepCSV_probbb;   //!
  TBranch        *b_slimmedJets_m_btag_BoostedDoubleSecondaryVertexAK8;   //!
  TBranch        *b_slimmedJets_m_btag_BoostedDoubleSecondaryVertexCA15;   //!
  TBranch        *b_slimmedJets_m_JEC_factor_raw;   //!
  TBranch        *b_slimmedJets_m_JER_factor_raw;   //!
  TBranch        *b_slimmedJets_m_JEC_L1factor_raw;   //!
  TBranch        *b_slimmedJets_m_genjet_index;   //!
  TBranch        *b_slimmedJets_m_hadronFlavor;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackMomentum;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackEta;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackEtaRel;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackDeltaR;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackSip3dVal;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackSip3dSig;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackSip2dVal;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackSip2dSig;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackDecayLenVal;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackChi2;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackNTotalHits;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackNPixelHits;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackPtRel;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackPPar;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackPtRatio;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackPParRatio;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackJetDistVal;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackJetDistSig;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackGhostTrackDistVal;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackGhostTrackDistSig;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackGhostTrackWeight;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_FlightDistance2dVal;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_FlightDistance2dSig;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_FlightDistance3dVal;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_FlightDistance3dSig;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexJetDeltaR;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_JetNSecondaryVertices;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexNTracks;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_SecondaryVertex;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_slimmedJets_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_slimmedJets_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_slimmedJets_m_lepton_keys;   //!
  TBranch        *b_slimmedJets_tags_tagdata;   //!
  TBranch        *b_slimmedMETs_m_pt;   //!
  TBranch        *b_slimmedMETs_m_phi;   //!
  TBranch        *b_slimmedMETs_m_mEtSig;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_JetEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_JetEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_JetResUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_JetResDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_UnclusteredEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_UnclusteredEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_ElectronEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_ElectronEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_TauEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_TauEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_MuonEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPx_MuonEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_JetEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_JetEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_JetResUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_JetResDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_UnclusteredEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_UnclusteredEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_ElectronEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_ElectronEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_TauEnUp;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_TauEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_MuonEnDown;   //!
  TBranch        *b_slimmedMETs_m_shiftedPy_MuonEnUp;   //!
  TBranch        *b_slimmedMETs_m_rawCHS_px;   //!
  TBranch        *b_slimmedMETs_m_rawCHS_py;   //!
  TBranch        *b_slimmedMETs_m_uncorr_pt;   //!
  TBranch        *b_slimmedMETs_m_uncorr_phi;   //!
  TBranch        *b_slimmedMuonsUSER_;   //!
  TBranch        *b_slimmedMuonsUSER_m_charge;   //!
  TBranch        *b_slimmedMuonsUSER_m_pt;   //!
  TBranch        *b_slimmedMuonsUSER_m_eta;   //!
  TBranch        *b_slimmedMuonsUSER_m_phi;   //!
  TBranch        *b_slimmedMuonsUSER_m_energy;   //!
  TBranch        *b_slimmedMuonsUSER_sel_bits;   //!
  TBranch        *b_slimmedMuonsUSER_m_dxy;   //!
  TBranch        *b_slimmedMuonsUSER_m_dxy_error;   //!
  TBranch        *b_slimmedMuonsUSER_m_dz;   //!
  TBranch        *b_slimmedMuonsUSER_m_dz_error;   //!
  TBranch        *b_slimmedMuonsUSER_m_globalTrack_normalizedChi2;   //!
  TBranch        *b_slimmedMuonsUSER_m_globalTrack_numberOfValidMuonHits;   //!
  TBranch        *b_slimmedMuonsUSER_m_numberOfMatchedStations;   //!
  TBranch        *b_slimmedMuonsUSER_m_innerTrack_trackerLayersWithMeasurement;   //!
  TBranch        *b_slimmedMuonsUSER_m_innerTrack_numberOfValidPixelHits;   //!
  TBranch        *b_slimmedMuonsUSER_m_innerTrack_validFraction;   //!
  TBranch        *b_slimmedMuonsUSER_m_combinedQuality_chi2LocalPosition;   //!
  TBranch        *b_slimmedMuonsUSER_m_combinedQuality_trkKink;   //!
  TBranch        *b_slimmedMuonsUSER_m_segmentCompatibility;   //!
  TBranch        *b_slimmedMuonsUSER_m_sumChargedHadronPt;   //!
  TBranch        *b_slimmedMuonsUSER_m_sumNeutralHadronEt;   //!
  TBranch        *b_slimmedMuonsUSER_m_sumPhotonEt;   //!
  TBranch        *b_slimmedMuonsUSER_m_sumPUPt;   //!
  TBranch        *b_slimmedMuonsUSER_m_pfMINIIso_CH;   //!
  TBranch        *b_slimmedMuonsUSER_m_pfMINIIso_NH;   //!
  TBranch        *b_slimmedMuonsUSER_m_pfMINIIso_Ph;   //!
  TBranch        *b_slimmedMuonsUSER_m_pfMINIIso_PU;   //!
  TBranch        *b_slimmedMuonsUSER_m_pfMINIIso_NH_pfwgt;   //!
  TBranch        *b_slimmedMuonsUSER_m_pfMINIIso_Ph_pfwgt;   //!
  TBranch        *b_slimmedMuonsUSER_m_simType;   //!
  TBranch        *b_slimmedMuonsUSER_m_simFlavor;   //!
  TBranch        *b_slimmedMuonsUSER_m_simPdgId;   //!
  TBranch        *b_slimmedMuonsUSER_m_simMotherPdgId;   //!
  TBranch        *b_slimmedMuonsUSER_m_simHeaviestMotherFlavor;   //!
  TBranch        *b_slimmedMuonsUSER_m_tunePTrackPt;   //!
  TBranch        *b_slimmedMuonsUSER_m_tunePTrackEta;   //!
  TBranch        *b_slimmedMuonsUSER_m_tunePTrackPhi;   //!
  TBranch        *b_slimmedMuonsUSER_m_tunePTrackType;   //!
  TBranch        *b_slimmedMuonsUSER_m_source_candidates;   //!
  TBranch        *b_slimmedMuonsUSER_tags_tagdata;   //!
  TBranch        *b_slimmedTaus_;   //!
  TBranch        *b_slimmedTaus_m_charge;   //!
  TBranch        *b_slimmedTaus_m_pt;   //!
  TBranch        *b_slimmedTaus_m_eta;   //!
  TBranch        *b_slimmedTaus_m_phi;   //!
  TBranch        *b_slimmedTaus_m_energy;   //!
  TBranch        *b_slimmedTaus_id_bits;   //!
  TBranch        *b_slimmedTaus_m_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
  TBranch        *b_slimmedTaus_m_byIsolationMVA3newDMwoLTraw;   //!
  TBranch        *b_slimmedTaus_m_byIsolationMVArun2v1DBnewDMwLTraw;   //!
  TBranch        *b_slimmedTaus_m_chargedIsoPtSum;   //!
  TBranch        *b_slimmedTaus_m_neutralIsoPtSum;   //!
  TBranch        *b_slimmedTaus_m_puCorrPtSum;   //!
  TBranch        *b_slimmedTaus_m_neutralIsoPtSumWeight;   //!
  TBranch        *b_slimmedTaus_m_footprintCorrection;   //!
  TBranch        *b_slimmedTaus_m_photonPtSumOutsideSignalCone;   //!
  TBranch        *b_slimmedTaus_m_decayMode;   //!
  TBranch        *b_slimmedTaus_tags_tagdata;   //!
  TBranch        *b_weight_pu;   //!
  TBranch        *b_weight_pu_down;   //!
  TBranch        *b_weight_pu_up;   //!

  AnalysisTree(TString name="uhh2.AnalysisModuleRunner.MC_DYJets.root",TString outname="uhh2.testUHHDY.root");
  virtual ~AnalysisTree();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalysisTree_cxx
AnalysisTree::AnalysisTree(TString name, TString outname_) : fChain(0), outname(outname_)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TTree *tree=0;
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(name);
    if (!f || !f->IsOpen()) {
      f = new TFile(name);
    }
    f->GetObject("AnalysisTree",tree);

  }
  Init(tree);
}

AnalysisTree::~AnalysisTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t AnalysisTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t AnalysisTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void AnalysisTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("GenParticles", &GenParticles_, &b_GenParticles_);
  fChain->SetBranchAddress("GenParticles.m_charge", GenParticles_m_charge, &b_GenParticles_m_charge);
  fChain->SetBranchAddress("GenParticles.m_pt", GenParticles_m_pt, &b_GenParticles_m_pt);
  fChain->SetBranchAddress("GenParticles.m_eta", GenParticles_m_eta, &b_GenParticles_m_eta);
  fChain->SetBranchAddress("GenParticles.m_phi", GenParticles_m_phi, &b_GenParticles_m_phi);
  fChain->SetBranchAddress("GenParticles.m_energy", GenParticles_m_energy, &b_GenParticles_m_energy);
  fChain->SetBranchAddress("GenParticles.m_pdgId", GenParticles_m_pdgId, &b_GenParticles_m_pdgId);
  fChain->SetBranchAddress("GenParticles.m_status", GenParticles_m_status, &b_GenParticles_m_status);
  fChain->SetBranchAddress("GenParticles.m_index", GenParticles_m_index, &b_GenParticles_m_index);
  fChain->SetBranchAddress("GenParticles.m_mother1", GenParticles_m_mother1, &b_GenParticles_m_mother1);
  fChain->SetBranchAddress("GenParticles.m_mother2", GenParticles_m_mother2, &b_GenParticles_m_mother2);
  fChain->SetBranchAddress("GenParticles.m_daughter1", GenParticles_m_daughter1, &b_GenParticles_m_daughter1);
  fChain->SetBranchAddress("GenParticles.m_daughter2", GenParticles_m_daughter2, &b_GenParticles_m_daughter2);
  fChain->SetBranchAddress("GenParticles.m_spin", GenParticles_m_spin, &b_GenParticles_m_spin);
  fChain->SetBranchAddress("beamspot_x0", &beamspot_x0, &b_beamspot_x0);
  fChain->SetBranchAddress("beamspot_y0", &beamspot_y0, &b_beamspot_y0);
  fChain->SetBranchAddress("beamspot_z0", &beamspot_z0, &b_beamspot_z0);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("m_binningValues", &m_binningValues, &b_genInfo_m_binningValues);
  fChain->SetBranchAddress("m_weights", &m_weights, &b_genInfo_m_weights);
  fChain->SetBranchAddress("m_systweights", &m_systweights, &b_genInfo_m_systweights);
  fChain->SetBranchAddress("m_originalXWGTUP", &m_originalXWGTUP, &b_genInfo_m_originalXWGTUP);
  fChain->SetBranchAddress("m_alphaQCD", &m_alphaQCD, &b_genInfo_m_alphaQCD);
  fChain->SetBranchAddress("m_alphaQED", &m_alphaQED, &b_genInfo_m_alphaQED);
  fChain->SetBranchAddress("m_qScale", &m_qScale, &b_genInfo_m_qScale);
  fChain->SetBranchAddress("m_pdf_id1", &m_pdf_id1, &b_genInfo_m_pdf_id1);
  fChain->SetBranchAddress("m_pdf_id2", &m_pdf_id2, &b_genInfo_m_pdf_id2);
  fChain->SetBranchAddress("m_pdf_x1", &m_pdf_x1, &b_genInfo_m_pdf_x1);
  fChain->SetBranchAddress("m_pdf_x2", &m_pdf_x2, &b_genInfo_m_pdf_x2);
  fChain->SetBranchAddress("m_pdf_xPDF1", &m_pdf_xPDF1, &b_genInfo_m_pdf_xPDF1);
  fChain->SetBranchAddress("m_pdf_xPDF2", &m_pdf_xPDF2, &b_genInfo_m_pdf_xPDF2);
  fChain->SetBranchAddress("m_pdf_scalePDF", &m_pdf_scalePDF, &b_genInfo_m_pdf_scalePDF);
  fChain->SetBranchAddress("m_pileup_NumInteractions_intime", &m_pileup_NumInteractions_intime, &b_genInfo_m_pileup_NumInteractions_intime);
  fChain->SetBranchAddress("m_pileup_NumInteractions_ootbefore", &m_pileup_NumInteractions_ootbefore, &b_genInfo_m_pileup_NumInteractions_ootbefore);
  fChain->SetBranchAddress("m_pileup_NumInteractions_ootafter", &m_pileup_NumInteractions_ootafter, &b_genInfo_m_pileup_NumInteractions_ootafter);
  fChain->SetBranchAddress("m_pileup_TrueNumInteractions", &m_pileup_TrueNumInteractions, &b_genInfo_m_pileup_TrueNumInteractions);
  fChain->SetBranchAddress("m_PU_pT_hat_max", &m_PU_pT_hat_max, &b_genInfo_m_PU_pT_hat_max);
  fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
  fChain->SetBranchAddress("lumi_weights", &lumi_weights, &b_lumi_weights);
  fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
  fChain->SetBranchAddress("my_weights", &my_weights, &b_my_weights);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices", &offlineSlimmedPrimaryVertices_, &b_offlineSlimmedPrimaryVertices_);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_x", offlineSlimmedPrimaryVertices_m_x, &b_offlineSlimmedPrimaryVertices_m_x);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_y", offlineSlimmedPrimaryVertices_m_y, &b_offlineSlimmedPrimaryVertices_m_y);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_z", offlineSlimmedPrimaryVertices_m_z, &b_offlineSlimmedPrimaryVertices_m_z);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_nTracks", offlineSlimmedPrimaryVertices_m_nTracks, &b_offlineSlimmedPrimaryVertices_m_nTracks);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_chi2", offlineSlimmedPrimaryVertices_m_chi2, &b_offlineSlimmedPrimaryVertices_m_chi2);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_ndof", offlineSlimmedPrimaryVertices_m_ndof, &b_offlineSlimmedPrimaryVertices_m_ndof);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS", &packedPatJetsAk8CHSJets_SoftDropCHS_, &b_packedPatJetsAk8CHSJets_SoftDropCHS_);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_charge", packedPatJetsAk8CHSJets_SoftDropCHS_m_charge, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_charge);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_pt", packedPatJetsAk8CHSJets_SoftDropCHS_m_pt, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_pt);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_eta", packedPatJetsAk8CHSJets_SoftDropCHS_m_eta, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_eta);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_phi", packedPatJetsAk8CHSJets_SoftDropCHS_m_phi, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_phi);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_energy", packedPatJetsAk8CHSJets_SoftDropCHS_m_energy, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_energy);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_pdgId", packedPatJetsAk8CHSJets_SoftDropCHS_m_pdgId, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_pdgId);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_jetArea", packedPatJetsAk8CHSJets_SoftDropCHS_m_jetArea, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_jetArea);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_numberOfDaughters", packedPatJetsAk8CHSJets_SoftDropCHS_m_numberOfDaughters, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_numberOfDaughters);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_neutralEmEnergyFraction", packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralEmEnergyFraction, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralEmEnergyFraction);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_neutralHadronEnergyFraction", packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronEnergyFraction, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronEnergyFraction);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_chargedEmEnergyFraction", packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedEmEnergyFraction, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedEmEnergyFraction);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_chargedHadronEnergyFraction", packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedHadronEnergyFraction, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedHadronEnergyFraction);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_muonEnergyFraction", packedPatJetsAk8CHSJets_SoftDropCHS_m_muonEnergyFraction, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_muonEnergyFraction);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_photonEnergyFraction", packedPatJetsAk8CHSJets_SoftDropCHS_m_photonEnergyFraction, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_photonEnergyFraction);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_chargedMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_chargedMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_neutralMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_muonMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_muonMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_muonMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_electronMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_electronMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_electronMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_photonMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_photonMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_photonMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_puppiMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_puppiMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_puppiMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_neutralPuppiMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralPuppiMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralPuppiMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_neutralHadronPuppiMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronPuppiMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_neutralHadronPuppiMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_photonPuppiMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_photonPuppiMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_photonPuppiMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_HFHadronPuppiMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_HFHadronPuppiMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_HFHadronPuppiMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_HFEMPuppiMultiplicity", packedPatJetsAk8CHSJets_SoftDropCHS_m_HFEMPuppiMultiplicity, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_HFEMPuppiMultiplicity);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btag_combinedSecondaryVertex", packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertex, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertex);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btag_combinedSecondaryVertexMVA", packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertexMVA, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_combinedSecondaryVertexMVA);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btag_DeepCSV_probb", packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probb, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probb);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btag_DeepCSV_probbb", packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probbb, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_DeepCSV_probbb);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btag_BoostedDoubleSecondaryVertexAK8", packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexAK8, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexAK8);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btag_BoostedDoubleSecondaryVertexCA15", packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexCA15, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btag_BoostedDoubleSecondaryVertexCA15);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_JEC_factor_raw", packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_factor_raw, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_factor_raw);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_JER_factor_raw", packedPatJetsAk8CHSJets_SoftDropCHS_m_JER_factor_raw, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_JER_factor_raw);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_JEC_L1factor_raw", packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_L1factor_raw, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_JEC_L1factor_raw);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_genjet_index", packedPatJetsAk8CHSJets_SoftDropCHS_m_genjet_index, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_genjet_index);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_hadronFlavor", packedPatJetsAk8CHSJets_SoftDropCHS_m_hadronFlavor, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_hadronFlavor);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackMomentum", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackMomentum, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackMomentum);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackEta", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEta, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEta);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackEtaRel", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEtaRel, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackEtaRel);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackDeltaR", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDeltaR, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDeltaR);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackSip3dVal", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dVal, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dVal);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackSip3dSig", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSig, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSig);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackSip2dVal", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dVal, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dVal);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackSip2dSig", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dSig, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip2dSig);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackDecayLenVal", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDecayLenVal, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackDecayLenVal);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackChi2", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackChi2, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackChi2);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackNTotalHits", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNTotalHits, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNTotalHits);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackNPixelHits", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNPixelHits, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackNPixelHits);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackPtRel", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRel, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRel);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackPPar", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPPar, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPPar);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackPtRatio", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRatio, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPtRatio);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackPParRatio", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPParRatio, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackPParRatio);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackJetDistVal", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistVal, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistVal);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackJetDistSig", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistSig, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackJetDistSig);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackGhostTrackDistVal", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistVal, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistVal);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackGhostTrackDistSig", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistSig, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackDistSig);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackGhostTrackWeight", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackWeight, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackGhostTrackWeight);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_FlightDistance2dVal", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dVal, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dVal);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_FlightDistance2dSig", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dSig, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance2dSig);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_FlightDistance3dVal", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dVal, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dVal);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_FlightDistance3dSig", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dSig, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_FlightDistance3dSig);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexJetDeltaR", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexJetDeltaR, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexJetDeltaR);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_JetNSecondaryVertices", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_JetNSecondaryVertices, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_JetNSecondaryVertices);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexNTracks", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNTracks, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNTracks);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_SecondaryVertex", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_SecondaryVertex, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_SecondaryVertex);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexChi2", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexChi2, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexChi2);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexNdof", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNdof, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNdof);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexNormalizedChi2", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNormalizedChi2, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNormalizedChi2);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexCategoryJTC", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexCategoryJTC, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexCategoryJTC);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexMassJTC", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexMassJTC, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexMassJTC);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_VertexEnergyRatioJTC", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexEnergyRatioJTC, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexEnergyRatioJTC);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.m_TrackSip3dSigAboveCharmJTC", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_btaginfo.tags.tagdata", packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_tags_tagdata, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_tags_tagdata);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_lepton_keys", packedPatJetsAk8CHSJets_SoftDropCHS_m_lepton_keys, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_lepton_keys);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.tags.tagdata", packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata, &b_packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_subjets", packedPatJetsAk8CHSJets_SoftDropCHS_m_subjets, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_subjets);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_qjets_volatility", packedPatJetsAk8CHSJets_SoftDropCHS_m_qjets_volatility, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_qjets_volatility);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau1", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau2", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau3", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau4", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau1_groomed", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1_groomed, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau1_groomed);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau2_groomed", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2_groomed, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau2_groomed);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau3_groomed", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3_groomed, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau3_groomed);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_tau4_groomed", packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4_groomed, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_tau4_groomed);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_ecfN2_beta1", packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta1, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta1);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_ecfN2_beta2", packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta2, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN2_beta2);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_ecfN3_beta1", packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta1, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta1);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_ecfN3_beta2", packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta2, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_ecfN3_beta2);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_mvahiggsdiscr", packedPatJetsAk8CHSJets_SoftDropCHS_m_mvahiggsdiscr, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_mvahiggsdiscr);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_prunedmass", packedPatJetsAk8CHSJets_SoftDropCHS_m_prunedmass, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_prunedmass);
  fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.m_softdropmass", packedPatJetsAk8CHSJets_SoftDropCHS_m_softdropmass, &b_packedPatJetsAk8CHSJets_SoftDropCHS_m_softdropmass);
  //    fChain->SetBranchAddress("packedPatJetsAk8CHSJets_SoftDropCHS.tags.tagdata", packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata, &b_packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata);
  fChain->SetBranchAddress("rho", &rho, &b_rho);
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("slimmedElectronsUSER", &slimmedElectronsUSER_, &b_slimmedElectronsUSER_);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_charge", slimmedElectronsUSER_m_charge, &b_slimmedElectronsUSER_m_charge);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_pt", slimmedElectronsUSER_m_pt, &b_slimmedElectronsUSER_m_pt);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_eta", slimmedElectronsUSER_m_eta, &b_slimmedElectronsUSER_m_eta);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_phi", slimmedElectronsUSER_m_phi, &b_slimmedElectronsUSER_m_phi);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_energy", slimmedElectronsUSER_m_energy, &b_slimmedElectronsUSER_m_energy);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_ptError", slimmedElectronsUSER_m_ptError, &b_slimmedElectronsUSER_m_ptError);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_etaError", slimmedElectronsUSER_m_etaError, &b_slimmedElectronsUSER_m_etaError);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_phiError", slimmedElectronsUSER_m_phiError, &b_slimmedElectronsUSER_m_phiError);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_supercluster_eta", slimmedElectronsUSER_m_supercluster_eta, &b_slimmedElectronsUSER_m_supercluster_eta);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_supercluster_phi", slimmedElectronsUSER_m_supercluster_phi, &b_slimmedElectronsUSER_m_supercluster_phi);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_dB", slimmedElectronsUSER_m_dB, &b_slimmedElectronsUSER_m_dB);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_neutralHadronIso", slimmedElectronsUSER_m_neutralHadronIso, &b_slimmedElectronsUSER_m_neutralHadronIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_chargedHadronIso", slimmedElectronsUSER_m_chargedHadronIso, &b_slimmedElectronsUSER_m_chargedHadronIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_photonIso", slimmedElectronsUSER_m_photonIso, &b_slimmedElectronsUSER_m_photonIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_trackIso", slimmedElectronsUSER_m_trackIso, &b_slimmedElectronsUSER_m_trackIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_puChargedHadronIso", slimmedElectronsUSER_m_puChargedHadronIso, &b_slimmedElectronsUSER_m_puChargedHadronIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_gsfTrack_trackerExpectedHitsInner_numberOfLostHits", slimmedElectronsUSER_m_gsfTrack_trackerExpectedHitsInner_numberOfLostHits, &b_slimmedElectronsUSER_m_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_gsfTrack_px", slimmedElectronsUSER_m_gsfTrack_px, &b_slimmedElectronsUSER_m_gsfTrack_px);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_gsfTrack_py", slimmedElectronsUSER_m_gsfTrack_py, &b_slimmedElectronsUSER_m_gsfTrack_py);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_gsfTrack_pz", slimmedElectronsUSER_m_gsfTrack_pz, &b_slimmedElectronsUSER_m_gsfTrack_pz);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_gsfTrack_vx", slimmedElectronsUSER_m_gsfTrack_vx, &b_slimmedElectronsUSER_m_gsfTrack_vx);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_gsfTrack_vy", slimmedElectronsUSER_m_gsfTrack_vy, &b_slimmedElectronsUSER_m_gsfTrack_vy);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_gsfTrack_vz", slimmedElectronsUSER_m_gsfTrack_vz, &b_slimmedElectronsUSER_m_gsfTrack_vz);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_passconversionveto", slimmedElectronsUSER_m_passconversionveto, &b_slimmedElectronsUSER_m_passconversionveto);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_dEtaIn", slimmedElectronsUSER_m_dEtaIn, &b_slimmedElectronsUSER_m_dEtaIn);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_dPhiIn", slimmedElectronsUSER_m_dPhiIn, &b_slimmedElectronsUSER_m_dPhiIn);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_sigmaIEtaIEta", slimmedElectronsUSER_m_sigmaIEtaIEta, &b_slimmedElectronsUSER_m_sigmaIEtaIEta);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_HoverE", slimmedElectronsUSER_m_HoverE, &b_slimmedElectronsUSER_m_HoverE);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_fbrem", slimmedElectronsUSER_m_fbrem, &b_slimmedElectronsUSER_m_fbrem);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_EoverPIn", slimmedElectronsUSER_m_EoverPIn, &b_slimmedElectronsUSER_m_EoverPIn);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_EcalEnergy", slimmedElectronsUSER_m_EcalEnergy, &b_slimmedElectronsUSER_m_EcalEnergy);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_hcalOverEcal", slimmedElectronsUSER_m_hcalOverEcal, &b_slimmedElectronsUSER_m_hcalOverEcal);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_ecalPFClusterIso", slimmedElectronsUSER_m_ecalPFClusterIso, &b_slimmedElectronsUSER_m_ecalPFClusterIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_hcalPFClusterIso", slimmedElectronsUSER_m_hcalPFClusterIso, &b_slimmedElectronsUSER_m_hcalPFClusterIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_dr03TkSumPt", slimmedElectronsUSER_m_dr03TkSumPt, &b_slimmedElectronsUSER_m_dr03TkSumPt);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_mvaIso", slimmedElectronsUSER_m_mvaIso, &b_slimmedElectronsUSER_m_mvaIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_mvaNoIso", slimmedElectronsUSER_m_mvaNoIso, &b_slimmedElectronsUSER_m_mvaNoIso);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_AEff", slimmedElectronsUSER_m_AEff, &b_slimmedElectronsUSER_m_AEff);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_pfMINIIso_CH", slimmedElectronsUSER_m_pfMINIIso_CH, &b_slimmedElectronsUSER_m_pfMINIIso_CH);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_pfMINIIso_NH", slimmedElectronsUSER_m_pfMINIIso_NH, &b_slimmedElectronsUSER_m_pfMINIIso_NH);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_pfMINIIso_Ph", slimmedElectronsUSER_m_pfMINIIso_Ph, &b_slimmedElectronsUSER_m_pfMINIIso_Ph);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_pfMINIIso_PU", slimmedElectronsUSER_m_pfMINIIso_PU, &b_slimmedElectronsUSER_m_pfMINIIso_PU);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_pfMINIIso_NH_pfwgt", slimmedElectronsUSER_m_pfMINIIso_NH_pfwgt, &b_slimmedElectronsUSER_m_pfMINIIso_NH_pfwgt);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_pfMINIIso_Ph_pfwgt", slimmedElectronsUSER_m_pfMINIIso_Ph_pfwgt, &b_slimmedElectronsUSER_m_pfMINIIso_Ph_pfwgt);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_source_candidates", slimmedElectronsUSER_m_source_candidates, &b_slimmedElectronsUSER_m_source_candidates);
  fChain->SetBranchAddress("slimmedElectronsUSER.tags.tagdata", slimmedElectronsUSER_tags_tagdata, &b_slimmedElectronsUSER_tags_tagdata);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_Nclusters", slimmedElectronsUSER_m_Nclusters, &b_slimmedElectronsUSER_m_Nclusters);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_Class", slimmedElectronsUSER_m_Class, &b_slimmedElectronsUSER_m_Class);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_isEcalDriven", slimmedElectronsUSER_m_isEcalDriven, &b_slimmedElectronsUSER_m_isEcalDriven);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_dxy", slimmedElectronsUSER_m_dxy, &b_slimmedElectronsUSER_m_dxy);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_dEtaInSeed", slimmedElectronsUSER_m_dEtaInSeed, &b_slimmedElectronsUSER_m_dEtaInSeed);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_full5x5_e1x5", slimmedElectronsUSER_m_full5x5_e1x5, &b_slimmedElectronsUSER_m_full5x5_e1x5);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_full5x5_e2x5Max", slimmedElectronsUSER_m_full5x5_e2x5Max, &b_slimmedElectronsUSER_m_full5x5_e2x5Max);
  fChain->SetBranchAddress("slimmedElectronsUSER.m_full5x5_e5x5", slimmedElectronsUSER_m_full5x5_e5x5, &b_slimmedElectronsUSER_m_full5x5_e5x5);
  fChain->SetBranchAddress("slimmedGenJets", &slimmedGenJets_, &b_slimmedGenJets_);
  fChain->SetBranchAddress("slimmedGenJets.m_charge", slimmedGenJets_m_charge, &b_slimmedGenJets_m_charge);
  fChain->SetBranchAddress("slimmedGenJets.m_pt", slimmedGenJets_m_pt, &b_slimmedGenJets_m_pt);
  fChain->SetBranchAddress("slimmedGenJets.m_eta", slimmedGenJets_m_eta, &b_slimmedGenJets_m_eta);
  fChain->SetBranchAddress("slimmedGenJets.m_phi", slimmedGenJets_m_phi, &b_slimmedGenJets_m_phi);
  fChain->SetBranchAddress("slimmedGenJets.m_energy", slimmedGenJets_m_energy, &b_slimmedGenJets_m_energy);
  fChain->SetBranchAddress("slimmedJets", &slimmedJets_, &b_slimmedJets_);
  fChain->SetBranchAddress("slimmedJets.m_charge", slimmedJets_m_charge, &b_slimmedJets_m_charge);
  fChain->SetBranchAddress("slimmedJets.m_pt", slimmedJets_m_pt, &b_slimmedJets_m_pt);
  fChain->SetBranchAddress("slimmedJets.m_eta", slimmedJets_m_eta, &b_slimmedJets_m_eta);
  fChain->SetBranchAddress("slimmedJets.m_phi", slimmedJets_m_phi, &b_slimmedJets_m_phi);
  fChain->SetBranchAddress("slimmedJets.m_energy", slimmedJets_m_energy, &b_slimmedJets_m_energy);
  fChain->SetBranchAddress("slimmedJets.m_pdgId", slimmedJets_m_pdgId, &b_slimmedJets_m_pdgId);
  fChain->SetBranchAddress("slimmedJets.m_jetArea", slimmedJets_m_jetArea, &b_slimmedJets_m_jetArea);
  fChain->SetBranchAddress("slimmedJets.m_numberOfDaughters", slimmedJets_m_numberOfDaughters, &b_slimmedJets_m_numberOfDaughters);
  fChain->SetBranchAddress("slimmedJets.m_neutralEmEnergyFraction", slimmedJets_m_neutralEmEnergyFraction, &b_slimmedJets_m_neutralEmEnergyFraction);
  fChain->SetBranchAddress("slimmedJets.m_neutralHadronEnergyFraction", slimmedJets_m_neutralHadronEnergyFraction, &b_slimmedJets_m_neutralHadronEnergyFraction);
  fChain->SetBranchAddress("slimmedJets.m_chargedEmEnergyFraction", slimmedJets_m_chargedEmEnergyFraction, &b_slimmedJets_m_chargedEmEnergyFraction);
  fChain->SetBranchAddress("slimmedJets.m_chargedHadronEnergyFraction", slimmedJets_m_chargedHadronEnergyFraction, &b_slimmedJets_m_chargedHadronEnergyFraction);
  fChain->SetBranchAddress("slimmedJets.m_muonEnergyFraction", slimmedJets_m_muonEnergyFraction, &b_slimmedJets_m_muonEnergyFraction);
  fChain->SetBranchAddress("slimmedJets.m_photonEnergyFraction", slimmedJets_m_photonEnergyFraction, &b_slimmedJets_m_photonEnergyFraction);
  fChain->SetBranchAddress("slimmedJets.m_chargedMultiplicity", slimmedJets_m_chargedMultiplicity, &b_slimmedJets_m_chargedMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_neutralMultiplicity", slimmedJets_m_neutralMultiplicity, &b_slimmedJets_m_neutralMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_muonMultiplicity", slimmedJets_m_muonMultiplicity, &b_slimmedJets_m_muonMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_electronMultiplicity", slimmedJets_m_electronMultiplicity, &b_slimmedJets_m_electronMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_photonMultiplicity", slimmedJets_m_photonMultiplicity, &b_slimmedJets_m_photonMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_puppiMultiplicity", slimmedJets_m_puppiMultiplicity, &b_slimmedJets_m_puppiMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_neutralPuppiMultiplicity", slimmedJets_m_neutralPuppiMultiplicity, &b_slimmedJets_m_neutralPuppiMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_neutralHadronPuppiMultiplicity", slimmedJets_m_neutralHadronPuppiMultiplicity, &b_slimmedJets_m_neutralHadronPuppiMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_photonPuppiMultiplicity", slimmedJets_m_photonPuppiMultiplicity, &b_slimmedJets_m_photonPuppiMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_HFHadronPuppiMultiplicity", slimmedJets_m_HFHadronPuppiMultiplicity, &b_slimmedJets_m_HFHadronPuppiMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_HFEMPuppiMultiplicity", slimmedJets_m_HFEMPuppiMultiplicity, &b_slimmedJets_m_HFEMPuppiMultiplicity);
  fChain->SetBranchAddress("slimmedJets.m_btag_combinedSecondaryVertex", slimmedJets_m_btag_combinedSecondaryVertex, &b_slimmedJets_m_btag_combinedSecondaryVertex);
  fChain->SetBranchAddress("slimmedJets.m_btag_combinedSecondaryVertexMVA", slimmedJets_m_btag_combinedSecondaryVertexMVA, &b_slimmedJets_m_btag_combinedSecondaryVertexMVA);
  fChain->SetBranchAddress("slimmedJets.m_btag_DeepCSV_probb", slimmedJets_m_btag_DeepCSV_probb, &b_slimmedJets_m_btag_DeepCSV_probb);
  fChain->SetBranchAddress("slimmedJets.m_btag_DeepCSV_probbb", slimmedJets_m_btag_DeepCSV_probbb, &b_slimmedJets_m_btag_DeepCSV_probbb);
  fChain->SetBranchAddress("slimmedJets.m_btag_BoostedDoubleSecondaryVertexAK8", slimmedJets_m_btag_BoostedDoubleSecondaryVertexAK8, &b_slimmedJets_m_btag_BoostedDoubleSecondaryVertexAK8);
  fChain->SetBranchAddress("slimmedJets.m_btag_BoostedDoubleSecondaryVertexCA15", slimmedJets_m_btag_BoostedDoubleSecondaryVertexCA15, &b_slimmedJets_m_btag_BoostedDoubleSecondaryVertexCA15);
  fChain->SetBranchAddress("slimmedJets.m_JEC_factor_raw", slimmedJets_m_JEC_factor_raw, &b_slimmedJets_m_JEC_factor_raw);
  fChain->SetBranchAddress("slimmedJets.m_JER_factor_raw", slimmedJets_m_JER_factor_raw, &b_slimmedJets_m_JER_factor_raw);
  fChain->SetBranchAddress("slimmedJets.m_JEC_L1factor_raw", slimmedJets_m_JEC_L1factor_raw, &b_slimmedJets_m_JEC_L1factor_raw);
  fChain->SetBranchAddress("slimmedJets.m_genjet_index", slimmedJets_m_genjet_index, &b_slimmedJets_m_genjet_index);
  fChain->SetBranchAddress("slimmedJets.m_hadronFlavor", slimmedJets_m_hadronFlavor, &b_slimmedJets_m_hadronFlavor);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackMomentum", slimmedJets_m_btaginfo_m_TrackMomentum, &b_slimmedJets_m_btaginfo_m_TrackMomentum);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackEta", slimmedJets_m_btaginfo_m_TrackEta, &b_slimmedJets_m_btaginfo_m_TrackEta);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackEtaRel", slimmedJets_m_btaginfo_m_TrackEtaRel, &b_slimmedJets_m_btaginfo_m_TrackEtaRel);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackDeltaR", slimmedJets_m_btaginfo_m_TrackDeltaR, &b_slimmedJets_m_btaginfo_m_TrackDeltaR);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackSip3dVal", slimmedJets_m_btaginfo_m_TrackSip3dVal, &b_slimmedJets_m_btaginfo_m_TrackSip3dVal);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackSip3dSig", slimmedJets_m_btaginfo_m_TrackSip3dSig, &b_slimmedJets_m_btaginfo_m_TrackSip3dSig);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackSip2dVal", slimmedJets_m_btaginfo_m_TrackSip2dVal, &b_slimmedJets_m_btaginfo_m_TrackSip2dVal);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackSip2dSig", slimmedJets_m_btaginfo_m_TrackSip2dSig, &b_slimmedJets_m_btaginfo_m_TrackSip2dSig);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackDecayLenVal", slimmedJets_m_btaginfo_m_TrackDecayLenVal, &b_slimmedJets_m_btaginfo_m_TrackDecayLenVal);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackChi2", slimmedJets_m_btaginfo_m_TrackChi2, &b_slimmedJets_m_btaginfo_m_TrackChi2);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackNTotalHits", slimmedJets_m_btaginfo_m_TrackNTotalHits, &b_slimmedJets_m_btaginfo_m_TrackNTotalHits);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackNPixelHits", slimmedJets_m_btaginfo_m_TrackNPixelHits, &b_slimmedJets_m_btaginfo_m_TrackNPixelHits);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackPtRel", slimmedJets_m_btaginfo_m_TrackPtRel, &b_slimmedJets_m_btaginfo_m_TrackPtRel);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackPPar", slimmedJets_m_btaginfo_m_TrackPPar, &b_slimmedJets_m_btaginfo_m_TrackPPar);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackPtRatio", slimmedJets_m_btaginfo_m_TrackPtRatio, &b_slimmedJets_m_btaginfo_m_TrackPtRatio);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackPParRatio", slimmedJets_m_btaginfo_m_TrackPParRatio, &b_slimmedJets_m_btaginfo_m_TrackPParRatio);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackJetDistVal", slimmedJets_m_btaginfo_m_TrackJetDistVal, &b_slimmedJets_m_btaginfo_m_TrackJetDistVal);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackJetDistSig", slimmedJets_m_btaginfo_m_TrackJetDistSig, &b_slimmedJets_m_btaginfo_m_TrackJetDistSig);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackGhostTrackDistVal", slimmedJets_m_btaginfo_m_TrackGhostTrackDistVal, &b_slimmedJets_m_btaginfo_m_TrackGhostTrackDistVal);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackGhostTrackDistSig", slimmedJets_m_btaginfo_m_TrackGhostTrackDistSig, &b_slimmedJets_m_btaginfo_m_TrackGhostTrackDistSig);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackGhostTrackWeight", slimmedJets_m_btaginfo_m_TrackGhostTrackWeight, &b_slimmedJets_m_btaginfo_m_TrackGhostTrackWeight);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_FlightDistance2dVal", slimmedJets_m_btaginfo_m_FlightDistance2dVal, &b_slimmedJets_m_btaginfo_m_FlightDistance2dVal);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_FlightDistance2dSig", slimmedJets_m_btaginfo_m_FlightDistance2dSig, &b_slimmedJets_m_btaginfo_m_FlightDistance2dSig);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_FlightDistance3dVal", slimmedJets_m_btaginfo_m_FlightDistance3dVal, &b_slimmedJets_m_btaginfo_m_FlightDistance3dVal);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_FlightDistance3dSig", slimmedJets_m_btaginfo_m_FlightDistance3dSig, &b_slimmedJets_m_btaginfo_m_FlightDistance3dSig);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexJetDeltaR", slimmedJets_m_btaginfo_m_VertexJetDeltaR, &b_slimmedJets_m_btaginfo_m_VertexJetDeltaR);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_JetNSecondaryVertices", slimmedJets_m_btaginfo_m_JetNSecondaryVertices, &b_slimmedJets_m_btaginfo_m_JetNSecondaryVertices);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexNTracks", slimmedJets_m_btaginfo_m_VertexNTracks, &b_slimmedJets_m_btaginfo_m_VertexNTracks);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_SecondaryVertex", slimmedJets_m_btaginfo_m_SecondaryVertex, &b_slimmedJets_m_btaginfo_m_SecondaryVertex);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexChi2", slimmedJets_m_btaginfo_m_VertexChi2, &b_slimmedJets_m_btaginfo_m_VertexChi2);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexNdof", slimmedJets_m_btaginfo_m_VertexNdof, &b_slimmedJets_m_btaginfo_m_VertexNdof);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexNormalizedChi2", slimmedJets_m_btaginfo_m_VertexNormalizedChi2, &b_slimmedJets_m_btaginfo_m_VertexNormalizedChi2);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexCategoryJTC", slimmedJets_m_btaginfo_m_VertexCategoryJTC, &b_slimmedJets_m_btaginfo_m_VertexCategoryJTC);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexMassJTC", slimmedJets_m_btaginfo_m_VertexMassJTC, &b_slimmedJets_m_btaginfo_m_VertexMassJTC);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_VertexEnergyRatioJTC", slimmedJets_m_btaginfo_m_VertexEnergyRatioJTC, &b_slimmedJets_m_btaginfo_m_VertexEnergyRatioJTC);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.m_TrackSip3dSigAboveCharmJTC", slimmedJets_m_btaginfo_m_TrackSip3dSigAboveCharmJTC, &b_slimmedJets_m_btaginfo_m_TrackSip3dSigAboveCharmJTC);
  fChain->SetBranchAddress("slimmedJets.m_btaginfo.tags.tagdata", slimmedJets_m_btaginfo_tags_tagdata, &b_slimmedJets_m_btaginfo_tags_tagdata);
  fChain->SetBranchAddress("slimmedJets.m_lepton_keys", slimmedJets_m_lepton_keys, &b_slimmedJets_m_lepton_keys);
  fChain->SetBranchAddress("slimmedJets.tags.tagdata", slimmedJets_tags_tagdata, &b_slimmedJets_tags_tagdata);
  fChain->SetBranchAddress("m_pt", &m_pt, &b_slimmedMETs_m_pt);
  fChain->SetBranchAddress("m_phi", &m_phi, &b_slimmedMETs_m_phi);
  fChain->SetBranchAddress("m_mEtSig", &m_mEtSig, &b_slimmedMETs_m_mEtSig);
  fChain->SetBranchAddress("m_shiftedPx_JetEnUp", &m_shiftedPx_JetEnUp, &b_slimmedMETs_m_shiftedPx_JetEnUp);
  fChain->SetBranchAddress("m_shiftedPx_JetEnDown", &m_shiftedPx_JetEnDown, &b_slimmedMETs_m_shiftedPx_JetEnDown);
  fChain->SetBranchAddress("m_shiftedPx_JetResUp", &m_shiftedPx_JetResUp, &b_slimmedMETs_m_shiftedPx_JetResUp);
  fChain->SetBranchAddress("m_shiftedPx_JetResDown", &m_shiftedPx_JetResDown, &b_slimmedMETs_m_shiftedPx_JetResDown);
  fChain->SetBranchAddress("m_shiftedPx_UnclusteredEnUp", &m_shiftedPx_UnclusteredEnUp, &b_slimmedMETs_m_shiftedPx_UnclusteredEnUp);
  fChain->SetBranchAddress("m_shiftedPx_UnclusteredEnDown", &m_shiftedPx_UnclusteredEnDown, &b_slimmedMETs_m_shiftedPx_UnclusteredEnDown);
  fChain->SetBranchAddress("m_shiftedPx_ElectronEnUp", &m_shiftedPx_ElectronEnUp, &b_slimmedMETs_m_shiftedPx_ElectronEnUp);
  fChain->SetBranchAddress("m_shiftedPx_ElectronEnDown", &m_shiftedPx_ElectronEnDown, &b_slimmedMETs_m_shiftedPx_ElectronEnDown);
  fChain->SetBranchAddress("m_shiftedPx_TauEnUp", &m_shiftedPx_TauEnUp, &b_slimmedMETs_m_shiftedPx_TauEnUp);
  fChain->SetBranchAddress("m_shiftedPx_TauEnDown", &m_shiftedPx_TauEnDown, &b_slimmedMETs_m_shiftedPx_TauEnDown);
  fChain->SetBranchAddress("m_shiftedPx_MuonEnDown", &m_shiftedPx_MuonEnDown, &b_slimmedMETs_m_shiftedPx_MuonEnDown);
  fChain->SetBranchAddress("m_shiftedPx_MuonEnUp", &m_shiftedPx_MuonEnUp, &b_slimmedMETs_m_shiftedPx_MuonEnUp);
  fChain->SetBranchAddress("m_shiftedPy_JetEnUp", &m_shiftedPy_JetEnUp, &b_slimmedMETs_m_shiftedPy_JetEnUp);
  fChain->SetBranchAddress("m_shiftedPy_JetEnDown", &m_shiftedPy_JetEnDown, &b_slimmedMETs_m_shiftedPy_JetEnDown);
  fChain->SetBranchAddress("m_shiftedPy_JetResUp", &m_shiftedPy_JetResUp, &b_slimmedMETs_m_shiftedPy_JetResUp);
  fChain->SetBranchAddress("m_shiftedPy_JetResDown", &m_shiftedPy_JetResDown, &b_slimmedMETs_m_shiftedPy_JetResDown);
  fChain->SetBranchAddress("m_shiftedPy_UnclusteredEnUp", &m_shiftedPy_UnclusteredEnUp, &b_slimmedMETs_m_shiftedPy_UnclusteredEnUp);
  fChain->SetBranchAddress("m_shiftedPy_UnclusteredEnDown", &m_shiftedPy_UnclusteredEnDown, &b_slimmedMETs_m_shiftedPy_UnclusteredEnDown);
  fChain->SetBranchAddress("m_shiftedPy_ElectronEnUp", &m_shiftedPy_ElectronEnUp, &b_slimmedMETs_m_shiftedPy_ElectronEnUp);
  fChain->SetBranchAddress("m_shiftedPy_ElectronEnDown", &m_shiftedPy_ElectronEnDown, &b_slimmedMETs_m_shiftedPy_ElectronEnDown);
  fChain->SetBranchAddress("m_shiftedPy_TauEnUp", &m_shiftedPy_TauEnUp, &b_slimmedMETs_m_shiftedPy_TauEnUp);
  fChain->SetBranchAddress("m_shiftedPy_TauEnDown", &m_shiftedPy_TauEnDown, &b_slimmedMETs_m_shiftedPy_TauEnDown);
  fChain->SetBranchAddress("m_shiftedPy_MuonEnDown", &m_shiftedPy_MuonEnDown, &b_slimmedMETs_m_shiftedPy_MuonEnDown);
  fChain->SetBranchAddress("m_shiftedPy_MuonEnUp", &m_shiftedPy_MuonEnUp, &b_slimmedMETs_m_shiftedPy_MuonEnUp);
  fChain->SetBranchAddress("m_rawCHS_px", &m_rawCHS_px, &b_slimmedMETs_m_rawCHS_px);
  fChain->SetBranchAddress("m_rawCHS_py", &m_rawCHS_py, &b_slimmedMETs_m_rawCHS_py);
  fChain->SetBranchAddress("m_uncorr_pt", &m_uncorr_pt, &b_slimmedMETs_m_uncorr_pt);
  fChain->SetBranchAddress("m_uncorr_phi", &m_uncorr_phi, &b_slimmedMETs_m_uncorr_phi);
  fChain->SetBranchAddress("slimmedMuonsUSER", &slimmedMuonsUSER_, &b_slimmedMuonsUSER_);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_charge", &slimmedMuonsUSER_m_charge, &b_slimmedMuonsUSER_m_charge);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pt", &slimmedMuonsUSER_m_pt, &b_slimmedMuonsUSER_m_pt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_eta", &slimmedMuonsUSER_m_eta, &b_slimmedMuonsUSER_m_eta);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_phi", &slimmedMuonsUSER_m_phi, &b_slimmedMuonsUSER_m_phi);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_energy", &slimmedMuonsUSER_m_energy, &b_slimmedMuonsUSER_m_energy);
  fChain->SetBranchAddress("slimmedMuonsUSER.sel_bits", &slimmedMuonsUSER_sel_bits, &b_slimmedMuonsUSER_sel_bits);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dxy", &slimmedMuonsUSER_m_dxy, &b_slimmedMuonsUSER_m_dxy);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dxy_error", &slimmedMuonsUSER_m_dxy_error, &b_slimmedMuonsUSER_m_dxy_error);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dz", &slimmedMuonsUSER_m_dz, &b_slimmedMuonsUSER_m_dz);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dz_error", &slimmedMuonsUSER_m_dz_error, &b_slimmedMuonsUSER_m_dz_error);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_globalTrack_normalizedChi2", &slimmedMuonsUSER_m_globalTrack_normalizedChi2, &b_slimmedMuonsUSER_m_globalTrack_normalizedChi2);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_globalTrack_numberOfValidMuonHits", &slimmedMuonsUSER_m_globalTrack_numberOfValidMuonHits, &b_slimmedMuonsUSER_m_globalTrack_numberOfValidMuonHits);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_numberOfMatchedStations", &slimmedMuonsUSER_m_numberOfMatchedStations, &b_slimmedMuonsUSER_m_numberOfMatchedStations);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_innerTrack_trackerLayersWithMeasurement", &slimmedMuonsUSER_m_innerTrack_trackerLayersWithMeasurement, &b_slimmedMuonsUSER_m_innerTrack_trackerLayersWithMeasurement);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_innerTrack_numberOfValidPixelHits", &slimmedMuonsUSER_m_innerTrack_numberOfValidPixelHits, &b_slimmedMuonsUSER_m_innerTrack_numberOfValidPixelHits);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_innerTrack_validFraction", &slimmedMuonsUSER_m_innerTrack_validFraction, &b_slimmedMuonsUSER_m_innerTrack_validFraction);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_combinedQuality_chi2LocalPosition", &slimmedMuonsUSER_m_combinedQuality_chi2LocalPosition, &b_slimmedMuonsUSER_m_combinedQuality_chi2LocalPosition);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_combinedQuality_trkKink", &slimmedMuonsUSER_m_combinedQuality_trkKink, &b_slimmedMuonsUSER_m_combinedQuality_trkKink);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_segmentCompatibility", &slimmedMuonsUSER_m_segmentCompatibility, &b_slimmedMuonsUSER_m_segmentCompatibility);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumChargedHadronPt", &slimmedMuonsUSER_m_sumChargedHadronPt, &b_slimmedMuonsUSER_m_sumChargedHadronPt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumNeutralHadronEt", &slimmedMuonsUSER_m_sumNeutralHadronEt, &b_slimmedMuonsUSER_m_sumNeutralHadronEt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumPhotonEt", &slimmedMuonsUSER_m_sumPhotonEt, &b_slimmedMuonsUSER_m_sumPhotonEt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumPUPt", &slimmedMuonsUSER_m_sumPUPt, &b_slimmedMuonsUSER_m_sumPUPt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_CH", &slimmedMuonsUSER_m_pfMINIIso_CH, &b_slimmedMuonsUSER_m_pfMINIIso_CH);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_NH", &slimmedMuonsUSER_m_pfMINIIso_NH, &b_slimmedMuonsUSER_m_pfMINIIso_NH);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_Ph", &slimmedMuonsUSER_m_pfMINIIso_Ph, &b_slimmedMuonsUSER_m_pfMINIIso_Ph);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_PU", &slimmedMuonsUSER_m_pfMINIIso_PU, &b_slimmedMuonsUSER_m_pfMINIIso_PU);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_NH_pfwgt", &slimmedMuonsUSER_m_pfMINIIso_NH_pfwgt, &b_slimmedMuonsUSER_m_pfMINIIso_NH_pfwgt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_Ph_pfwgt", &slimmedMuonsUSER_m_pfMINIIso_Ph_pfwgt, &b_slimmedMuonsUSER_m_pfMINIIso_Ph_pfwgt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simType", &slimmedMuonsUSER_m_simType, &b_slimmedMuonsUSER_m_simType);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simFlavor", &slimmedMuonsUSER_m_simFlavor, &b_slimmedMuonsUSER_m_simFlavor);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simPdgId", &slimmedMuonsUSER_m_simPdgId, &b_slimmedMuonsUSER_m_simPdgId);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simMotherPdgId", &slimmedMuonsUSER_m_simMotherPdgId, &b_slimmedMuonsUSER_m_simMotherPdgId);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simHeaviestMotherFlavor", &slimmedMuonsUSER_m_simHeaviestMotherFlavor, &b_slimmedMuonsUSER_m_simHeaviestMotherFlavor);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackPt", &slimmedMuonsUSER_m_tunePTrackPt, &b_slimmedMuonsUSER_m_tunePTrackPt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackEta", &slimmedMuonsUSER_m_tunePTrackEta, &b_slimmedMuonsUSER_m_tunePTrackEta);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackPhi", &slimmedMuonsUSER_m_tunePTrackPhi, &b_slimmedMuonsUSER_m_tunePTrackPhi);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackType", &slimmedMuonsUSER_m_tunePTrackType, &b_slimmedMuonsUSER_m_tunePTrackType);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_source_candidates", &slimmedMuonsUSER_m_source_candidates, &b_slimmedMuonsUSER_m_source_candidates);
  fChain->SetBranchAddress("slimmedMuonsUSER.tags.tagdata", &slimmedMuonsUSER_tags_tagdata, &b_slimmedMuonsUSER_tags_tagdata);
  fChain->SetBranchAddress("slimmedTaus", &slimmedTaus_, &b_slimmedTaus_);
  fChain->SetBranchAddress("slimmedTaus.m_charge", slimmedTaus_m_charge, &b_slimmedTaus_m_charge);
  fChain->SetBranchAddress("slimmedTaus.m_pt", slimmedTaus_m_pt, &b_slimmedTaus_m_pt);
  fChain->SetBranchAddress("slimmedTaus.m_eta", slimmedTaus_m_eta, &b_slimmedTaus_m_eta);
  fChain->SetBranchAddress("slimmedTaus.m_phi", slimmedTaus_m_phi, &b_slimmedTaus_m_phi);
  fChain->SetBranchAddress("slimmedTaus.m_energy", slimmedTaus_m_energy, &b_slimmedTaus_m_energy);
  fChain->SetBranchAddress("slimmedTaus.id_bits", slimmedTaus_id_bits, &b_slimmedTaus_id_bits);
  fChain->SetBranchAddress("slimmedTaus.m_byCombinedIsolationDeltaBetaCorrRaw3Hits", slimmedTaus_m_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_slimmedTaus_m_byCombinedIsolationDeltaBetaCorrRaw3Hits);
  fChain->SetBranchAddress("slimmedTaus.m_byIsolationMVA3newDMwoLTraw", slimmedTaus_m_byIsolationMVA3newDMwoLTraw, &b_slimmedTaus_m_byIsolationMVA3newDMwoLTraw);
  fChain->SetBranchAddress("slimmedTaus.m_byIsolationMVArun2v1DBnewDMwLTraw", slimmedTaus_m_byIsolationMVArun2v1DBnewDMwLTraw, &b_slimmedTaus_m_byIsolationMVArun2v1DBnewDMwLTraw);
  fChain->SetBranchAddress("slimmedTaus.m_chargedIsoPtSum", slimmedTaus_m_chargedIsoPtSum, &b_slimmedTaus_m_chargedIsoPtSum);
  fChain->SetBranchAddress("slimmedTaus.m_neutralIsoPtSum", slimmedTaus_m_neutralIsoPtSum, &b_slimmedTaus_m_neutralIsoPtSum);
  fChain->SetBranchAddress("slimmedTaus.m_puCorrPtSum", slimmedTaus_m_puCorrPtSum, &b_slimmedTaus_m_puCorrPtSum);
  fChain->SetBranchAddress("slimmedTaus.m_neutralIsoPtSumWeight", slimmedTaus_m_neutralIsoPtSumWeight, &b_slimmedTaus_m_neutralIsoPtSumWeight);
  fChain->SetBranchAddress("slimmedTaus.m_footprintCorrection", slimmedTaus_m_footprintCorrection, &b_slimmedTaus_m_footprintCorrection);
  fChain->SetBranchAddress("slimmedTaus.m_photonPtSumOutsideSignalCone", slimmedTaus_m_photonPtSumOutsideSignalCone, &b_slimmedTaus_m_photonPtSumOutsideSignalCone);
  fChain->SetBranchAddress("slimmedTaus.m_decayMode", slimmedTaus_m_decayMode, &b_slimmedTaus_m_decayMode);
  fChain->SetBranchAddress("slimmedTaus.tags.tagdata", slimmedTaus_tags_tagdata, &b_slimmedTaus_tags_tagdata);
  fChain->SetBranchAddress("weight_pu", &weight_pu, &b_weight_pu);
  fChain->SetBranchAddress("weight_pu_down", &weight_pu_down, &b_weight_pu_down);
  fChain->SetBranchAddress("weight_pu_up", &weight_pu_up, &b_weight_pu_up);
  Notify();
}

Bool_t AnalysisTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void AnalysisTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t AnalysisTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef AnalysisTree_cxx
