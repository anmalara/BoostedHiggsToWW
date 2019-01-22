//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jan 20 12:34:42 2019 by ROOT version 6.10/09
// from TTree AnalysisTree/AnalysisTree
// found on file: origin_DY.root
//////////////////////////////////////////////////////////

#ifndef AnalysisTree_h
#define AnalysisTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "UHH2/core/include/Particle.h"
#include "vector"
#include "UHH2/core/include/GenTopJet.h"
#include "UHH2/core/include/MET.h"
#include "vector"
#include "UHH2/core/include/PrimaryVertex.h"
#include "UHH2/core/include/GenInfo.h"
#include "vector"
#include "UHH2/core/include/FlavorParticle.h"
#include "UHH2/core/include/GenParticle.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "UHH2/core/include/RecParticle.h"
#include "UHH2/core/include/Electron.h"
#include "UHH2/core/include/Tags.h"
#include "vector"
#include "UHH2/core/include/Muon.h"
#include "vector"
#include "UHH2/core/include/Tau.h"
#include "vector"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/JetBTagInfo.h"
#include "vector"
#include "UHH2/core/include/TopJet.h"

class AnalysisTree {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  TString outname;

  // Fixed size dimensions of array or collections stored in the TTree if any.
  static constexpr Int_t kMaxslimmedGenJets = 25;
  static constexpr Int_t kMaxslimmedGenJetsAK8 = 8;
  static constexpr Int_t kMaxak8GenJetsSoftDrop = 6;
  static constexpr Int_t kMaxofflineSlimmedPrimaryVertices = 202;
  static constexpr Int_t kMaxGenParticles = 9;
  static constexpr Int_t kMaxslimmedElectronsUSER = 13;
  static constexpr Int_t kMaxslimmedMuonsUSER = 29;
  static constexpr Int_t kMaxslimmedTaus = 11;
  static constexpr Int_t kMaxslimmedJets = 154;
  static constexpr Int_t kMaxupdatedPatJetsSlimmedJetsPuppi = 20;
  static constexpr Int_t kMaxupdatedPatJetsPatJetsAK8PFPUPPI = 21;
  static constexpr Int_t kMaxpatJetsAK8PFCHS = 49;
  static constexpr Int_t kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi = 6;
  static constexpr Int_t kMaxpackedPatJetsAk8CHSJets_SoftDropCHS = 15;
  static constexpr Int_t kMaxpatJetsHepTopTagCHSPacked_daughters = 5;

  // Declaration of leaf types
  Int_t           run;
  Long64_t        event;
  Int_t           luminosityBlock;
  Bool_t          isRealData;
  Float_t         rho;
  Float_t         beamspot_x0;
  Float_t         beamspot_y0;
  Float_t         beamspot_z0;
  Int_t           slimmedGenJets_;
  Short_t         slimmedGenJets_m_charge[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_pt[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_eta[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_phi[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Float_t         slimmedGenJets_m_energy[kMaxslimmedGenJets];   //[slimmedGenJets_]
  Int_t           slimmedGenJetsAK8_;
  Short_t         slimmedGenJetsAK8_m_charge[kMaxslimmedGenJetsAK8];   //[slimmedGenJetsAK8_]
  Float_t         slimmedGenJetsAK8_m_pt[kMaxslimmedGenJetsAK8];   //[slimmedGenJetsAK8_]
  Float_t         slimmedGenJetsAK8_m_eta[kMaxslimmedGenJetsAK8];   //[slimmedGenJetsAK8_]
  Float_t         slimmedGenJetsAK8_m_phi[kMaxslimmedGenJetsAK8];   //[slimmedGenJetsAK8_]
  Float_t         slimmedGenJetsAK8_m_energy[kMaxslimmedGenJetsAK8];   //[slimmedGenJetsAK8_]
  Int_t           ak8GenJetsSoftDrop_;
  Short_t         ak8GenJetsSoftDrop_m_charge[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Float_t         ak8GenJetsSoftDrop_m_pt[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Float_t         ak8GenJetsSoftDrop_m_eta[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Float_t         ak8GenJetsSoftDrop_m_phi[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Float_t         ak8GenJetsSoftDrop_m_energy[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  vector<Particle> ak8GenJetsSoftDrop_m_subjets[kMaxak8GenJetsSoftDrop];
  Double_t        ak8GenJetsSoftDrop_m_tau1[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Double_t        ak8GenJetsSoftDrop_m_tau2[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Double_t        ak8GenJetsSoftDrop_m_tau3[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Double_t        ak8GenJetsSoftDrop_m_chf[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Double_t        ak8GenJetsSoftDrop_m_cef[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Double_t        ak8GenJetsSoftDrop_m_nhf[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
  Double_t        ak8GenJetsSoftDrop_m_nef[kMaxak8GenJetsSoftDrop];   //[ak8GenJetsSoftDrop_]
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
  //MET             *slimmedMETsPuppi;
  // Float_t         m_pt;
  // Float_t         m_phi;
  // Float_t         m_mEtSig;
  // Float_t         m_shiftedPx_JetEnUp;
  // Float_t         m_shiftedPx_JetEnDown;
  // Float_t         m_shiftedPx_JetResUp;
  // Float_t         m_shiftedPx_JetResDown;
  // Float_t         m_shiftedPx_UnclusteredEnUp;
  // Float_t         m_shiftedPx_UnclusteredEnDown;
  // Float_t         m_shiftedPx_ElectronEnUp;
  // Float_t         m_shiftedPx_ElectronEnDown;
  // Float_t         m_shiftedPx_TauEnUp;
  // Float_t         m_shiftedPx_TauEnDown;
  // Float_t         m_shiftedPx_MuonEnDown;
  // Float_t         m_shiftedPx_MuonEnUp;
  // Float_t         m_shiftedPy_JetEnUp;
  // Float_t         m_shiftedPy_JetEnDown;
  // Float_t         m_shiftedPy_JetResUp;
  // Float_t         m_shiftedPy_JetResDown;
  // Float_t         m_shiftedPy_UnclusteredEnUp;
  // Float_t         m_shiftedPy_UnclusteredEnDown;
  // Float_t         m_shiftedPy_ElectronEnUp;
  // Float_t         m_shiftedPy_ElectronEnDown;
  // Float_t         m_shiftedPy_TauEnUp;
  // Float_t         m_shiftedPy_TauEnDown;
  // Float_t         m_shiftedPy_MuonEnDown;
  // Float_t         m_shiftedPy_MuonEnUp;
  // Float_t         m_rawCHS_px;
  // Float_t         m_rawCHS_py;
  // Float_t         m_uncorr_pt;
  // Float_t         m_uncorr_phi;
  Int_t           offlineSlimmedPrimaryVertices_;
  Float_t         offlineSlimmedPrimaryVertices_m_x[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_y[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_z[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  UInt_t          offlineSlimmedPrimaryVertices_m_nTracks[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_chi2[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
  Float_t         offlineSlimmedPrimaryVertices_m_ndof[kMaxofflineSlimmedPrimaryVertices];   //[offlineSlimmedPrimaryVertices_]
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
  vector<string>  *triggerNames;
  vector<bool>    *triggerResults;
  vector<int>     *triggerPrescales;
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
  //vector<TLorentzVector> slimmedJets_m_btaginfo_m_SecondaryVertex[kMaxslimmedJets];
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
  Int_t           updatedPatJetsSlimmedJetsPuppi_;
  Short_t         updatedPatJetsSlimmedJetsPuppi_m_charge[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_pt[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_eta[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_phi[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_energy[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_pdgId[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_jetArea[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_numberOfDaughters[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_neutralEmEnergyFraction[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_neutralHadronEnergyFraction[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_chargedEmEnergyFraction[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_chargedHadronEnergyFraction[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_muonEnergyFraction[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_photonEnergyFraction[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_chargedMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_neutralMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_muonMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_electronMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_photonMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_puppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_neutralPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_neutralHadronPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_photonPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_HFHadronPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_HFEMPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertex[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertexMVA[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probb[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probbb[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexAK8[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexCA15[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_JEC_factor_raw[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_JEC_L1factor_raw[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_genjet_index[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_hadronFlavor[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackMomentum[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEta[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEtaRel[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDeltaR[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dVal[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSig[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dVal[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dSig[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDecayLenVal[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackChi2[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNTotalHits[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNPixelHits[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRel[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPPar[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRatio[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPParRatio[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistVal[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistSig[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistVal[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistSig[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackWeight[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dVal[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dSig[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dVal[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dSig[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexJetDeltaR[kMaxupdatedPatJetsSlimmedJetsPuppi];
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_JetNSecondaryVertices[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  vector<int>     updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNTracks[kMaxupdatedPatJetsSlimmedJetsPuppi];
  //vector<TLorentzVector> updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_SecondaryVertex[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexChi2[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNdof[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<float>   updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNormalizedChi2[kMaxupdatedPatJetsSlimmedJetsPuppi];
  Int_t           updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexCategoryJTC[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexMassJTC[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexEnergyRatioJTC[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  Float_t         updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC[kMaxupdatedPatJetsSlimmedJetsPuppi];   //[updatedPatJetsSlimmedJetsPuppi_]
  map<int,float>  updatedPatJetsSlimmedJetsPuppi_m_btaginfo_tags_tagdata[kMaxupdatedPatJetsSlimmedJetsPuppi];
  vector<long>    updatedPatJetsSlimmedJetsPuppi_m_lepton_keys[kMaxupdatedPatJetsSlimmedJetsPuppi];
  map<int,float>  updatedPatJetsSlimmedJetsPuppi_tags_tagdata[kMaxupdatedPatJetsSlimmedJetsPuppi];
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_;
  Short_t         updatedPatJetsPatJetsAK8PFPUPPI_m_charge[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_pt[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_eta[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_phi[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_energy[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_pdgId[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_jetArea[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_numberOfDaughters[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_neutralEmEnergyFraction[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronEnergyFraction[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_chargedEmEnergyFraction[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_chargedHadronEnergyFraction[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_muonEnergyFraction[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_photonEnergyFraction[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_chargedMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_neutralMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_muonMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_electronMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_photonMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_puppiMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_neutralPuppiMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronPuppiMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_photonPuppiMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_HFHadronPuppiMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_HFEMPuppiMultiplicity[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertex[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertexMVA[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probb[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probbb[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexAK8[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexCA15[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_factor_raw[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_L1factor_raw[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_genjet_index[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_hadronFlavor[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackMomentum[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEta[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEtaRel[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDeltaR[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dVal[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSig[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dVal[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dSig[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDecayLenVal[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackChi2[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNTotalHits[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNPixelHits[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRel[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPPar[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRatio[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPParRatio[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistVal[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistSig[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistVal[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistSig[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackWeight[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dVal[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dSig[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dVal[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dSig[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexJetDeltaR[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_JetNSecondaryVertices[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  vector<int>     updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNTracks[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  //vector<TLorentzVector> updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_SecondaryVertex[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexChi2[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNdof[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<float>   updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNormalizedChi2[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  Int_t           updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexCategoryJTC[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexMassJTC[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexEnergyRatioJTC[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  Float_t         updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSigAboveCharmJTC[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];   //[updatedPatJetsPatJetsAK8PFPUPPI_]
  map<int,float>  updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_tags_tagdata[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  vector<long>    updatedPatJetsPatJetsAK8PFPUPPI_m_lepton_keys[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  map<int,float>  updatedPatJetsPatJetsAK8PFPUPPI_tags_tagdata[kMaxupdatedPatJetsPatJetsAK8PFPUPPI];
  Int_t           patJetsAK8PFCHS_;
  Short_t         patJetsAK8PFCHS_m_charge[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_pt[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_eta[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_phi[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_energy[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_pdgId[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_jetArea[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_numberOfDaughters[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_neutralEmEnergyFraction[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_neutralHadronEnergyFraction[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_chargedEmEnergyFraction[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_chargedHadronEnergyFraction[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_muonEnergyFraction[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_photonEnergyFraction[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_chargedMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_neutralMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_muonMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_electronMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_photonMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_puppiMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_neutralPuppiMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_neutralHadronPuppiMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_photonPuppiMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_HFHadronPuppiMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_HFEMPuppiMultiplicity[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btag_combinedSecondaryVertex[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btag_combinedSecondaryVertexMVA[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btag_DeepCSV_probb[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btag_DeepCSV_probbb[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexAK8[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexCA15[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_JEC_factor_raw[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_JEC_L1factor_raw[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_genjet_index[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Int_t           patJetsAK8PFCHS_m_hadronFlavor[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackMomentum[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackEta[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackEtaRel[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackDeltaR[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dVal[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSig[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dVal[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dSig[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackDecayLenVal[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackChi2[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackNTotalHits[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackNPixelHits[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackPtRel[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackPPar[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackPtRatio[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackPParRatio[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistVal[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistSig[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistVal[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistSig[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackWeight[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dVal[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dSig[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dVal[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dSig[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_VertexJetDeltaR[kMaxpatJetsAK8PFCHS];
  Int_t           patJetsAK8PFCHS_m_btaginfo_m_JetNSecondaryVertices[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  vector<int>     patJetsAK8PFCHS_m_btaginfo_m_VertexNTracks[kMaxpatJetsAK8PFCHS];
  //vector<TLorentzVector> patJetsAK8PFCHS_m_btaginfo_m_SecondaryVertex[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_VertexChi2[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_VertexNdof[kMaxpatJetsAK8PFCHS];
  vector<float>   patJetsAK8PFCHS_m_btaginfo_m_VertexNormalizedChi2[kMaxpatJetsAK8PFCHS];
  Int_t           patJetsAK8PFCHS_m_btaginfo_m_VertexCategoryJTC[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btaginfo_m_VertexMassJTC[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btaginfo_m_VertexEnergyRatioJTC[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  Float_t         patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC[kMaxpatJetsAK8PFCHS];   //[patJetsAK8PFCHS_]
  map<int,float>  patJetsAK8PFCHS_m_btaginfo_tags_tagdata[kMaxpatJetsAK8PFCHS];
  vector<long>    patJetsAK8PFCHS_m_lepton_keys[kMaxpatJetsAK8PFCHS];
  map<int,float>  patJetsAK8PFCHS_tags_tagdata[kMaxpatJetsAK8PFCHS];
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_;
  Short_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_charge[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pt[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_eta[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_phi[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_energy[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pdgId[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_jetArea[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_numberOfDaughters[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralEmEnergyFraction[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronEnergyFraction[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedEmEnergyFraction[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedHadronEnergyFraction[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonEnergyFraction[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonEnergyFraction[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_electronMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_puppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFHadronPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFEMPuppiMultiplicity[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertex[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertexMVA[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probb[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probbb[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexAK8[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexCA15[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_factor_raw[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_L1factor_raw[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_genjet_index[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_hadronFlavor[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackMomentum[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEta[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEtaRel[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDeltaR[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dVal[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSig[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dVal[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dSig[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDecayLenVal[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackChi2[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNTotalHits[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNPixelHits[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRel[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPPar[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRatio[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPParRatio[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistVal[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistSig[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistVal[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistSig[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackWeight[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dVal[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dSig[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dVal[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dSig[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexJetDeltaR[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_JetNSecondaryVertices[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  vector<int>     updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNTracks[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  //vector<TLorentzVector> updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_SecondaryVertex[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexChi2[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNdof[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<float>   updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNormalizedChi2[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  Int_t           updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexCategoryJTC[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexMassJTC[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexEnergyRatioJTC[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  map<int,float>  updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_tags_tagdata[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<long>    updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_lepton_keys[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  map<int,float>  updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  vector<Jet>     updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_subjets[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_qjets_volatility[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1_groomed[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2_groomed[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3_groomed[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4_groomed[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta1[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta2[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta1[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta2[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_mvahiggsdiscr[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_prunedmass[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  Float_t         updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_softdropmass[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];   //[updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_]
  // map<int,float>  updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata[kMaxupdatedPatJetsSlimmedJetsAK8_SoftDropPuppi];
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
  //vector<TLorentzVector> packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_SecondaryVertex[kMaxpackedPatJetsAk8CHSJets_SoftDropCHS];
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
  Int_t           patJetsHepTopTagCHSPacked_daughters_;
  Short_t         patJetsHepTopTagCHSPacked_daughters_m_charge[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_pt[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_eta[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_phi[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_energy[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_pdgId[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_jetArea[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_numberOfDaughters[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_neutralEmEnergyFraction[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_neutralHadronEnergyFraction[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_chargedEmEnergyFraction[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_chargedHadronEnergyFraction[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_muonEnergyFraction[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_photonEnergyFraction[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_chargedMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_neutralMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_muonMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_electronMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_photonMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_puppiMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_neutralPuppiMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_neutralHadronPuppiMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_photonPuppiMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_HFHadronPuppiMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_HFEMPuppiMultiplicity[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertex[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertexMVA[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probb[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probbb[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexAK8[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexCA15[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_JEC_factor_raw[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_JEC_L1factor_raw[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_genjet_index[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_hadronFlavor[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackMomentum[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEta[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEtaRel[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDeltaR[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dVal[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSig[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dVal[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dSig[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDecayLenVal[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackChi2[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNTotalHits[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNPixelHits[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRel[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPPar[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRatio[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPParRatio[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistVal[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistSig[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistVal[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistSig[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackWeight[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dVal[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dSig[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dVal[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dSig[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexJetDeltaR[kMaxpatJetsHepTopTagCHSPacked_daughters];
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_JetNSecondaryVertices[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  vector<int>     patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNTracks[kMaxpatJetsHepTopTagCHSPacked_daughters];
  //vector<TLorentzVector> patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_SecondaryVertex[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexChi2[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNdof[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<float>   patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNormalizedChi2[kMaxpatJetsHepTopTagCHSPacked_daughters];
  Int_t           patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexCategoryJTC[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexMassJTC[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexEnergyRatioJTC[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSigAboveCharmJTC[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  map<int,float>  patJetsHepTopTagCHSPacked_daughters_m_btaginfo_tags_tagdata[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<long>    patJetsHepTopTagCHSPacked_daughters_m_lepton_keys[kMaxpatJetsHepTopTagCHSPacked_daughters];
  map<int,float>  patJetsHepTopTagCHSPacked_daughters_tags_tagdata[kMaxpatJetsHepTopTagCHSPacked_daughters];
  vector<Jet>     patJetsHepTopTagCHSPacked_daughters_m_subjets[kMaxpatJetsHepTopTagCHSPacked_daughters];
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_qjets_volatility[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau1[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau2[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau3[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau4[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau1_groomed[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau2_groomed[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau3_groomed[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_tau4_groomed[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta1[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta2[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta1[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta2[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_mvahiggsdiscr[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_prunedmass[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  Float_t         patJetsHepTopTagCHSPacked_daughters_m_softdropmass[kMaxpatJetsHepTopTagCHSPacked_daughters];   //[patJetsHepTopTagCHSPacked_daughters_]
  // map<int,float>  patJetsHepTopTagCHSPacked_daughters_tags_tagdata[kMaxpatJetsHepTopTagCHSPacked_daughters];

  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_event;   //!
  TBranch        *b_luminosityBlock;   //!
  TBranch        *b_isRealData;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_beamspot_x0;   //!
  TBranch        *b_beamspot_y0;   //!
  TBranch        *b_beamspot_z0;   //!
  TBranch        *b_slimmedGenJets_;   //!
  TBranch        *b_slimmedGenJets_m_charge;   //!
  TBranch        *b_slimmedGenJets_m_pt;   //!
  TBranch        *b_slimmedGenJets_m_eta;   //!
  TBranch        *b_slimmedGenJets_m_phi;   //!
  TBranch        *b_slimmedGenJets_m_energy;   //!
  TBranch        *b_slimmedGenJetsAK8_;   //!
  TBranch        *b_slimmedGenJetsAK8_m_charge;   //!
  TBranch        *b_slimmedGenJetsAK8_m_pt;   //!
  TBranch        *b_slimmedGenJetsAK8_m_eta;   //!
  TBranch        *b_slimmedGenJetsAK8_m_phi;   //!
  TBranch        *b_slimmedGenJetsAK8_m_energy;   //!
  TBranch        *b_ak8GenJetsSoftDrop_;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_charge;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_pt;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_eta;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_phi;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_energy;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_subjets;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_tau1;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_tau2;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_tau3;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_chf;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_cef;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_nhf;   //!
  TBranch        *b_ak8GenJetsSoftDrop_m_nef;   //!
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
  TBranch        *b_slimmedMETsPuppi_m_pt;   //!
  TBranch        *b_slimmedMETsPuppi_m_phi;   //!
  TBranch        *b_slimmedMETsPuppi_m_mEtSig;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_JetEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_JetEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_JetResUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_JetResDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_UnclusteredEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_UnclusteredEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_ElectronEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_ElectronEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_TauEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_TauEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_MuonEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPx_MuonEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_JetEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_JetEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_JetResUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_JetResDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_UnclusteredEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_UnclusteredEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_ElectronEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_ElectronEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_TauEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_TauEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_MuonEnDown;   //!
  TBranch        *b_slimmedMETsPuppi_m_shiftedPy_MuonEnUp;   //!
  TBranch        *b_slimmedMETsPuppi_m_rawCHS_px;   //!
  TBranch        *b_slimmedMETsPuppi_m_rawCHS_py;   //!
  TBranch        *b_slimmedMETsPuppi_m_uncorr_pt;   //!
  TBranch        *b_slimmedMETsPuppi_m_uncorr_phi;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_x;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_y;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_z;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_nTracks;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_chi2;   //!
  TBranch        *b_offlineSlimmedPrimaryVertices_m_ndof;   //!
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
  TBranch        *b_triggerNames;   //!
  TBranch        *b_triggerResults;   //!
  TBranch        *b_triggerPrescales;   //!
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
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_charge;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_pt;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_eta;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_phi;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_energy;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_pdgId;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_jetArea;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_numberOfDaughters;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_neutralEmEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_neutralHadronEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_chargedEmEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_chargedHadronEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_muonEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_photonEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_chargedMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_neutralMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_muonMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_electronMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_photonMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_puppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_neutralPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_neutralHadronPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_photonPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_HFHadronPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_HFEMPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertex;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertexMVA;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probb;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probbb;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexAK8;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexCA15;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_JEC_factor_raw;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_JEC_L1factor_raw;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_genjet_index;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_hadronFlavor;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackMomentum;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEta;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEtaRel;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDeltaR;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDecayLenVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackChi2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNTotalHits;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNPixelHits;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRel;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPPar;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRatio;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPParRatio;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackWeight;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexJetDeltaR;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_JetNSecondaryVertices;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNTracks;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_m_lepton_keys;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsPuppi_tags_tagdata;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_charge;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_pt;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_eta;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_phi;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_energy;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_pdgId;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_jetArea;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_numberOfDaughters;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralEmEnergyFraction;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronEnergyFraction;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_chargedEmEnergyFraction;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_chargedHadronEnergyFraction;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_muonEnergyFraction;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_photonEnergyFraction;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_chargedMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_muonMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_electronMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_photonMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_puppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_photonPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_HFHadronPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_HFEMPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertex;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertexMVA;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probb;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probbb;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexAK8;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexCA15;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_factor_raw;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_L1factor_raw;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_genjet_index;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_hadronFlavor;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackMomentum;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEta;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEtaRel;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDeltaR;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dVal;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSig;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dVal;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dSig;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDecayLenVal;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackChi2;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNTotalHits;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNPixelHits;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRel;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPPar;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRatio;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPParRatio;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistVal;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistSig;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistVal;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistSig;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackWeight;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dVal;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dSig;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dVal;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dSig;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexJetDeltaR;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_JetNSecondaryVertices;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNTracks;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_m_lepton_keys;   //!
  TBranch        *b_updatedPatJetsPatJetsAK8PFPUPPI_tags_tagdata;   //!
  TBranch        *b_patJetsAK8PFCHS_;   //!
  TBranch        *b_patJetsAK8PFCHS_m_charge;   //!
  TBranch        *b_patJetsAK8PFCHS_m_pt;   //!
  TBranch        *b_patJetsAK8PFCHS_m_eta;   //!
  TBranch        *b_patJetsAK8PFCHS_m_phi;   //!
  TBranch        *b_patJetsAK8PFCHS_m_energy;   //!
  TBranch        *b_patJetsAK8PFCHS_m_pdgId;   //!
  TBranch        *b_patJetsAK8PFCHS_m_jetArea;   //!
  TBranch        *b_patJetsAK8PFCHS_m_numberOfDaughters;   //!
  TBranch        *b_patJetsAK8PFCHS_m_neutralEmEnergyFraction;   //!
  TBranch        *b_patJetsAK8PFCHS_m_neutralHadronEnergyFraction;   //!
  TBranch        *b_patJetsAK8PFCHS_m_chargedEmEnergyFraction;   //!
  TBranch        *b_patJetsAK8PFCHS_m_chargedHadronEnergyFraction;   //!
  TBranch        *b_patJetsAK8PFCHS_m_muonEnergyFraction;   //!
  TBranch        *b_patJetsAK8PFCHS_m_photonEnergyFraction;   //!
  TBranch        *b_patJetsAK8PFCHS_m_chargedMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_neutralMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_muonMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_electronMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_photonMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_puppiMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_neutralPuppiMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_neutralHadronPuppiMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_photonPuppiMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_HFHadronPuppiMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_HFEMPuppiMultiplicity;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btag_combinedSecondaryVertex;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btag_combinedSecondaryVertexMVA;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btag_DeepCSV_probb;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btag_DeepCSV_probbb;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexAK8;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexCA15;   //!
  TBranch        *b_patJetsAK8PFCHS_m_JEC_factor_raw;   //!
  TBranch        *b_patJetsAK8PFCHS_m_JEC_L1factor_raw;   //!
  TBranch        *b_patJetsAK8PFCHS_m_genjet_index;   //!
  TBranch        *b_patJetsAK8PFCHS_m_hadronFlavor;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackMomentum;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackEta;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackEtaRel;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackDeltaR;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dVal;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSig;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dVal;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dSig;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackDecayLenVal;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackChi2;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackNTotalHits;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackNPixelHits;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackPtRel;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackPPar;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackPtRatio;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackPParRatio;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistVal;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistSig;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistVal;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistSig;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackWeight;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dVal;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dSig;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dVal;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dSig;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexJetDeltaR;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_JetNSecondaryVertices;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexNTracks;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_patJetsAK8PFCHS_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_patJetsAK8PFCHS_m_lepton_keys;   //!
  TBranch        *b_patJetsAK8PFCHS_tags_tagdata;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_charge;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pt;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_eta;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_phi;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_energy;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pdgId;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_jetArea;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_numberOfDaughters;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralEmEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedEmEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedHadronEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonEnergyFraction;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_electronMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_puppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFHadronPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFEMPuppiMultiplicity;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertex;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertexMVA;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probb;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probbb;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexAK8;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexCA15;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_factor_raw;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_L1factor_raw;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_genjet_index;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_hadronFlavor;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackMomentum;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEta;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEtaRel;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDeltaR;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDecayLenVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackChi2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNTotalHits;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNPixelHits;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRel;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPPar;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRatio;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPParRatio;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackWeight;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dVal;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dSig;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexJetDeltaR;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_JetNSecondaryVertices;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNTracks;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_lepton_keys;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_subjets;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_qjets_volatility;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1_groomed;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2_groomed;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3_groomed;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4_groomed;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta1;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta1;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta2;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_mvahiggsdiscr;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_prunedmass;   //!
  TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_softdropmass;   //!
  // TBranch        *b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata;   //!
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
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_m_lepton_keys;   //!
  TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata;   //!
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
  // TBranch        *b_packedPatJetsAk8CHSJets_SoftDropCHS_tags_tagdata;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_charge;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_pt;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_eta;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_phi;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_energy;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_pdgId;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_jetArea;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_numberOfDaughters;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_neutralEmEnergyFraction;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_neutralHadronEnergyFraction;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_chargedEmEnergyFraction;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_chargedHadronEnergyFraction;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_muonEnergyFraction;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_photonEnergyFraction;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_chargedMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_neutralMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_muonMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_electronMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_photonMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_puppiMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_neutralPuppiMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_neutralHadronPuppiMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_photonPuppiMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_HFHadronPuppiMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_HFEMPuppiMultiplicity;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertex;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertexMVA;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probb;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probbb;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexAK8;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexCA15;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_JEC_factor_raw;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_JEC_L1factor_raw;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_genjet_index;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_hadronFlavor;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackMomentum;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEta;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEtaRel;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDeltaR;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dVal;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSig;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dVal;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dSig;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDecayLenVal;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackChi2;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNTotalHits;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNPixelHits;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRel;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPPar;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRatio;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPParRatio;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistVal;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistSig;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistVal;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistSig;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackWeight;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dVal;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dSig;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dVal;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dSig;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexJetDeltaR;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_JetNSecondaryVertices;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNTracks;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexChi2;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNdof;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNormalizedChi2;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexCategoryJTC;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexMassJTC;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexEnergyRatioJTC;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSigAboveCharmJTC;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_tags_tagdata;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_lepton_keys;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_tags_tagdata;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_subjets;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_qjets_volatility;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau1;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau2;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau3;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau4;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau1_groomed;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau2_groomed;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau3_groomed;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_tau4_groomed;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta1;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta2;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta1;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta2;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_mvahiggsdiscr;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_prunedmass;   //!
  TBranch        *b_patJetsHepTopTagCHSPacked_daughters_m_softdropmass;   //!
  // TBranch        *b_patJetsHepTopTagCHSPacked_daughters_tags_tagdata;   //!

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
    f->GetObject("boostedAK8/events",tree);

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

  // Set object pointer
  triggerNames = 0;
  triggerResults = 0;
  triggerPrescales = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
  fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
  fChain->SetBranchAddress("rho", &rho, &b_rho);
  fChain->SetBranchAddress("beamspot_x0", &beamspot_x0, &b_beamspot_x0);
  fChain->SetBranchAddress("beamspot_y0", &beamspot_y0, &b_beamspot_y0);
  fChain->SetBranchAddress("beamspot_z0", &beamspot_z0, &b_beamspot_z0);
  fChain->SetBranchAddress("slimmedGenJets", &slimmedGenJets_, &b_slimmedGenJets_);
  fChain->SetBranchAddress("slimmedGenJets.m_charge", slimmedGenJets_m_charge, &b_slimmedGenJets_m_charge);
  fChain->SetBranchAddress("slimmedGenJets.m_pt", slimmedGenJets_m_pt, &b_slimmedGenJets_m_pt);
  fChain->SetBranchAddress("slimmedGenJets.m_eta", slimmedGenJets_m_eta, &b_slimmedGenJets_m_eta);
  fChain->SetBranchAddress("slimmedGenJets.m_phi", slimmedGenJets_m_phi, &b_slimmedGenJets_m_phi);
  fChain->SetBranchAddress("slimmedGenJets.m_energy", slimmedGenJets_m_energy, &b_slimmedGenJets_m_energy);
  fChain->SetBranchAddress("slimmedGenJetsAK8", &slimmedGenJetsAK8_, &b_slimmedGenJetsAK8_);
  fChain->SetBranchAddress("slimmedGenJetsAK8.m_charge", slimmedGenJetsAK8_m_charge, &b_slimmedGenJetsAK8_m_charge);
  fChain->SetBranchAddress("slimmedGenJetsAK8.m_pt", slimmedGenJetsAK8_m_pt, &b_slimmedGenJetsAK8_m_pt);
  fChain->SetBranchAddress("slimmedGenJetsAK8.m_eta", slimmedGenJetsAK8_m_eta, &b_slimmedGenJetsAK8_m_eta);
  fChain->SetBranchAddress("slimmedGenJetsAK8.m_phi", slimmedGenJetsAK8_m_phi, &b_slimmedGenJetsAK8_m_phi);
  fChain->SetBranchAddress("slimmedGenJetsAK8.m_energy", slimmedGenJetsAK8_m_energy, &b_slimmedGenJetsAK8_m_energy);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop", &ak8GenJetsSoftDrop_, &b_ak8GenJetsSoftDrop_);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_charge", ak8GenJetsSoftDrop_m_charge, &b_ak8GenJetsSoftDrop_m_charge);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_pt", ak8GenJetsSoftDrop_m_pt, &b_ak8GenJetsSoftDrop_m_pt);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_eta", ak8GenJetsSoftDrop_m_eta, &b_ak8GenJetsSoftDrop_m_eta);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_phi", ak8GenJetsSoftDrop_m_phi, &b_ak8GenJetsSoftDrop_m_phi);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_energy", ak8GenJetsSoftDrop_m_energy, &b_ak8GenJetsSoftDrop_m_energy);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_subjets", ak8GenJetsSoftDrop_m_subjets, &b_ak8GenJetsSoftDrop_m_subjets);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_tau1", ak8GenJetsSoftDrop_m_tau1, &b_ak8GenJetsSoftDrop_m_tau1);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_tau2", ak8GenJetsSoftDrop_m_tau2, &b_ak8GenJetsSoftDrop_m_tau2);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_tau3", ak8GenJetsSoftDrop_m_tau3, &b_ak8GenJetsSoftDrop_m_tau3);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_chf", ak8GenJetsSoftDrop_m_chf, &b_ak8GenJetsSoftDrop_m_chf);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_cef", ak8GenJetsSoftDrop_m_cef, &b_ak8GenJetsSoftDrop_m_cef);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_nhf", ak8GenJetsSoftDrop_m_nhf, &b_ak8GenJetsSoftDrop_m_nhf);
  fChain->SetBranchAddress("ak8GenJetsSoftDrop.m_nef", ak8GenJetsSoftDrop_m_nef, &b_ak8GenJetsSoftDrop_m_nef);
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
  //    fChain->SetBranchAddress("m_pt", &m_pt, &b_slimmedMETsPuppi_m_pt);
  //    fChain->SetBranchAddress("m_phi", &m_phi, &b_slimmedMETsPuppi_m_phi);
  //    fChain->SetBranchAddress("m_mEtSig", &m_mEtSig, &b_slimmedMETsPuppi_m_mEtSig);
  //    fChain->SetBranchAddress("m_shiftedPx_JetEnUp", &m_shiftedPx_JetEnUp, &b_slimmedMETsPuppi_m_shiftedPx_JetEnUp);
  //    fChain->SetBranchAddress("m_shiftedPx_JetEnDown", &m_shiftedPx_JetEnDown, &b_slimmedMETsPuppi_m_shiftedPx_JetEnDown);
  //    fChain->SetBranchAddress("m_shiftedPx_JetResUp", &m_shiftedPx_JetResUp, &b_slimmedMETsPuppi_m_shiftedPx_JetResUp);
  //    fChain->SetBranchAddress("m_shiftedPx_JetResDown", &m_shiftedPx_JetResDown, &b_slimmedMETsPuppi_m_shiftedPx_JetResDown);
  //    fChain->SetBranchAddress("m_shiftedPx_UnclusteredEnUp", &m_shiftedPx_UnclusteredEnUp, &b_slimmedMETsPuppi_m_shiftedPx_UnclusteredEnUp);
  //    fChain->SetBranchAddress("m_shiftedPx_UnclusteredEnDown", &m_shiftedPx_UnclusteredEnDown, &b_slimmedMETsPuppi_m_shiftedPx_UnclusteredEnDown);
  //    fChain->SetBranchAddress("m_shiftedPx_ElectronEnUp", &m_shiftedPx_ElectronEnUp, &b_slimmedMETsPuppi_m_shiftedPx_ElectronEnUp);
  //    fChain->SetBranchAddress("m_shiftedPx_ElectronEnDown", &m_shiftedPx_ElectronEnDown, &b_slimmedMETsPuppi_m_shiftedPx_ElectronEnDown);
  //    fChain->SetBranchAddress("m_shiftedPx_TauEnUp", &m_shiftedPx_TauEnUp, &b_slimmedMETsPuppi_m_shiftedPx_TauEnUp);
  //    fChain->SetBranchAddress("m_shiftedPx_TauEnDown", &m_shiftedPx_TauEnDown, &b_slimmedMETsPuppi_m_shiftedPx_TauEnDown);
  //    fChain->SetBranchAddress("m_shiftedPx_MuonEnDown", &m_shiftedPx_MuonEnDown, &b_slimmedMETsPuppi_m_shiftedPx_MuonEnDown);
  //    fChain->SetBranchAddress("m_shiftedPx_MuonEnUp", &m_shiftedPx_MuonEnUp, &b_slimmedMETsPuppi_m_shiftedPx_MuonEnUp);
  //    fChain->SetBranchAddress("m_shiftedPy_JetEnUp", &m_shiftedPy_JetEnUp, &b_slimmedMETsPuppi_m_shiftedPy_JetEnUp);
  //    fChain->SetBranchAddress("m_shiftedPy_JetEnDown", &m_shiftedPy_JetEnDown, &b_slimmedMETsPuppi_m_shiftedPy_JetEnDown);
  //    fChain->SetBranchAddress("m_shiftedPy_JetResUp", &m_shiftedPy_JetResUp, &b_slimmedMETsPuppi_m_shiftedPy_JetResUp);
  //    fChain->SetBranchAddress("m_shiftedPy_JetResDown", &m_shiftedPy_JetResDown, &b_slimmedMETsPuppi_m_shiftedPy_JetResDown);
  //    fChain->SetBranchAddress("m_shiftedPy_UnclusteredEnUp", &m_shiftedPy_UnclusteredEnUp, &b_slimmedMETsPuppi_m_shiftedPy_UnclusteredEnUp);
  //    fChain->SetBranchAddress("m_shiftedPy_UnclusteredEnDown", &m_shiftedPy_UnclusteredEnDown, &b_slimmedMETsPuppi_m_shiftedPy_UnclusteredEnDown);
  //    fChain->SetBranchAddress("m_shiftedPy_ElectronEnUp", &m_shiftedPy_ElectronEnUp, &b_slimmedMETsPuppi_m_shiftedPy_ElectronEnUp);
  //    fChain->SetBranchAddress("m_shiftedPy_ElectronEnDown", &m_shiftedPy_ElectronEnDown, &b_slimmedMETsPuppi_m_shiftedPy_ElectronEnDown);
  //    fChain->SetBranchAddress("m_shiftedPy_TauEnUp", &m_shiftedPy_TauEnUp, &b_slimmedMETsPuppi_m_shiftedPy_TauEnUp);
  //    fChain->SetBranchAddress("m_shiftedPy_TauEnDown", &m_shiftedPy_TauEnDown, &b_slimmedMETsPuppi_m_shiftedPy_TauEnDown);
  //    fChain->SetBranchAddress("m_shiftedPy_MuonEnDown", &m_shiftedPy_MuonEnDown, &b_slimmedMETsPuppi_m_shiftedPy_MuonEnDown);
  //    fChain->SetBranchAddress("m_shiftedPy_MuonEnUp", &m_shiftedPy_MuonEnUp, &b_slimmedMETsPuppi_m_shiftedPy_MuonEnUp);
  //    fChain->SetBranchAddress("m_rawCHS_px", &m_rawCHS_px, &b_slimmedMETsPuppi_m_rawCHS_px);
  //    fChain->SetBranchAddress("m_rawCHS_py", &m_rawCHS_py, &b_slimmedMETsPuppi_m_rawCHS_py);
  //    fChain->SetBranchAddress("m_uncorr_pt", &m_uncorr_pt, &b_slimmedMETsPuppi_m_uncorr_pt);
  //    fChain->SetBranchAddress("m_uncorr_phi", &m_uncorr_phi, &b_slimmedMETsPuppi_m_uncorr_phi);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices", &offlineSlimmedPrimaryVertices_, &b_offlineSlimmedPrimaryVertices_);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_x", offlineSlimmedPrimaryVertices_m_x, &b_offlineSlimmedPrimaryVertices_m_x);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_y", offlineSlimmedPrimaryVertices_m_y, &b_offlineSlimmedPrimaryVertices_m_y);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_z", offlineSlimmedPrimaryVertices_m_z, &b_offlineSlimmedPrimaryVertices_m_z);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_nTracks", offlineSlimmedPrimaryVertices_m_nTracks, &b_offlineSlimmedPrimaryVertices_m_nTracks);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_chi2", offlineSlimmedPrimaryVertices_m_chi2, &b_offlineSlimmedPrimaryVertices_m_chi2);
  fChain->SetBranchAddress("offlineSlimmedPrimaryVertices.m_ndof", offlineSlimmedPrimaryVertices_m_ndof, &b_offlineSlimmedPrimaryVertices_m_ndof);
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
  fChain->SetBranchAddress("triggerNames", &triggerNames, &b_triggerNames);
  fChain->SetBranchAddress("triggerResults", &triggerResults, &b_triggerResults);
  fChain->SetBranchAddress("triggerPrescales", &triggerPrescales, &b_triggerPrescales);
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
  fChain->SetBranchAddress("slimmedMuonsUSER", &slimmedMuonsUSER_, &b_slimmedMuonsUSER_);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_charge", slimmedMuonsUSER_m_charge, &b_slimmedMuonsUSER_m_charge);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pt", slimmedMuonsUSER_m_pt, &b_slimmedMuonsUSER_m_pt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_eta", slimmedMuonsUSER_m_eta, &b_slimmedMuonsUSER_m_eta);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_phi", slimmedMuonsUSER_m_phi, &b_slimmedMuonsUSER_m_phi);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_energy", slimmedMuonsUSER_m_energy, &b_slimmedMuonsUSER_m_energy);
  fChain->SetBranchAddress("slimmedMuonsUSER.sel_bits", slimmedMuonsUSER_sel_bits, &b_slimmedMuonsUSER_sel_bits);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dxy", slimmedMuonsUSER_m_dxy, &b_slimmedMuonsUSER_m_dxy);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dxy_error", slimmedMuonsUSER_m_dxy_error, &b_slimmedMuonsUSER_m_dxy_error);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dz", slimmedMuonsUSER_m_dz, &b_slimmedMuonsUSER_m_dz);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_dz_error", slimmedMuonsUSER_m_dz_error, &b_slimmedMuonsUSER_m_dz_error);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_globalTrack_normalizedChi2", slimmedMuonsUSER_m_globalTrack_normalizedChi2, &b_slimmedMuonsUSER_m_globalTrack_normalizedChi2);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_globalTrack_numberOfValidMuonHits", slimmedMuonsUSER_m_globalTrack_numberOfValidMuonHits, &b_slimmedMuonsUSER_m_globalTrack_numberOfValidMuonHits);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_numberOfMatchedStations", slimmedMuonsUSER_m_numberOfMatchedStations, &b_slimmedMuonsUSER_m_numberOfMatchedStations);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_innerTrack_trackerLayersWithMeasurement", slimmedMuonsUSER_m_innerTrack_trackerLayersWithMeasurement, &b_slimmedMuonsUSER_m_innerTrack_trackerLayersWithMeasurement);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_innerTrack_numberOfValidPixelHits", slimmedMuonsUSER_m_innerTrack_numberOfValidPixelHits, &b_slimmedMuonsUSER_m_innerTrack_numberOfValidPixelHits);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_innerTrack_validFraction", slimmedMuonsUSER_m_innerTrack_validFraction, &b_slimmedMuonsUSER_m_innerTrack_validFraction);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_combinedQuality_chi2LocalPosition", slimmedMuonsUSER_m_combinedQuality_chi2LocalPosition, &b_slimmedMuonsUSER_m_combinedQuality_chi2LocalPosition);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_combinedQuality_trkKink", slimmedMuonsUSER_m_combinedQuality_trkKink, &b_slimmedMuonsUSER_m_combinedQuality_trkKink);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_segmentCompatibility", slimmedMuonsUSER_m_segmentCompatibility, &b_slimmedMuonsUSER_m_segmentCompatibility);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumChargedHadronPt", slimmedMuonsUSER_m_sumChargedHadronPt, &b_slimmedMuonsUSER_m_sumChargedHadronPt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumNeutralHadronEt", slimmedMuonsUSER_m_sumNeutralHadronEt, &b_slimmedMuonsUSER_m_sumNeutralHadronEt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumPhotonEt", slimmedMuonsUSER_m_sumPhotonEt, &b_slimmedMuonsUSER_m_sumPhotonEt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_sumPUPt", slimmedMuonsUSER_m_sumPUPt, &b_slimmedMuonsUSER_m_sumPUPt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_CH", slimmedMuonsUSER_m_pfMINIIso_CH, &b_slimmedMuonsUSER_m_pfMINIIso_CH);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_NH", slimmedMuonsUSER_m_pfMINIIso_NH, &b_slimmedMuonsUSER_m_pfMINIIso_NH);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_Ph", slimmedMuonsUSER_m_pfMINIIso_Ph, &b_slimmedMuonsUSER_m_pfMINIIso_Ph);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_PU", slimmedMuonsUSER_m_pfMINIIso_PU, &b_slimmedMuonsUSER_m_pfMINIIso_PU);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_NH_pfwgt", slimmedMuonsUSER_m_pfMINIIso_NH_pfwgt, &b_slimmedMuonsUSER_m_pfMINIIso_NH_pfwgt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_pfMINIIso_Ph_pfwgt", slimmedMuonsUSER_m_pfMINIIso_Ph_pfwgt, &b_slimmedMuonsUSER_m_pfMINIIso_Ph_pfwgt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simType", slimmedMuonsUSER_m_simType, &b_slimmedMuonsUSER_m_simType);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simFlavor", slimmedMuonsUSER_m_simFlavor, &b_slimmedMuonsUSER_m_simFlavor);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simPdgId", slimmedMuonsUSER_m_simPdgId, &b_slimmedMuonsUSER_m_simPdgId);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simMotherPdgId", slimmedMuonsUSER_m_simMotherPdgId, &b_slimmedMuonsUSER_m_simMotherPdgId);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_simHeaviestMotherFlavor", slimmedMuonsUSER_m_simHeaviestMotherFlavor, &b_slimmedMuonsUSER_m_simHeaviestMotherFlavor);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackPt", slimmedMuonsUSER_m_tunePTrackPt, &b_slimmedMuonsUSER_m_tunePTrackPt);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackEta", slimmedMuonsUSER_m_tunePTrackEta, &b_slimmedMuonsUSER_m_tunePTrackEta);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackPhi", slimmedMuonsUSER_m_tunePTrackPhi, &b_slimmedMuonsUSER_m_tunePTrackPhi);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_tunePTrackType", slimmedMuonsUSER_m_tunePTrackType, &b_slimmedMuonsUSER_m_tunePTrackType);
  fChain->SetBranchAddress("slimmedMuonsUSER.m_source_candidates", slimmedMuonsUSER_m_source_candidates, &b_slimmedMuonsUSER_m_source_candidates);
  fChain->SetBranchAddress("slimmedMuonsUSER.tags.tagdata", slimmedMuonsUSER_tags_tagdata, &b_slimmedMuonsUSER_tags_tagdata);
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
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi", &updatedPatJetsSlimmedJetsPuppi_, &b_updatedPatJetsSlimmedJetsPuppi_);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_charge", updatedPatJetsSlimmedJetsPuppi_m_charge, &b_updatedPatJetsSlimmedJetsPuppi_m_charge);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_pt", updatedPatJetsSlimmedJetsPuppi_m_pt, &b_updatedPatJetsSlimmedJetsPuppi_m_pt);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_eta", updatedPatJetsSlimmedJetsPuppi_m_eta, &b_updatedPatJetsSlimmedJetsPuppi_m_eta);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_phi", updatedPatJetsSlimmedJetsPuppi_m_phi, &b_updatedPatJetsSlimmedJetsPuppi_m_phi);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_energy", updatedPatJetsSlimmedJetsPuppi_m_energy, &b_updatedPatJetsSlimmedJetsPuppi_m_energy);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_pdgId", updatedPatJetsSlimmedJetsPuppi_m_pdgId, &b_updatedPatJetsSlimmedJetsPuppi_m_pdgId);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_jetArea", updatedPatJetsSlimmedJetsPuppi_m_jetArea, &b_updatedPatJetsSlimmedJetsPuppi_m_jetArea);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_numberOfDaughters", updatedPatJetsSlimmedJetsPuppi_m_numberOfDaughters, &b_updatedPatJetsSlimmedJetsPuppi_m_numberOfDaughters);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_neutralEmEnergyFraction", updatedPatJetsSlimmedJetsPuppi_m_neutralEmEnergyFraction, &b_updatedPatJetsSlimmedJetsPuppi_m_neutralEmEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_neutralHadronEnergyFraction", updatedPatJetsSlimmedJetsPuppi_m_neutralHadronEnergyFraction, &b_updatedPatJetsSlimmedJetsPuppi_m_neutralHadronEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_chargedEmEnergyFraction", updatedPatJetsSlimmedJetsPuppi_m_chargedEmEnergyFraction, &b_updatedPatJetsSlimmedJetsPuppi_m_chargedEmEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_chargedHadronEnergyFraction", updatedPatJetsSlimmedJetsPuppi_m_chargedHadronEnergyFraction, &b_updatedPatJetsSlimmedJetsPuppi_m_chargedHadronEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_muonEnergyFraction", updatedPatJetsSlimmedJetsPuppi_m_muonEnergyFraction, &b_updatedPatJetsSlimmedJetsPuppi_m_muonEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_photonEnergyFraction", updatedPatJetsSlimmedJetsPuppi_m_photonEnergyFraction, &b_updatedPatJetsSlimmedJetsPuppi_m_photonEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_chargedMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_chargedMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_chargedMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_neutralMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_neutralMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_neutralMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_muonMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_muonMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_muonMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_electronMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_electronMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_electronMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_photonMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_photonMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_photonMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_puppiMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_puppiMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_puppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_neutralPuppiMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_neutralPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_neutralPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_neutralHadronPuppiMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_neutralHadronPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_neutralHadronPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_photonPuppiMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_photonPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_photonPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_HFHadronPuppiMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_HFHadronPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_HFHadronPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_HFEMPuppiMultiplicity", updatedPatJetsSlimmedJetsPuppi_m_HFEMPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsPuppi_m_HFEMPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btag_combinedSecondaryVertex", updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertex, &b_updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertex);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btag_combinedSecondaryVertexMVA", updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertexMVA, &b_updatedPatJetsSlimmedJetsPuppi_m_btag_combinedSecondaryVertexMVA);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btag_DeepCSV_probb", updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probb, &b_updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probb);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btag_DeepCSV_probbb", updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probbb, &b_updatedPatJetsSlimmedJetsPuppi_m_btag_DeepCSV_probbb);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btag_BoostedDoubleSecondaryVertexAK8", updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexAK8, &b_updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexAK8);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btag_BoostedDoubleSecondaryVertexCA15", updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexCA15, &b_updatedPatJetsSlimmedJetsPuppi_m_btag_BoostedDoubleSecondaryVertexCA15);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_JEC_factor_raw", updatedPatJetsSlimmedJetsPuppi_m_JEC_factor_raw, &b_updatedPatJetsSlimmedJetsPuppi_m_JEC_factor_raw);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_JEC_L1factor_raw", updatedPatJetsSlimmedJetsPuppi_m_JEC_L1factor_raw, &b_updatedPatJetsSlimmedJetsPuppi_m_JEC_L1factor_raw);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_genjet_index", updatedPatJetsSlimmedJetsPuppi_m_genjet_index, &b_updatedPatJetsSlimmedJetsPuppi_m_genjet_index);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_hadronFlavor", updatedPatJetsSlimmedJetsPuppi_m_hadronFlavor, &b_updatedPatJetsSlimmedJetsPuppi_m_hadronFlavor);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackMomentum", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackMomentum, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackMomentum);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackEta", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEta, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEta);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackEtaRel", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEtaRel, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackEtaRel);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackDeltaR", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDeltaR, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDeltaR);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackSip3dVal", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dVal, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackSip3dSig", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSig, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackSip2dVal", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dVal, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackSip2dSig", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dSig, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip2dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackDecayLenVal", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDecayLenVal, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackDecayLenVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackChi2", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackChi2, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackChi2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackNTotalHits", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNTotalHits, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNTotalHits);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackNPixelHits", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNPixelHits, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackNPixelHits);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackPtRel", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRel, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRel);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackPPar", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPPar, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPPar);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackPtRatio", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRatio, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPtRatio);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackPParRatio", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPParRatio, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackPParRatio);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackJetDistVal", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistVal, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackJetDistSig", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistSig, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackJetDistSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackGhostTrackDistVal", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistVal, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackGhostTrackDistSig", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistSig, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackDistSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackGhostTrackWeight", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackWeight, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackGhostTrackWeight);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_FlightDistance2dVal", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dVal, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_FlightDistance2dSig", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dSig, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance2dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_FlightDistance3dVal", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dVal, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_FlightDistance3dSig", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dSig, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_FlightDistance3dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexJetDeltaR", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexJetDeltaR, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexJetDeltaR);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_JetNSecondaryVertices", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_JetNSecondaryVertices, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_JetNSecondaryVertices);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexNTracks", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNTracks, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNTracks);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexChi2", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexChi2, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexChi2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexNdof", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNdof, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNdof);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexNormalizedChi2", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNormalizedChi2, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexNormalizedChi2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexCategoryJTC", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexCategoryJTC, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexCategoryJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexMassJTC", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexMassJTC, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexMassJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_VertexEnergyRatioJTC", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexEnergyRatioJTC, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_VertexEnergyRatioJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.m_TrackSip3dSigAboveCharmJTC", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_btaginfo.tags.tagdata", updatedPatJetsSlimmedJetsPuppi_m_btaginfo_tags_tagdata, &b_updatedPatJetsSlimmedJetsPuppi_m_btaginfo_tags_tagdata);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.m_lepton_keys", updatedPatJetsSlimmedJetsPuppi_m_lepton_keys, &b_updatedPatJetsSlimmedJetsPuppi_m_lepton_keys);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsPuppi.tags.tagdata", updatedPatJetsSlimmedJetsPuppi_tags_tagdata, &b_updatedPatJetsSlimmedJetsPuppi_tags_tagdata);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI", &updatedPatJetsPatJetsAK8PFPUPPI_, &b_updatedPatJetsPatJetsAK8PFPUPPI_);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_charge", updatedPatJetsPatJetsAK8PFPUPPI_m_charge, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_charge);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_pt", updatedPatJetsPatJetsAK8PFPUPPI_m_pt, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_pt);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_eta", updatedPatJetsPatJetsAK8PFPUPPI_m_eta, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_eta);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_phi", updatedPatJetsPatJetsAK8PFPUPPI_m_phi, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_phi);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_energy", updatedPatJetsPatJetsAK8PFPUPPI_m_energy, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_energy);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_pdgId", updatedPatJetsPatJetsAK8PFPUPPI_m_pdgId, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_pdgId);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_jetArea", updatedPatJetsPatJetsAK8PFPUPPI_m_jetArea, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_jetArea);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_numberOfDaughters", updatedPatJetsPatJetsAK8PFPUPPI_m_numberOfDaughters, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_numberOfDaughters);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_neutralEmEnergyFraction", updatedPatJetsPatJetsAK8PFPUPPI_m_neutralEmEnergyFraction, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralEmEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_neutralHadronEnergyFraction", updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronEnergyFraction, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_chargedEmEnergyFraction", updatedPatJetsPatJetsAK8PFPUPPI_m_chargedEmEnergyFraction, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_chargedEmEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_chargedHadronEnergyFraction", updatedPatJetsPatJetsAK8PFPUPPI_m_chargedHadronEnergyFraction, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_chargedHadronEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_muonEnergyFraction", updatedPatJetsPatJetsAK8PFPUPPI_m_muonEnergyFraction, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_muonEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_photonEnergyFraction", updatedPatJetsPatJetsAK8PFPUPPI_m_photonEnergyFraction, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_photonEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_chargedMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_chargedMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_chargedMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_neutralMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_neutralMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_muonMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_muonMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_muonMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_electronMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_electronMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_electronMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_photonMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_photonMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_photonMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_puppiMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_puppiMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_puppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_neutralPuppiMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_neutralPuppiMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_neutralHadronPuppiMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronPuppiMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_neutralHadronPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_photonPuppiMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_photonPuppiMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_photonPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_HFHadronPuppiMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_HFHadronPuppiMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_HFHadronPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_HFEMPuppiMultiplicity", updatedPatJetsPatJetsAK8PFPUPPI_m_HFEMPuppiMultiplicity, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_HFEMPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btag_combinedSecondaryVertex", updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertex, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertex);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btag_combinedSecondaryVertexMVA", updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertexMVA, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_combinedSecondaryVertexMVA);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btag_DeepCSV_probb", updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probb, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probb);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btag_DeepCSV_probbb", updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probbb, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_DeepCSV_probbb);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btag_BoostedDoubleSecondaryVertexAK8", updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexAK8, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexAK8);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btag_BoostedDoubleSecondaryVertexCA15", updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexCA15, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btag_BoostedDoubleSecondaryVertexCA15);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_JEC_factor_raw", updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_factor_raw, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_factor_raw);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_JEC_L1factor_raw", updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_L1factor_raw, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_JEC_L1factor_raw);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_genjet_index", updatedPatJetsPatJetsAK8PFPUPPI_m_genjet_index, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_genjet_index);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_hadronFlavor", updatedPatJetsPatJetsAK8PFPUPPI_m_hadronFlavor, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_hadronFlavor);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackMomentum", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackMomentum, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackMomentum);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackEta", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEta, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEta);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackEtaRel", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEtaRel, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackEtaRel);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackDeltaR", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDeltaR, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDeltaR);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackSip3dVal", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dVal, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dVal);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackSip3dSig", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSig, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSig);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackSip2dVal", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dVal, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dVal);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackSip2dSig", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dSig, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip2dSig);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackDecayLenVal", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDecayLenVal, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackDecayLenVal);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackChi2", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackChi2, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackChi2);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackNTotalHits", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNTotalHits, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNTotalHits);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackNPixelHits", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNPixelHits, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackNPixelHits);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackPtRel", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRel, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRel);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackPPar", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPPar, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPPar);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackPtRatio", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRatio, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPtRatio);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackPParRatio", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPParRatio, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackPParRatio);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackJetDistVal", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistVal, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistVal);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackJetDistSig", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistSig, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackJetDistSig);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackGhostTrackDistVal", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistVal, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistVal);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackGhostTrackDistSig", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistSig, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackDistSig);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackGhostTrackWeight", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackWeight, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackGhostTrackWeight);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_FlightDistance2dVal", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dVal, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dVal);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_FlightDistance2dSig", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dSig, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance2dSig);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_FlightDistance3dVal", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dVal, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dVal);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_FlightDistance3dSig", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dSig, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_FlightDistance3dSig);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexJetDeltaR", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexJetDeltaR, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexJetDeltaR);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_JetNSecondaryVertices", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_JetNSecondaryVertices, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_JetNSecondaryVertices);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexNTracks", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNTracks, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNTracks);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexChi2", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexChi2, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexChi2);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexNdof", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNdof, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNdof);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexNormalizedChi2", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNormalizedChi2, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexNormalizedChi2);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexCategoryJTC", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexCategoryJTC, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexCategoryJTC);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexMassJTC", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexMassJTC, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexMassJTC);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_VertexEnergyRatioJTC", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexEnergyRatioJTC, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_VertexEnergyRatioJTC);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.m_TrackSip3dSigAboveCharmJTC", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSigAboveCharmJTC, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_m_TrackSip3dSigAboveCharmJTC);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_btaginfo.tags.tagdata", updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_tags_tagdata, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_btaginfo_tags_tagdata);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.m_lepton_keys", updatedPatJetsPatJetsAK8PFPUPPI_m_lepton_keys, &b_updatedPatJetsPatJetsAK8PFPUPPI_m_lepton_keys);
  fChain->SetBranchAddress("updatedPatJetsPatJetsAK8PFPUPPI.tags.tagdata", updatedPatJetsPatJetsAK8PFPUPPI_tags_tagdata, &b_updatedPatJetsPatJetsAK8PFPUPPI_tags_tagdata);
  fChain->SetBranchAddress("patJetsAK8PFCHS", &patJetsAK8PFCHS_, &b_patJetsAK8PFCHS_);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_charge", patJetsAK8PFCHS_m_charge, &b_patJetsAK8PFCHS_m_charge);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_pt", patJetsAK8PFCHS_m_pt, &b_patJetsAK8PFCHS_m_pt);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_eta", patJetsAK8PFCHS_m_eta, &b_patJetsAK8PFCHS_m_eta);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_phi", patJetsAK8PFCHS_m_phi, &b_patJetsAK8PFCHS_m_phi);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_energy", patJetsAK8PFCHS_m_energy, &b_patJetsAK8PFCHS_m_energy);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_pdgId", patJetsAK8PFCHS_m_pdgId, &b_patJetsAK8PFCHS_m_pdgId);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_jetArea", patJetsAK8PFCHS_m_jetArea, &b_patJetsAK8PFCHS_m_jetArea);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_numberOfDaughters", patJetsAK8PFCHS_m_numberOfDaughters, &b_patJetsAK8PFCHS_m_numberOfDaughters);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_neutralEmEnergyFraction", patJetsAK8PFCHS_m_neutralEmEnergyFraction, &b_patJetsAK8PFCHS_m_neutralEmEnergyFraction);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_neutralHadronEnergyFraction", patJetsAK8PFCHS_m_neutralHadronEnergyFraction, &b_patJetsAK8PFCHS_m_neutralHadronEnergyFraction);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_chargedEmEnergyFraction", patJetsAK8PFCHS_m_chargedEmEnergyFraction, &b_patJetsAK8PFCHS_m_chargedEmEnergyFraction);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_chargedHadronEnergyFraction", patJetsAK8PFCHS_m_chargedHadronEnergyFraction, &b_patJetsAK8PFCHS_m_chargedHadronEnergyFraction);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_muonEnergyFraction", patJetsAK8PFCHS_m_muonEnergyFraction, &b_patJetsAK8PFCHS_m_muonEnergyFraction);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_photonEnergyFraction", patJetsAK8PFCHS_m_photonEnergyFraction, &b_patJetsAK8PFCHS_m_photonEnergyFraction);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_chargedMultiplicity", patJetsAK8PFCHS_m_chargedMultiplicity, &b_patJetsAK8PFCHS_m_chargedMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_neutralMultiplicity", patJetsAK8PFCHS_m_neutralMultiplicity, &b_patJetsAK8PFCHS_m_neutralMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_muonMultiplicity", patJetsAK8PFCHS_m_muonMultiplicity, &b_patJetsAK8PFCHS_m_muonMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_electronMultiplicity", patJetsAK8PFCHS_m_electronMultiplicity, &b_patJetsAK8PFCHS_m_electronMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_photonMultiplicity", patJetsAK8PFCHS_m_photonMultiplicity, &b_patJetsAK8PFCHS_m_photonMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_puppiMultiplicity", patJetsAK8PFCHS_m_puppiMultiplicity, &b_patJetsAK8PFCHS_m_puppiMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_neutralPuppiMultiplicity", patJetsAK8PFCHS_m_neutralPuppiMultiplicity, &b_patJetsAK8PFCHS_m_neutralPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_neutralHadronPuppiMultiplicity", patJetsAK8PFCHS_m_neutralHadronPuppiMultiplicity, &b_patJetsAK8PFCHS_m_neutralHadronPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_photonPuppiMultiplicity", patJetsAK8PFCHS_m_photonPuppiMultiplicity, &b_patJetsAK8PFCHS_m_photonPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_HFHadronPuppiMultiplicity", patJetsAK8PFCHS_m_HFHadronPuppiMultiplicity, &b_patJetsAK8PFCHS_m_HFHadronPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_HFEMPuppiMultiplicity", patJetsAK8PFCHS_m_HFEMPuppiMultiplicity, &b_patJetsAK8PFCHS_m_HFEMPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btag_combinedSecondaryVertex", patJetsAK8PFCHS_m_btag_combinedSecondaryVertex, &b_patJetsAK8PFCHS_m_btag_combinedSecondaryVertex);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btag_combinedSecondaryVertexMVA", patJetsAK8PFCHS_m_btag_combinedSecondaryVertexMVA, &b_patJetsAK8PFCHS_m_btag_combinedSecondaryVertexMVA);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btag_DeepCSV_probb", patJetsAK8PFCHS_m_btag_DeepCSV_probb, &b_patJetsAK8PFCHS_m_btag_DeepCSV_probb);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btag_DeepCSV_probbb", patJetsAK8PFCHS_m_btag_DeepCSV_probbb, &b_patJetsAK8PFCHS_m_btag_DeepCSV_probbb);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btag_BoostedDoubleSecondaryVertexAK8", patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexAK8, &b_patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexAK8);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btag_BoostedDoubleSecondaryVertexCA15", patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexCA15, &b_patJetsAK8PFCHS_m_btag_BoostedDoubleSecondaryVertexCA15);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_JEC_factor_raw", patJetsAK8PFCHS_m_JEC_factor_raw, &b_patJetsAK8PFCHS_m_JEC_factor_raw);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_JEC_L1factor_raw", patJetsAK8PFCHS_m_JEC_L1factor_raw, &b_patJetsAK8PFCHS_m_JEC_L1factor_raw);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_genjet_index", patJetsAK8PFCHS_m_genjet_index, &b_patJetsAK8PFCHS_m_genjet_index);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_hadronFlavor", patJetsAK8PFCHS_m_hadronFlavor, &b_patJetsAK8PFCHS_m_hadronFlavor);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackMomentum", patJetsAK8PFCHS_m_btaginfo_m_TrackMomentum, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackMomentum);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackEta", patJetsAK8PFCHS_m_btaginfo_m_TrackEta, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackEta);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackEtaRel", patJetsAK8PFCHS_m_btaginfo_m_TrackEtaRel, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackEtaRel);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackDeltaR", patJetsAK8PFCHS_m_btaginfo_m_TrackDeltaR, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackDeltaR);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackSip3dVal", patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dVal, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dVal);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackSip3dSig", patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSig, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSig);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackSip2dVal", patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dVal, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dVal);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackSip2dSig", patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dSig, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip2dSig);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackDecayLenVal", patJetsAK8PFCHS_m_btaginfo_m_TrackDecayLenVal, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackDecayLenVal);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackChi2", patJetsAK8PFCHS_m_btaginfo_m_TrackChi2, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackChi2);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackNTotalHits", patJetsAK8PFCHS_m_btaginfo_m_TrackNTotalHits, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackNTotalHits);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackNPixelHits", patJetsAK8PFCHS_m_btaginfo_m_TrackNPixelHits, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackNPixelHits);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackPtRel", patJetsAK8PFCHS_m_btaginfo_m_TrackPtRel, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackPtRel);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackPPar", patJetsAK8PFCHS_m_btaginfo_m_TrackPPar, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackPPar);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackPtRatio", patJetsAK8PFCHS_m_btaginfo_m_TrackPtRatio, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackPtRatio);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackPParRatio", patJetsAK8PFCHS_m_btaginfo_m_TrackPParRatio, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackPParRatio);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackJetDistVal", patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistVal, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistVal);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackJetDistSig", patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistSig, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackJetDistSig);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackGhostTrackDistVal", patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistVal, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistVal);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackGhostTrackDistSig", patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistSig, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackDistSig);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackGhostTrackWeight", patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackWeight, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackGhostTrackWeight);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_FlightDistance2dVal", patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dVal, &b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dVal);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_FlightDistance2dSig", patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dSig, &b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance2dSig);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_FlightDistance3dVal", patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dVal, &b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dVal);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_FlightDistance3dSig", patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dSig, &b_patJetsAK8PFCHS_m_btaginfo_m_FlightDistance3dSig);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexJetDeltaR", patJetsAK8PFCHS_m_btaginfo_m_VertexJetDeltaR, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexJetDeltaR);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_JetNSecondaryVertices", patJetsAK8PFCHS_m_btaginfo_m_JetNSecondaryVertices, &b_patJetsAK8PFCHS_m_btaginfo_m_JetNSecondaryVertices);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexNTracks", patJetsAK8PFCHS_m_btaginfo_m_VertexNTracks, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexNTracks);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexChi2", patJetsAK8PFCHS_m_btaginfo_m_VertexChi2, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexChi2);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexNdof", patJetsAK8PFCHS_m_btaginfo_m_VertexNdof, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexNdof);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexNormalizedChi2", patJetsAK8PFCHS_m_btaginfo_m_VertexNormalizedChi2, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexNormalizedChi2);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexCategoryJTC", patJetsAK8PFCHS_m_btaginfo_m_VertexCategoryJTC, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexCategoryJTC);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexMassJTC", patJetsAK8PFCHS_m_btaginfo_m_VertexMassJTC, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexMassJTC);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_VertexEnergyRatioJTC", patJetsAK8PFCHS_m_btaginfo_m_VertexEnergyRatioJTC, &b_patJetsAK8PFCHS_m_btaginfo_m_VertexEnergyRatioJTC);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.m_TrackSip3dSigAboveCharmJTC", patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC, &b_patJetsAK8PFCHS_m_btaginfo_m_TrackSip3dSigAboveCharmJTC);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_btaginfo.tags.tagdata", patJetsAK8PFCHS_m_btaginfo_tags_tagdata, &b_patJetsAK8PFCHS_m_btaginfo_tags_tagdata);
  fChain->SetBranchAddress("patJetsAK8PFCHS.m_lepton_keys", patJetsAK8PFCHS_m_lepton_keys, &b_patJetsAK8PFCHS_m_lepton_keys);
  fChain->SetBranchAddress("patJetsAK8PFCHS.tags.tagdata", patJetsAK8PFCHS_tags_tagdata, &b_patJetsAK8PFCHS_tags_tagdata);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi", &updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_charge", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_charge, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_charge);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_pt", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pt, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pt);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_eta", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_eta, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_eta);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_phi", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_phi, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_phi);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_energy", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_energy, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_energy);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_pdgId", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pdgId, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_pdgId);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_jetArea", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_jetArea, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_jetArea);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_numberOfDaughters", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_numberOfDaughters, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_numberOfDaughters);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_neutralEmEnergyFraction", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralEmEnergyFraction, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralEmEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_neutralHadronEnergyFraction", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronEnergyFraction, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_chargedEmEnergyFraction", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedEmEnergyFraction, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedEmEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_chargedHadronEnergyFraction", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedHadronEnergyFraction, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedHadronEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_muonEnergyFraction", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonEnergyFraction, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_photonEnergyFraction", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonEnergyFraction, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonEnergyFraction);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_chargedMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_chargedMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_neutralMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_muonMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_muonMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_electronMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_electronMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_electronMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_photonMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_puppiMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_puppiMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_puppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_neutralPuppiMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_neutralHadronPuppiMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_neutralHadronPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_photonPuppiMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_photonPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_HFHadronPuppiMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFHadronPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFHadronPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_HFEMPuppiMultiplicity", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFEMPuppiMultiplicity, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_HFEMPuppiMultiplicity);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btag_combinedSecondaryVertex", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertex, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertex);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btag_combinedSecondaryVertexMVA", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertexMVA, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_combinedSecondaryVertexMVA);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btag_DeepCSV_probb", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probb, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probb);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btag_DeepCSV_probbb", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probbb, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_DeepCSV_probbb);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btag_BoostedDoubleSecondaryVertexAK8", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexAK8, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexAK8);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btag_BoostedDoubleSecondaryVertexCA15", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexCA15, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btag_BoostedDoubleSecondaryVertexCA15);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_JEC_factor_raw", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_factor_raw, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_factor_raw);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_JEC_L1factor_raw", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_L1factor_raw, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_JEC_L1factor_raw);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_genjet_index", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_genjet_index, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_genjet_index);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_hadronFlavor", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_hadronFlavor, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_hadronFlavor);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackMomentum", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackMomentum, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackMomentum);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackEta", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEta, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEta);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackEtaRel", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEtaRel, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackEtaRel);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackDeltaR", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDeltaR, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDeltaR);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackSip3dVal", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dVal, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackSip3dSig", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSig, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackSip2dVal", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dVal, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackSip2dSig", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dSig, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip2dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackDecayLenVal", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDecayLenVal, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackDecayLenVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackChi2", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackChi2, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackChi2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackNTotalHits", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNTotalHits, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNTotalHits);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackNPixelHits", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNPixelHits, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackNPixelHits);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackPtRel", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRel, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRel);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackPPar", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPPar, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPPar);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackPtRatio", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRatio, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPtRatio);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackPParRatio", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPParRatio, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackPParRatio);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackJetDistVal", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistVal, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackJetDistSig", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistSig, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackJetDistSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackGhostTrackDistVal", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistVal, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackGhostTrackDistSig", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistSig, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackDistSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackGhostTrackWeight", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackWeight, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackGhostTrackWeight);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_FlightDistance2dVal", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dVal, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_FlightDistance2dSig", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dSig, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance2dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_FlightDistance3dVal", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dVal, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dVal);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_FlightDistance3dSig", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dSig, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_FlightDistance3dSig);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexJetDeltaR", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexJetDeltaR, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexJetDeltaR);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_JetNSecondaryVertices", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_JetNSecondaryVertices, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_JetNSecondaryVertices);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexNTracks", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNTracks, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNTracks);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexChi2", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexChi2, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexChi2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexNdof", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNdof, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNdof);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexNormalizedChi2", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNormalizedChi2, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexNormalizedChi2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexCategoryJTC", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexCategoryJTC, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexCategoryJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexMassJTC", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexMassJTC, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexMassJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_VertexEnergyRatioJTC", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexEnergyRatioJTC, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_VertexEnergyRatioJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.m_TrackSip3dSigAboveCharmJTC", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_m_TrackSip3dSigAboveCharmJTC);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_btaginfo.tags.tagdata", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_tags_tagdata, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_btaginfo_tags_tagdata);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_lepton_keys", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_lepton_keys, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_lepton_keys);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.tags.tagdata", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_subjets", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_subjets, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_subjets);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_qjets_volatility", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_qjets_volatility, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_qjets_volatility);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau1", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau2", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau3", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau4", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau1_groomed", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1_groomed, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau1_groomed);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau2_groomed", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2_groomed, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau2_groomed);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau3_groomed", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3_groomed, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau3_groomed);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_tau4_groomed", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4_groomed, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_tau4_groomed);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_ecfN2_beta1", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta1, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta1);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_ecfN2_beta2", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta2, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN2_beta2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_ecfN3_beta1", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta1, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta1);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_ecfN3_beta2", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta2, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_ecfN3_beta2);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_mvahiggsdiscr", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_mvahiggsdiscr, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_mvahiggsdiscr);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_prunedmass", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_prunedmass, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_prunedmass);
  fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.m_softdropmass", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_softdropmass, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_m_softdropmass);
  //    fChain->SetBranchAddress("updatedPatJetsSlimmedJetsAK8_SoftDropPuppi.tags.tagdata", updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata, &b_updatedPatJetsSlimmedJetsAK8_SoftDropPuppi_tags_tagdata);
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
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters", &patJetsHepTopTagCHSPacked_daughters_, &b_patJetsHepTopTagCHSPacked_daughters_);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_charge", patJetsHepTopTagCHSPacked_daughters_m_charge, &b_patJetsHepTopTagCHSPacked_daughters_m_charge);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_pt", patJetsHepTopTagCHSPacked_daughters_m_pt, &b_patJetsHepTopTagCHSPacked_daughters_m_pt);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_eta", patJetsHepTopTagCHSPacked_daughters_m_eta, &b_patJetsHepTopTagCHSPacked_daughters_m_eta);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_phi", patJetsHepTopTagCHSPacked_daughters_m_phi, &b_patJetsHepTopTagCHSPacked_daughters_m_phi);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_energy", patJetsHepTopTagCHSPacked_daughters_m_energy, &b_patJetsHepTopTagCHSPacked_daughters_m_energy);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_pdgId", patJetsHepTopTagCHSPacked_daughters_m_pdgId, &b_patJetsHepTopTagCHSPacked_daughters_m_pdgId);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_jetArea", patJetsHepTopTagCHSPacked_daughters_m_jetArea, &b_patJetsHepTopTagCHSPacked_daughters_m_jetArea);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_numberOfDaughters", patJetsHepTopTagCHSPacked_daughters_m_numberOfDaughters, &b_patJetsHepTopTagCHSPacked_daughters_m_numberOfDaughters);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_neutralEmEnergyFraction", patJetsHepTopTagCHSPacked_daughters_m_neutralEmEnergyFraction, &b_patJetsHepTopTagCHSPacked_daughters_m_neutralEmEnergyFraction);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_neutralHadronEnergyFraction", patJetsHepTopTagCHSPacked_daughters_m_neutralHadronEnergyFraction, &b_patJetsHepTopTagCHSPacked_daughters_m_neutralHadronEnergyFraction);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_chargedEmEnergyFraction", patJetsHepTopTagCHSPacked_daughters_m_chargedEmEnergyFraction, &b_patJetsHepTopTagCHSPacked_daughters_m_chargedEmEnergyFraction);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_chargedHadronEnergyFraction", patJetsHepTopTagCHSPacked_daughters_m_chargedHadronEnergyFraction, &b_patJetsHepTopTagCHSPacked_daughters_m_chargedHadronEnergyFraction);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_muonEnergyFraction", patJetsHepTopTagCHSPacked_daughters_m_muonEnergyFraction, &b_patJetsHepTopTagCHSPacked_daughters_m_muonEnergyFraction);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_photonEnergyFraction", patJetsHepTopTagCHSPacked_daughters_m_photonEnergyFraction, &b_patJetsHepTopTagCHSPacked_daughters_m_photonEnergyFraction);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_chargedMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_chargedMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_chargedMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_neutralMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_neutralMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_neutralMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_muonMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_muonMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_muonMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_electronMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_electronMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_electronMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_photonMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_photonMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_photonMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_puppiMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_puppiMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_puppiMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_neutralPuppiMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_neutralPuppiMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_neutralPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_neutralHadronPuppiMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_neutralHadronPuppiMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_neutralHadronPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_photonPuppiMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_photonPuppiMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_photonPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_HFHadronPuppiMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_HFHadronPuppiMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_HFHadronPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_HFEMPuppiMultiplicity", patJetsHepTopTagCHSPacked_daughters_m_HFEMPuppiMultiplicity, &b_patJetsHepTopTagCHSPacked_daughters_m_HFEMPuppiMultiplicity);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btag_combinedSecondaryVertex", patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertex, &b_patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertex);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btag_combinedSecondaryVertexMVA", patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertexMVA, &b_patJetsHepTopTagCHSPacked_daughters_m_btag_combinedSecondaryVertexMVA);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btag_DeepCSV_probb", patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probb, &b_patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probb);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btag_DeepCSV_probbb", patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probbb, &b_patJetsHepTopTagCHSPacked_daughters_m_btag_DeepCSV_probbb);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btag_BoostedDoubleSecondaryVertexAK8", patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexAK8, &b_patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexAK8);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btag_BoostedDoubleSecondaryVertexCA15", patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexCA15, &b_patJetsHepTopTagCHSPacked_daughters_m_btag_BoostedDoubleSecondaryVertexCA15);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_JEC_factor_raw", patJetsHepTopTagCHSPacked_daughters_m_JEC_factor_raw, &b_patJetsHepTopTagCHSPacked_daughters_m_JEC_factor_raw);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_JEC_L1factor_raw", patJetsHepTopTagCHSPacked_daughters_m_JEC_L1factor_raw, &b_patJetsHepTopTagCHSPacked_daughters_m_JEC_L1factor_raw);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_genjet_index", patJetsHepTopTagCHSPacked_daughters_m_genjet_index, &b_patJetsHepTopTagCHSPacked_daughters_m_genjet_index);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_hadronFlavor", patJetsHepTopTagCHSPacked_daughters_m_hadronFlavor, &b_patJetsHepTopTagCHSPacked_daughters_m_hadronFlavor);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackMomentum", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackMomentum, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackMomentum);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackEta", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEta, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEta);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackEtaRel", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEtaRel, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackEtaRel);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackDeltaR", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDeltaR, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDeltaR);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackSip3dVal", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dVal, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dVal);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackSip3dSig", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSig, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSig);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackSip2dVal", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dVal, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dVal);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackSip2dSig", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dSig, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip2dSig);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackDecayLenVal", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDecayLenVal, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackDecayLenVal);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackChi2", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackChi2, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackChi2);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackNTotalHits", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNTotalHits, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNTotalHits);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackNPixelHits", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNPixelHits, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackNPixelHits);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackPtRel", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRel, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRel);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackPPar", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPPar, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPPar);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackPtRatio", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRatio, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPtRatio);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackPParRatio", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPParRatio, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackPParRatio);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackJetDistVal", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistVal, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistVal);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackJetDistSig", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistSig, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackJetDistSig);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackGhostTrackDistVal", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistVal, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistVal);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackGhostTrackDistSig", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistSig, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackDistSig);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackGhostTrackWeight", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackWeight, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackGhostTrackWeight);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_FlightDistance2dVal", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dVal, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dVal);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_FlightDistance2dSig", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dSig, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance2dSig);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_FlightDistance3dVal", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dVal, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dVal);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_FlightDistance3dSig", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dSig, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_FlightDistance3dSig);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexJetDeltaR", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexJetDeltaR, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexJetDeltaR);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_JetNSecondaryVertices", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_JetNSecondaryVertices, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_JetNSecondaryVertices);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexNTracks", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNTracks, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNTracks);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexChi2", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexChi2, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexChi2);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexNdof", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNdof, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNdof);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexNormalizedChi2", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNormalizedChi2, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexNormalizedChi2);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexCategoryJTC", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexCategoryJTC, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexCategoryJTC);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexMassJTC", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexMassJTC, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexMassJTC);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_VertexEnergyRatioJTC", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexEnergyRatioJTC, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_VertexEnergyRatioJTC);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.m_TrackSip3dSigAboveCharmJTC", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSigAboveCharmJTC, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_m_TrackSip3dSigAboveCharmJTC);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_btaginfo.tags.tagdata", patJetsHepTopTagCHSPacked_daughters_m_btaginfo_tags_tagdata, &b_patJetsHepTopTagCHSPacked_daughters_m_btaginfo_tags_tagdata);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_lepton_keys", patJetsHepTopTagCHSPacked_daughters_m_lepton_keys, &b_patJetsHepTopTagCHSPacked_daughters_m_lepton_keys);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.tags.tagdata", patJetsHepTopTagCHSPacked_daughters_tags_tagdata, &b_patJetsHepTopTagCHSPacked_daughters_tags_tagdata);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_subjets", patJetsHepTopTagCHSPacked_daughters_m_subjets, &b_patJetsHepTopTagCHSPacked_daughters_m_subjets);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_qjets_volatility", patJetsHepTopTagCHSPacked_daughters_m_qjets_volatility, &b_patJetsHepTopTagCHSPacked_daughters_m_qjets_volatility);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau1", patJetsHepTopTagCHSPacked_daughters_m_tau1, &b_patJetsHepTopTagCHSPacked_daughters_m_tau1);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau2", patJetsHepTopTagCHSPacked_daughters_m_tau2, &b_patJetsHepTopTagCHSPacked_daughters_m_tau2);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau3", patJetsHepTopTagCHSPacked_daughters_m_tau3, &b_patJetsHepTopTagCHSPacked_daughters_m_tau3);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau4", patJetsHepTopTagCHSPacked_daughters_m_tau4, &b_patJetsHepTopTagCHSPacked_daughters_m_tau4);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau1_groomed", patJetsHepTopTagCHSPacked_daughters_m_tau1_groomed, &b_patJetsHepTopTagCHSPacked_daughters_m_tau1_groomed);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau2_groomed", patJetsHepTopTagCHSPacked_daughters_m_tau2_groomed, &b_patJetsHepTopTagCHSPacked_daughters_m_tau2_groomed);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau3_groomed", patJetsHepTopTagCHSPacked_daughters_m_tau3_groomed, &b_patJetsHepTopTagCHSPacked_daughters_m_tau3_groomed);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_tau4_groomed", patJetsHepTopTagCHSPacked_daughters_m_tau4_groomed, &b_patJetsHepTopTagCHSPacked_daughters_m_tau4_groomed);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_ecfN2_beta1", patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta1, &b_patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta1);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_ecfN2_beta2", patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta2, &b_patJetsHepTopTagCHSPacked_daughters_m_ecfN2_beta2);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_ecfN3_beta1", patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta1, &b_patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta1);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_ecfN3_beta2", patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta2, &b_patJetsHepTopTagCHSPacked_daughters_m_ecfN3_beta2);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_mvahiggsdiscr", patJetsHepTopTagCHSPacked_daughters_m_mvahiggsdiscr, &b_patJetsHepTopTagCHSPacked_daughters_m_mvahiggsdiscr);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_prunedmass", patJetsHepTopTagCHSPacked_daughters_m_prunedmass, &b_patJetsHepTopTagCHSPacked_daughters_m_prunedmass);
  fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.m_softdropmass", patJetsHepTopTagCHSPacked_daughters_m_softdropmass, &b_patJetsHepTopTagCHSPacked_daughters_m_softdropmass);
  //    fChain->SetBranchAddress("patJetsHepTopTagCHSPacked_daughters.tags.tagdata", patJetsHepTopTagCHSPacked_daughters_tags_tagdata, &b_patJetsHepTopTagCHSPacked_daughters_tags_tagdata);
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
