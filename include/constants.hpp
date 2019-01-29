#pragma once

#ifndef  CONSTANTS_H
#define  CONSTANTS_H

// DEFINE GLOBAL VARIABLES

// taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis
constexpr static int s_runnr_A = 297019; //up to this one, including this one
constexpr static int s_runnr_B = 299329; //up to this one, including this one
constexpr static int s_runnr_C = 302029; //up to this one, including this one
constexpr static int s_runnr_D = 303434; //up to this one, including this one
constexpr static int s_runnr_E = 304826; //up to this one, including this one
constexpr static int s_runnr_F = 306462; //up to this one, including this one

const float ZMASS  = 91.0;
const float ZWIDTH = 10.0;

const int min_leptons = 2;
const float min_pt_lepton = 30.0;
const float max_eta_lepton = 2.4;
const float iso_lepton = 0.15;
const float min_delta_phi = 2.7;

#define  JETVARIABLE topjets

enum taggers {NN_IsHiggs=1000, NN_IsQCD, NN_IsTop};

enum Conditions { ZMatch, HMatch, WWfullLep, WWfullHad, WWsemiLep };
enum Decay {nodecay, leptonic, semileptonic, hadronic };

#endif
