#pragma once


#define MAKE_JEC(jecv,jetLabel)         					                   \
if(JEC_Version == #jecv){                                            \
  JEC_corr_B_##jetLabel  = JERFiles::jecv##_B_L123_##jetLabel##_DATA;\
  JEC_corr_C_##jetLabel  = JERFiles::jecv##_C_L123_##jetLabel##_DATA;\
  JEC_corr_D_##jetLabel  = JERFiles::jecv##_D_L123_##jetLabel##_DATA;\
  JEC_corr_E_##jetLabel  = JERFiles::jecv##_E_L123_##jetLabel##_DATA;\
  JEC_corr_F_##jetLabel  = JERFiles::jecv##_F_L123_##jetLabel##_DATA;\
  JEC_corr_MC_##jetLabel = JERFiles::jecv##_L123_##jetLabel##_MC;    \
}                                                                    \


#define SET_MET_FILTERS()                                                                                                                         \
/* https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2*/                                                                    \
metfilters_selection.reset(new AndSelection(ctx, "metfilters"));                                                                                  \
metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");                                                           \
metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");                                                     \
metfilters_selection->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");                             \
metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");                     \
/* if (!is_mc) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");  TODO Not recommended for MC, but do check */  \
metfilters_selection->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter");                                       \
metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");                                                           \
metfilters_selection->add<TriggerSelection>("goodVertices", "Flag_goodVertices");                                                                 \
metfilters_selection->add<TriggerSelection>("ecalBadCalibFilter", "Flag_ecalBadCalibFilter");                                                     \
/* if(pvfilter) metfilters_selection->add<NPVSelection>("1 good PV",1,-1,pvid); TODO */                                                           \


#define SET_JETLEPTON_CLEANER(class_name, jetCollection, jet_type)                                                  \
if(is_mc)	class_name##_MC.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_corr_MC_##jetCollection, #jet_type )); \
else {                                                                                                              \
  class_name##_B.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_corr_B_##jetCollection, #jet_type));            \
  class_name##_C.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_corr_C_##jetCollection, #jet_type));            \
  class_name##_D.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_corr_D_##jetCollection, #jet_type));            \
  class_name##_E.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_corr_E_##jetCollection, #jet_type));            \
  class_name##_F.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_corr_F_##jetCollection, #jet_type));            \
}                                                                                                                   \

// #define  value


#define SET_JET_CORRECTION(class_name, jetCollection, jet_type)               \
jet_type##_corrector_B.reset(new class_name(ctx, JEC_corr_B_##jetCollection));\
jet_type##_corrector_C.reset(new class_name(ctx, JEC_corr_C_##jetCollection));\
jet_type##_corrector_D.reset(new class_name(ctx, JEC_corr_D_##jetCollection));\
jet_type##_corrector_E.reset(new class_name(ctx, JEC_corr_E_##jetCollection));\
jet_type##_corrector_F.reset(new class_name(ctx, JEC_corr_F_##jetCollection));\


#define APPLY_PROCESS(process,method)             \
if(is_mc) process##_MC->method(event);            \
else{                                             \
  if     (apply_B){ process##_B->method(event); } \
  else if(apply_C){ process##_C->method(event); } \
  else if(apply_D){ process##_D->method(event); } \
  else if(apply_E){ process##_E->method(event); } \
  else if(apply_F){ process##_F->method(event); } \
}                                                 \



/*
█ ███████ ███    ██ ██████      ███    ███  █████   ██████ ██████   ██████  ███████
█ ██      ████   ██ ██   ██     ████  ████ ██   ██ ██      ██   ██ ██    ██ ██
█ █████   ██ ██  ██ ██   ██     ██ ████ ██ ███████ ██      ██████  ██    ██ ███████
█ ██      ██  ██ ██ ██   ██     ██  ██  ██ ██   ██ ██      ██   ██ ██    ██      ██
█ ███████ ██   ████ ██████      ██      ██ ██   ██  ██████ ██   ██  ██████  ███████
*/
