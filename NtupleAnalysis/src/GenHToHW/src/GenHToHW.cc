// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

// User
#include "Auxiliary/interface/Table.h"
#include "Auxiliary/interface/Tools.h"
#include "Tools/interface/MCTools.h"
#include "Tools/interface/DirectionalCut.h"
#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"

// ROOT
#include "TDirectory.h"
#include "Math/VectorUtil.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

struct PtComparator
{
  bool operator() (const genParticle p1, const genParticle p2) const { return ( p1.pt() > p2.pt() ); }
  bool operator() (const math::XYZTLorentzVector p1, const math::XYZTLorentzVector p2) const { return ( p1.pt() > p2.pt() ); }
};


class GenHToHW: public BaseSelector {
public:
  explicit GenHToHW(const ParameterSet& config, const TH1* skimCounters);
  virtual ~GenHToHW() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;
  virtual vector<genParticle> GetAllCopies(const vector<genParticle> genParticles, genParticle genP);
  virtual vector<GenJet> GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCut, std::vector<float> etaCut);
  virtual vector<GenJet> GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCut, std::vector<float> etaCut, vector<genParticle> genParticlesToMatch);
  virtual vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy=true, const bool hasNoDaughters=false);
  virtual double GetMt(const math::XYVector tau1, const math::XYVector tau2, const math::XYVector muon, const math::XYVector& met);
private:
  // Input parameters
  const double cfg_Verbose;
  const ParameterSet PSet_ElectronSelection;
  const double cfg_ElectronPtCut;
  const double cfg_ElectronEtaCut;
  const ParameterSet PSet_MuonSelection;
  const double cfg_MuonPtCut;
  const double cfg_MuonEtaCut;
  const ParameterSet PSet_TauSelection;
  const double cfg_TauPtCut;
  const double cfg_TauEtaCut;
  const ParameterSet PSet_JetSelection;
  const std::vector<float> cfg_JetPtCuts;
  const std::vector<float> cfg_JetEtaCuts;
  const DirectionalCut<int> cfg_JetNumberCut;
  const DirectionalCut<float> cfg_HtCut;
  const ParameterSet PSet_BJetSelection;
  const std::vector<float> cfg_BJetPtCuts;
  const std::vector<float> cfg_BJetEtaCuts;
  const DirectionalCut<int> cfg_BJetNumberCut;
  const ParameterSet PSet_METSelection;
  const DirectionalCut<float> cfg_METCut;
  const HistogramSettings cfg_PtBinSetting;
  const HistogramSettings cfg_EtaBinSetting;
  const HistogramSettings cfg_PhiBinSetting;
  const HistogramSettings cfg_MassBinSetting;
  const HistogramSettings cfg_DeltaEtaBinSetting;
  const HistogramSettings cfg_DeltaPhiBinSetting;
  const HistogramSettings cfg_DeltaRBinSetting;
  
  Tools auxTools;
  
  // Event Counters
  Count cAllEvents;
  Count cAllEventsAC;
  Count cTrigger;
  Count cElectronVeto;
  Count cMuonVeto;
  Count cTauVeto;
  Count cJetSelection;
  Count cBJetSelection;
  Count cMETSelection;
  // Count cTopologySelection;
  Count cTopSelection;
  Count cSelected;
  // BR Counters
  Count cInclusive;
  Count cbHt_HPlus;
  Count cbHt_HBoson;
  Count cbHt_WBoson;
  Count cbHt_tbW_BQuark;
  Count cHWh_Wqq_Quark;
  Count cHWh_Wqq_AntiQuark;
  Count cHWh_Wqq_Leptons;
  Count ctwb_Wqq_Quark;
  Count ctwb_Wqq_AntiQuark;
  Count ctwb_Wqq_Leptons;
  Count cHWh_hDiTau_tau;
  //Quark Counters
  Count clquark1;
  Count clquark2;
  Count clquark3;
  Count clquark4;
  Count csquark1;
  Count csquark2;
  Count csquark3;
  Count csquark4;

  // GenParticles
  WrappedTH1 *hParticle_Pt;
  WrappedTH1 *h_bHt_TQuark_Pt;
  WrappedTH1 *h_bHt_tbW_WBoson_Pt;
  WrappedTH1 *h_bHt_tbW_BQuark_Pt;
  WrappedTH1 *h_tbH_HPlus_Pt;
  WrappedTH1 *h_tbH_HWh_HBoson_Pt;
  WrappedTH1 *h_tbH_HWh_WBoson_Pt;
  WrappedTH1 *h_tbH_HWh_hDiTau_Tau_Pt;
  WrappedTH1 *h_hDiTau_TauDiMu_DiMu_Pt;
  WrappedTH1 *h_hDiTau_TauDiMu_Pion_Pt;
  WrappedTH1 *h_Wqq_Quark_Pt;
  WrappedTH1 *h_Wqq_AntiQuark_Pt;
  WrappedTH1 *h_HWh_Wmn_Muon_Pt;
  WrappedTH1 *h_tbH_BQuark_Pt;
  WrappedTH1 *h_leadBQuark_pt;
  WrappedTH1 *h_Wqq_Quarks_Pt;
  WrappedTH1 *h_Wln_Lepton_Pt;
  //
  WrappedTH1 *h_tbH_BQuark_dPhiTH_Cut1_Pt;
  WrappedTH1 *h_tbH_BQuark_dPhiTH_Cut2_Pt;
  WrappedTH2 *h_bHt_TQuark_HWh_HBoson_dPhi_Vs_Pt;
  //
  WrappedTH1 *h_bHt_TQuark_Eta;
  WrappedTH1 *h_bHt_tbW_WBoson_Eta;
  WrappedTH1 *h_bHt_tbW_BQuark_Eta;
  WrappedTH1 *h_tbH_HPlus_Eta;
  WrappedTH1 *h_tbH_HWh_HBoson_Eta;
  WrappedTH1 *h_tbH_HWh_WBoson_Eta;
  WrappedTH1 *h_tbH_HWh_hDiTau_Tau_Eta;
  WrappedTH1 *h_hDiTau_TauDiMu_DiMu_Eta;
  WrappedTH1 *h_hDiTau_TauDiMu_Pion_Eta;
  WrappedTH1 *h_Wqq_Quark_Eta;
  WrappedTH1 *h_Wqq_AntiQuark_Eta;
  WrappedTH1 *h_HWh_Wmn_Muon_Eta;
  WrappedTH1 *h_tbH_BQuark_Eta;
  WrappedTH1 *h_Wqq_Quarks_Eta;
  WrappedTH1 *h_Wln_Lepton_Eta;
  //
  WrappedTH1 *h_bHt_TQuark_Rap;
  WrappedTH1 *h_bHt_tbW_WBoson_Rap;
  WrappedTH1 *h_bHt_tbW_BQuark_Rap;
  WrappedTH1 *h_tbH_HPlus_Rap;
  WrappedTH1 *h_tbH_HWh_HBoson_Rap;
  WrappedTH1 *h_tbH_HWh_WBoson_Rap;
  WrappedTH1 *h_tbH_HWh_hDiTau_Tau_Rap;
  WrappedTH1 *h_hDiTau_TauDiMu_DiMu_Rap;
  WrappedTH1 *h_hDiTau_TauDiMu_Pion_Rap;
  WrappedTH1 *h_Wqq_Quark_Rap;
  WrappedTH1 *h_Wqq_AntiQuark_Rap;
  WrappedTH1 *h_HWh_Wmn_Muon_Rap;
  WrappedTH1 *h_tbH_BQuark_Rap;
  WrappedTH1 *h_Wqq_Quarks_Rap;
  WrappedTH1 *h_Wln_Lepton_Rap;
  //
  WrappedTH1 *h_tbH_HWh_WBoson_HBoson_dR;
  WrappedTH1 *h_tbH_HPlus_top_dR;
  WrappedTH1 *h_bHt_TQuark_HWh_HBoson_dR;
  WrappedTH1 *h_twb_WBoson_HWh_HBoson_dR;
  WrappedTH1 *h_twb_BQuark_HWh_HBoson_dR;
  WrappedTH1 *h_twb_TQuark_HWh_WBoson_dR;
  WrappedTH1 *h_twb_BQuark_HWh_WBoson_dR;
  WrappedTH1 *h_twb_WBoson_HWh_WBoson_dR;
  WrappedTH1 *h_twb_WBoson_twb_BQuark_dR;
  WrappedTH1 *h_Wqq_Quark_HWh_HBoson_dR;
  WrappedTH1 *h_Wqq_AntiQuark_HWh_HBoson_dR;
  WrappedTH1 *h_twb_Wqq_Quark_HWh_HBoson_dR;
  WrappedTH1 *h_twb_Wqq_AntiQuark_HWh_HBoson_dR;
  WrappedTH1 *h_twb_Wqq_Quark_AntiQuark_dR;
  WrappedTH1 *h_HWh_Wmn_Muon_HWh_HBoson_dR;
  WrappedTH1 *h_tbH_BQuark_HWh_HPlus_dR;
  WrappedTH1 *h_tbH_BQuark_twb_TQuark_dR;
  WrappedTH1 *h_Wqq_leadQuark_tbqq_Quark_dR;
  WrappedTH1 *h_tbH_BQuark_tbqq_leadQuark_dR;
  WrappedTH1 *h_tbH_BQuark_tbqq_subleadQuark_dR;
  WrappedTH1 *h_bHt_tbW_BQuark_tbH_BQuark_dR;
  WrappedTH1 *h_HWh_hDiTau_Tau_Tau_dR;
  WrappedTH1 *h_HWh_Wmuon_Muon_HWh_h_Tau_dR;
  //
  WrappedTH1 *h_tbH_HWh_WBoson_HBoson_dPhi;
  WrappedTH1 *h_tbH_HPlus_top_dPhi;
  WrappedTH1 *h_bHt_TQuark_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_WBoson_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_BQuark_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_TQuark_HWh_WBoson_dPhi;
  WrappedTH1 *h_twb_BQuark_HWh_WBoson_dPhi;
  WrappedTH1 *h_twb_WBoson_HWh_WBoson_dPhi;
  WrappedTH1 *h_twb_WBoson_twb_BQuark_dPhi;
  WrappedTH1 *h_Wqq_Quark_HWh_HBoson_dPhi;
  WrappedTH1 *h_Wqq_AntiQuark_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_Wqq_Quark_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_Wqq_AntiQuark_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_Wqq_Quark_AntiQuark_dPhi;
  WrappedTH1 *h_HWh_Wmn_Muon_HWh_HBoson_dPhi;
  WrappedTH1 *h_tbH_BQuark_HWh_HPlus_dPhi;
  WrappedTH1 *h_tbH_BQuark_twb_TQuark_dPhi;
  WrappedTH1 *h_Wqq_leadQuark_tbqq_Quark_dPhi;
  WrappedTH1 *h_HWh_WBoson_neut_HWh_HBoson_neut_Dphi;
  WrappedTH1 *h_HWh_HBoson_neut_Met_dPhi;
  WrappedTH1 *h_HWh_hDiTau_Tau_Tau_dPhi;
  WrappedTH1 *h_HWh_Wmuon_Muon_HWh_h_Tau_dPhi;
  //
  WrappedTH1 *h_tbH_HWh_WBoson_HBoson_dEta;
  WrappedTH1 *h_tbH_HPlus_top_dEta;
  WrappedTH1 *h_bHt_TQuark_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_WBoson_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_BQuark_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_TQuark_HWh_WBoson_dEta;
  WrappedTH1 *h_twb_BQuark_HWh_WBoson_dEta;
  WrappedTH1 *h_twb_WBoson_HWh_WBoson_dEta;
  WrappedTH1 *h_twb_WBoson_twb_BQuark_dEta;
  WrappedTH1 *h_Wqq_Quark_HWh_HBoson_dEta;
  WrappedTH1 *h_Wqq_AntiQuark_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_Wqq_Quark_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_Wqq_AntiQuark_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_Wqq_Quark_AntiQuark_dEta;
  WrappedTH1 *h_HWh_Wmn_Muon_HWh_HBoson_dEta;
  WrappedTH1 *h_tbH_BQuark_HWh_HPlus_dEta;
  WrappedTH1 *h_tbH_BQuark_twb_TQuark_dEta;
  WrappedTH1 *h_Wqq_leadQuark_tbqq_Quark_dEta;
  //
  WrappedTH1 *h_tbH_WBoson_HWh_HBoson_dRap;
  //
  WrappedTH1 *h_muonPtTauVsPt;
  WrappedTH1 *h_MET_Pt;
  WrappedTH1 *h_MET_Et;
  //
  WrappedTH2 *h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_tbH_HPlus_top_dEta_Vs_dPhi;
  WrappedTH2 *h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_WBoson_twb_BQuark_dEta_Vs_dPhi;
  WrappedTH2 *h_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_Wqq_Quark_AntiQuark_dEta_Vs_dPhi;
  WrappedTH2 *h_HWh_Wmn_Muon_HWh_HBoson_dEta_Vs_dPhi;
  //
  WrappedTH2 *h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_tbH_HPlus_top_dR_Vs_Pt;
  WrappedTH2 *h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_WBoson_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_BQuark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_TQuark_HWh_WBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_BQuark_HWh_WBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_WBoson_HWh_WBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_WBoson_twb_BQuark_dR_Vs_Pt;
  WrappedTH2 *h_Wqq_Quark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_Wqq_Quark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_Wqq_Quark_AntiQuark_dR_Vs_Pt;
  WrappedTH2 *HWh_Wmn_Muon_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_tbH_BQuark_HWh_HPlus_dR_Vs_Pt;
  WrappedTH2 *h_tbH_BQuark_twb_TQuark_dR_Vs_Pt;
  WrappedTH2 *h_Wqq_leadQuark_tbqq_Quark_dR_Vs_Pt;
  WrappedTH2 *h_tbH_BQuark_tbqq_leadQuark_dR_Vs_Pt;
  WrappedTH2 *h_tbH_BQuark_tbqq_subleadQuark_dR_Vs_Pt;
  WrappedTH2 *h_HWh_hDiTau_Tau_Tau_dR_Vs_Pt;
  //
  WrappedTH1 *h_HWh_hDiTau_Tau_MET_dPhi;
  WrappedTH1 *h_twb_Wqq_Quark_MET_dPhi;
  WrappedTH1 *h_HWh_Wmn_Muon_MET_dPhi;
  WrappedTH1 *h_tbH_BQuark_MET_dPhi;
  WrappedTH1 *h_twb_BQuark_MET_dPhi;
  WrappedTH1 *h_HWh_DiTau_TauNeut_dPhi;
  //
  WrappedTH1 *h_tbH_HPlus_Charge;
  //
  WrappedTH1 *h_InitialEOverFinalE;
  //
  WrappedTH1 *h_mT_ChargeHiggs_M;
  WrappedTH1 *h_mT_ChargeHiggsMeTTau_M;
  WrappedTH1 *h_mT_ChargeHiggsMeTTau_Cut1_M;
  WrappedTH1 *h_mT_ChargeHiggsMeTTau_Cut2_M;
  WrappedTH1 *h_mT_ChargeHiggs_assW_hadr_M;
  WrappedTH1 *h_mT_ChargeHiggs_assW_lept_M;
  WrappedTH1 *h_mT_ChargeHiggs_assMuon_M;
  WrappedTH1 *h_mT_ChargeHiggs_assqq_M;
  WrappedTH1 *h_mT_ChargeHiggs_taudPhils_M;
  WrappedTH1 *h_mT_ChargeHiggs_taudPhigr_M;
  WrappedTH1 *h_mT_ChargeHiggs_OneAnyMuon_M;
  WrappedTH1 *h_mT_ChargeHiggs_hptMuon_M;
  //
  WrappedTH2 *h_HWh_HBoson_neut_Met_dPhi_Vs_HWh_WBoson_neut_HWh_HBoson_neut_dPhi;
  //
  WrappedTH1 *h_leadquarknumber;
  WrappedTH1 *h_subleadquarknumber;
  WrappedTH1 *h_eventcount;
};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(GenHToHW);

GenHToHW::GenHToHW(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    cfg_Verbose(config.getParameter<bool>("verbose")),
    PSet_ElectronSelection(config.getParameter<ParameterSet>("ElectronSelection")),
    cfg_ElectronPtCut(config.getParameter<float>("ElectronSelection.electronPtCut")),  
    cfg_ElectronEtaCut(config.getParameter<float>("ElectronSelection.electronEtaCut")),
    PSet_MuonSelection(config.getParameter<ParameterSet>("MuonSelection")),
    cfg_MuonPtCut(config.getParameter<float>("MuonSelection.muonPtCut")),
    cfg_MuonEtaCut(config.getParameter<float>("MuonSelection.muonEtaCut")),
    PSet_TauSelection(config.getParameter<ParameterSet>("TauSelection")),
    cfg_TauPtCut(config.getParameter<float>("TauSelection.tauPtCut")),
    cfg_TauEtaCut(config.getParameter<float>("TauSelection.tauEtaCut")),
    PSet_JetSelection(config.getParameter<ParameterSet>("JetSelection")),
    cfg_JetPtCuts(config.getParameter<std::vector<float>>("JetSelection.jetPtCuts")),
    cfg_JetEtaCuts(config.getParameter<std::vector<float>>("JetSelection.jetEtaCuts")),
    cfg_JetNumberCut(config, "JetSelection.numberOfJetsCut"),
    cfg_HtCut(config, "JetSelection.HTCut"),
    PSet_BJetSelection(config.getParameter<ParameterSet>("BJetSelection")),
    cfg_BJetPtCuts(config.getParameter<std::vector<float>>("BJetSelection.jetPtCuts")),
    cfg_BJetEtaCuts(config.getParameter<std::vector<float>>("BJetSelection.jetEtaCuts")),
    cfg_BJetNumberCut(config, "BJetSelection.numberOfBJetsCut"),
    PSet_METSelection(config.getParameter<ParameterSet>("METSelection")),
    cfg_METCut(config, "METSelection.METCut"),
    // PSet_TopologySelection(config.getParameter<ParameterSet>("TopologySelection")),
    // cfg_SphericityCut(config, "TopologySelection.SphericityCut"),
    // cfg_AplanarityCut(config, "TopologySelection.AplanarityCut"),
    // cfg_PlanarityCut(config, "TopologySelection.PlanarityCut"),
    // cfg_CircularityCut(config, "TopologySelection.CircularityCut"),
    // cfg_Y23Cut(config, "TopologySelection.Y23Cut"),
    // cfg_CparameterCut(config, "TopologySelection.CparameterCut"),
    // cfg_DparameterCut(config, "TopologySelection.DparameterCut"),
    // cfg_FoxWolframMomentCut(config, "TopologySelection.FoxWolframMomentCut"),
    // cfg_AlphaTCut(config, "TopologySelection.AlphaTCut"),
    // cfg_CentralityCut(config, "TopologySelection.CentralityCut"),
    // PSet_TopSelection(config.getParameter<ParameterSet>("TopSelection")),
    cfg_PtBinSetting(config.getParameter<ParameterSet>("CommonPlots.ptBins")),
    cfg_EtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.etaBins")),
    cfg_PhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.phiBins")),
    cfg_MassBinSetting(config.getParameter<ParameterSet>("CommonPlots.invMassBins")),
    cfg_DeltaEtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaEtaBins")),
    cfg_DeltaPhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaPhiBins")),
    cfg_DeltaRBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaRBins")),
    cAllEvents(fEventCounter.addCounter("All events")),
    cAllEventsAC(fEventCounter.addCounter("All events AC")),
    cTrigger(fEventCounter.addCounter("Trigger")),
    cElectronVeto(fEventCounter.addCounter("e-veto")),
    cMuonVeto(fEventCounter.addCounter("#mu-veto")),
    cTauVeto(fEventCounter.addCounter("#tau-veto")),
    cJetSelection(fEventCounter.addCounter("Jets + H_{T}")),
    cBJetSelection(fEventCounter.addCounter("b-jets")),
    cMETSelection(fEventCounter.addCounter("MET")),
    // cTopologySelection(fEventCounter.addCounter("Topology")),
    cTopSelection(fEventCounter.addCounter("Top")),
    cSelected(fEventCounter.addCounter("All Selections")),
    cInclusive(fEventCounter.addSubCounter("Branching", "All events")),
    cbHt_HPlus(fEventCounter.addSubCounter("Branching", "H+")),
    cbHt_HBoson(fEventCounter.addSubCounter("Branching", "H+->HW, H")),
    cbHt_WBoson(fEventCounter.addSubCounter("Branching", "H+->HW, W")),
    cbHt_tbW_BQuark(fEventCounter.addSubCounter("Branching", "t->bW, b")),
    cHWh_Wqq_Quark(fEventCounter.addSubCounter("Branching", "H+->HW, W->qq, q")),
    cHWh_Wqq_AntiQuark(fEventCounter.addSubCounter("Branching", "H+->HW, W->qq, qbar")),
    cHWh_Wqq_Leptons(fEventCounter.addSubCounter("Branching", "H+->HW, W->l v")),
    ctwb_Wqq_Quark(fEventCounter.addSubCounter("Branching", "t->bW, W->qq, q")),
    ctwb_Wqq_AntiQuark(fEventCounter.addSubCounter("Branching", "t->bW, W->qq, qbar")),
    ctwb_Wqq_Leptons(fEventCounter.addSubCounter("Branching", "t->bW, W->l v")),
    cHWh_hDiTau_tau(fEventCounter.addSubCounter("Branching", "H+->HW, H->tau, tau")),
    clquark1(fEventCounter.addSubCounter("Branching", "lead BQuark (t->bW)")),
    clquark2(fEventCounter.addSubCounter("Branching", "lead Quark (W->qq)")),
    clquark3(fEventCounter.addSubCounter("Branching", "lead AntiQuark (W->qq)")),
    clquark4(fEventCounter.addSubCounter("Branching", "lead BQuark (gg->tbH)")),
    csquark1(fEventCounter.addSubCounter("Branching", "sub-lead BQuark (t->bW)")),
    csquark2(fEventCounter.addSubCounter("Branching", "sub-lead Quark (W->qq)")),
    csquark3(fEventCounter.addSubCounter("Branching", "sub-lead AntiQuark (W->qq)")),
    csquark4(fEventCounter.addSubCounter("Branching", "sub-lead BQuark (gg->tbH)"))
{ }

void GenHToHW::book(TDirectory *dir) {

  // Fixed binning
  const int nBinsPt   = 4*cfg_PtBinSetting.bins();
  const double minPt  = cfg_PtBinSetting.min();
  const double maxPt  = 4*cfg_PtBinSetting.max();

  const int nBinsEta  = 2*cfg_EtaBinSetting.bins();
  const double minEta = cfg_EtaBinSetting.min();
  const double maxEta = 2*cfg_EtaBinSetting.max();

  const int nBinsRap  = cfg_EtaBinSetting.bins();
  const double minRap = cfg_EtaBinSetting.min();
  const double maxRap = cfg_EtaBinSetting.max();

  // const int nBinsPhi  = cfg_PhiBinSetting.bins();
  // const double minPhi = cfg_PhiBinSetting.min();
  // const double maxPhi = cfg_PhiBinSetting.max();

  const int nBinsM  = cfg_MassBinSetting.bins();
  const double minM = cfg_MassBinSetting.min();
  const double maxM = cfg_MassBinSetting.max();
  
  const int nBinsdEta  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindEta = cfg_DeltaEtaBinSetting.min();
  const double maxdEta = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdRap  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindRap = cfg_DeltaEtaBinSetting.min();
  const double maxdRap = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdPhi  = cfg_DeltaPhiBinSetting.bins();
  const double mindPhi = cfg_DeltaPhiBinSetting.min();
  const double maxdPhi = 1.2*cfg_DeltaPhiBinSetting.max(); 

  const int nBinsdR  = 2*cfg_DeltaRBinSetting.bins();
  const double mindR = cfg_DeltaRBinSetting.min();
  const double maxdR = cfg_DeltaRBinSetting.max();
  
  const int nBinsToOne  = 0.2*2*cfg_DeltaEtaBinSetting.bins();
  const double minToOne = 0.0;
  const double maxToOne = 1.0;

  const int nBinsCha =  5;
  const double minCha   = -2.5;
  const double maxCha   =  2.5;

  TDirectory* th1 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH1");
  TDirectory* th2 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH2");
    
  // GenParticles  
  hParticle_Pt                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "AllParticle_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_TQuark_Pt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Pt"            , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tbW_WBoson_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_WBoson_Pt"        , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tbW_BQuark_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_BQuark_Pt"        , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HPlus_Pt                = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HWh_HBoson_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_HBoson_Pt"      , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HWh_WBoson_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_Pt"      , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HWh_hDiTau_Tau_Pt       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_hDiTau_Tau_Pt"  , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_hDiTau_TauDiMu_DiMu_Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_DiMu_Pt" , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_hDiTau_TauDiMu_Pion_Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_Pion_Pt" , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt); 
  h_Wqq_Quark_Pt                = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quark_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_Wqq_AntiQuark_Pt            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_AntiQuark_Pt"       , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_HWh_Wmn_Muon_Pt             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmn_Muon_Pt"        , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_BQuark_Pt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_BQuark_dPhiTH_Cut1_Pt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_dPhiTH_Cut1_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_BQuark_dPhiTH_Cut2_Pt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_dPhiTH_Cut2_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_leadBQuark_pt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_leadBQuark_pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_Wqq_Quarks_Pt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quarks_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_Wln_Lepton_Pt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wln_Lepton_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  //
  h_bHt_TQuark_Eta              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Eta"           , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tbW_WBoson_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_WBoson_Eta"       , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tbW_BQuark_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_BQuark_Eta"       , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HPlus_Eta               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_Eta"          , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HWh_HBoson_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_HBoson_Eta"     , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HWh_WBoson_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_Eta"     , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HWh_hDiTau_Tau_Eta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_hDiTau_Tau_Eta" , ";#eta", nBinsEta, minEta, maxEta);
  h_hDiTau_TauDiMu_DiMu_Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_DiMu_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_hDiTau_TauDiMu_Pion_Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_Pion_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_Wqq_Quark_Eta               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quark_Eta"          , ";#eta", nBinsEta, minEta, maxEta);
  h_Wqq_AntiQuark_Eta           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_AntiQuark_Eta"      , ";#eta", nBinsEta, minEta, maxEta);
  h_HWh_Wmn_Muon_Eta            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmn_Muon_Eta"       , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_BQuark_Eta              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_Eta"         , ";#eta", nBinsEta, minEta, maxEta);
  h_Wqq_Quarks_Eta              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quarks_Eta"         , ";#eta", nBinsEta, minEta, maxEta);
  h_Wln_Lepton_Eta              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wln_Lepton_Eta"         , ";#eta", nBinsEta, minEta, maxEta);
  //
  h_bHt_TQuark_Rap              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Rap"           , ";#omega", nBinsRap, minRap, maxRap);
  h_bHt_tbW_WBoson_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_WBoson_Rap"       , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_tbW_BQuark_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_BQuark_Rap"       , ";#omega", nBinsEta, minRap, maxRap);
  h_tbH_HPlus_Rap               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_Rap"          , ";#omega", nBinsRap, minRap, maxRap);
  h_tbH_HWh_HBoson_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_HBoson_Rap"     , ";#omega", nBinsRap, minRap, maxRap);  
  h_tbH_HWh_WBoson_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_Rap"     , ";#omega", nBinsRap, minRap, maxRap);	
  h_tbH_HWh_hDiTau_Tau_Rap      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_hDiTau_Tau_Rap" , ";#omega", nBinsRap, minRap, maxRap);
  h_hDiTau_TauDiMu_DiMu_Rap     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_DiMu_Rap", ";#omega", nBinsRap, minRap, maxRap);
  h_hDiTau_TauDiMu_Pion_Rap     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_Pion_Rap", ";#omega", nBinsRap, minRap, maxRap);
  h_Wqq_Quark_Rap               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quark_Rap"          , ";#omega", nBinsRap, minRap, maxRap);
  h_Wqq_AntiQuark_Rap           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_AntiQuark_Rap"      , ";#omega", nBinsRap, minRap, maxRap);
  h_HWh_Wmn_Muon_Rap            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmn_Muon_Rap"       , ";#omega", nBinsRap, minRap, maxRap);
  h_tbH_BQuark_Rap              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_Rap"         , ";#omega", nBinsRap, minRap, maxRap);
  h_Wqq_Quarks_Rap              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quarks_Rap"         , ";#omega", nBinsRap, minRap, maxRap);
  h_Wln_Lepton_Rap              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wln_Lepton_Rap"         , ";#omega", nBinsRap, minRap, maxRap);
  //
  h_tbH_HWh_WBoson_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_tbH_HPlus_top_dR            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_top_dR"              , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_TQuark_HWh_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_TQuark_HWh_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_WBoson_HWh_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_BQuark_HWh_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_TQuark_HWh_WBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_TQuark_HWh_WBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_BQuark_HWh_WBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_WBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_WBoson_HWh_WBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_WBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_WBoson_twb_BQuark_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_twb_BQuark_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_Wqq_Quark_HWh_HBoson_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quark_HWh_HBoson_dR"       , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_Wqq_AntiQuark_HWh_HBoson_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_AntiQuark_HWh_HBoson_dR"   , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_Wqq_Quark_HWh_HBoson_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_Quark_HWh_HBoson_dR"   , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_Wqq_AntiQuark_HWh_HBoson_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_AntiQuark_HWh_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_Wqq_Quark_AntiQuark_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_Quark_AntiQuark_dR"    , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_HWh_Wmn_Muon_HWh_HBoson_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmn_Muon_HWh_HBoson_dR"    , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_tbH_BQuark_HWh_HPlus_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_HWh_HPlus_dR"       , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_tbH_BQuark_twb_TQuark_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_twb_TQuark_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_Wqq_leadQuark_tbqq_Quark_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_leadQuark_tbqq_Quark_dR"   , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_tbH_BQuark_tbqq_leadQuark_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_tbqq_leadQuark_dR" , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_tbH_BQuark_tbqq_subleadQuark_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_tbqq_subleadQuark_dR", ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_tbW_BQuark_tbH_BQuark_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_tbW_BQuark_tbH_BQuark_dR" , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_HWh_hDiTau_Tau_Tau_dR       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_hDiTau_Tau_Tau_dR"         , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_HWh_Wmuon_Muon_HWh_h_Tau_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmuon_Muon_HWh_h_Tau_dR"   , ";#DeltaR", nBinsdR, mindR, maxdR);
  //
  h_tbH_HWh_WBoson_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_tbH_HPlus_top_dPhi          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_top_dPhi"         , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_TQuark_HWh_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_TQuark_HWh_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_TQuark_HWh_WBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_TQuark_HWh_WBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_WBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_WBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_WBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_WBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_twb_BQuark_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_twb_BQuark_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_Wqq_Quark_HWh_HBoson_dPhi   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quark_HWh_HBoson_dPhi"  , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi); 
  h_Wqq_AntiQuark_HWh_HBoson_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_AntiQuark_HWh_HBoson_dPhi" , ";#Delta#phi"  , nBinsdPhi, mindPhi, maxdPhi); 
  h_twb_Wqq_Quark_HWh_HBoson_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_Quark_HWh_HBoson_dPhi" , ";#Delta#phi"  , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_Wqq_AntiQuark_HWh_HBoson_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_AntiQuark_HWh_HBoson_dPhi" , ";#Delta#phi"     , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_Wqq_Quark_AntiQuark_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_Quark_AntiQuark_dPhi"   , ";#Delta#phi"  , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_Wmn_Muon_HWh_HBoson_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmn_Muon_HWh_HBoson_dPhi"   , ";#Delta#phi"  , nBinsdPhi, mindPhi, maxdPhi);
  h_tbH_BQuark_HWh_HPlus_dPhi   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_HWh_HPlus_dPhi" , ";#Delta#phi"        , nBinsdPhi, mindPhi, maxdPhi);
  h_tbH_BQuark_twb_TQuark_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_twb_TQuark_dPhi", ";#Delta#phi"        , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_hDiTau_Tau_MET_dPhi     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_hDiTau_Tau_MET_dPhi"   , ";#Delta#phi"        , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_Wqq_Quark_MET_dPhi      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_Quark_MET_dPhi"    , ";#Delta#phi"        , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_Wmn_Muon_MET_dPhi       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmn_Muon_MET_dPhi"     , ";#Delta#phi"        , nBinsdPhi, mindPhi, maxdPhi);
  h_tbH_BQuark_MET_dPhi         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_MET_dPhi"       , ";#Delta#phi"        , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_MET_dPhi         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_MET_dPhi"       , ";#Delta#phi"        , nBinsdPhi, mindPhi, maxdPhi);
  h_Wqq_leadQuark_tbqq_Quark_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_leadQuark_tbqq_Quark_dPhi", ";#Delta#phi"   , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_WBoson_neut_HWh_HBoson_neut_Dphi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_WBoson_neut_HWh_HBoson_neut_Dphi" , ";#Delta#phi" , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_HBoson_neut_Met_dPhi    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_HBoson_neut_Met_dPhi"       , ";#Delta#phi"   , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_DiTau_TauNeut_dPhi      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_DiTau_TauNeut_dPhi"       , ";#Delta#phi"     , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_hDiTau_Tau_Tau_dPhi     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_hDiTau_Tau_Tau_dPhi"       , ";#Delta#phi"    , nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_Wmuon_Muon_HWh_h_Tau_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmuon_Muon_HWh_h_Tau_dPhi" , ";#Delta#phi"    , nBinsdPhi, mindPhi, maxdPhi);
  //
  h_tbH_HWh_WBoson_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_tbH_HPlus_top_dEta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_top_dEta"         , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_bHt_TQuark_HWh_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_TQuark_HWh_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_WBoson_HWh_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_BQuark_HWh_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_TQuark_HWh_WBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_TQuark_HWh_WBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_BQuark_HWh_WBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_WBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_WBoson_HWh_WBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_WBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_WBoson_twb_BQuark_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_twb_BQuark_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_Wqq_Quark_HWh_HBoson_dEta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_Quark_HWh_HBoson_dEta"  , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_Wqq_AntiQuark_HWh_HBoson_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_AntiQuark_HWh_HBoson_dEta" , ";#Delta#eta"  , nBinsdEta, mindEta, maxdEta);
  h_twb_Wqq_Quark_HWh_HBoson_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_Quark_HWh_HBoson_dEta" , ";#Delta#eta"  , nBinsdEta, mindEta, maxdEta);
  h_twb_Wqq_AntiQuark_HWh_HBoson_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_AntiQuark_HWh_HBoson_dEta" , ";#Delta#eta"      , nBinsdEta, mindEta, maxdEta);
  h_twb_Wqq_Quark_AntiQuark_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_Wqq_Quark_AntiQuark_dEta"   , ";#Delta#eta"  , nBinsdEta, mindEta, maxdEta);
  h_HWh_Wmn_Muon_HWh_HBoson_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HWh_Wmn_Muon_HWh_HBoson_dEta"   , ";#Delta#eta"  , nBinsdEta, mindEta, maxdEta);
  h_tbH_BQuark_HWh_HPlus_dEta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_HWh_HPlus_dEta"  , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_tbH_BQuark_twb_TQuark_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_BQuark_twb_TQuark_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_Wqq_leadQuark_tbqq_Quark_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_Wqq_leadQuark_tbqq_Quark_dEta", ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
    //
  h_tbH_WBoson_HWh_HBoson_dRap  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_WBoson_HWh_HBoson_dRap" , ";#Delta#omega"     , nBinsdRap, mindRap, maxdRap);
  //
  h_muonPtTauVsPt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_muonPtTauVsPt" , ";p_{T}^{#tau} (GeV/c)" , nBinsToOne, minToOne, maxToOne);
  h_MET_Pt                      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_MET_Pt"                     , ";p_{T} (GeV/c)"    , nBinsPt, minPt, maxPt);
  h_MET_Et                      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_MET_Et"                     , ";E_{T} (GeV/c)"    , nBinsPt, minPt, maxPt);
  //
  h_tbH_HPlus_Charge            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_Charge", ""                             , nBinsCha, minCha, maxCha);          
  //
  h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_tbH_HPlus_top_dEta_Vs_dPhi         = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_HPlus_top_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)"        , nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)" , nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_twb_BQuark_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_twb_BQuark_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)" , nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_Wqq_Quark_AntiQuark_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_Wqq_Quark_AntiQuark_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_HWh_Wmn_Muon_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_HWh_Wmn_Muon_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  //
  h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_tbH_HPlus_top_dR_Vs_Pt             = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_HPlus_top_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_WBoson_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_BQuark_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_TQuark_HWh_WBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_TQuark_HWh_WBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_BQuark_HWh_WBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_WBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_WBoson_HWh_WBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_WBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_WBoson_twb_BQuark_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_twb_BQuark_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_Wqq_Quark_HWh_HBoson_dR_Vs_Pt      = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_Wqq_Quark_HWh_HBoson_dR_Vs_Pt" , ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_Wqq_Quark_HWh_HBoson_dR_Vs_Pt  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_Wqq_Quark_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_Wqq_Quark_AntiQuark_dR_Vs_Pt   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_Wqq_Quark_AntiQuark_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  HWh_Wmn_Muon_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "HWh_Wmn_Muon_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_tbH_BQuark_HWh_HPlus_dR_Vs_Pt      = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_BQuark_HWh_HPlus_dR_Vs_Pt" , ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_tbH_BQuark_twb_TQuark_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_BQuark_twb_TQuark_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_Wqq_leadQuark_tbqq_Quark_dR_Vs_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_Wqq_leadQuark_tbqq_Quark_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_tbH_BQuark_tbqq_leadQuark_dR_Vs_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_BQuark_tbqq_leadQuark_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_tbH_BQuark_tbqq_subleadQuark_dR_Vs_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_BQuark_tbqq_subleadQuark_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_HWh_hDiTau_Tau_Tau_dR_Vs_Pt        = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_HWh_hDiTau_Tau_Tau_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  //
  h_bHt_TQuark_HWh_HBoson_dPhi_Vs_Pt   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_bHt_TQuark_HWh_HBoson_dPhi_Vs_Pt", ";#Delta#phi;p_{T} (GeV/c)", nBinsdPhi, mindPhi, maxdPhi, nBinsPt, minPt, maxPt);
  //
  h_InitialEOverFinalE                 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_InitialEOverFinalE" , ";E_{i}/E_{f}" , nBinsToOne, minToOne, maxToOne);
  //
  h_mT_ChargeHiggs_M                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_M"             , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggsMeTTau_M            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggsMeTTau_M"       , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggsMeTTau_Cut1_M       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggsMeTTau_Cut1_M"  , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggsMeTTau_Cut2_M       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggsMeTTau_Cut2_M"  , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_assW_hadr_M        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_assW_hadr_M"   , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_assW_lept_M        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_assW_lept_M"   , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_assMuon_M          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_assMuon_M"     , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_assqq_M            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_assqq_M"       , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_taudPhils_M        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_taudPhils_M"   , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_taudPhigr_M        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_taudPhigr_M"   , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_OneAnyMuon_M       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_OneAnyMuon_M"  , ";M_{T} (GeV)" , 200, minM, 2000.);
  h_mT_ChargeHiggs_hptMuon_M          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_mT_ChargeHiggs_hptMuon_M"     , ";M_{T} (GeV)" , 200, minM, 2000.);
  //
  h_HWh_HBoson_neut_Met_dPhi_Vs_HWh_WBoson_neut_HWh_HBoson_neut_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_HWh_HBoson_neut_Met_dPhi_Vs_HWh_WBoson_neut_HWh_HBoson_neut_dPhi", ";#Delta#phi_{#nu_{#mu}-#nu_{#tau}};#Delta#phi_{MeT-#nu_{#tau}}", nBinsdPhi, mindPhi, maxdPhi,nBinsdPhi, mindPhi, maxdPhi);
  //
  h_leadquarknumber                    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_leadquarknumber"  ,  ";LeadQuark" , 6, -0.5, 5.5);
  h_subleadquarknumber                 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_subleadquarknumber"  ,  ";SubLeadQuark" , 6, -0.5, 5.5);
  h_eventcount                         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_eventcount"  ,  ";SubLeadQuark" , 4, -0.5, 3.5);
  return;
}

void GenHToHW::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
}




void GenHToHW::process(Long64_t entry) {

  if ( !fEvent.isMC() ) return;
  
  // Create MCools object
  MCTools mcTools(fEvent);
  
  // Increment Counter
  cAllEvents.increment();

  //================================================================================================
  // 1) Apply trigger
  //================================================================================================
  // if (cfg_Verbose) std::cout << "=== Trigger" << std::endl;
  // if ( !(fEvent.passTriggerDecision()) ) return;
  // cTrigger.increment();


  //================================================================================================
  // 2) MET filters (to remove events with spurious sources of fake MET)       
  //================================================================================================


  //================================================================================================
  // 3) Primarty Vertex (Check that a PV exists)
  //================================================================================================


  //================================================================================================
  // 4) Electron veto (fully hadronic + orthogonality)  
  //================================================================================================
  //  if (cfg_Verbose) std::cout << "=== Electron veto" << std::endl;
  //vector<genParticle> selectedElectrons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_ElectronPtCut, cfg_ElectronEtaCut, 11, true, false);
  //if (0)
  // {
  //  std::cout << "\nnElectrons = " << selectedElectrons.size() << std::endl;
  //  for (auto& p: selectedElectrons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
  //}
  //if ( selectedElectrons.size() > 0 ) return;
  //cElectronVeto.increment();


  //================================================================================================
  // 5) Muon selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Muon selection" << std::endl;
  vector<genParticle> selectedMuons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_MuonPtCut, cfg_MuonEtaCut, 13, true, false);
  if (0)
    {
      std::cout << "nMuons = " << selectedMuons.size() << std::endl;
      for (auto& p: selectedMuons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    } 
  if ( selectedMuons.size() < 1 ) return;
  cMuonVeto.increment();

  //================================================================================================
  // 6) Tau selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Tau selection" << std::endl;
  vector<genParticle> selectedTaus = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_TauPtCut, cfg_TauEtaCut, 15, true, false);
  if (0) 
    {
      std::cout << "nTaus = " << selectedTaus.size() << std::endl;
      for (auto& p: selectedTaus) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }
  if ( selectedTaus.size() < 2 ) return;
  cTauVeto.increment();

  /*
  //================================================================================================
  // 7) Jet Selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Jet Selection" << std::endl;
  vector<GenJet> selectedJets = GetGenJets(fEvent.genjets(), cfg_JetPtCuts, cfg_JetEtaCuts);
  if (0)
    {
      std::cout << "nJets = " << selectedJets.size() << std::endl;
      for (auto& p: selectedJets) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
    }
  if (!cfg_JetNumberCut.passedCut(selectedJets.size())) return;

  // HT Selection
  double genJ_HT = 0.0;
  std::vector<math::XYZTLorentzVector> selJets_p4;
  math::XYZTLorentzVector jet_p4;
  for(auto jet: selectedJets) 
    {
        jet_p4 = jet.p4();
	genJ_HT += jet.pt();
	selJets_p4.push_back( jet_p4 );
    }

  if ( !cfg_HtCut.passedCut(genJ_HT) ) return;
  cJetSelection.increment();
  */


  /*
  //================================================================================================
  // 8) BJet Selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== BJet Selection" << std::endl;
  vector<genParticle> selectedBQuarks = GetGenParticles(fEvent.genparticles().getGenParticles(), 10, 3, 5, true, false);
  std::sort( selectedBQuarks.begin(), selectedBQuarks.end(), PtComparator() );
  if (0) 
    {
      std::cout << "nBQuarks = " << selectedBQuarks.size() << std::endl;
      // for (auto& p: selectedBQuarks) mcTools.PrintGenParticle(p);
      for (auto& p: selectedBQuarks) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }

  // Match b-quarks with GenJets
  vector<GenJet> selectedBJets = GetGenJets(selectedJets, cfg_BJetPtCuts, cfg_BJetEtaCuts, selectedBQuarks);
  if (0) 
    {
      std::cout << "nBJets = " << selectedBJets.size() << std::endl;
      for (auto& p: selectedBJets) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
      std::cout << "" << std::endl;
    }

  // Get selected jets excluding the matched bjets
  bool isBJet = false;
  std::vector<math::XYZTLorentzVector> selJets_NoBJets_p4;
  // For-loop: Selected jets
  for (auto& jet: selectedJets) 
    {
      isBJet = false;
	    
      // For-loop: Selected bjets
      for (auto& bjet: selectedBJets) 
	{
	  double dR = ROOT::Math::VectorUtil::DeltaR(jet.p4(), bjet.p4());
	  if (dR < 0.01) isBJet = true;
	}
      if (isBJet) continue;
      jet_p4 = jet.p4();
      selJets_NoBJets_p4.push_back(jet_p4);
    }

  if (!cfg_BJetNumberCut.passedCut(selectedBJets.size())) return;
  cBJetSelection.increment();
  */

  //================================================================================================
  // 9) BJet SF
  //================================================================================================
  
  
  //================================================================================================
  // 10) MET selection 
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== MET Selection" << std::endl;
  if (!cfg_METCut.passedCut(fEvent.genMET().et())) return;
  if (0) std::cout << "=== MET = " << fEvent.genMET().et() << std::endl;      
  cMETSelection.increment();
  
  
  //================================================================================================
  // 11) Topology selection 
  //================================================================================================
  
  /*
  //================================================================================================
  // 12) Top selection 
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Top Selection" << std::endl;
  cTopSelection.increment();
  */

  //================================================================================================
  // All Selections
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();


  bool bSkipEvent = false;
  // 4-momenta                                                                                                                                                               
  math::XYZTLorentzVector Htb_HPlus_p4;
  math::XYZTLorentzVector HWh_HBoson_p4;
  math::XYZTLorentzVector HWh_WBoson_p4;
  math::XYZTLorentzVector HWh_hDiTau_p4;
  math::XYZTLorentzVector bHt_TQuark_p4;
  math::XYZTLorentzVector bHt_tbW_WBoson_p4;
  math::XYZTLorentzVector bHt_tbW_BQuark_p4;
  std::vector<math::XYZTLorentzVector> v_ditau_htt_p4;
  std::vector<math::XYZTLorentzVector> v_dimu_htt_p4;
  std::vector<math::XYZTLorentzVector> v_pion_htt_p4;
  std::vector<math::XYZTLorentzVector> v_HWh_Wqq_Quark_p4;
  std::vector<math::XYZTLorentzVector> v_HWh_Wqq_AntiQuark_p4;
  std::vector<math::XYZTLorentzVector> v_Wqq_Quark_p4;
  std::vector<math::XYZTLorentzVector> v_Wqq_AntiQuark_p4;
  std::vector<math::XYZTLorentzVector> v_twb_Wqq_Quark_p4;
  std::vector<math::XYZTLorentzVector> v_twb_Wqq_AntiQuark_p4;
  math::XYZTLorentzVector HWh_hDiTau_HadrTau_p4;
  math::XYZTLorentzVector HWh_hDiTau_TauDiMu_p4;
  math::XYZTLorentzVector Wqq_Quark_p4;
  math::XYZTLorentzVector Wqq_AnriQuark_p4;
  math::XYZTLorentzVector Wmuon_p4;
  math::XYZTLorentzVector TopWmuon_p4;
  math::XYZTLorentzVector Q1_p4;
  math::XYZTLorentzVector l1_p4;
  std::vector<math::XYZTLorentzVector> v_Q1_p4;
  std::vector<math::XYZTLorentzVector> v_l1_p4;
  std::vector<math::XYZTLorentzVector> v_twb_Wmuon_Muon_p4;
  std::vector<math::XYZTLorentzVector> v_HWh_Wmuon_Muon_p4;
  math::XYZTLorentzVector neutrino_p4;
  std::vector<math::XYZTLorentzVector> v_neutrino_p4;
  math::XYZTLorentzVector assoBQuark_p4;
  std::vector<math::XYZTLorentzVector> v_assoBQuark_p4;
  //math::XYZTLorentzVector visibleTau_p4;
  std::vector<math::XYZTLorentzVector> v_visibleTau_p4;
  std::vector<math::XYVector> v_visibleTau_p2;
  math::XYVector visibleTau_p2;
  std::vector<float> v_muonPtVsTauPt;
  //math::XYZTLorentzVector tauvish_p4;
  std::vector<math::XYZTLorentzVector> v_HiggsMom_p4;
  math::XYZTLorentzVector HiggsMom_p4;
  std::vector<math::XYZTLorentzVector>  v_twb_topproducts;
  math::XYZTLorentzVector neut_from_W_p4;
  math::XYZTLorentzVector neut_from_tau_p4;
  std::vector<math::XYZTLorentzVector> v_neut_from_tau_p4;
  std::vector<genParticle> g_ditau_htt;
  std::vector<genParticle> g_Wmuon;
  std::vector<genParticle> g_assWmuon;
  std::vector<genParticle> g_neut_from_tau;
  std::vector<genParticle> g_assWQuark;
  std::vector<genParticle> g_assWAntiQuark;
  std::vector<genParticle> g_muon_p4;
  math::XYZTLorentzVector muon_p4;
  std::vector<math::XYZTLorentzVector> v_muon_p4;

  bool assWHadr = false;
  // Define the table                                                                                                                                                                                     
  Table table("Evt | Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | D0 (mm) | Lxy (mm) | Mom | Daughters", "Text"); //LaTeX or Text
 
  int row = 0;
  
  int qourkc = 0;
  // For-loop: GenParticles
  for (auto& p: fEvent.genparticles().getGenParticles()) {
        
    hParticle_Pt -> Fill(p.pt());
    //int genMom_index = -1; 
    //double genMom_pdgId = -999.99; 
    // Particle properties
    int genP_pdgId       = p.pdgId();
    
    short genP_index     = p.index();
    int genP_status      = p.status();
    double genP_pt       = p.pt();
    double genP_eta      = p.eta();
    double genP_phi      = p.phi();
    double genP_energy   = p.e(); 
    int genP_charge      = p.charge();
    
    // Associated genParticles
    std::vector<genParticle> genP_daughters;
    for (unsigned int i=0; i < p.daughters().size(); i++) genP_daughters.push_back(fEvent.genparticles().getGenParticles()[p.daughters().at(i)]);
    std::vector<genParticle> genP_mothers;
    for (unsigned int i=0; i < p.mothers().size(); i++) genP_mothers.push_back(fEvent.genparticles().getGenParticles()[p.mothers().at(i)]);
    std::vector<genParticle> genP_grandMothers;
    std::vector<genParticle> genP_grandDaughters;
    std::vector<genParticle> genP_sisters;
    std::vector<genParticle> genP_allCopies = GetAllCopies(fEvent.genparticles().getGenParticles(), p);
    std::vector<unsigned int> genGMoms_index;
    std::vector<unsigned int> genGMoms_pdgId;
    std::vector<unsigned int> genMoms_index;
    std::vector<unsigned int> genMoms_pdgId;
    std::vector<unsigned int> genDaus_index;
    std::vector<unsigned int> genDaus_pdgId;
    

    // Removed real vertex from Tree (to save size)                                                                                                                                                       
    ROOT::Math::XYZPoint vtxIdeal;
    vtxIdeal.SetXYZ(0, 0, 0);
    double genP_vtxX = vtxIdeal.X(); // p.vtxX()*10; // in mm                                                                                                                                             
    double genP_vtxY = vtxIdeal.Y(); // p.vtxY()*10; // in mm                                                                                                                                              
    double genP_vtxZ = vtxIdeal.Z(); // p.vtxZ()*10; // in mm


    // Daughter, Mom and Grand-mom properties                                                                                                                              
   
    genParticle m;
    genParticle g;
    genParticle d;
    genParticle s;
    genParticle firstMom;
    genParticle lastMom;
    // genParticle mopart;

    if (genP_daughters.size() > 0) d = genP_daughters.at(0);
    if (p.mothers().size() > 0)
      {
        m = genP_mothers.at(0);                                                                                                                                    
	for (unsigned int i=0; i < m.mothers().size(); i++) genP_grandMothers.push_back(fEvent.genparticles().getGenParticles()[m.mothers().at(i)]);
        if (m.mothers().size() > 0) g = genP_grandMothers.at(0);
	for (unsigned int i=0; i < m.daughters().size(); i++) genP_sisters.push_back(fEvent.genparticles().getGenParticles()[m.daughters().at(i)]);
	if (m.daughters().size() > 0) s = genP_sisters.at(0);
      }
  	// For convenience, save the pdgIds in vectors                                                                                                                          

    for (unsigned int i=0; i < genP_grandMothers.size(); i++)
      {
        if (genP_grandMothers.at(i).index() < 0) continue;
        genGMoms_index.push_back(genP_grandMothers.at(i).index());
        genGMoms_pdgId.push_back(genP_grandMothers.at(i).pdgId());
      }
    for (unsigned int i=0; i < genP_mothers.size(); i++)
      {
        if (genP_mothers.at(i).index() < 0) continue;
        genMoms_index.push_back(genP_mothers.at(i).index());
        genMoms_pdgId.push_back(genP_mothers.at(i).pdgId());
      }
    for (unsigned int i=0; i < genP_daughters.size(); i++)
      {
        if (genP_daughters.at(i).index() < 0) continue;
        genDaus_index.push_back(genP_daughters.at(i).index());
        genDaus_pdgId.push_back(genP_daughters.at(i).pdgId());
      }
    
    if (genP_mothers.size() > 0){                                                                                                                                                                          
      firstMom = genP_mothers.at(0);                                                                                                                                                          
      lastMom = genP_mothers.at(genP_mothers.size()-1);                                                                                                                                       
    }  

    // Properties that need to be calculated                                                                                                                               
    bool bIsLastCopy  = std::find(genDaus_pdgId.begin(), genDaus_pdgId.end(), genP_pdgId) == genDaus_pdgId.end();
    double genP_Lxy  = 0.0;
    double genP_d0   = 0.0;
    if (genP_daughters.size() > 0 && genP_mothers.size() > 0)
      {
        genP_d0  = mcTools.GetD0 (p, m, d, vtxIdeal); // in mm                                                                                                                                            
        genP_Lxy = mcTools.GetLxy(p, m, d, vtxIdeal); // in mm                                                                                                                                             
      }




    // Print genParticle properties or decay tree ?                                                                                                                                                       
    if (0)
      {
        mcTools.PrintGenParticle(p);
        mcTools.PrintGenDaughters(p);
      }
    if (0)
      {
	// Add table rows
	table.AddRowColumn(row, auxTools.ToString(entry)           );
	table.AddRowColumn(row, auxTools.ToString(genP_index)      );
	table.AddRowColumn(row, auxTools.ToString(genP_pdgId)      );
	table.AddRowColumn(row, auxTools.ToString(genP_status)     );
	table.AddRowColumn(row, auxTools.ToString(genP_charge)     );
	table.AddRowColumn(row, auxTools.ToString(genP_pt , 3)     );
	table.AddRowColumn(row, auxTools.ToString(genP_eta, 4)     );
	table.AddRowColumn(row, auxTools.ToString(genP_phi, 3)     );
	table.AddRowColumn(row, auxTools.ToString(genP_energy, 3)  );
	table.AddRowColumn(row, "(" + auxTools.ToString(genP_vtxX, 3) + ", " + auxTools.ToString(genP_vtxY, 3)  + ", " + auxTools.ToString(genP_vtxZ, 3) + ")" );
	table.AddRowColumn(row, auxTools.ToString(genP_d0 , 3)     );
	table.AddRowColumn(row, auxTools.ToString(genP_Lxy, 3)     );
	if (genMoms_index.size() < 6)
	  {
	    table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genMoms_index) );
	    // table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genMoms_pdgId) );                                                                                                                    
	  }
	else table.AddRowColumn(row, ".. Too many .." );
	if (genDaus_index.size() < 6)
	  {
	    table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genDaus_index) );
	    // table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genDaus_pdgId) );                                                                                                                 
	  }
	else table.AddRowColumn(row, ".. Too many .." );
	row++;
      }
    

    if (0) std::cout << ".. Collection of Charge Higgs .." << std::endl;
    if (abs(genP_pdgId) == 37) 
      {
	if (!bIsLastCopy) continue;
	int HPMomId = mcTools.RecursivelyLookForMotherId(p);
        if (0) std::cout << "HPMomId= " << HPMomId << std::endl;
	cbHt_HPlus.increment();
	Htb_HPlus_p4 = p.p4();
	int cHc = genP_charge;
	h_tbH_HPlus_Charge -> Fill(cHc);
       	// For-loop: All daughters 
	for (auto& m: genP_mothers)
	  {
	    HiggsMom_p4 = m.p4();
	    v_HiggsMom_p4.push_back(HiggsMom_p4);
	  }

	for (auto& d: genP_daughters)
          {
            if ( abs(d.pdgId()) == 25) // h0 = H0_{1}                                                                                                                        
              {
                cbHt_HBoson.increment();
		HWh_HBoson_p4 = d.p4();
              }
            else if ( abs(d.pdgId()) == 35) // H0 = H0_{2}                                                                                                                   
              {
                cbHt_HBoson.increment();
		HWh_HBoson_p4 = d.p4();
              }
            else if ( abs(d.pdgId()) == 24)// W+ or W-                                                                                                                       
              {		
                cbHt_WBoson.increment();
                HWh_WBoson_p4 = d.p4();
	      }
            else
              {
		std::cout << "*** KinematicsHToHW::process() H- daughter is " << d.pdgId() << " (H- bIsLastCopy = " << bIsLastCopy << ")" << std::endl;
                
              }

          }// for (auto& d: genP_daughters)  
      } // if (abs(genP_pdgId) == 37)

    if (0) std::cout << "=== Associated top" << std::endl;
    if ( abs(genP_pdgId) == 6)
      {
        if (!bIsLastCopy) continue;
	int TopMomId = mcTools.RecursivelyLookForMotherId(p);
        if (0) std::cout << "TopMomId= " << TopMomId << std::endl;
        // Fill histos                                                                                                                                                                                     
        if (0) std::cout << "TopIndx=" << p.index() << "---->" << "TopId=" << genP_pdgId  << std::endl;
        bHt_TQuark_p4 = p.p4();
        // For-loop: All daughters                                                                                                                                                                         
        for (auto& d: genP_daughters)
          {
            if ( abs(d.pdgId()) == 5)
              {
                if (0) std::cout << "TopIndx2=" << p.index() << "---->"  << "TopId2=" << genP_pdgId << "---->" << d.pdgId() << std::endl;
                bHt_tbW_BQuark_p4 = d.p4();
		v_twb_topproducts.push_back(bHt_tbW_BQuark_p4);
		/*
		if ( d.pt()  < cfg_MuonPtCut )
		  {
		    bSkipEvent = true;
		  }
		if ( abs(d.eta()) > cfg_MuonEtaCut) 
		  {
		    bSkipEvent = true;
		  }
		*/
	      }
            else if ( abs(d.pdgId()) == 24)
              {
                if (0) std::cout << "TopIndx3=" << p.index() << "---->" << "TopId3=" << genP_pdgId << "---->" << d.pdgId() << std::endl;
                // NOTE: t->Wb dominant, but t->Ws and t->Wd also possible!                                                                                                                                
                bHt_tbW_WBoson_p4 = d.p4();
              }
            //else bSkipEvent = true;
          } // for (auto& d: genP_daughters)                                                                                                                                                               
      }// if (abs(genP_pdgId) == 6)  

    if ( abs(genP_pdgId) == 24)
      {
        if (!bIsLastCopy) continue;
	for (auto& d: genP_daughters)
	  {
	    if (mcTools.IsQuark(d.pdgId()))
	      {
		Q1_p4 = d.p4();
		v_Q1_p4.push_back(Q1_p4);
	      }
	    else if ( mcTools.IsLepton(d.pdgId()))
	      {
		l1_p4 = d.p4();
		v_l1_p4.push_back(l1_p4);
	      }
	  }
      }
    
    if (0) std::cout << "=== W+ or W-" << std::endl;
    if ( abs(genP_pdgId) == 24)
      {  	
	if (!bIsLastCopy) continue;
	int mompdgId = mcTools.RecursivelyLookForMotherId(p);
	if (0) std::cout << "m Id= " << mompdgId  <<  std::endl;
       	if (abs(mompdgId) == 6) 
	  {	
	    if (0) std::cout << "genP_daughters size from top= " << genP_daughters.size()  <<  std::endl;
	    for (auto& d: genP_daughters)
	      {
		int dId = d.pdgId();
		if (d.pdgId() == 10 )
		  {
		    dId = d.pdgId()/10;
		  }
		else if (d.pdgId() > 19) 
		  {
		    dId = d.pdgId()/10;
		  }
		if (0) std::cout << "Quark from top Id= " << dId  <<  std::endl;
		if (mcTools.IsQuark(dId)) // t->Wb, W+->qq                                                                                                                                
		  {
		    assWHadr = true;
		    if (dId > 0)
		      {
			// Quarks
			Wqq_Quark_p4 = d.p4();
			v_Wqq_Quark_p4.push_back(Wqq_Quark_p4);
			v_twb_topproducts.push_back(Wqq_Quark_p4);
			ctwb_Wqq_Quark.increment();
			v_twb_Wqq_Quark_p4.push_back(Wqq_Quark_p4);
			g_assWQuark.push_back(d);
			/*
			if ( d.pt()  < cfg_MuonPtCut ) 
			  {
			    bSkipEvent = true;
			  }			
			if ( abs(d.eta()) > cfg_MuonEtaCut)
			  {
			    bSkipEvent = true;
			  }
			*/
		      }
		    else
		      {
			// Anti-Quarks
			Wqq_AnriQuark_p4 = d.p4();
			v_Wqq_AntiQuark_p4.push_back(Wqq_AnriQuark_p4);
			v_twb_topproducts.push_back(Wqq_AnriQuark_p4);
			ctwb_Wqq_AntiQuark.increment();
			v_twb_Wqq_AntiQuark_p4.push_back(Wqq_AnriQuark_p4);
			g_assWAntiQuark.push_back(d);
			/*
			if ( d.pt()  < cfg_MuonPtCut )
                          {
                            bSkipEvent = true;
                          }
                        if ( abs(d.eta()) > cfg_MuonEtaCut)
                          {
                            bSkipEvent = true;
                          }
			*/
		      }
		  }//if ( mcTools.IsQuark(d.pdgId()) )
		else if ( mcTools.IsLepton(dId) ) // H+->tb, t->bW+, W+->l v                                                                                                                       
		  {
		    ctwb_Wqq_Leptons.increment();
		    if (abs(dId) == 13)
                      {
			g_assWmuon.push_back(d);
			TopWmuon_p4 = d.p4();
			v_twb_Wmuon_Muon_p4.push_back(TopWmuon_p4);
		      }
		  }//else if ( mcTools.IsLepton(d.pdgId() ) ) H+->Wh, W+->l v
	      }//for (auto& d: genP_daughters)
	  }
	else if (abs(mompdgId) == 37)
	  {
	    if (0) std::cout << "genP_daughters size from H+= " << genP_daughters.size()  <<  std::endl;
	    for (auto& d: genP_daughters)
              {
                int dId = d.pdgId();
                if (d.pdgId() == 10 )
                  {
                    dId = d.pdgId()/10;
                  }
                else if (d.pdgId() > 19)
                  {
                    dId = d.pdgId()/10;
                  }
		if (0) std::cout << "Quark from H+ Id= " << dId <<  std::endl;
                if (mcTools.IsQuark(dId)) // t->Wb, W+->qq                                                                                                                                                 
                  {
		    //bSkipEvent = true;
                    if (dId > 0)
                      {
                        // Quarks                                                                                                                                                                          
                        Wqq_Quark_p4 = d.p4();
                        v_Wqq_Quark_p4.push_back(Wqq_Quark_p4);
                        cHWh_Wqq_Quark.increment();
			v_HWh_Wqq_Quark_p4.push_back(Wqq_Quark_p4);
                      }
		    else
		      {
			// Anti-Quarks                                                                                                                                                                     
			Wqq_AnriQuark_p4 = d.p4();
			v_Wqq_AntiQuark_p4.push_back(Wqq_AnriQuark_p4);
			cHWh_Wqq_AntiQuark.increment();
			v_HWh_Wqq_AntiQuark_p4.push_back(Wqq_AnriQuark_p4);
		      }
		  }//if ( mcTools.IsQuark(d.pdgId()) )                                                                                                                                                     
		else if ( mcTools.IsLepton(dId ) ) // H+->tb, t->bW+, W+->l v                                                                                                                              
                  {
		    cHWh_Wqq_Leptons.increment();
		    if (abs(dId) == 13)
		      {
			if ( d.pt()  < 5.0 ) continue;
			if ( abs(d.eta()) > cfg_MuonEtaCut) continue;
			Wmuon_p4 = d.p4();
			v_HWh_Wmuon_Muon_p4.push_back(Wmuon_p4);
			g_Wmuon.push_back(d);
		      }
		  }//else if ( mcTools.IsLepton(d.pdgId() ) ) H+->Wh, W+->l v                                                                                                                              
		if (mcTools.IsNeutrino(dId)) 
		  {
		    neut_from_W_p4 = d.p4();
		    if (0) std::cout << "neut from W" << std::endl;
		  }
	      }//for (auto& d: genP_daughters)                                                                                                                                                             
          }
      }
	

    math::XYZTLorentzVector visibleTau_p4;
    if ( abs(genP_pdgId) == 15)
      {
	if (!bIsLastCopy) continue;
	if (p.pt() < cfg_TauPtCut ) continue;
	if (abs(p.eta()) > cfg_TauEtaCut) continue;
	int momTaupdgId = mcTools.RecursivelyLookForMotherId(p);
        if ( abs(momTaupdgId) == 25 || abs(momTaupdgId) == 35) 
	  { 
	    bool hadrtau = true;
	    if (p.pt() < cfg_TauPtCut ) continue;
	    if (abs(p.eta()) > cfg_TauEtaCut) continue;
	    cHWh_hDiTau_tau.increment();
	    for (auto& d: genP_daughters)                                                                                                                                                                  
	      {   
		math::XYZTLorentzVector tauvish_p4;
		math::XYVector tauvish_p2;
		if ( abs(d.pdgId()) == 211)                                                                                                                                                                
		  {                                                                                                                                                                                        
		    HWh_hDiTau_HadrTau_p4 = d.p4();
		    // visibleTau_p4 += HWh_hDiTau_HadrTau_p4; 
		    v_pion_htt_p4.push_back(HWh_hDiTau_HadrTau_p4);                                                                                                                                       
		  }                                                                                                                                                                                        
		if ( abs(d.pdgId()) == 13)                                                                                                                                                            
		  {
		    hadrtau = false;
		    bSkipEvent = true;
		    if ( d.pt()  < 5.0 ) continue; 
		    if ( abs(d.eta()) > cfg_MuonEtaCut) continue;
		    HWh_hDiTau_TauDiMu_p4 = d.p4();
		    double muonPt = d.pt();
		    double TauPtmuonPt = muonPt / p.pt();
		    if (0) std::cout << muonPt  << std::endl;
		    v_muonPtVsTauPt.push_back(TauPtmuonPt);
		    // visibleTau_p4 += HWh_hDiTau_TauDiMu_p4;
		    v_dimu_htt_p4.push_back(HWh_hDiTau_TauDiMu_p4);                                                                                                                                      
		  }
		
		if (mcTools.IsNeutrino(d.pdgId()))
		  {
		    if (0) std::cout <<"neut from tau" << std::endl;
		    neut_from_tau_p4 = d.p4();
		    v_neut_from_tau_p4.push_back(neut_from_tau_p4);
		    g_neut_from_tau.push_back(d);
		  }
		if (!mcTools.IsNeutrino(d.pdgId())) 
		  {
		    if (0) std::cout <<"tau d id= " << d.pdgId() <<  std::endl;
		    tauvish_p4 = d.p4();
		    tauvish_p2 = d.p2();
		    visibleTau_p2 += tauvish_p2; 
		    visibleTau_p4 += tauvish_p4;
		  }
	      }
	    if (hadrtau)
	      //if ( abs(d.pdgId()) != 13)
	      {
		if (0) std::cout << "hadrtau " << hadrtau << std::endl;
		HWh_hDiTau_p4 = p.p4();
		v_ditau_htt_p4.push_back(HWh_hDiTau_p4);
		g_ditau_htt.push_back(p);
		v_visibleTau_p4.push_back(visibleTau_p4);
		v_visibleTau_p2.push_back(visibleTau_p2);
	      }
	  }
      }
    
    //muon 
    if (0) std::cout << ".. Collection of muon" << std::endl;
    if ( abs(genP_pdgId) == 13)
      {
	if (!bIsLastCopy) continue;
	g_muon_p4.push_back(p);
	muon_p4 = p.p4();
	v_muon_p4.push_back(muon_p4);
      }


    //associated bottom quark
    if ( abs(genP_pdgId) == 5)
      { 
	bool asquark = false;
	int BQmompdgId = mcTools.RecursivelyLookForMotherId(p);
	if (0) std::cout << "BQmompdgId= " << BQmompdgId << std::endl;
	for (auto& s: genP_sisters)
	  {
	    if (0) std::cout << "=== sisters =" << s.pdgId() << std::endl; 
	    if(abs(s.pdgId()) == 37 )
	      {
		asquark = true;
	      }
	  }
	if (asquark)
	  {
	    /*
	    if (p.pt() < 20.0)
	      {
		bSkipEvent = true;
	      }
	    if (abs(p.eta()) > 2.4)
	      {
		bSkipEvent = true;
	      }
	    */
	    assoBQuark_p4 = p.p4();
	    v_assoBQuark_p4.push_back(assoBQuark_p4);
	  }
      }
     

    // Missing Et     
    if (mcTools.IsNeutrino(genP_pdgId))
      { 
	neutrino_p4 = p.p4();
	v_neutrino_p4.push_back(neutrino_p4);			
      }




  } //for (auto& p: fEvent.genparticles().getGenParticles())
  int allevent = 1 ;
  h_eventcount -> Fill(allevent);
  if ( bSkipEvent ) return;  
  cAllEventsAC.increment();
  //met calculation
  math::XYZTLorentzVector Met_p4;
  
  if (v_neutrino_p4.size() > 0)
    {
      for (size_t i=0; i < v_neutrino_p4.size(); i++)
        {
          Met_p4 += v_neutrino_p4.at(i);
        }
    }
  
  if (fEvent.genMET().et() >= -30 )
    {
      h_MET_Et -> Fill(fEvent.genMET().et());
      h_MET_Pt -> Fill(Met_p4.pt());
    }
  //com

  // loop over tau from H^{0}
  math::XYZTLorentzVector tau_p4;
  if (v_ditau_htt_p4.size() > 0 ) {
    //for-loop: All taus from H^{0}
    for (size_t i=0; i < v_ditau_htt_p4.size(); i++)
      {
	double tau_Pt  = v_ditau_htt_p4.at(i).pt();
	double tau_Eta = v_ditau_htt_p4.at(i).eta();
	double tau_Rap = mcTools.GetRapidity(v_ditau_htt_p4.at(i));
	h_tbH_HWh_hDiTau_Tau_Pt  -> Fill(tau_Pt);
	h_tbH_HWh_hDiTau_Tau_Eta -> Fill(tau_Eta);
	h_tbH_HWh_hDiTau_Tau_Rap -> Fill(tau_Rap);
	tau_p4 += v_ditau_htt_p4.at(i);
      } // for (size_t i=0; i < v_ditau_htt_p4.size(); i++)
  } //if (v_ditau_htt_p4.size() > 1 )
  
 
  if (v_muonPtVsTauPt.size() > 0)
    {
      for (size_t i=0; i < v_muonPtVsTauPt.size(); i++)
	{
	  h_muonPtTauVsPt -> Fill(v_muonPtVsTauPt.at(i));
	}
    }
  math::XYZTLorentzVector InEnergy;
  math::XYZTLorentzVector FiEnergy;
  FiEnergy = Htb_HPlus_p4 + bHt_tbW_BQuark_p4 + assoBQuark_p4;
  if (v_HiggsMom_p4.size() > 0)
    {
      for (size_t i=0; i < v_HiggsMom_p4.size(); i++)
	{
	  InEnergy += v_HiggsMom_p4.at(i);
	}
      h_InitialEOverFinalE -> Fill(FiEnergy.e()/InEnergy.e()); 
    }
  
  
  // DR H0-boson - W-boson
  // double bQuarkPair_dEta = std::abs(bQuarks_p4.at(deltaRMin_i).eta() - bQuarks_p4.at(deltaRMin_j).eta());
  //double bQuarkPair_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bQuarks_p4.at(deltaRMin_i), bQuarks_p4.at(deltaRMin_j)));
 
  //Correlations
  double tbH_WBoson_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(HWh_HBoson_p4, HWh_WBoson_p4);
  double tbH_WBoson_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(HWh_HBoson_p4, HWh_WBoson_p4)); 
  double tbH_WBoson_HWh_HBoson_dEta = std::abs(HWh_HBoson_p4.eta() - HWh_WBoson_p4.eta());
  double tbH_WBoson_HWh_HBoson_dY   = std::abs(mcTools.GetRapidity(HWh_HBoson_p4) - mcTools.GetRapidity(HWh_WBoson_p4));

  double tbH_HPlus_top_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_TQuark_p4, Htb_HPlus_p4);
  double tbH_HPlus_top_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_TQuark_p4, Htb_HPlus_p4));
  double tbH_HPlus_top_dEta = std::abs(bHt_TQuark_p4.eta() - Htb_HPlus_p4.eta());

  double bHt_TQuark_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_TQuark_p4, HWh_HBoson_p4);
  double bHt_TQuark_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_TQuark_p4, HWh_HBoson_p4));
  double bHt_TQuark_HWh_HBoson_dEta = std::abs(bHt_TQuark_p4.eta() - HWh_HBoson_p4.eta());
  if (v_assoBQuark_p4.size() > 0)
    {
      h_bHt_TQuark_HWh_HBoson_dPhi_Vs_Pt -> Fill(bHt_TQuark_HWh_HBoson_dPhi , assoBQuark_p4.pt());
      if (bHt_TQuark_HWh_HBoson_dPhi > 2.5) 
	{
	  h_tbH_BQuark_dPhiTH_Cut1_Pt    -> Fill(assoBQuark_p4.pt());
	}
      if (bHt_TQuark_HWh_HBoson_dPhi < 2)
	{
          h_tbH_BQuark_dPhiTH_Cut2_Pt    -> Fill(assoBQuark_p4.pt());
        }
    }
  

  double twb_WBoson_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_tbW_WBoson_p4, HWh_HBoson_p4);
  double twb_WBoson_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_tbW_WBoson_p4, HWh_HBoson_p4));
  double twb_WBoson_HWh_HBoson_dEta = std::abs(bHt_tbW_WBoson_p4.eta() - HWh_HBoson_p4.eta());

  double twb_BQuark_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_tbW_BQuark_p4, HWh_HBoson_p4);
  double twb_BQuark_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_tbW_BQuark_p4, HWh_HBoson_p4));
  double twb_BQuark_HWh_HBoson_dEta = std::abs(bHt_tbW_BQuark_p4.eta() - HWh_HBoson_p4.eta());
  
  double twb_TQuark_HWh_WBoson_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_TQuark_p4, HWh_WBoson_p4);
  double twb_TQuark_HWh_WBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_TQuark_p4, HWh_WBoson_p4));
  double twb_TQuark_HWh_WBoson_dEta = std::abs(bHt_TQuark_p4.eta() - HWh_WBoson_p4.eta());

  double twb_BQuark_HWh_WBoson_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_tbW_BQuark_p4, HWh_WBoson_p4);
  double twb_BQuark_HWh_WBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_tbW_BQuark_p4, HWh_WBoson_p4));
  double twb_BQuark_HWh_WBoson_dEta = std::abs(bHt_tbW_BQuark_p4.eta() - HWh_WBoson_p4.eta());
  
  double twb_WBoson_HWh_WBoson_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_tbW_WBoson_p4, HWh_WBoson_p4);
  double twb_WBoson_HWh_WBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_tbW_WBoson_p4, HWh_WBoson_p4));
  double twb_WBoson_HWh_WBoson_dEta = std::abs(bHt_tbW_WBoson_p4.eta() - HWh_WBoson_p4.eta());

  double twb_WBoson_twb_BQuark_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_tbW_WBoson_p4, bHt_tbW_BQuark_p4);
  double twb_WBoson_twb_BQuark_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_tbW_WBoson_p4, bHt_tbW_BQuark_p4));
  double twb_WBoson_twb_BQuark_dEta = std::abs(bHt_tbW_WBoson_p4.eta() - bHt_tbW_BQuark_p4.eta());

  // Fill Gen Histos 
  if (0) std::cout << "===qourkc=" << " " << qourkc << "HWh_Wqq_Quark_p4 size = " << v_HWh_Wqq_Quark_p4.size() << std::endl;
  // g->tt (Associated top)
  h_bHt_TQuark_Pt ->Fill( bHt_TQuark_p4.pt()  );
  h_bHt_TQuark_Eta->Fill( bHt_TQuark_p4.eta() );
  h_bHt_TQuark_Rap->Fill( mcTools.GetRapidity(bHt_TQuark_p4) );
  // tb->H+, H+                                                                    
  h_tbH_HPlus_Pt ->Fill( Htb_HPlus_p4.pt()  ); 
  h_tbH_HPlus_Eta->Fill( Htb_HPlus_p4.eta() );
  h_tbH_HPlus_Rap->Fill( mcTools.GetRapidity(Htb_HPlus_p4) );
  // g->tt, t->bW, W-boson (NOTE: t->Wb dominant, but t->Ws and t->Wd also possible!)                                                                                       
  h_bHt_tbW_WBoson_Pt ->Fill( bHt_tbW_WBoson_p4.pt()  );
  h_bHt_tbW_WBoson_Eta->Fill( bHt_tbW_WBoson_p4.eta() );
  h_bHt_tbW_WBoson_Rap->Fill( mcTools.GetRapidity(bHt_tbW_WBoson_p4) );
  // g->tt, t->bW, b (NOTE: t->Wb dominant, but t->Ws and t->Wd also possible!)                                                                                             
  h_bHt_tbW_BQuark_Pt ->Fill( bHt_tbW_BQuark_p4.pt()  );
  h_bHt_tbW_BQuark_Eta->Fill( bHt_tbW_BQuark_p4.eta() );
  h_bHt_tbW_BQuark_Rap->Fill( mcTools.GetRapidity(bHt_tbW_BQuark_p4) );
  //H+->HW, H-boson
  h_tbH_HWh_HBoson_Pt  -> Fill(HWh_HBoson_p4.pt());
  h_tbH_HWh_HBoson_Eta -> Fill(HWh_HBoson_p4.eta());
  h_tbH_HWh_HBoson_Rap -> Fill( mcTools.GetRapidity(HWh_HBoson_p4) );
  // H+-> HW, W-boson
  h_tbH_HWh_WBoson_Pt  -> Fill(HWh_WBoson_p4.pt());
  h_tbH_HWh_WBoson_Eta -> Fill(HWh_WBoson_p4.eta());
  h_tbH_HWh_WBoson_Rap -> Fill( mcTools.GetRapidity(HWh_WBoson_p4) );
  
  // H0-> DiTau -> DiMu
  if (v_dimu_htt_p4.size() > 0 ) 
    {
      for (size_t i=0; i < v_dimu_htt_p4.size(); i++)
	{
	  h_hDiTau_TauDiMu_DiMu_Pt  -> Fill(v_dimu_htt_p4.at(i).pt());
	  h_hDiTau_TauDiMu_DiMu_Eta -> Fill(v_dimu_htt_p4.at(i).eta());
	  h_hDiTau_TauDiMu_DiMu_Rap -> Fill( mcTools.GetRapidity(v_dimu_htt_p4.at(i)) );
	  
	}
    }
  // H0-> DiTau -> Pion
  if (v_pion_htt_p4.size() > 0 )
    {  
      for (size_t i=0; i < v_pion_htt_p4.size(); i++)
	{
	  h_hDiTau_TauDiMu_Pion_Pt  -> Fill(v_pion_htt_p4.at(i).pt());
	  h_hDiTau_TauDiMu_Pion_Eta -> Fill(v_pion_htt_p4.at(i).eta());
	  h_hDiTau_TauDiMu_Pion_Rap -> Fill( mcTools.GetRapidity(v_pion_htt_p4.at(i)) );
	} 
    }
  //t->bW, W -> qq
  if (v_Wqq_Quark_p4.size() > 0)
    {
      for (size_t i=0; i < v_Wqq_Quark_p4.size(); i++)
	{
	  h_Wqq_Quark_Pt        -> Fill(v_Wqq_Quark_p4.at(i).pt());
	  h_Wqq_Quark_Eta       -> Fill(v_Wqq_Quark_p4.at(i).eta());
	  h_Wqq_Quark_Rap       -> Fill(mcTools.GetRapidity(v_Wqq_Quark_p4.at(i)));
	  double Wqq_Quark_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(HWh_HBoson_p4, v_Wqq_Quark_p4.at(i));
          double Wqq_Quark_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(HWh_HBoson_p4, v_Wqq_Quark_p4.at(i)));
          double Wqq_Quark_HWh_HBoson_dEta = std::abs(HWh_HBoson_p4.eta() - v_Wqq_Quark_p4.at(i).eta());
          h_Wqq_Quark_HWh_HBoson_dR   -> Fill(Wqq_Quark_HWh_HBoson_dR);
          h_Wqq_Quark_HWh_HBoson_dPhi -> Fill(Wqq_Quark_HWh_HBoson_dPhi);
          h_Wqq_Quark_HWh_HBoson_dEta -> Fill(Wqq_Quark_HWh_HBoson_dEta);
          h_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi -> Fill(Wqq_Quark_HWh_HBoson_dEta , Wqq_Quark_HWh_HBoson_dPhi);
          h_Wqq_Quark_HWh_HBoson_dR_Vs_Pt     -> Fill(Wqq_Quark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
	}
    }
  if (v_Wqq_AntiQuark_p4.size() > 0)
    {
      for (size_t i=0; i < v_Wqq_AntiQuark_p4.size(); i++)
	{
	  h_Wqq_AntiQuark_Pt    -> Fill(v_Wqq_AntiQuark_p4.at(i).pt());
	  h_Wqq_AntiQuark_Eta   -> Fill(v_Wqq_AntiQuark_p4.at(i).eta());
	  h_Wqq_AntiQuark_Rap   -> Fill(mcTools.GetRapidity(v_Wqq_AntiQuark_p4.at(i)));
	  double Wqq_AntiQuark_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(HWh_HBoson_p4, v_Wqq_AntiQuark_p4.at(i));
          double Wqq_AntiQuark_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(HWh_HBoson_p4, v_Wqq_AntiQuark_p4.at(i)));
          double Wqq_AntiQuark_HWh_HBoson_dEta = std::abs(HWh_HBoson_p4.eta() - v_Wqq_AntiQuark_p4.at(i).eta());
          h_Wqq_AntiQuark_HWh_HBoson_dR   -> Fill(Wqq_AntiQuark_HWh_HBoson_dR);
          h_Wqq_AntiQuark_HWh_HBoson_dPhi -> Fill(Wqq_AntiQuark_HWh_HBoson_dPhi);
          h_Wqq_AntiQuark_HWh_HBoson_dEta -> Fill(Wqq_AntiQuark_HWh_HBoson_dEta);
          h_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi -> Fill(Wqq_AntiQuark_HWh_HBoson_dEta , Wqq_AntiQuark_HWh_HBoson_dPhi);
          h_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt     -> Fill(Wqq_AntiQuark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
	}
    }
  
  if (v_twb_Wqq_Quark_p4.size() > 0)
    {
      for (size_t i=0; i < v_twb_Wqq_Quark_p4.size(); i++)
        {
	  double twb_Wqq_Quark_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(HWh_HBoson_p4, v_twb_Wqq_Quark_p4.at(i));
          double twb_Wqq_Quark_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(HWh_HBoson_p4, v_twb_Wqq_Quark_p4.at(i)));
          double twb_Wqq_Quark_HWh_HBoson_dEta = std::abs(HWh_HBoson_p4.eta() - v_twb_Wqq_Quark_p4.at(i).eta());
          h_twb_Wqq_Quark_HWh_HBoson_dR   -> Fill(twb_Wqq_Quark_HWh_HBoson_dR);
          h_twb_Wqq_Quark_HWh_HBoson_dPhi -> Fill(twb_Wqq_Quark_HWh_HBoson_dPhi);
          h_twb_Wqq_Quark_HWh_HBoson_dEta -> Fill(twb_Wqq_Quark_HWh_HBoson_dEta);
	  h_twb_Wqq_Quark_HWh_HBoson_dEta_Vs_dPhi -> Fill(twb_Wqq_Quark_HWh_HBoson_dEta , twb_Wqq_Quark_HWh_HBoson_dPhi);
          h_twb_Wqq_Quark_HWh_HBoson_dR_Vs_Pt     -> Fill(twb_Wqq_Quark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
        }
    }
  
  if (v_twb_Wqq_AntiQuark_p4.size() > 0)
    {
      for (size_t i=0; i < v_twb_Wqq_AntiQuark_p4.size(); i++)
        {
          double twb_Wqq_AntiQuark_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(HWh_HBoson_p4, v_twb_Wqq_AntiQuark_p4.at(i));
          double twb_Wqq_AntiQuark_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(HWh_HBoson_p4, v_twb_Wqq_AntiQuark_p4.at(i)));
          double twb_Wqq_AntiQuark_HWh_HBoson_dEta = std::abs(HWh_HBoson_p4.eta() - v_twb_Wqq_AntiQuark_p4.at(i).eta());
          h_twb_Wqq_AntiQuark_HWh_HBoson_dR   -> Fill(twb_Wqq_AntiQuark_HWh_HBoson_dR);
          h_twb_Wqq_AntiQuark_HWh_HBoson_dPhi -> Fill(twb_Wqq_AntiQuark_HWh_HBoson_dPhi);
          h_twb_Wqq_AntiQuark_HWh_HBoson_dEta -> Fill(twb_Wqq_AntiQuark_HWh_HBoson_dEta);
	  h_twb_Wqq_AntiQuark_HWh_HBoson_dEta_Vs_dPhi -> Fill(twb_Wqq_AntiQuark_HWh_HBoson_dEta , twb_Wqq_AntiQuark_HWh_HBoson_dPhi);
          h_twb_Wqq_AntiQuark_HWh_HBoson_dR_Vs_Pt     -> Fill(twb_Wqq_AntiQuark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
        }
    }
  
  math::XYZTLorentzVector Quark_p4;
  if (v_twb_Wqq_Quark_p4.size() > 0)
    {
      for (size_t i=0; i < v_twb_Wqq_Quark_p4.size(); i++)
	{
	  if (v_twb_Wqq_AntiQuark_p4.size() == 0) continue;
	  for (size_t j=0; j < v_twb_Wqq_AntiQuark_p4.size(); j++)
	    {
	      double twb_Wqq_Quark_AntiQuark_dR   = ROOT::Math::VectorUtil::DeltaR(v_twb_Wqq_Quark_p4.at(i), v_twb_Wqq_AntiQuark_p4.at(j));
	      double twb_Wqq_Quark_AntiQuark_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(v_twb_Wqq_Quark_p4.at(i), v_twb_Wqq_AntiQuark_p4.at(j)));
	      double twb_Wqq_Quark_AntiQuark_dEta = std::abs(v_twb_Wqq_Quark_p4.at(i).eta() - v_twb_Wqq_AntiQuark_p4.at(j).eta());
	      h_twb_Wqq_Quark_AntiQuark_dR       -> Fill(twb_Wqq_Quark_AntiQuark_dR);
	      h_twb_Wqq_Quark_AntiQuark_dPhi     -> Fill(twb_Wqq_Quark_AntiQuark_dPhi);
	      h_twb_Wqq_Quark_AntiQuark_dEta     -> Fill(twb_Wqq_Quark_AntiQuark_dEta);
	      h_twb_Wqq_Quark_AntiQuark_dEta_Vs_dPhi -> Fill(twb_Wqq_Quark_AntiQuark_dEta , twb_Wqq_Quark_AntiQuark_dPhi);
	      h_twb_Wqq_Quark_AntiQuark_dR_Vs_Pt     -> Fill(twb_Wqq_Quark_AntiQuark_dR , bHt_tbW_BQuark_p4.pt());
	      Quark_p4 += v_twb_Wqq_Quark_p4.at(i) + v_twb_Wqq_AntiQuark_p4.at(j);
	    }
	}
      double twb_Wqq_Quark_MET_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(Quark_p4, Met_p4));
      h_twb_Wqq_Quark_MET_dPhi -> Fill(twb_Wqq_Quark_MET_dPhi);
    }
  
  if (v_Q1_p4.size() > 0)
    {
      for (size_t i=0; i < v_Q1_p4.size(); i++)
	{
	  h_Wqq_Quarks_Pt  -> Fill(v_Q1_p4.at(i).pt());
	  h_Wqq_Quarks_Eta -> Fill(v_Q1_p4.at(i).eta());
	  h_Wqq_Quarks_Rap -> Fill(mcTools.GetRapidity(v_Q1_p4.at(i)));
	}
    }
  if (v_l1_p4.size() > 0)
    {
      for (size_t i=0; i < v_l1_p4.size(); i++)
        {
          h_Wln_Lepton_Pt  -> Fill(v_l1_p4.at(i).pt());
          h_Wln_Lepton_Eta -> Fill(v_l1_p4.at(i).eta());
          h_Wln_Lepton_Rap -> Fill(mcTools.GetRapidity(v_l1_p4.at(i)));
        }
    }
  

  double dphineu = 2323.0;
  math::XYZTLorentzVector DiNeur_from_tau_p4;
  if (v_neut_from_tau_p4.size() > 0)
    {
      for (size_t i=0; i < v_neut_from_tau_p4.size(); i++)
        {
          DiNeur_from_tau_p4 += v_neut_from_tau_p4.at(i);
        }
      if (v_HWh_Wmuon_Muon_p4.size() > 0 )
        {
          double HWh_WBoson_neut_HWh_HBoson_neut_Dphi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(DiNeur_from_tau_p4 , neut_from_W_p4));
          double HWh_HBoson_neut_Met_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(DiNeur_from_tau_p4 , Met_p4));
	  dphineu = HWh_WBoson_neut_HWh_HBoson_neut_Dphi;
	  h_HWh_WBoson_neut_HWh_HBoson_neut_Dphi -> Fill(HWh_WBoson_neut_HWh_HBoson_neut_Dphi);
          h_HWh_HBoson_neut_Met_dPhi -> Fill(HWh_HBoson_neut_Met_dPhi);
          h_HWh_HBoson_neut_Met_dPhi_Vs_HWh_WBoson_neut_HWh_HBoson_neut_dPhi -> Fill(HWh_WBoson_neut_HWh_HBoson_neut_Dphi , HWh_HBoson_neut_Met_dPhi);
        }
    }

  math::XYZTLorentzVector divistau_p4;
  if (v_visibleTau_p4.size() > 0 )
    {
      if (v_HWh_Wmuon_Muon_p4.size() > 0 )
	{
	  for (size_t i=0; i < v_visibleTau_p4.size(); i++)
	    {
	      divistau_p4 += v_visibleTau_p4.at(i);
	      double HWh_Wmuon_Muon_HWh_h_Tau_dR = ROOT::Math::VectorUtil::DeltaR(Wmuon_p4 , v_visibleTau_p4.at(i));
	      h_HWh_Wmuon_Muon_HWh_h_Tau_dR -> Fill(HWh_Wmuon_Muon_HWh_h_Tau_dR);
	    }
	  double HWh_Wmuon_Muon_HWh_h_Tau_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(divistau_p4, Wmuon_p4 ));
	  h_HWh_Wmuon_Muon_HWh_h_Tau_dPhi -> Fill(HWh_Wmuon_Muon_HWh_h_Tau_dPhi);
	}
    }
  double taudPhi = 0.;
  if (v_visibleTau_p4.size() > 1)
    {  
      for (size_t i=0; i < v_visibleTau_p4.size()-1; i++)
        {
	  for (size_t j= i+1; j < v_visibleTau_p4.size(); j++)
	    {
	      double HWh_hDiTau_Tau_Tau_dR = ROOT::Math::VectorUtil::DeltaR(v_visibleTau_p4.at(i), v_visibleTau_p4.at(j));
	      double HWh_hDiTau_Tau_Tau_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(v_visibleTau_p4.at(i), v_visibleTau_p4.at(j)));
	      taudPhi = HWh_hDiTau_Tau_Tau_dPhi;
	      h_HWh_hDiTau_Tau_Tau_dR -> Fill(HWh_hDiTau_Tau_Tau_dR);
	      h_HWh_hDiTau_Tau_Tau_dPhi -> Fill(HWh_hDiTau_Tau_Tau_dPhi);
	      if (v_visibleTau_p4.at(i).pt() > v_visibleTau_p4.at(j).pt())
		{
		  h_HWh_hDiTau_Tau_Tau_dR_Vs_Pt     -> Fill(HWh_hDiTau_Tau_Tau_dR , v_visibleTau_p4.at(i).pt());
		}
	      else
		{
		  h_HWh_hDiTau_Tau_Tau_dR_Vs_Pt     -> Fill(HWh_hDiTau_Tau_Tau_dR , v_visibleTau_p4.at(j).pt());
		}
	    }
	}
    }
  
  if (v_visibleTau_p4.size() > 0 ) 
    {
      for (size_t i=0; i < v_visibleTau_p4.size(); i++)
	{
	  double HWh_hDiTau_Tau_MET_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(v_visibleTau_p4.at(i), Met_p4));
	  double HWh_DiTau_TauNeut_dPhi  = std::abs(ROOT::Math::VectorUtil::DeltaPhi(v_visibleTau_p4.at(i), DiNeur_from_tau_p4));
	  h_HWh_hDiTau_Tau_MET_dPhi     -> Fill(HWh_hDiTau_Tau_MET_dPhi);
	  h_HWh_DiTau_TauNeut_dPhi      -> Fill(HWh_DiTau_TauNeut_dPhi);
	}
    }
  
  
  //Transverse Mass of Charge Higgs
  
  // Definitions
  //if (g_ditau_htt.size() > 1)
  if (v_visibleTau_p2.size() > 1 )
    {
      const math::XYVector met   = (const math::XYVector) fEvent.genMET().p2();
      if (g_Wmuon.size() > 0)
	{
	  const math::XYVector muon  = g_Wmuon.at(0).p2();
	  //Transverse Mass
	  double mT_genP = GetMt(v_visibleTau_p2.at(0), v_visibleTau_p2.at(1), muon, met);
	  h_mT_ChargeHiggs_M -> Fill(mT_genP);
	  if (taudPhi < 2.5 ) {
	    h_mT_ChargeHiggs_taudPhils_M -> Fill(mT_genP);
	  }
	  if (taudPhi > 2. ) {
            h_mT_ChargeHiggs_taudPhigr_M -> Fill(mT_genP);
          }
	  if (assWHadr)
	    {
	      h_mT_ChargeHiggs_assW_hadr_M -> Fill(mT_genP);
	    }
	  else
	    {
	      h_mT_ChargeHiggs_assW_lept_M -> Fill(mT_genP);
	    }
	  if (g_neut_from_tau.size() > 1)
	    {
	      const math::XYVector metfromtau = g_neut_from_tau.at(0).p2() + g_neut_from_tau.at(1).p2();
	      double mT_CHoTauMet = GetMt(v_visibleTau_p2.at(0), v_visibleTau_p2.at(1), muon, metfromtau);
	      h_mT_ChargeHiggsMeTTau_M -> Fill(mT_CHoTauMet);
	      if (dphineu < 1.047)
		{
		  h_mT_ChargeHiggsMeTTau_Cut1_M -> Fill(mT_CHoTauMet);
		}
	      else
		{
		  h_mT_ChargeHiggsMeTTau_Cut2_M -> Fill(mT_CHoTauMet);
		}
	    }
	  if (g_assWQuark.size() > 0)
	    {
	      const math::XYVector assq = g_assWQuark.at(0).p2();
	      const math::XYVector assaq = g_assWAntiQuark.at(0).p2();
	      double mT_assqq = GetMt(assq, assaq, muon, met);
	      h_mT_ChargeHiggs_assqq_M -> Fill(mT_assqq);
	    }
	} //if (g_Wmuon.size() > 0)
      if (g_muon_p4.size() == 1)
	{
	  const math::XYVector muonany = g_muon_p4.at(0).p2();
	  double mT_genPanymuon = GetMt(v_visibleTau_p2.at(0), v_visibleTau_p2.at(1), muonany, met);
	  h_mT_ChargeHiggs_OneAnyMuon_M -> Fill(mT_genPanymuon);
	}
      if (g_muon_p4.size() > 1)
        {
	  double hpt = -20.;
	  genParticle g_hpt_muon;
	  for (size_t i=0; i < g_muon_p4.size(); i++)
	    {
	      double mpt = g_muon_p4.at(i).pt();
	      if (mpt > hpt)
		{
		  g_hpt_muon = g_muon_p4.at(i);
		}
	    }
	  const math::XYVector muonhpt = g_hpt_muon.p2();
	  double mT_genPanymuon = GetMt(v_visibleTau_p2.at(0), v_visibleTau_p2.at(1), muonhpt, met);
	  h_mT_ChargeHiggs_hptMuon_M -> Fill(mT_genPanymuon);
	}
      if (g_assWmuon.size() > 0)
	{
	  const math::XYVector assmuon = g_assWmuon.at(0).p2();
	  double mT_assWgenP = GetMt(v_visibleTau_p2.at(0), v_visibleTau_p2.at(1), assmuon, met);
	  h_mT_ChargeHiggs_assMuon_M -> Fill(mT_assWgenP);
	}
    }
  
  
  // H^+ -> WH, W-> muon,neutrino
  
  if (v_HWh_Wmuon_Muon_p4.size() > 0)
    {
      h_HWh_Wmn_Muon_Pt    -> Fill(Wmuon_p4.pt());
      h_HWh_Wmn_Muon_Eta   -> Fill(Wmuon_p4.eta());
      h_HWh_Wmn_Muon_Rap   -> Fill(mcTools.GetRapidity(Wmuon_p4));
      double HWh_Wmn_Muon_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(HWh_HBoson_p4, Wmuon_p4);
      double HWh_Wmn_Muon_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(HWh_HBoson_p4, Wmuon_p4));
      double HWh_Wmn_Muon_HWh_HBoson_dEta = std::abs(HWh_HBoson_p4.eta() - Wmuon_p4.eta());
      h_HWh_Wmn_Muon_HWh_HBoson_dR        -> Fill(HWh_Wmn_Muon_HWh_HBoson_dR);
      h_HWh_Wmn_Muon_HWh_HBoson_dPhi      -> Fill(HWh_Wmn_Muon_HWh_HBoson_dPhi);
      h_HWh_Wmn_Muon_HWh_HBoson_dEta      -> Fill(HWh_Wmn_Muon_HWh_HBoson_dEta);
      h_HWh_Wmn_Muon_HWh_HBoson_dEta_Vs_dPhi -> Fill(HWh_Wmn_Muon_HWh_HBoson_dEta , HWh_Wmn_Muon_HWh_HBoson_dPhi);
      HWh_Wmn_Muon_HWh_HBoson_dR_Vs_Pt       -> Fill(HWh_Wmn_Muon_HWh_HBoson_dR , HWh_HBoson_p4.pt());
      
      double HWh_Wmn_Muon_MET_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(Wmuon_p4 , Met_p4));
      h_HWh_Wmn_Muon_MET_dPhi -> Fill(HWh_Wmn_Muon_MET_dPhi);
    }

  //associated bottom quark
  
  math::XYZTLorentzVector tbH_BQuark;
  if (v_assoBQuark_p4.size() > 0) 
    {
      h_tbH_BQuark_Pt    -> Fill(assoBQuark_p4.pt());
      h_tbH_BQuark_Eta   -> Fill(assoBQuark_p4.eta());
      h_tbH_BQuark_Rap   -> Fill(mcTools.GetRapidity(assoBQuark_p4));
      double tbH_BQuark_HWh_HPlus_dR   = ROOT::Math::VectorUtil::DeltaR(Htb_HPlus_p4, assoBQuark_p4);
      double tbH_BQuark_HWh_HPlus_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(Htb_HPlus_p4, assoBQuark_p4));
      double tbH_BQuark_HWh_HPlus_dEta = std::abs(Htb_HPlus_p4.eta() - assoBQuark_p4.eta());
      double tbH_BQuark_twb_TQuark_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_TQuark_p4, assoBQuark_p4);
      double tbH_BQuark_twb_TQuark_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_TQuark_p4, assoBQuark_p4));
      double tbH_BQuark_twb_TQuark_dEta = std::abs(bHt_TQuark_p4.eta() - assoBQuark_p4.eta());
      h_tbH_BQuark_HWh_HPlus_dR        -> Fill(tbH_BQuark_HWh_HPlus_dR);
      h_tbH_BQuark_HWh_HPlus_dPhi      -> Fill(tbH_BQuark_HWh_HPlus_dPhi);
      h_tbH_BQuark_HWh_HPlus_dEta      -> Fill(tbH_BQuark_HWh_HPlus_dEta);
      h_tbH_BQuark_HWh_HPlus_dR_Vs_Pt  -> Fill(tbH_BQuark_HWh_HPlus_dR , assoBQuark_p4.pt());
      h_tbH_BQuark_twb_TQuark_dR       -> Fill(tbH_BQuark_twb_TQuark_dR);
      h_tbH_BQuark_twb_TQuark_dPhi     -> Fill(tbH_BQuark_twb_TQuark_dPhi);
      h_tbH_BQuark_twb_TQuark_dEta     -> Fill(tbH_BQuark_twb_TQuark_dEta);
      h_tbH_BQuark_twb_TQuark_dR_Vs_Pt -> Fill(tbH_BQuark_twb_TQuark_dR , assoBQuark_p4.pt());
      tbH_BQuark = assoBQuark_p4;
      
      double tbH_BQuark_MET_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(tbH_BQuark, Met_p4));
      h_tbH_BQuark_MET_dPhi -> Fill(tbH_BQuark_MET_dPhi);
    }
  double tbH_BQuark_MET_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_tbW_BQuark_p4,Met_p4));
  h_twb_BQuark_MET_dPhi -> Fill(tbH_BQuark_MET_dPhi);


  math::XYZTLorentzVector leadBQuark_p4;
  //b-b properties
  if (v_assoBQuark_p4.size() > 0) 
    {
      double leadBQuark_pt = bHt_tbW_BQuark_p4.pt();
      if (leadBQuark_pt > assoBQuark_p4.pt())
	{
	  leadBQuark_p4 = bHt_tbW_BQuark_p4;
	  h_leadBQuark_pt -> Fill(leadBQuark_pt);
	}
      else
	{
	  leadBQuark_p4 = assoBQuark_p4;
	  h_leadBQuark_pt -> Fill(assoBQuark_p4.pt());
	}
      
      double bHt_tbW_BQuark_tbH_BQuark_dR = ROOT::Math::VectorUtil::DeltaR(bHt_tbW_BQuark_p4 , assoBQuark_p4);
      h_bHt_tbW_BQuark_tbH_BQuark_dR -> Fill(bHt_tbW_BQuark_tbH_BQuark_dR);
    }

  math::XYZTLorentzVector Wqq_leadQuark_p4;
  math::XYZTLorentzVector Wqq_sublead_p4;
  math::XYZTLorentzVector twb_topproducts_p4;
  int leadquarknumber = 1;
  int subleadquarknumber = 1;
  if (v_twb_Wqq_Quark_p4.size() > 0)
    {
      twb_topproducts_p4 = bHt_tbW_BQuark_p4 + Wqq_Quark_p4 + Wqq_AnriQuark_p4;
      Wqq_sublead_p4  = bHt_tbW_BQuark_p4;
      Wqq_leadQuark_p4 = bHt_tbW_BQuark_p4;
      if (Wqq_Quark_p4.pt() > Wqq_leadQuark_p4.pt()) 
	{
	  Wqq_leadQuark_p4 = Wqq_Quark_p4;
	  leadquarknumber = 2;
	}
      else 
	{
	  Wqq_sublead_p4 = Wqq_Quark_p4;
	  subleadquarknumber = 2;
	}
      if (Wqq_AnriQuark_p4.pt() > Wqq_leadQuark_p4.pt())
	{
	Wqq_leadQuark_p4 = Wqq_AnriQuark_p4;
	leadquarknumber = 3;
	}
      else if (Wqq_AnriQuark_p4.pt() > Wqq_sublead_p4.pt())
	{
	  Wqq_sublead_p4 = Wqq_AnriQuark_p4;
	  subleadquarknumber = 3;
	}
      if (assoBQuark_p4.pt() > Wqq_leadQuark_p4.pt())
	{
	  Wqq_leadQuark_p4 = assoBQuark_p4;
	  leadquarknumber = 4;
	}
      else if (assoBQuark_p4.pt() > Wqq_sublead_p4.pt())
	{
	  Wqq_sublead_p4 = assoBQuark_p4;
	  subleadquarknumber = 4;
	}
      int eventpass = 2;
      h_eventcount -> Fill(eventpass);
      h_leadquarknumber -> Fill(leadquarknumber);
      h_subleadquarknumber -> Fill(subleadquarknumber);
      if (leadquarknumber == 1)
	{
	  clquark1.increment();
	}
      else if (leadquarknumber == 2)
	{
	  clquark2.increment();
	}
      else if (leadquarknumber == 3)
	{
          clquark3.increment();
        }
      else
	{
	  clquark4.increment();
	}
      if (subleadquarknumber == 1)
	{
          csquark1.increment();
	}
      else if (subleadquarknumber == 2)
	{
          csquark2.increment();
        }
      else if (subleadquarknumber == 3)
        {
          csquark3.increment();
        }
      else
	{
          csquark4.increment();
	}


      for (size_t i=0; i < v_twb_topproducts.size(); i++)
	{
	  double Wqq_leadQuark_tbqq_Quark_dR = ROOT::Math::VectorUtil::DeltaR(Wqq_leadQuark_p4 , v_twb_topproducts.at(i));
	  double Wqq_leadQuark_tbqq_Quark_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(Wqq_leadQuark_p4 , v_twb_topproducts.at(i)));
	  double Wqq_leadQuark_tbqq_Quark_dEta = std::abs(Wqq_leadQuark_p4.eta() - v_twb_topproducts.at(i).eta());
	  if (Wqq_leadQuark_tbqq_Quark_dR > 0.0 )
	    {
	      h_Wqq_leadQuark_tbqq_Quark_dR -> Fill(Wqq_leadQuark_tbqq_Quark_dR);
	      h_Wqq_leadQuark_tbqq_Quark_dR_Vs_Pt -> Fill( Wqq_leadQuark_tbqq_Quark_dR , Wqq_leadQuark_p4.pt());
	      h_Wqq_leadQuark_tbqq_Quark_dPhi -> Fill(Wqq_leadQuark_tbqq_Quark_dPhi);
	      h_Wqq_leadQuark_tbqq_Quark_dEta -> Fill(Wqq_leadQuark_tbqq_Quark_dEta);
	    }
	}
      if (v_assoBQuark_p4.size() > 0)
	{
	  double tbH_BQuark_tbqq_leadQuark_dR = ROOT::Math::VectorUtil::DeltaR(Wqq_leadQuark_p4 , assoBQuark_p4);
	  double tbH_BQuark_tbqq_subleadQuark_dR = ROOT::Math::VectorUtil::DeltaR(Wqq_sublead_p4 , assoBQuark_p4);
	  h_tbH_BQuark_tbqq_leadQuark_dR -> Fill(tbH_BQuark_tbqq_leadQuark_dR);
	  h_tbH_BQuark_tbqq_leadQuark_dR_Vs_Pt -> Fill(tbH_BQuark_tbqq_leadQuark_dR , Wqq_leadQuark_p4.pt());
	  h_tbH_BQuark_tbqq_subleadQuark_dR -> Fill(tbH_BQuark_tbqq_subleadQuark_dR);
	  h_tbH_BQuark_tbqq_subleadQuark_dR_Vs_Pt -> Fill(tbH_BQuark_tbqq_subleadQuark_dR , Wqq_sublead_p4.pt());
	}
    }
  

  // Fill dR Histos
  h_tbH_HWh_WBoson_HBoson_dR   -> Fill(tbH_WBoson_HWh_HBoson_dR);
  h_tbH_HPlus_top_dR           -> Fill(tbH_HPlus_top_dR);
  h_bHt_TQuark_HWh_HBoson_dR   -> Fill(bHt_TQuark_HWh_HBoson_dR);
  h_twb_WBoson_HWh_HBoson_dR   -> Fill(twb_WBoson_HWh_HBoson_dR);
  h_twb_BQuark_HWh_HBoson_dR   -> Fill(twb_BQuark_HWh_HBoson_dR);
  h_twb_TQuark_HWh_WBoson_dR   -> Fill(twb_TQuark_HWh_WBoson_dR);
  h_twb_BQuark_HWh_WBoson_dR   -> Fill(twb_BQuark_HWh_WBoson_dR);
  h_twb_WBoson_HWh_WBoson_dR   -> Fill(twb_WBoson_HWh_WBoson_dR);
  h_twb_WBoson_twb_BQuark_dR   -> Fill(twb_WBoson_twb_BQuark_dR);
  // Fill dPhi Histos
  h_tbH_HWh_WBoson_HBoson_dPhi -> Fill(tbH_WBoson_HWh_HBoson_dPhi);
  h_tbH_HPlus_top_dPhi         -> Fill(tbH_HPlus_top_dPhi);
  h_bHt_TQuark_HWh_HBoson_dPhi -> Fill(bHt_TQuark_HWh_HBoson_dPhi);
  h_twb_WBoson_HWh_HBoson_dPhi -> Fill(twb_WBoson_HWh_HBoson_dPhi);
  h_twb_BQuark_HWh_HBoson_dPhi -> Fill(twb_BQuark_HWh_HBoson_dPhi);
  h_twb_TQuark_HWh_WBoson_dPhi -> Fill(twb_TQuark_HWh_WBoson_dPhi);
  h_twb_BQuark_HWh_WBoson_dPhi -> Fill(twb_BQuark_HWh_WBoson_dPhi);
  h_twb_WBoson_HWh_WBoson_dPhi -> Fill(twb_WBoson_HWh_WBoson_dPhi);
  h_twb_WBoson_twb_BQuark_dPhi -> Fill(twb_WBoson_twb_BQuark_dPhi);
  // Fill dEta Histos
  h_tbH_HWh_WBoson_HBoson_dEta -> Fill(tbH_WBoson_HWh_HBoson_dEta);
  h_tbH_HPlus_top_dEta         -> Fill(tbH_HPlus_top_dEta);
  h_bHt_TQuark_HWh_HBoson_dEta -> Fill(bHt_TQuark_HWh_HBoson_dEta);
  h_twb_WBoson_HWh_HBoson_dEta -> Fill(twb_WBoson_HWh_HBoson_dEta);
  h_twb_BQuark_HWh_HBoson_dEta -> Fill(twb_BQuark_HWh_HBoson_dEta);
  h_twb_TQuark_HWh_WBoson_dEta -> Fill(twb_TQuark_HWh_WBoson_dEta);
  h_twb_BQuark_HWh_WBoson_dEta -> Fill(twb_BQuark_HWh_WBoson_dEta);
  h_twb_WBoson_HWh_WBoson_dEta -> Fill(twb_WBoson_HWh_WBoson_dEta);
  h_twb_WBoson_twb_BQuark_dEta -> Fill(twb_WBoson_twb_BQuark_dEta);
  // Fill dY Histos
  h_tbH_WBoson_HWh_HBoson_dRap   -> Fill(tbH_WBoson_HWh_HBoson_dY);
  // dEta_Vs_dPhi
  h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi -> Fill(tbH_WBoson_HWh_HBoson_dEta , tbH_WBoson_HWh_HBoson_dPhi);
  h_tbH_HPlus_top_dEta_Vs_dPhi         -> Fill(tbH_HPlus_top_dEta , tbH_HPlus_top_dPhi);
  h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi  -> Fill(bHt_TQuark_HWh_HBoson_dEta , bHt_TQuark_HWh_HBoson_dPhi);
  h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi -> Fill(twb_WBoson_HWh_HBoson_dEta , twb_WBoson_HWh_HBoson_dPhi);
  h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi -> Fill(twb_BQuark_HWh_HBoson_dEta , twb_BQuark_HWh_HBoson_dPhi);
  h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi -> Fill(twb_TQuark_HWh_WBoson_dEta , twb_TQuark_HWh_WBoson_dPhi);
  h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi -> Fill(twb_BQuark_HWh_WBoson_dEta , twb_BQuark_HWh_WBoson_dPhi);
  h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi -> Fill(twb_WBoson_HWh_WBoson_dEta , twb_WBoson_HWh_WBoson_dPhi);
  h_twb_WBoson_twb_BQuark_dEta_Vs_dPhi -> Fill(twb_WBoson_twb_BQuark_dEta , twb_WBoson_twb_BQuark_dPhi);
  //dR_Vs_Pt
  h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt -> Fill(tbH_WBoson_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_tbH_HPlus_top_dR_Vs_Pt         -> Fill(tbH_HPlus_top_dR , Htb_HPlus_p4.pt());
  h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt -> Fill(bHt_TQuark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_twb_WBoson_HWh_HBoson_dR_Vs_Pt -> Fill(twb_WBoson_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_twb_BQuark_HWh_HBoson_dR_Vs_Pt -> Fill(twb_BQuark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_twb_TQuark_HWh_WBoson_dR_Vs_Pt -> Fill(twb_TQuark_HWh_WBoson_dR , HWh_WBoson_p4.pt());
  h_twb_BQuark_HWh_WBoson_dR_Vs_Pt -> Fill(twb_BQuark_HWh_WBoson_dR , HWh_WBoson_p4.pt());
  h_twb_WBoson_HWh_WBoson_dR_Vs_Pt -> Fill(twb_WBoson_HWh_WBoson_dR , HWh_WBoson_p4.pt());
  h_twb_WBoson_twb_BQuark_dR_Vs_Pt -> Fill(twb_WBoson_twb_BQuark_dR , bHt_tbW_BQuark_p4.pt());
  
  


  return;
}

double GenHToHW::GetMt(const math::XYVector tau1, const math::XYVector tau2, const math::XYVector muon, const math::XYVector& met) {
  // Use scalar sums to get the transverse mass
  double metEt  = met.R();
  double tau1Et = tau1.r();
  double tau2Et = tau2.r();
  double muonEt = muon.r();
  double mT     =-999.9;
  double mTSq   =   0.0;
  double EtSq   = (metEt + tau1Et + tau2Et + muonEt) * (metEt + tau1Et + tau2Et + muonEt);
  double EtXSq  = (tau1.x() + tau2.x() + muon.x() + met.x()) * (tau1.x() + tau2.x() + muon.x() + met.x());
  double EtYSq  = (tau1.y() + tau2.y() + muon.y() + met.y()) * (tau1.y() + tau2.y() + muon.y() + met.y());
  mTSq = EtSq - (EtXSq + EtYSq);

  if (mTSq >= 0) mT = std::sqrt(mTSq);
  return mT;
}

vector<GenJet> GenHToHW::GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCuts, std::vector<float> etaCuts, vector<genParticle> genParticlesToMatch)
{
  /*
    Jet-Flavour Definitions (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools)

    Algorithmic definition: (NOTE: Algorithmic definition is used by default for all b-tagging purposes)
    - Try to find the parton that most likely determines the properties of the jet and assign that flavour as the true flavour
    - Here, the ¡°final state¡± partons (after showering, radiation) are analyzed (within ¦¤R < 0.3 of the reconstructed jet axis). 
    Partons selected for the algorithmic definition are those partons that don't have other partons as daughters, 
    without any explicit requirement on their status (for Pythia6 these are status=2 partons).
    - Jets from radiation are matched with full efficiency
    -If there is a b/c within the jet cone: label as b/c
    -Otherwise: assign flavour of the hardest parton
   */

  // Definitions
  std::vector<GenJet> jets;
  unsigned int jet_index   = -1;
  unsigned int ptCut_index  = 0;
  unsigned int etaCut_index = 0;

  // For-loop: All genParticles
  for (auto& p: genParticlesToMatch) 
    {
      
      // Comparison variables
      double dR   = 1e6;
      double dPt  = 1e6;
      double dEta = 1e6;
      double dPhi = 1e6;
	  
      // For-loop: All Generated Jets
      for (auto jet: genJets) 
	{
	  
	  // Jet index (for pT and eta cuts)
	  jet_index++;
	  
	  dPt  = jet.pt() - p.pt();
	  dEta = jet.eta() - p.eta();
	  dPhi = jet.phi() - p.phi();
       	  dR   = ROOT::Math::VectorUtil::DeltaR(jet.p4(), p.p4());
      
	  // Fail to match
	  if (dR > 0.3) continue;
	  
	  // Apply cuts after the matching
	  const float ptCut  = ptCuts.at(ptCut_index);
	  const float etaCut = etaCuts.at(etaCut_index);
	  if (jet.pt() < ptCut) continue;
	  if (std::abs(jet.eta()) > etaCut) continue;
	  
	  // Save this particle
	  jets.push_back(jet);
	  if (0) std::cout << "dR = " << dR << ": dPt = " << dPt << ", dEta = " << dEta << ", dPhi = " << dPhi << std::endl;

	  // Increment cut index only. Cannot be bigger than the size of the cut list provided
	  if (ptCut_index  < ptCuts.size()-1  ) ptCut_index++;
	  if (etaCut_index < etaCuts.size()-1  ) etaCut_index++;
	  break;
	}
    }
  if (0) std::cout << "bjets.size() = " << jets.size() << std::endl;
  return jets;
}

vector<GenJet> GenHToHW::GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts)
{
  std::vector<GenJet> jets;

  // Definitions
  unsigned int jet_index   = -1;
  unsigned int ptCut_index  = 0;
  unsigned int etaCut_index = 0;

  // For-loop: All Generated Jets
  for (auto jet: genJets) 
    {

      // Jet index (for pT and eta cuts)
      jet_index++;
      
      // Apply cuts
      const float ptCut  = ptCuts.at(ptCut_index);
      const float etaCut = etaCuts.at(etaCut_index);
      if (jet.pt() < ptCut) continue;
      if (std::abs(jet.eta()) > etaCut) continue;

      // Save this particle
      jets.push_back(jet);

      // Increment cut index only. Cannot be bigger than the size of the cut list provided
      if (ptCut_index  < ptCuts.size()-1  ) ptCut_index++;
      if (etaCut_index < etaCuts.size()-1  ) etaCut_index++;
    }
  return jets;
}

vector<genParticle> GenHToHW::GetAllCopies(const vector<genParticle> genParticles, genParticle genP)
{

  std::vector<genParticle> copies;
  short genP_maxIndex = genP.index();

  // For-loop: All genParticles
  for (size_t i = genP_maxIndex; i>=0; i--)
    {
      
      genParticle p = genParticles.at(i);

      if (p.mothers().size() < 1) break;

      genParticle m = genParticles.at(i);
      if ( m.pdgId() == genP.pdgId() )
	{
	  copies.push_back(p);
	  i = m.index();
	}
    }

  // Sort vector with ascending index order
  std::reverse(copies.begin(),copies.end());
  // for (auto& p: copies) 
  //   {
  //     std::cout << "p.pdgId() = " << p.pdgId() << ", p.index() = " << p.index()  << std::endl;
  //   }
  // std::cout << "" << std::endl;

  return copies;
}

vector<genParticle> GenHToHW::GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy, const bool hasNoDaughters)
  {
  
  std::vector<genParticle> particles;

  // For-loop: All genParticles
  for (auto& p: genParticles) 
    {

      // std::cout << "p.pdgId() = " << p.pdgId() << std::endl;

      // Find last copy of a given particle
      // if (isLastCopy) if (!p.isLastCopy()) continue; // crashes
      if (isLastCopy)
	{
	  // fixme: iro
	  if (abs(p.pdgId()) == 5){ if (p.status() != 23) continue;} // Consider only status=23 (outgoing) particles
	  else{ if (p.status() != 1 and p.status() != 2) continue;}
	}
      
      // Commonly enables for parton-based jet flavour definition
      if (hasNoDaughters) if (p.daughters().size() > 0) continue;

      // Consider only particles
      if (std::abs(p.pdgId()) != pdgId) continue;

      // Apply cuts
      if ( p.pt() < ptCut ) continue;
      if (std::abs(p.eta()) > etaCut) continue;
      
      // Save this particle
      particles.push_back(p);
    }

  // std::cout << "Returning " << particles.size() << " particles." << std::endl;
  return particles;
  }
