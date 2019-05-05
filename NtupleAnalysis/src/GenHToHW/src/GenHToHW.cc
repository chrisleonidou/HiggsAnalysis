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
  Count cbHt_tbW_WBoson;
  Count cbHt_tbW_Wqq_Quark;
  Count cbHt_tbW_Wqq_AntiQuark;
  Count cbHt_tbW_Wqq_Leptons;

  
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
  //
  WrappedTH1 *h_tbH_HWh_WBoson_HBoson_dR;
  WrappedTH1 *h_tbH_HPlus_top_dR;
  WrappedTH1 *h_bHt_TQuark_HWh_HBoson_dR;
  WrappedTH1 *h_twb_WBoson_HWh_HBoson_dR;
  WrappedTH1 *h_twb_BQuark_HWh_HBoson_dR;
  WrappedTH1 *h_twb_TQuark_HWh_WBoson_dR;
  WrappedTH1 *h_twb_BQuark_HWh_WBoson_dR;
  WrappedTH1 *h_twb_WBoson_HWh_WBoson_dR;
  //
  WrappedTH1 *h_tbH_HWh_WBoson_HBoson_dPhi;
  WrappedTH1 *h_tbH_HPlus_top_dPhi;
  WrappedTH1 *h_bHt_TQuark_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_WBoson_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_BQuark_HWh_HBoson_dPhi;
  WrappedTH1 *h_twb_TQuark_HWh_WBoson_dPhi;
  WrappedTH1 *h_twb_BQuark_HWh_WBoson_dPhi;
  WrappedTH1 *h_twb_WBoson_HWh_WBoson_dPhi;
  //
  WrappedTH1 *h_tbH_HWh_WBoson_HBoson_dEta;
  WrappedTH1 *h_tbH_HPlus_top_dEta;
  WrappedTH1 *h_bHt_TQuark_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_WBoson_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_BQuark_HWh_HBoson_dEta;
  WrappedTH1 *h_twb_TQuark_HWh_WBoson_dEta;
  WrappedTH1 *h_twb_BQuark_HWh_WBoson_dEta;
  WrappedTH1 *h_twb_WBoson_HWh_WBoson_dEta;
  //
  WrappedTH2 *h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_tbH_HPlus_top_dEta_Vs_dPhi;
  WrappedTH2 *h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi;
  WrappedTH2 *h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi;
  //
  WrappedTH2 *h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_tbH_HPlus_top_dR_Vs_Pt;
  WrappedTH2 *h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_WBoson_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_BQuark_HWh_HBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_TQuark_HWh_WBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_BQuark_HWh_WBoson_dR_Vs_Pt;
  WrappedTH2 *h_twb_WBoson_HWh_WBoson_dR_Vs_Pt;
  
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
    cbHt_tbW_BQuark(fEventCounter.addSubCounter("Branching", "H+->tb, t->bW, b")),
    cbHt_tbW_WBoson(fEventCounter.addSubCounter("Branching", "H+->tb, t->bW, W+")),
    cbHt_tbW_Wqq_Quark(fEventCounter.addSubCounter("Branching", "H+->tb, t->bW, W->qq, q")),
    cbHt_tbW_Wqq_AntiQuark(fEventCounter.addSubCounter("Branching", "H+->tb, t->bW, W->qq, qbar")),
    cbHt_tbW_Wqq_Leptons(fEventCounter.addSubCounter("Branching", "H+->tb, t->bW, W->l v"))
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

  // const int nBinsM  = cfg_MassBinSetting.bins();
  // const double minM = cfg_MassBinSetting.min();
  // const double maxM = cfg_MassBinSetting.max();
  
  const int nBinsdEta  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindEta = cfg_DeltaEtaBinSetting.min();
  const double maxdEta = 2*cfg_DeltaEtaBinSetting.max();

  // const int nBinsdRap  = 2*cfg_DeltaEtaBinSetting.bins();
  // const double mindRap = cfg_DeltaEtaBinSetting.min();
  // const double maxdRap = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdPhi  = cfg_DeltaPhiBinSetting.bins();
  const double mindPhi = cfg_DeltaPhiBinSetting.min();
  const double maxdPhi = 1.2*cfg_DeltaPhiBinSetting.max(); //desp

  const int nBinsdR  = 2*cfg_DeltaRBinSetting.bins();
  const double mindR = cfg_DeltaRBinSetting.min();
  const double maxdR = cfg_DeltaRBinSetting.max();

  TDirectory* th1 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH1");
  TDirectory* th2 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH2");
    
  // GenParticles  
  hParticle_Pt                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "AllParticle_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_TQuark_Pt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tbW_WBoson_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_WBoson_Pt"       , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tbW_BQuark_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_BQuark_Pt"       , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HPlus_Pt                = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HWh_HBoson_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_HBoson_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HWh_WBoson_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_tbH_HWh_hDiTau_Tau_Pt       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_hDiTau_Tau_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_hDiTau_TauDiMu_DiMu_Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_DiMu_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_hDiTau_TauDiMu_Pion_Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_Pion_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt); 
  //
  h_bHt_TQuark_Eta              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Eta"           , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tbW_WBoson_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_WBoson_Eta"       , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tbW_BQuark_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_BQuark_Eta"       , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HPlus_Eta               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HWh_HBoson_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_HBoson_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HWh_WBoson_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  h_tbH_HWh_hDiTau_Tau_Eta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_hDiTau_Tau_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  h_hDiTau_TauDiMu_DiMu_Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_DiMu_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  h_hDiTau_TauDiMu_Pion_Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_Pion_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  //
  h_bHt_TQuark_Rap              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Rap"           , ";#omega", nBinsRap, minRap, maxRap);
  h_bHt_tbW_WBoson_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_WBoson_Rap"       , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_tbW_BQuark_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tbW_BQuark_Rap"       , ";#omega", nBinsEta, minRap, maxRap);
  h_tbH_HPlus_Rap               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_Rap"            , ";#omega", nBinsRap, minRap, maxRap);
  h_tbH_HWh_HBoson_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_HBoson_Rap"            , ";#omega", nBinsRap, minRap, maxRap);  
  h_tbH_HWh_WBoson_Rap          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_Rap"            , ";#omega", nBinsRap, minRap, maxRap);	
  h_tbH_HWh_hDiTau_Tau_Rap      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_hDiTau_Tau_Rap"            , ";#omega", nBinsRap, minRap, maxRap);
  h_hDiTau_TauDiMu_DiMu_Rap     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_DiMu_Rap"            , ";#omega", nBinsRap, minRap, maxRap);
  h_hDiTau_TauDiMu_Pion_Rap     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_hDiTau_TauDiMu_Pion_Rap"            , ";#omega", nBinsRap, minRap, maxRap);
  //
  h_tbH_HWh_WBoson_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_tbH_HPlus_top_dR            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_top_dR"              , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_TQuark_HWh_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_TQuark_HWh_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_WBoson_HWh_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_BQuark_HWh_HBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_HBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_TQuark_HWh_WBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_TQuark_HWh_WBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_BQuark_HWh_WBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_WBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_twb_WBoson_HWh_WBoson_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_WBoson_dR"      , ";#DeltaR", nBinsdR, mindR, maxdR);
  //
  h_tbH_HWh_WBoson_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_tbH_HPlus_top_dPhi          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_top_dPhi"         , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_TQuark_HWh_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_TQuark_HWh_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_HBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_HBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_TQuark_HWh_WBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_TQuark_HWh_WBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_WBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_WBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_WBoson_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_WBoson_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);
  //
  h_tbH_HWh_WBoson_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HWh_WBoson_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_tbH_HPlus_top_dEta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_tbH_HPlus_top_dEta"         , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_bHt_TQuark_HWh_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_TQuark_HWh_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_WBoson_HWh_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_BQuark_HWh_HBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_HBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_TQuark_HWh_WBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_TQuark_HWh_WBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_BQuark_HWh_WBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_BQuark_HWh_WBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_twb_WBoson_HWh_WBoson_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_twb_WBoson_HWh_WBoson_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  //
  h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_tbH_HPlus_top_dEta_Vs_dPhi         = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_HPlus_top_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  //
  h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_tbH_HPlus_top_dR_Vs_Pt             = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_tbH_HPlus_top_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_WBoson_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_BQuark_HWh_HBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_HBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_TQuark_HWh_WBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_TQuark_HWh_WBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_BQuark_HWh_WBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_BQuark_HWh_WBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);
  h_twb_WBoson_HWh_WBoson_dR_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "h_twb_WBoson_HWh_WBoson_dR_Vs_Pt", ";#DeltaR;p_{T} (GeV/c)", nBinsdR, mindR, maxdR, nBinsPt, minPt, maxPt);

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
  if (cfg_Verbose) std::cout << "=== Electron veto" << std::endl;
  vector<genParticle> selectedElectrons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_ElectronPtCut, cfg_ElectronEtaCut, 11, true, false);
  if (0)
    {
      std::cout << "\nnElectrons = " << selectedElectrons.size() << std::endl;
      for (auto& p: selectedElectrons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }
  if ( selectedElectrons.size() > 0 ) return;
  cElectronVeto.increment();


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

  //================================================================================================
  // 12) Top selection 
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Top Selection" << std::endl;
  cTopSelection.increment();


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
  math::XYZTLorentzVector HWh_hDiTau_HadrTau_p4;
  math::XYZTLorentzVector HWh_hDiTau_TauDiMu_p4;


  // For-loop: GenParticles
  for (auto& p: fEvent.genparticles().getGenParticles()) {
        
    hParticle_Pt -> Fill(p.pt());
    //int genMom_index = -1; 
    //double genMom_pdgId = -999.99; 
    // Particle properties
    int genP_pdgId       = p.pdgId();
    /*
    short genP_index     = p.index();
    int genP_status      = p.status();
    double genP_pt       = p.pt();
    double genP_eta      = p.eta();
    double genP_phi      = p.phi();
    double genP_energy   = p.e(); 
    int genP_charge      = p.charge();
    */
    // Associated genParticles
    std::vector<genParticle> genP_daughters;
    for (unsigned int i=0; i < p.daughters().size(); i++) genP_daughters.push_back(fEvent.genparticles().getGenParticles()[p.daughters().at(i)]);
    std::vector<genParticle> genP_mothers;
    for (unsigned int i=0; i < p.mothers().size(); i++) genP_mothers.push_back(fEvent.genparticles().getGenParticles()[p.mothers().at(i)]);
    std::vector<genParticle> genP_grandMothers;
    std::vector<genParticle> genP_allCopies = GetAllCopies(fEvent.genparticles().getGenParticles(), p);
    std::vector<unsigned int> genGMoms_index;
    std::vector<unsigned int> genGMoms_pdgId;
    std::vector<unsigned int> genMoms_index;
    std::vector<unsigned int> genMoms_pdgId;
    std::vector<unsigned int> genDaus_index;
    std::vector<unsigned int> genDaus_pdgId; 

    // Daughter, Mom and Grand-mom properties                                                                                                                              
 
    genParticle m;
    genParticle g;
    genParticle d;
    genParticle firstMom;
    genParticle lastMom;

    if (genP_daughters.size() > 0) d = genP_daughters.at(0);
    if (p.mothers().size() > 0)
      {
        m = genP_mothers.at(0);                                                                                                                                    
	for (unsigned int i=0; i < m.mothers().size(); i++) genP_grandMothers.push_back(fEvent.genparticles().getGenParticles()[m.mothers().at(i)]);
        if (m.mothers().size() > 0) g = genP_grandMothers.at(0);                                                                                                   
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
    
    if (0) std::cout << ".. Collection of Charge Higgs .." << std::endl;
    if (abs(genP_pdgId) == 37) // chris1
      {
	if (!bIsLastCopy) continue;
	
	cbHt_HPlus.increment();
	Htb_HPlus_p4 = p.p4();
      
	// For-loop: All daughters 
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
    
    if ( abs(genP_pdgId) == 15)
      {
	if (!bIsLastCopy) continue;
	if (p.pt() < cfg_TauPtCut ) continue;                                                                                                                                                  
	if (abs(p.eta()) > cfg_TauEtaCut) continue;
	if (abs(firstMom.pdgId()) == 25 || abs(firstMom.pdgId()) == 35) 
	  {
	    HWh_hDiTau_p4 = p.p4();
	    v_ditau_htt_p4.push_back (HWh_hDiTau_p4);
	    for (auto& d: genP_daughters)
	      {
		if ( abs(d.pdgId()) == 211)
		  {
		    HWh_hDiTau_HadrTau_p4 = d.p4();
		    v_pion_htt_p4.push_back (HWh_hDiTau_HadrTau_p4);
		  }
		else if ( abs(d.pdgId()) == 13)
		  {
		    HWh_hDiTau_TauDiMu_p4 = d.p4();
		    v_dimu_htt_p4.push_back (HWh_hDiTau_TauDiMu_p4);
		    
		  }
	      }
	  }
      }
    
    

    if (0) std::cout << "=== Associated top" << std::endl;
    if ( abs(genP_pdgId) == 6)
      {
        if (!bIsLastCopy) continue;
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
              }
            else if ( abs(d.pdgId()) == 24)
              {
		if (0) std::cout << "TopIndx3=" << p.index() << "---->" << "TopId3=" << genP_pdgId << "---->" << d.pdgId() << std::endl;
                // NOTE: t->Wb dominant, but t->Ws and t->Wd also possible!                                                                                                 
                bHt_tbW_WBoson_p4 = d.p4();
              }
            else bSkipEvent = true;
          } // for (auto& d: genP_daughters)                                                                                                                                
      }// if (abs(genP_pdgId) == 6) 




  } //for (auto& p: fEvent.genparticles().getGenParticles())
  
  
  if ( bSkipEvent ) return; 
  //com

  // loop over tau from H^{0}
  
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
	
      } // for (size_t i=0; i < v_ditau_htt_p4.size(); i++)
    
  } //if (v_ditau_htt_p4.size() > 1 )
  

  // DR H0-boson - W-boson
  // double bQuarkPair_dEta = std::abs(bQuarks_p4.at(deltaRMin_i).eta() - bQuarks_p4.at(deltaRMin_j).eta());
  //double bQuarkPair_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bQuarks_p4.at(deltaRMin_i), bQuarks_p4.at(deltaRMin_j)));
 
  //Correlations
  double tbH_WBoson_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(HWh_HBoson_p4, HWh_WBoson_p4);
  double tbH_WBoson_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(HWh_HBoson_p4, HWh_WBoson_p4)); 
  double tbH_WBoson_HWh_HBoson_dEta = std::abs(HWh_HBoson_p4.eta() - HWh_WBoson_p4.eta());
  
  if (tbH_WBoson_HWh_HBoson_dPhi > 3.15) 
    {
      if (0) std::cout << "dphi>pi" << tbH_WBoson_HWh_HBoson_dPhi <<  std::endl;
    }
  double tbH_HPlus_top_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_TQuark_p4, Htb_HPlus_p4);
  double tbH_HPlus_top_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_TQuark_p4, Htb_HPlus_p4));
  double tbH_HPlus_top_dEta = std::abs(bHt_TQuark_p4.eta() - Htb_HPlus_p4.eta());

  double bHt_TQuark_HWh_HBoson_dR   = ROOT::Math::VectorUtil::DeltaR(bHt_TQuark_p4, HWh_HBoson_p4);
  double bHt_TQuark_HWh_HBoson_dPhi = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bHt_TQuark_p4, HWh_HBoson_p4));
  double bHt_TQuark_HWh_HBoson_dEta = std::abs(bHt_TQuark_p4.eta() - HWh_HBoson_p4.eta());

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

  // Fill Gen Histos chris
  
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
      h_hDiTau_TauDiMu_DiMu_Pt  -> Fill(HWh_hDiTau_TauDiMu_p4.pt());
      h_hDiTau_TauDiMu_DiMu_Eta -> Fill(HWh_hDiTau_TauDiMu_p4.eta());
      h_hDiTau_TauDiMu_DiMu_Rap -> Fill( mcTools.GetRapidity(HWh_hDiTau_TauDiMu_p4) );
    }  
  // H0-> DiTau -> Pion
  if (v_pion_htt_p4.size() > 0 )
    {  
      h_hDiTau_TauDiMu_Pion_Pt  -> Fill(HWh_hDiTau_HadrTau_p4.pt());
      h_hDiTau_TauDiMu_Pion_Eta -> Fill(HWh_hDiTau_HadrTau_p4.eta());
      h_hDiTau_TauDiMu_Pion_Rap -> Fill( mcTools.GetRapidity(HWh_hDiTau_HadrTau_p4) );
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
  // Fill dPhi Histos
  h_tbH_HWh_WBoson_HBoson_dPhi -> Fill(tbH_WBoson_HWh_HBoson_dPhi);
  h_tbH_HPlus_top_dPhi         -> Fill(tbH_HPlus_top_dPhi);
  h_bHt_TQuark_HWh_HBoson_dPhi -> Fill(bHt_TQuark_HWh_HBoson_dPhi);
  h_twb_WBoson_HWh_HBoson_dPhi -> Fill(twb_WBoson_HWh_HBoson_dPhi);
  h_twb_BQuark_HWh_HBoson_dPhi -> Fill(twb_BQuark_HWh_HBoson_dPhi);
  h_twb_TQuark_HWh_WBoson_dPhi -> Fill(twb_TQuark_HWh_WBoson_dPhi);
  h_twb_BQuark_HWh_WBoson_dPhi -> Fill(twb_BQuark_HWh_WBoson_dPhi);
  h_twb_WBoson_HWh_WBoson_dPhi -> Fill(twb_WBoson_HWh_WBoson_dPhi);
  // Fill dEta Histos
  h_tbH_HWh_WBoson_HBoson_dEta -> Fill(tbH_WBoson_HWh_HBoson_dEta);
  h_tbH_HPlus_top_dEta         -> Fill(tbH_HPlus_top_dEta);
  h_bHt_TQuark_HWh_HBoson_dEta -> Fill(bHt_TQuark_HWh_HBoson_dEta);
  h_twb_WBoson_HWh_HBoson_dEta -> Fill(twb_WBoson_HWh_HBoson_dEta);
  h_twb_BQuark_HWh_HBoson_dEta -> Fill(twb_BQuark_HWh_HBoson_dEta);
  h_twb_TQuark_HWh_WBoson_dEta -> Fill(twb_TQuark_HWh_WBoson_dEta);
  h_twb_BQuark_HWh_WBoson_dEta -> Fill(twb_BQuark_HWh_WBoson_dEta);
  h_twb_WBoson_HWh_WBoson_dEta -> Fill(twb_WBoson_HWh_WBoson_dEta);
  // dEta_Vs_dPhi
  h_tbH_WBoson_HWh_HBoson_dEta_Vs_dPhi -> Fill(tbH_WBoson_HWh_HBoson_dEta , tbH_WBoson_HWh_HBoson_dPhi);
  h_tbH_HPlus_top_dEta_Vs_dPhi         -> Fill(tbH_HPlus_top_dEta , tbH_HPlus_top_dPhi);
  h_bHt_TQuark_HWh_HBoso_dEta_Vs_dPhi  -> Fill(bHt_TQuark_HWh_HBoson_dEta , bHt_TQuark_HWh_HBoson_dPhi);
  h_twb_WBoson_HWh_HBoson_dEta_Vs_dPhi -> Fill(twb_WBoson_HWh_HBoson_dEta , twb_WBoson_HWh_HBoson_dPhi);
  h_twb_BQuark_HWh_HBoson_dEta_Vs_dPhi -> Fill(twb_BQuark_HWh_HBoson_dEta , twb_BQuark_HWh_HBoson_dPhi);
  h_twb_TQuark_HWh_WBoson_dEta_Vs_dPhi -> Fill(twb_TQuark_HWh_WBoson_dEta , twb_TQuark_HWh_WBoson_dPhi);
  h_twb_BQuark_HWh_WBoson_dEta_Vs_dPhi -> Fill(twb_BQuark_HWh_WBoson_dEta , twb_BQuark_HWh_WBoson_dPhi);
  h_twb_WBoson_HWh_WBoson_dEta_Vs_dPhi -> Fill(twb_WBoson_HWh_WBoson_dEta , twb_WBoson_HWh_WBoson_dPhi);
  //dR_Vs_Pt
  h_tbH_WBoson_HWh_HBoson_dR_Vs_Pt -> Fill(tbH_WBoson_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_tbH_HPlus_top_dR_Vs_Pt         -> Fill(tbH_HPlus_top_dR , Htb_HPlus_p4.pt());
  h_bHt_TQuark_HWh_HBoson_dR_Vs_Pt -> Fill(bHt_TQuark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_twb_WBoson_HWh_HBoson_dR_Vs_Pt -> Fill(twb_WBoson_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_twb_BQuark_HWh_HBoson_dR_Vs_Pt -> Fill(twb_BQuark_HWh_HBoson_dR , HWh_HBoson_p4.pt());
  h_twb_TQuark_HWh_WBoson_dR_Vs_Pt -> Fill(twb_TQuark_HWh_WBoson_dR , HWh_WBoson_p4.pt());
  h_twb_BQuark_HWh_WBoson_dR_Vs_Pt -> Fill(twb_BQuark_HWh_WBoson_dR , HWh_WBoson_p4.pt());
  h_twb_WBoson_HWh_WBoson_dR_Vs_Pt -> Fill(twb_WBoson_HWh_WBoson_dR , HWh_WBoson_p4.pt());
  
  return;
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
