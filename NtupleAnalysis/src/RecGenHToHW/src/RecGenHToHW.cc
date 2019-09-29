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


class RecGenHToHW: public BaseSelector {
public:
  explicit RecGenHToHW(const ParameterSet& config, const TH1* skimCounters);
  virtual ~RecGenHToHW() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;
  virtual vector<genParticle> GetAllCopies(const vector<genParticle> genParticles, genParticle genP);
  virtual vector<GenJet> GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCut, std::vector<float> etaCut);
  Bool_t IsBjet(math::XYZTLorentzVector jet,vector<GenJet> selectedBJets);
  virtual vector<GenJet> GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCut, std::vector<float> etaCut, vector<genParticle> genParticlesToMatch);
  virtual vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy=true, const bool hasNoDaughters=false);
  virtual double GetMt(const math::XYVector tau1, const math::XYVector tau2, const math::XYVector muon, const math::XYVector& met);
  virtual double GetWMt(const math::XYVector muon, const math::XYVector& met);
  virtual double GetWMt2(const math::XYZTLorentzVector muon, const math::XYZTLorentzVector neut);
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
  Count cMETSelection;
  // Count cTopologySelection;
  Count cSelected;
  Count cFatTau;
  // BR Counters
  Count cInclusive;
  Count cbHt_HPlus;
  Count cbHt_HBoson;
  Count cbHt_WBoson;


  // GenParticles                                                                                                                                                                                           
  WrappedTH1 *hParticle_Pt;
  WrappedTH1 *h_htau1_htau2_dR;
  WrappedTH1 *h_HPus_dRCut_Pt;
  WrappedTH1 *h_bHt_tbW_BQuark_Pt;
  WrappedTH1 *h_bHt_tbW_BQuark_Eta;
  WrappedTH1 *h_bHt_tWb_Wqq_Pt;
  WrappedTH1 *h_bHt_tWb_Wqq_Eta;
  WrappedTH1 *h_bHt_tWb_Wqq_qq_dR;
  WrappedTH1 *h_bHt_tWb_Wqq_qq_dEta;

  //For Efficiency plots
  WrappedTH1Triplet *h_TopQuarkPt_isGenuineTop;
  //WrappedTH1 *h_TopQuarkPt_isGenuineTop;
  
  //Genuine Top
  WrappedTH1 *h_GenuineTop_ldgJet_Pt;
  WrappedTH1 *h_GenuineTop_subldgJet_Pt;
  WrappedTH1 *h_GenuineTop_Bjet_Pt;
  WrappedTH1 *h_GenuineTop_Dijet_dR;
  WrappedTH1 *h_GenuineTop_ldgJet_Bjet_dR;
  WrappedTH1 *h_GenuineTop_subJet_Bjet_dR;
  WrappedTH1 *h_GenuineTop_Dijet_Bjet_dR;
  WrappedTH1 *h_GenuineTop_Dijet_Pt;
  WrappedTH1 *h_GenuineTop_Dijet_M;
  WrappedTH1 *h_GenuineTop_Trijet_Pt;
  WrappedTH1 *h_GenuineTop_Trijet_M;


  // GenJets                                                                                                                                
  WrappedTH1 *h_GenJets_N;
  WrappedTH1 *h_GenJet1_Pt;
  WrappedTH1 *h_GenJet2_Pt;
  WrappedTH1 *h_GenJet3_Pt;
  WrappedTH1 *h_GenJet4_Pt;
  WrappedTH1 *h_GenJet5_Pt;
  WrappedTH1 *h_GenJet6_Pt;
  //                                                                                                                                        
  WrappedTH1 *h_GenJet1_Eta;
  WrappedTH1 *h_GenJet2_Eta;
  WrappedTH1 *h_GenJet3_Eta;
  WrappedTH1 *h_GenJet4_Eta;
  WrappedTH1 *h_GenJet5_Eta;
  WrappedTH1 *h_GenJet6_Eta;
  // GenJets: Trijet with largest pT                                                                                                        
  WrappedTH1 *h_MaxTriJetPt_Pt;
  WrappedTH1 *h_MaxTriJetPt_Eta;
  WrappedTH1 *h_MaxTriJetPt_Rap;
  WrappedTH1 *h_MaxTriJetPt_Mass;
  WrappedTH1 *h_MaxTriJetPt_dEtaMax;
  WrappedTH1 *h_MaxTriJetPt_dPhiMax;
  WrappedTH1 *h_MaxTriJetPt_dRMax;
  WrappedTH1 *h_MaxTriJetPt_dEtaMin;
  WrappedTH1 *h_MaxTriJetPt_dPhiMin;
  WrappedTH1 *h_MaxTriJetPt_dRMin;
  WrappedTH1 *h_MaxTriJetPt_dEtaAverage;
  WrappedTH1 *h_MaxTriJetPt_dPhiAverage;
  WrappedTH1 *h_MaxTriJetPt_dRAverage;
  WrappedTH1 *h_TriJetMaxPt_bHt_TQuark_dR;
  WrappedTH1 *h_TriJetMaxPt_Over_bHt_TQuarkPt;
  WrappedTH1 *h_TriJetMaxPt_Top_M;
  WrappedTH1 *h_TriJetMaxPt_Top_OverCut_M;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(RecGenHToHW);

RecGenHToHW::RecGenHToHW(const ParameterSet& config, const TH1* skimCounters)
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
    cMETSelection(fEventCounter.addCounter("MET")),
    // cTopologySelection(fEventCounter.addCounter("Topology")),
    cSelected(fEventCounter.addCounter("All Selections")),
    cFatTau(fEventCounter.addCounter("FatTau")),
    cInclusive(fEventCounter.addSubCounter("Branching", "All events")),
    cbHt_HPlus(fEventCounter.addSubCounter("Branching", "H+")),
    cbHt_HBoson(fEventCounter.addSubCounter("Branching", "H+->HW, H")),
    cbHt_WBoson(fEventCounter.addSubCounter("Branching", "H+->HW, W"))
{ }

void RecGenHToHW::book(TDirectory *dir) {
  
  // Fixed binning
  const int nBinsPt   = 4*cfg_PtBinSetting.bins();
  const double minPt  = cfg_PtBinSetting.min();
  const double maxPt  = 4*cfg_PtBinSetting.max();
  
  const int nBinsdR  = 2*cfg_DeltaRBinSetting.bins();                                                                                                                                                        const double mindR = cfg_DeltaRBinSetting.min();                                                                                                                                                           const double maxdR = cfg_DeltaRBinSetting.max();  
  
 
  const int nBinsEta  = 2*cfg_EtaBinSetting.bins();
  const double minEta = cfg_EtaBinSetting.min();
  const double maxEta = 2*cfg_EtaBinSetting.max();
  
  const int nBinsRap  = cfg_EtaBinSetting.bins();
  const double minRap = cfg_EtaBinSetting.min();
  const double maxRap = cfg_EtaBinSetting.max();

  // const int nBinsPhi  = cfg_PhiBinSetting.bins();
  // const double minPhi = cfg_PhiBinSetting.min();
  // const double maxPhi = cfg_PhiBinSetting.max();
  const double minM = cfg_MassBinSetting.min();
  
  const int nBinsdEta  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindEta = cfg_DeltaEtaBinSetting.min();
  const double maxdEta = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdRap  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindRap = cfg_DeltaEtaBinSetting.min();
  const double maxdRap = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdPhi  = cfg_DeltaPhiBinSetting.bins();
  const double mindPhi = cfg_DeltaPhiBinSetting.min();
  const double maxdPhi = 1.2*cfg_DeltaPhiBinSetting.max(); 

  
  
  TDirectory* th1 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH1");
  TDirectory* th2 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH2");
  std::string myInclusiveLabel  = "AnalysisTriplets";
  std::string myFakeLabel       = myInclusiveLabel+"False";
  std::string myGenuineLabel    = myInclusiveLabel+"True";
  TDirectory* myInclusiveDir         = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myInclusiveLabel);
  TDirectory* myFakeDir = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myFakeLabel);
  TDirectory* myGenuineDir = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myGenuineLabel);
  std::vector<TDirectory*> myDirs = {myInclusiveDir, myFakeDir, myGenuineDir};
  //  TDirectory* myDirs = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, "myGenuineDir");

  // GenParticles  
  hParticle_Pt                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "AllParticle_Pt"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_htau1_htau2_dR              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_htau1_htau2_dR"         , ";#DeltaR"      , nBinsdR, mindR, maxdR);
  h_HPus_dRCut_Pt               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_HPus_dRCut_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tbW_BQuark_Pt           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_tbW_BQuark_Pt"      , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tbW_BQuark_Eta          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_tbW_BQuark_Eta"     , ";#eta"      , nBinsEta, minEta, maxEta);
  h_bHt_tWb_Wqq_Pt              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_tWb_Wqq_Pt"         , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tWb_Wqq_Eta             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_tWb_Wqq_Eta"        , ";#eta"      , nBinsEta, minEta, maxEta);
  h_bHt_tWb_Wqq_qq_dR           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_tWb_Wqq_qq_dR"      , ";#DeltaR"      , nBinsdR, mindR, maxdR);
  h_bHt_tWb_Wqq_qq_dEta         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_bHt_tWb_Wqq_qq_dEta"    , ";#Delta#eta"   , nBinsdEta, mindEta, maxdEta);
  
  //For Efficiency plots
  h_TopQuarkPt_isGenuineTop = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "h_TopQuarkPt_isGenuineTop",  ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  //h_TopQuarkPt_isGenuineTop = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_TopQuarkPt_isGenuineTop"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  
  //Genuine Top
  h_GenuineTop_ldgJet_Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_ldgJet_Pt"     , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineTop_subldgJet_Pt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_subldgJet_Pt"  , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineTop_Bjet_Pt        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_Bjet_Pt"       , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineTop_Dijet_Pt       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_Dijet_Pt"      , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineTop_Trijet_Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_Trijet_Pt"     , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineTop_Dijet_dR       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_Dijet_dR"      , ";#DeltaR"      , nBinsdR, mindR, maxdR);
  h_GenuineTop_ldgJet_Bjet_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_ldgJet_Bjet_dR", ";#DeltaR"      , nBinsdR, mindR, maxdR);
  h_GenuineTop_subJet_Bjet_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_subJet_Bjet_dR", ";#DeltaR"      , nBinsdR, mindR, maxdR);
  h_GenuineTop_Dijet_Bjet_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_Dijet_Bjet_dR" , ";#DeltaR"      , nBinsdR, mindR, maxdR);
  h_GenuineTop_Dijet_M        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_Dijet_M"       , ";M (GeV/c^{2})", 100    , minM , 200.);
  h_GenuineTop_Trijet_M       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineTop_Trijet_M"      , ";M (GeV/c^{2})", 100    , minM , 500.);

  // GenJets                                                                                                                                
  h_GenJets_N   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJets_N" , ";genJet multiplicity", 30, 0.0, 30.0);
  h_GenJet1_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet1_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet2_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet2_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet3_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet3_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet4_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet4_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet5_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet5_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet6_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet6_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  //                                                                                                                                        
  h_GenJet1_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet1_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet2_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet2_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet3_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet3_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet4_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet4_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet5_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet5_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet6_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet6_Eta", ";#eta", nBinsEta, minEta, maxEta);
  // GenJets: Trijet with largest pT                                                                                                        
  h_MaxTriJetPt_Pt       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_Pt"     , ";p_{T} (GeV/c)", nBinsPt  , minPt  , maxPt );
  h_MaxTriJetPt_Eta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_Eta"    , ";#eta"         , nBinsEta , minEta , maxEta);
  h_MaxTriJetPt_Rap      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_Rap"    , ";#omega"       , nBinsRap , minRap , maxRap);
  h_MaxTriJetPt_Mass     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_Mass"   , ";M (GeV/c^{2})", 200   , minM   , 1000.  );
  h_MaxTriJetPt_dEtaMax  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dEtaMax", ";#Delta#eta"   , nBinsdEta, mindEta, maxdEta);
  h_MaxTriJetPt_dPhiMax  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dPhiMax", ";#Delta#phi"   , nBinsdPhi, mindPhi, maxdPhi);
  h_MaxTriJetPt_dRMax    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dRMax"  , ";#DeltaR"      , nBinsdR  , mindR  , maxdR  );
  h_MaxTriJetPt_dEtaMin  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dEtaMin", ";#Delta#eta"   , nBinsdEta, mindEta, maxdEta);
  h_MaxTriJetPt_dPhiMin  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dPhiMin", ";#Delta#phi"   , nBinsdPhi, mindPhi, maxdPhi);
  h_MaxTriJetPt_dRMin    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dRMin"  , ";#DeltaR"      , nBinsdR  , mindR  , maxdR  );
  h_MaxTriJetPt_dEtaAverage  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dEtaAverage", ";#Delta#eta"  , nBinsdEta, mindEta, maxdEta);
  h_MaxTriJetPt_dPhiAverage  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dPhiAverage", ";#Delta#phi"  , nBinsdPhi, mindPhi, maxdPhi);
  h_MaxTriJetPt_dRAverage    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "MaxTriJetPt_dRAverage"  , ";#Delta#omega", nBinsdRap, mindRap, maxdRap);
  h_TriJetMaxPt_bHt_TQuark_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "TriJetMaxPt_bHt_TQuark_dR",";#DeltaR"      , nBinsdR  , mindR  , maxdR  );
  h_TriJetMaxPt_Over_bHt_TQuarkPt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "h_TriJetMaxPt_Over_bHt_TQuarkPt", "" ,nBinsdR  , mindR  , maxdR  );
  h_TriJetMaxPt_Top_M    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "h_TriJetMaxPt_Top_M", ";M (GeV/c^{2})", 100   , minM   , 500.  );
  h_TriJetMaxPt_Top_OverCut_M = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "h_TriJetMaxPt_Top_OverCut_M", ";M (GeV/c^{2})", 100   , minM   , 500.  );


  return;
}

void RecGenHToHW::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
}




void RecGenHToHW::process(Long64_t entry) {

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
  

  //================================================================================================
  // 7) Jet Selection
  //==selectedJets==============================================================================================
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
  vector<genParticle> selectedBQuarks = GetGenParticles(fEvent.genparticles().getGenParticles(), 10, 5, true, false);
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
  
  std::vector<math::XYZTLorentzVector> selBJets_p4;
  math::XYZTLorentzVector Bjet_p4;
  for(auto Bjet: selectedBJets)
    {
      Bjet_p4 = Bjet.p4();
      selBJets_p4.push_back( Bjet_p4 );
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
  
  // Fill Histograms                                                                                                                        
  h_GenJets_N->Fill(selectedJets.size());
  
  bool bSkipEvent = false;
  // 4-momenta                                                                                                                                                               
  math::XYZTLorentzVector Htb_HPlus_p4;
  math::XYZTLorentzVector bHt_TQuark_p4;
  math::XYZTLorentzVector bHt_tWb_Wqq_p4;
  math::XYZTLorentzVector bHt_tbW_BQuark_p4;
  math::XYZTLorentzVector bHt_tWb_WBoson_p4;
  math::XYZTLorentzVector bHt_HWh_WBoson_p4;
  math::XYZTLorentzVector neutFromTau_p4;
  math::XYZTLorentzVector bHt_HWh_h_vistau_p4;
  math::XYZTLorentzVector htau1_p4;
  math::XYZTLorentzVector htau2_p4;
  std::vector<math::XYZTLorentzVector> v_bHt_HWh_h_vistau_p4;
  std::vector<math::XYZTLorentzVector> v_bHt_tWb_Wqq_p4;
  
  //bool assWHadr = false;
  // Define the table                                                                                                                                                                                     
  Table table("Evt | Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | D0 (mm) | Lxy (mm) | Mom | Daughters", "Text"); //LaTeX or Text
 
  int row = 0;
  
  //int qourkc = 0;
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
	Htb_HPlus_p4 = p.p4();
      }
    
    if (0) std::cout << "=== Associated top" << std::endl;
    if ( abs(genP_pdgId) == 6)
      {
        if (!bIsLastCopy) continue;
	bHt_TQuark_p4 = p.p4();
	for (auto& d: genP_daughters)
	  {
	    if ( abs(d.pdgId()) == 5)
	      {
		bHt_tbW_BQuark_p4 = d.p4();
	      }
	  }
      }
    
    if (0) std::cout << "=== W+ or W-" << std::endl;
    if ( abs(genP_pdgId) == 24)
      {
        if (!bIsLastCopy) continue;
	int mompdgId = mcTools.RecursivelyLookForMotherId(p);
	if (abs(mompdgId) == 6)
          {
	    bHt_tWb_WBoson_p4 = p.p4();
	    for (auto& d: genP_daughters)
              {
		if ( mcTools.IsLepton(d.pdgId()))
		  {
		    bSkipEvent = true;
		  }
		else if (abs(d.pdgId()) == 15)
		  {
		    bSkipEvent = true;
		  }
		else 
		  {
		    bHt_tWb_Wqq_p4 = d.p4();
		    v_bHt_tWb_Wqq_p4.push_back(bHt_tWb_Wqq_p4);
		  }
	      }//for (auto& d: genP_daughters)
	  }//if (abs(mompdgId) == 6)
	else if (abs(mompdgId) == 37)
          {
	    bHt_HWh_WBoson_p4 = p.p4();
	    for (auto& d: genP_daughters)
              {
		if (mcTools.IsQuark(d.pdgId()) || abs(d.pdgId()) == 11)
		  {
		    bSkipEvent = true;
		  }
	      }
	  }
      }//if ( abs(genP_pdgId) == 24)
    
    if ( abs(genP_pdgId) == 15)
      {
        if (!bIsLastCopy) continue;
	int momTaupdgId = mcTools.RecursivelyLookForMotherId(p);
	if ( abs(momTaupdgId) == 25 || abs(momTaupdgId) == 35)
          {
	    for (auto& d: genP_daughters)                                                                                                  \
              {
		if ( mcTools.IsLepton(d.pdgId()))
                  {
                    bSkipEvent = true;
                  }
		if (abs(d.pdgId()) == 16 )
		  {
		    neutFromTau_p4 = d.p4();
		  }
	      }//for (auto& d: genP_daughters)
	    bHt_HWh_h_vistau_p4 = p.p4() - neutFromTau_p4;
	    v_bHt_HWh_h_vistau_p4.push_back(bHt_HWh_h_vistau_p4);
	  }//if ( abs(momTaupdgId) == 25 || abs(momTaupdgId) == 35)
      }//if ( abs(genP_pdgId) == 15)
    


  } //for (auto& p: fEvent.genparticles().getGenParticles()) 
  //h_eventcount -> Fill(allevent);
  if ( bSkipEvent ) return;  
  
  if (v_bHt_HWh_h_vistau_p4.size() > 0 )
    {
      for (size_t i=0; i < v_bHt_HWh_h_vistau_p4.size()-1; i++)
	{
	  for (size_t j=0; j < v_bHt_HWh_h_vistau_p4.size(); j++)
	    {
	      htau1_p4 = v_bHt_HWh_h_vistau_p4.at(j);
	      htau2_p4 = v_bHt_HWh_h_vistau_p4.at(i);
	      if (v_bHt_HWh_h_vistau_p4.at(i).pt() > v_bHt_HWh_h_vistau_p4.at(j).pt())
		{
		  htau1_p4 = v_bHt_HWh_h_vistau_p4.at(i);
		  htau2_p4 = v_bHt_HWh_h_vistau_p4.at(j);
		}
	    }
	}
    }
  double htau1_htau2_dR = ROOT::Math::VectorUtil::DeltaR(htau1_p4, htau2_p4);
  h_htau1_htau2_dR -> Fill(htau1_htau2_dR);
  if (htau1_htau2_dR <= 0.8)
    {
      h_HPus_dRCut_Pt ->Fill(Htb_HPlus_p4.pt());
      cFatTau.increment();
    }
  // GenPa histo

  h_bHt_tbW_BQuark_Pt -> Fill(bHt_tbW_BQuark_p4.pt());
  h_bHt_tbW_BQuark_Eta -> Fill(bHt_tbW_BQuark_p4.eta());
  
  for (size_t i=0; i < v_bHt_tWb_Wqq_p4.size(); i++)
    {
      h_bHt_tWb_Wqq_Pt -> Fill(v_bHt_tWb_Wqq_p4.at(i).Pt());
      h_bHt_tWb_Wqq_Eta -> Fill(v_bHt_tWb_Wqq_p4.at(i).eta());
    }
  math::XYZTLorentzVector GenT_LdgQuark_p4, GenT_SubldgQuark_p4;
  if (v_bHt_tWb_Wqq_p4.size() > 1)
    {
      if (v_bHt_tWb_Wqq_p4.at(0).pt() > v_bHt_tWb_Wqq_p4.at(1).pt() )
	{
	  GenT_LdgQuark_p4     = v_bHt_tWb_Wqq_p4.at(0);
	  GenT_SubldgQuark_p4 = v_bHt_tWb_Wqq_p4.at(1);
	}
      else
	{
	  GenT_LdgQuark_p4     = v_bHt_tWb_Wqq_p4.at(1);
          GenT_SubldgQuark_p4 = v_bHt_tWb_Wqq_p4.at(0);
	}
      double bHt_tWb_Wqq_qq_dR = ROOT::Math::VectorUtil::DeltaR(v_bHt_tWb_Wqq_p4.at(0), v_bHt_tWb_Wqq_p4.at(1));
      h_bHt_tWb_Wqq_qq_dR -> Fill(bHt_tWb_Wqq_qq_dR);
      double bHt_tWb_Wqq_qq_dEta = std::abs(v_bHt_tWb_Wqq_p4.at(0).eta() - v_bHt_tWb_Wqq_p4.at(1).eta());
      h_bHt_tWb_Wqq_qq_dEta -> Fill(bHt_tWb_Wqq_qq_dEta);
    }
  if(0) std::cout << "=== bag1" << std::endl;
  const double twoSigmaDpt = 0.32;
  const double dRcut    = 0.4;
  double dRmin  = 9999.9;
  // Do matching
  math::XYZTLorentzVector mcMatchedT_BJet(0,0,0,0), mcMatchedT_LdgJet(0,0,0,0), mcMatchedT_SubldgJet(0,0,0,0);
  
  if(0) std::cout << "=== bag2" << std::endl;
  for(size_t i=0; i < selBJets_p4.size(); i++)
    {
      double dR = ROOT::Math::VectorUtil::DeltaR(selBJets_p4.at(i) , bHt_tbW_BQuark_p4 );
      if(0) std::cout << "=== dR= " << dR << std::endl;
      double dPtOverPt = std::abs((selBJets_p4.at(i).pt() - bHt_tbW_BQuark_p4.pt())/bHt_tbW_BQuark_p4.pt());
      if (dR <= dRmin )
	{
	  if(0) std::cout << "=== dPtOverPt = " << dPtOverPt << std::endl;
	  if (dPtOverPt < twoSigmaDpt) 
	    {
	      if(0) std::cout << "=== dRmin = dR" << std::endl;
	      dRmin = dR;
	      mcMatchedT_BJet = selBJets_p4.at(i);
	    }
	}
    }
  //h_GenuineTop_Bjet_Pt           -> Fill(mcMatchedT_BJet.pt());
  double dR1min, dR2min, dPtOverPt1min, dPtOverPt2min;
  dR1min = dR2min = dPtOverPt1min = dPtOverPt2min = 99999.9;
  if(0) std::cout << "=== bag3" << std::endl;
  for(size_t i=0; i < selJets_p4.size(); i++)
    {
      // Skip the jets that are matched with bquarks
      double dRn = 9999.9;
      if (dRmin < dRcut && mcMatchedT_BJet.pt() > 0.)
	{
	  dRn = ROOT::Math::VectorUtil::DeltaR(selJets_p4.at(i) , mcMatchedT_BJet);
	  
	}
      if( dRn <= 0.8) continue;
      
      // Find dR for the two jets in top-decay dijet
      
      double dR1 = ROOT::Math::VectorUtil::DeltaR(selJets_p4.at(i), GenT_LdgQuark_p4);
      double dR2 = ROOT::Math::VectorUtil::DeltaR(selJets_p4.at(i), GenT_SubldgQuark_p4);
      
      // Require both jets to be within dR <= dRcut 
      // Calculate dPtOverPt for each jet in top-decay dijet   
      
      double dPtOverPt1 = std::abs((selJets_p4.at(i).pt() - GenT_LdgQuark_p4.pt())/GenT_LdgQuark_p4.pt());
      double dPtOverPt2 = std::abs((selJets_p4.at(i).pt() - GenT_SubldgQuark_p4.pt())/GenT_SubldgQuark_p4.pt());
      if(0) std::cout << "=== bag4" << std::endl;
      // Find which of the two is the correct match                                                                                 
      if (dR1 < dR2)
	{
	  // Is Jet1 closer in eta-phi AND has smaller pT difference?                                                               
	  if (dR1 < dR1min)
	    {
	      if (dPtOverPt1 < twoSigmaDpt)
		{
		  dR1min = dR1;
		  dPtOverPt1min= dPtOverPt1;
		  mcMatchedT_LdgJet = selJets_p4.at(i);
		}
	    }
	  
	  // Is Jet2 closer in eta-phi AND has smaller pT difference?                                                           
	  //else if (dR2 <= dRcut && dR2 < dR2min)  at least two matched jets                                                       
	  else if (dR2 < dR2min)                  //at least two matched jets                                                       
	    {
	      if (dPtOverPt2 < twoSigmaDpt)
		{
		  dR2min  = dR2;
		  dPtOverPt2min = dPtOverPt2;
		  mcMatchedT_SubldgJet = selJets_p4.at(i);
		}
	    }
	}
      else
	{
	  // Is Jet2 closer in eta-phi AND has smaller pT difference?                                                              
	  if (dR2 < dR2min)
	    {
	      if (dPtOverPt2 < twoSigmaDpt)
		{
		  dR2min  = dR2;
		  dPtOverPt2min = dPtOverPt2;
		  mcMatchedT_SubldgJet = selJets_p4.at(i);
		}
	    }

	  // Is Jet2 closer in eta-phi AND has smaller pT difference? 
	  
	  else if (dR1 < dR1min)                //at least two matched jets                                                         
	    {
	      if  (dPtOverPt1 < twoSigmaDpt)
		{
		  dR1min  = dR1;
		  dPtOverPt1min = dPtOverPt1;
		  mcMatchedT_LdgJet = selJets_p4.at(i);
		}
	    }
	}
    }//for(auto Bjet: selectedBJets)
  if(0) std::cout << "=== bag5" << std::endl;
  // Check if TOP is genuine                                                                                                       
  if (dR1min<= dRcut && dR2min <= dRcut && dRmin <= dRcut)
    {
      if(0) std::cout << "=== dR1min= " << dR1min << std::endl;
      if(0) std::cout << "=== dR2min= " << dR2min << std::endl;
      if(0) std::cout << "=== dRmin= " << dRmin << std::endl;
    }
  bool isGenuine = false;
  if  (mcMatchedT_LdgJet.pt() > 0. && mcMatchedT_SubldgJet.pt() > 0. && mcMatchedT_BJet.pt() > 0. )
    {
      isGenuine = (dR1min<= dRcut && dR2min <= dRcut && dRmin <= dRcut);
      bool twoJMatched = isGenuine;  //at least two matched jets two top decay products (used for plots)
      if (!twoJMatched) twoJMatched = max(dR1min, dR2min) <= dRcut || max(dR1min, dRmin) <= dRcut || max(dR2min, dRmin) <= dRcut;
    }
  if(0) std::cout << "=== bag5" << std::endl;
  h_TopQuarkPt_isGenuineTop -> Fill(isGenuine, bHt_TQuark_p4.pt());
  if (isGenuine)
    {
      if(0) std::cout << "=== bag51" << std::endl;
      math::XYZTLorentzVector ldgJet_p4;
      math::XYZTLorentzVector subldgJet_p4;
      math::XYZTLorentzVector Tbjet_p4, dijet_p4, trijet_p4;
      if(0) std::cout << "=== bag511" << std::endl;
      ldgJet_p4    = mcMatchedT_LdgJet;
      if(0) std::cout << "=== bag512" << std::endl;
      subldgJet_p4 = mcMatchedT_SubldgJet;
      if(0) std::cout << "=== bag513" << std::endl;
      Tbjet_p4 = mcMatchedT_BJet;
      if(0) std::cout << "=== bag58" << std::endl;
      dijet_p4  = ldgJet_p4 + subldgJet_p4;
      trijet_p4 = ldgJet_p4 + subldgJet_p4 + Tbjet_p4;
      if(0) std::cout << "=== bag52" << std::endl;
      double dR_j1j2 = ROOT::Math::VectorUtil::DeltaR(ldgJet_p4, subldgJet_p4);
      double dR_j1b  = ROOT::Math::VectorUtil::DeltaR(ldgJet_p4, Tbjet_p4);
      double dR_j2b  = ROOT::Math::VectorUtil::DeltaR(subldgJet_p4, Tbjet_p4);
      
      // bool mergedDijet            = dR_j1j2 < 0.8 && dR_j1b > 0.8 && dR_j2b > 0.8 && max(ldgJet.bjetDiscriminator(), subldgJet.bjetDiscriminator()) < 0.5426;
      // bool mergedTrijet_untaggedW = max(dR_j1j2, max(dR_j1b, dR_j2b)) < 0.8 && max(ldgJet.bjetDiscriminator(), subldgJet.bjetDiscriminator()) < 0.5426;
      /*
      bool mergedDijet            = dR_j1j2 < 0.8 && dR_j1b > 0.8 && dR_j2b > 0.8;
      bool mergedTrijet           = max(dR_j1j2, max(dR_j1b, dR_j2b)) < 0.8;
      bool mergedJB               = dR_j1j2 > 0.8 && max(dR_j1b, dR_j2b) > 0.8 && min(dR_j1b, dR_j2b) < 0.8;
      */
      double dR_dijet_bjet = ROOT::Math::VectorUtil::DeltaR(dijet_p4 , Tbjet_p4);
      if(0) std::cout << "=== bag53" << std::endl;
      h_GenuineTop_ldgJet_Pt         -> Fill(ldgJet_p4.pt());
      h_GenuineTop_subldgJet_Pt      -> Fill(subldgJet_p4.pt());
      h_GenuineTop_Bjet_Pt           -> Fill(Tbjet_p4.pt());
      h_GenuineTop_Dijet_dR          -> Fill(dR_j1j2);
      h_GenuineTop_ldgJet_Bjet_dR    -> Fill(dR_j1b);
      h_GenuineTop_subJet_Bjet_dR    -> Fill(dR_j2b); 
      h_GenuineTop_Dijet_Bjet_dR     -> Fill(dR_dijet_bjet);
      h_GenuineTop_Dijet_Pt          -> Fill(dijet_p4.pt());
      h_GenuineTop_Dijet_M           -> Fill(dijet_p4.M());
      h_GenuineTop_Trijet_Pt         -> Fill(trijet_p4.pt());
      h_GenuineTop_Trijet_M          -> Fill(trijet_p4.M());
    }



  //==================================================================================================================
  //jet rec
  if(0) std::cout << "=== bag6" << std::endl;
  int iJet = 0;
  if (selJets_p4.size() > 1) {
    // For-loop: All selected jets                                                                                                          
    for (size_t i=0; i < selJets_p4.size(); i++)
      {
        iJet++;
	double genJ_Pt  = selJets_p4.at(i).pt();
        double genJ_Eta = selJets_p4.at(i).eta();
        // double genJ_Rap = mcTools.GetRapidity(selJets_p4.at(i));                                                                         

	if (iJet==1)
          {
	    
            h_GenJet1_Pt -> Fill( genJ_Pt  );
            h_GenJet1_Eta-> Fill( genJ_Eta );
	    
          }
        else if (iJet==2)
          {
	    
            h_GenJet2_Pt -> Fill( genJ_Pt  );
            h_GenJet2_Eta-> Fill( genJ_Eta );

          }
        else if (iJet==3)
          {

            h_GenJet3_Pt -> Fill( genJ_Pt  );
            h_GenJet3_Eta-> Fill( genJ_Eta );
	  }
        else if (iJet==4)
          {

            h_GenJet4_Pt -> Fill( genJ_Pt  );
            h_GenJet4_Eta-> Fill( genJ_Eta );

          }
        else if (iJet==5)
          {

            h_GenJet5_Pt -> Fill( genJ_Pt  );
            h_GenJet5_Eta-> Fill( genJ_Eta );

          }
	else if (iJet==6)
          {

            h_GenJet6_Pt -> Fill( genJ_Pt  );
            h_GenJet6_Eta-> Fill( genJ_Eta );
	    
          }
	else{}
      } // for (size_t i=0; i < selJets_p4.size(); i++)
  }  





  //////////////////////////////////////////////////////////////////////////////////////////////////////                                    
  // GenJets: Trijet with largest pT                                                                                                        
  //////////////////////////////////////////////////////////////////////////////////////////////////////                                    
  if(0) std::cout << "=== Trijet (Max Pt)" << std::endl;
  int nTriJetCands = 0;
  double dR_ij     = 0.0;
  double dR_ik     = 0.0;
  double dR_jk     = 0.0;
  double dEta_ij   = 0.0;
  double dEta_ik   = 0.0;
  double dEta_jk   = 0.0;
  double dPhi_ij   = 0.0;
  double dPhi_ik   = 0.0;
  double dPhi_jk   = 0.0;
  double dRSum     = 0.0;
  double dEtaSum   = 0.0;
  double dPhiSum   = 0.0;
  std::vector<double> dEta_ijk;
  std::vector<double> dPhi_ijk;
  std::vector<double> dR_ijk;

  math::XYZTLorentzVector TriJetMaxPt_p4(0,0,0,0), SubldgTrijet_p4(0,0,0,0), j1(0,0,0,0), j2(0,0,0,0), j3(0,0,0,0);
  vector<math::XYZTLorentzVector> LdgTrijet_JETS, SubldgTrijet_JETS;

  if (selJets_p4.size() > 2)
    {

      // For-loop: All Selected Jets                                                                                                        
      for (auto i = selJets_p4.begin(); i != selJets_p4.end()-2; ++i) {

        // Initialise values                                                                                                                
	dR_ij     = 0.0;
	dR_ik     = 0.0;
        dR_jk     = 0.0;
        dEta_ij   = 0.0;
        dEta_ik   = 0.0;
        dEta_jk   = 0.0;
        dPhi_ij   = 0.0;
	dPhi_ik   = 0.0;
        dPhi_jk   = 0.0;

        // For-loop: All Selected Jets (nested)                                                                                             
        for (auto j = i+1; j != selJets_p4.end()-1; ++j) {
	  
	  // For-loop: All Selected Jets (doubly-nested)                                                                                      
          for (auto k = j+1; k != selJets_p4.end(); ++k) {
            nTriJetCands++;

            // Calculations                                                                                                                 
            dR_ij = ROOT::Math::VectorUtil::DeltaR(*i, *j);
            dR_ijk.push_back(dR_ij);
            dR_ik = ROOT::Math::VectorUtil::DeltaR(*i, *k);
            dR_ijk.push_back(dR_ik);
            dR_jk = ROOT::Math::VectorUtil::DeltaR(*j, *k);
            dR_ijk.push_back(dR_jk);

            dEta_ij = std::abs( i->Eta() - j->Eta() );
            dEta_ijk.push_back(dEta_ij);
            dEta_ik = std::abs( i->Eta() - k->Eta() );
            dEta_ijk.push_back(dEta_ik);
            dEta_jk = std::abs( j->Eta() - k->Eta() );
            dEta_ijk.push_back(dEta_jk);

            dPhi_ij = std::abs( ROOT::Math::VectorUtil::DeltaPhi(*i, *j) );
            dPhi_ijk.push_back(dEta_ij);
            dPhi_ik = std::abs( ROOT::Math::VectorUtil::DeltaPhi(*i, *k) );
            dPhi_ijk.push_back(dEta_ik);
            dPhi_jk = std::abs( ROOT::Math::VectorUtil::DeltaPhi(*j, *k) );
            dPhi_ijk.push_back(dEta_jk);

	    math::XYZTLorentzVector p4 = *i + *j + *k;
            Bool_t BjetInTrijet = IsBjet(*i,selectedBJets) || IsBjet(*j,selectedBJets) || IsBjet(*k,selectedBJets);

	    if ( p4.Pt() > TriJetMaxPt_p4.Pt() && BjetInTrijet){
	      SubldgTrijet_p4 = TriJetMaxPt_p4;
              SubldgTrijet_JETS.clear();
              SubldgTrijet_JETS.push_back(j1); SubldgTrijet_JETS.push_back(j2); SubldgTrijet_JETS.push_back(j3);
              TriJetMaxPt_p4 = p4;
              j1 = *i;
              j2 = *j;
	      j3 = *k;
              LdgTrijet_JETS.clear();
              LdgTrijet_JETS.push_back(*i); LdgTrijet_JETS.push_back(*j); LdgTrijet_JETS.push_back(*k);
            }

	    // Calculate                                                                                                                    
            dRSum   += (dR_ij   + dR_ik   + dR_jk  )/3;
            dEtaSum += (dEta_ij + dEta_ik + dEta_jk)/3;
            dPhiSum += (dPhi_ij + dPhi_ik + dPhi_jk)/3;
          }
        }
      }
    }

  // Fill Histos                                                                                                                            
  h_MaxTriJetPt_Pt          ->Fill( TriJetMaxPt_p4.Pt()  );
  h_MaxTriJetPt_Eta         ->Fill( TriJetMaxPt_p4.Eta() );
  h_MaxTriJetPt_Rap         ->Fill( mcTools.GetRapidity(TriJetMaxPt_p4));
  h_MaxTriJetPt_Mass        ->Fill( TriJetMaxPt_p4.M() );
  h_MaxTriJetPt_dEtaMax     ->Fill( *max_element(dEta_ijk.begin(), dEta_ijk.end() ) );
  h_MaxTriJetPt_dPhiMax     ->Fill( *max_element(dPhi_ijk.begin(), dPhi_ijk.end() ) );
  h_MaxTriJetPt_dRMax       ->Fill( *max_element(dR_ijk.begin()  , dR_ijk.end()   ) );
  h_MaxTriJetPt_dEtaMin     ->Fill( *min_element(dEta_ijk.begin(), dEta_ijk.end() ) );
  h_MaxTriJetPt_dPhiMin     ->Fill( *min_element(dPhi_ijk.begin(), dPhi_ijk.end() ) );
  h_MaxTriJetPt_dRMin       ->Fill( *min_element(dR_ijk.begin()  , dR_ijk.end()   ) );
  h_MaxTriJetPt_dEtaAverage ->Fill( dEtaSum/nTriJetCands );
  h_MaxTriJetPt_dPhiAverage ->Fill( dPhiSum/nTriJetCands );
  h_MaxTriJetPt_dRAverage   ->Fill( dRSum/nTriJetCands   );
  
  double TriJetMaxPt_bHt_TQuark_dR = ROOT::Math::VectorUtil::DeltaR(TriJetMaxPt_p4, bHt_TQuark_p4);
  double dPtOverPt1 = std::abs((TriJetMaxPt_p4.pt() - bHt_TQuark_p4.pt())/bHt_TQuark_p4.pt());
  if (TriJetMaxPt_bHt_TQuark_dR < 0.1 )
    {
      h_TriJetMaxPt_Top_M -> Fill(TriJetMaxPt_p4.M());
    }
  h_TriJetMaxPt_bHt_TQuark_dR -> Fill(TriJetMaxPt_bHt_TQuark_dR);
  h_TriJetMaxPt_Over_bHt_TQuarkPt -> Fill(dPtOverPt1);
  if (dPtOverPt1 < 0.32)
    {
      h_TriJetMaxPt_Top_OverCut_M ->Fill(TriJetMaxPt_p4.M());
    }
  
  
  return;
}



//calculation of W transverse mass
double RecGenHToHW::GetWMt(const math::XYVector muon, const math::XYVector& met){
  double metEt  = met.R();
  double muonEt = muon.r();
  double mT     = -999.9;
  double mTSq   =   0.0;
  double EtSq   = (metEt + muonEt) * (metEt + muonEt);
  double EtXSq  = (muon.x() + met.x()) * (muon.x() + met.x());
  double EtYSq  = (muon.y() + met.y()) * (muon.y() + met.y());
  mTSq = EtSq - (EtXSq + EtYSq);
  if (mTSq >= 0) mT = std::sqrt(mTSq);
  return mT;
}


//calculation of W transverse mass second approach
double RecGenHToHW::GetWMt2(const math::XYZTLorentzVector muon, const math::XYZTLorentzVector neut){
  double muonPt = muon.pt();
  double neutPt = neut.pt();
  double dPhi   = std::abs(ROOT::Math::VectorUtil::DeltaPhi(muon, neut));
  double mT     = -999.9;
  double mTSq   = 0.0;
  mTSq          = 2. * muonPt * neutPt * (1. - std::cos(dPhi));
  if (mTSq >= 0) mT = std::sqrt(mTSq);
  return mT;
}



double RecGenHToHW::GetMt(const math::XYVector tau1, const math::XYVector tau2, const math::XYVector muon, const math::XYVector& met) {
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



vector<GenJet> RecGenHToHW::GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCuts, std::vector<float> etaCuts, vector<genParticle> genParticlesToMatch)
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


Bool_t RecGenHToHW::IsBjet(math::XYZTLorentzVector jet, vector<GenJet> selectedBJets)
{
  for (auto& bjet: selectedBJets)
    {
      double dR = ROOT::Math::VectorUtil::DeltaR(jet, bjet.p4());
      if (dR < 0.01) return true;
    }
  return false;
}





vector<GenJet> RecGenHToHW::GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts)
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

vector<genParticle> RecGenHToHW::GetAllCopies(const vector<genParticle> genParticles, genParticle genP)
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

vector<genParticle> RecGenHToHW::GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy, const bool hasNoDaughters)
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
