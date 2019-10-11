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
  const genParticle GetLastCopy(const vector<genParticle> genParticles, const genParticle &p);
  virtual vector<genParticle> GetGenParticlesRec(const vector<genParticle> genParticles, const int pdgId);
  Bool_t AreHadroDau(const vector<genParticle> genParticles, const vector<genParticle> genP);
  Bool_t IshadronicTop(const vector<genParticle> genParticles, const vector<genParticle> genP);


private:
  CommonPlots fCommonPlots;
  TauSelection fTauSelection;
  JetSelection fJetSelection;
  //ElectronSelection fElectronSelection;
  //MuonSelection fMuonSelection;

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
  //Count cElectronVeto;
  Count cMuonVeto;
  Count cTauVeto;
  Count cJetSelection;
  //Count cMETSelection;
  // Count cTopologySelection;
  Count cSelected;
  Count cgenTop;
  //Count cFatTau;
  // BR Counters
  //Count cInclusive;
  //Count cbHt_HPlus;
  //Count cbHt_HBoson;
  //Count cbHt_WBoson;


  // GenParticles                                                                                                                                                                                           
  WrappedTH1 *hParticle_Pt;
  WrappedTH1 *h_htau1_htau2_dR;
  WrappedTH1 *h_HPus_dRCut_Pt;

  //For Efficiency plots
  WrappedTH1Triplet *h_TopQuarkPt_isGenuineTop;
  //WrappedTH1 *h_TopQuarkPt_isGenuineTop;
  //Genuine W from Top
  WrappedTH1 *h_GenuineWfromTop_Dijet_M;
  WrappedTH1 *h_GenuineWfromTop_ldgJet_Pt;
  WrappedTH1 *h_GenuineWfromTop_subldgJet_Pt;
  WrappedTH1 *h_GenuineWfromTop_Dijet_Pt;
  WrappedTH1 *h_GenuineWfromTop_Dijet_dR;
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
  // second reco top
  WrappedTH1 *h_SecApro_TQuarkCan_M;
  WrappedTH1 *h_SecApro_WBosonCan_M;
  WrappedTH1 *h_SecApro_TQuark_Topjet_dR;
  WrappedTH1 *h_SecApro_TQuark_M;
  WrappedTH1 *h_SecApro_WBoson_M;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(RecGenHToHW);

RecGenHToHW::RecGenHToHW(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2HWAnalysis, fHistoWrapper),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    //fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    //fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
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
    //cElectronVeto(fEventCounter.addCounter("e-veto")),
    cMuonVeto(fEventCounter.addCounter("#mu-veto")),
    cTauVeto(fEventCounter.addCounter("#tau-veto")),
    cJetSelection(fEventCounter.addCounter("Jets + H_{T}")),
    //cMETSelection(fEventCounter.addCounter("MET")),
    // cTopologySelection(fEventCounter.addCounter("Topology")),
    cSelected(fEventCounter.addCounter("All Selections")),
    cgenTop(fEventCounter.addSubCounter("Branching", "Genuine Top")) 
    //cFatTau(fEventCounter.addCounter("FatTau")),
    //cInclusive(fEventCounter.addSubCounter("Branching", "All events")),
    //cbHt_HPlus(fEventCounter.addSubCounter("Branching", "H+")),
    //cbHt_HBoson(fEventCounter.addSubCounter("Branching", "H+->HW, H")),
    //cbHt_WBoson(fEventCounter.addSubCounter("Branching", "H+->HW, W"))
{ }

void RecGenHToHW::book(TDirectory *dir) {
  
  // Book common plots histograms
  fCommonPlots.book(dir, isData());
  // Book histograms in event selection classes
  fTauSelection.bookHistograms(dir);
  fJetSelection.bookHistograms(dir);

  // Fixed binning
  const int nBinsPt   = 4*cfg_PtBinSetting.bins();
  const double minPt  = cfg_PtBinSetting.min();
  const double maxPt  = 4*cfg_PtBinSetting.max();
  
  const int nBinsdR  = 2*cfg_DeltaRBinSetting.bins();                                                                                                                                                        const double mindR = cfg_DeltaRBinSetting.min();                                                                                                                                                           const double maxdR = cfg_DeltaRBinSetting.max();  
  
  /*
  const int nBinsEta  = 2*cfg_EtaBinSetting.bins();
  const double minEta = cfg_EtaBinSetting.min();
  const double maxEta = 2*cfg_EtaBinSetting.max();
  
  const int nBinsRap  = cfg_EtaBinSetting.bins();
  const double minRap = cfg_EtaBinSetting.min();
  const double maxRap = cfg_EtaBinSetting.max();
  */
  // const int nBinsPhi  = cfg_PhiBinSetting.bins();
  // const double minPhi = cfg_PhiBinSetting.min();
  // const double maxPhi = cfg_PhiBinSetting.max();
  const double minM = cfg_MassBinSetting.min();
  /*
  const int nBinsdEta  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindEta = cfg_DeltaEtaBinSetting.min();
  const double maxdEta = 2*cfg_DeltaEtaBinSetting.max();
  
  const int nBinsdRap  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindRap = cfg_DeltaEtaBinSetting.min();
  const double maxdRap = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdPhi  = cfg_DeltaPhiBinSetting.bins();
  const double mindPhi = cfg_DeltaPhiBinSetting.min();
  const double maxdPhi = 1.2*cfg_DeltaPhiBinSetting.max(); 
  */
  
  
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
  
  //For Efficiency plots
  h_TopQuarkPt_isGenuineTop = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "h_TopQuarkPt_isGenuineTop",  ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  //h_TopQuarkPt_isGenuineTop = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_TopQuarkPt_isGenuineTop"           , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  //Genuine W from Top 
  h_GenuineWfromTop_Dijet_M      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineWfromTop_Dijet_M"     , ";M (GeV/c^{2})", 100    , minM , 200. );
  h_GenuineWfromTop_ldgJet_Pt    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineWfromTop_ldgJet_Pt"   , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineWfromTop_subldgJet_Pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineWfromTop_subldgJet_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineWfromTop_Dijet_Pt     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineWfromTop_Dijet_Pt"    , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenuineWfromTop_Dijet_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_GenuineWfromTop_Dijet_dR"    , ";#DeltaR"      , nBinsdR, mindR, maxdR);

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
  //sec top rec
  h_SecApro_TQuarkCan_M       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_SecApro_TQuarkCan_M"      , ";M (GeV/c^{2})", 100    , minM , 500.);
  h_SecApro_WBosonCan_M       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_SecApro_WBosonCan_M"      , ";M (GeV/c^{2})", 100    , minM , 200.);
  h_SecApro_TQuark_M          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_SecApro_TQuark_M"         , ";M (GeV/c^{2})", 100    , minM , 500.);
  h_SecApro_WBoson_M          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_SecApro_WBoson_M"         , ";M (GeV/c^{2})", 100    , minM , 200.);
  h_SecApro_TQuark_Topjet_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "h_SecApro_TQuark_Topjet_dR" , ";#DeltaR"      , nBinsdR, mindR, maxdR);
  

  return;
}

void RecGenHToHW::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
}




void RecGenHToHW::process(Long64_t entry) {

  //if ( !fEvent.isMC() ) return;
  
  // Create MCools object
  MCTools mcTools(fEvent);
  
  // Increment Counter
  fCommonPlots.initialize();
  fCommonPlots.setFactorisationBinForEvent(std::vector<float> {});
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
  //const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  //if (eData.hasIdentifiedElectrons())
  //return;

  
  /*
  if (cfg_Verbose) std::cout << "=== Electron veto" << std::endl;
  vector<genParticle> selectedElectrons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_ElectronPtCut, cfg_ElectronEtaCut, 11, true, false);
  //if (0)
  // {
  //  std::cout << "\nnElectrons = " << selectedElectrons.size() << std::endl;
  //  for (auto& p: selectedElectrons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
  //}
  if ( selectedElectrons.size() > 0 ) return;
  cElectronVeto.increment();
  */
  
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
  
  //const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  //if (!muData.hasIdentifiedMuons()) return;
  cMuonVeto.increment();
  /*
  //================================================================================================
  // asso. top selection
  //================================================================================================  
  vector<genParticle> selectedTop = GetGenParticles(fEvent.genparticles().getGenParticles(),0., 10., 6, true, false);
  bool hadronicTop = IshadronicTop(fEvent.genparticles().getGenParticles(), selectedTop);
  if (!hadronicTop) return;
  */

  
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  if (!tauData.hasIdentifiedTaus()) return;
  
  
  //================================================================================================
  // 6) Tau selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Tau selection" << std::endl;
  vector<genParticle> selectedTaus = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_TauPtCut, cfg_TauEtaCut, 15, true, false);
  if ( selectedTaus.size() < 2 ) return;
  cTauVeto.increment();
  
  


  /*
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
  */
  
  
  //================================================================================================
  // 7) Jet selection
  //=============================================================================================== 
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, tauData.getSelectedTau());
  if (!jetData.passedSelection()) return;
  cJetSelection.increment();
  std::vector<math::XYZTLorentzVector> selJets_p4;                                                                                                                                                        
  math::XYZTLorentzVector jet_p4; 
  for(auto jet: jetData.getSelectedJets())
    {
      jet_p4 = jet.p4();
      //genJ_HT += jet.pt();
      selJets_p4.push_back( jet_p4 );
    }
  

  /*
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
  //  cMETSelection.increment();
  
  
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
  std::vector<math::XYZTLorentzVector> v_bHt_TQuark_p4;
  std::vector<math::XYZTLorentzVector> v_bHt_tbW_BQuark_p4;
  std::vector<genParticle> GenTops;
  //std::vector<genParticle> GenTops_BQuark;
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
	GenTops.push_back(p); //chris
	v_bHt_TQuark_p4.push_back(bHt_TQuark_p4);
	for (auto& d: genP_daughters)
	  {
	    if ( abs(d.pdgId()) == 5)
	      {
		bHt_tbW_BQuark_p4 = d.p4();
		v_bHt_tbW_BQuark_p4.push_back(bHt_TQuark_p4);
		//genParticle b = GetLastCopy(fEvent.genparticles().getGenParticles(), d);
		//GenTops_BQuark.push_back(d);
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
  // if ( bSkipEvent ) return;  

  std::vector<genParticle> GenTops_BQuark;
  std::vector<genParticle> GenTops_SubldgQuark;
  std::vector<genParticle> GenTops_LdgQuark;
  int cbquarkbag = 0;
  if(0) std::cout << "bag1 "  << std::endl;
  //vector<genParticle> GenTops = GetGenParticlesRec(fEvent.genparticles().getGenParticles(), 6);
  if(0) std::cout << "GenTops size = " << GenTops.size() << std::endl;
  // For-loop: All top quarks                                                                                                               
  
  for (auto& top: GenTops){
    genParticle  Gen_Wa;
    bool FoundBQuark = false;
    std::vector<genParticle> quarks;
    genParticle bquark;
    // For-loop: Top quark daughters (Nested)                                                                                               
    for (size_t i=0; i<top.daughters().size(); i++)
      {
        int dau_index = top.daughters().at(i);
        genParticle dau = fEvent.genparticles().getGenParticles()[dau_index];
	if(0) std::cout << "dau.pdgId()) == " << dau.pdgId() << std::endl;
        // B-Quark                                                                                                                          
        if (std::abs(dau.pdgId()) ==  5)
          {
            bquark = dau;
	    if(0) std::cout << "GenBQuark pt = " << bquark.pt() << std::endl;
            FoundBQuark = true;
	    cbquarkbag++;
	  }
	
        // W-Boson   
	if (std::abs(dau.pdgId()) == 24)
          {
            // Get the last copy                                                                                                            
            genParticle W = GetLastCopy(fEvent.genparticles().getGenParticles(), dau);
	    Gen_Wa = W;
	    
            // For-loop: W-boson daughters                                                                                                  
            for (size_t idau=0; idau<W.daughters().size(); idau++)
              {
                // Find the decay products of W-boson 
		
                int Wdau_index = W.daughters().at(idau);
                genParticle Wdau = fEvent.genparticles().getGenParticles()[Wdau_index];
		
                // Consider only quarks as decaying products                                                                                
                if (std::abs(Wdau.pdgId()) > 5) continue;
		
                // Save daughter                                                                                                            
                quarks.push_back(Wdau);
              }//W-boson daughters		
	    
          }//W-boson
      }//Top-quark daughters                                                                                                                
    
    // Skip top if b-quark is not found (i.e. top decays to W and c)                                                                        
    if (!FoundBQuark) continue;
    // Skip top if it decays leptonically (the "quarks" vector will be empty causing errors)                                                
    if(0) std::cout << "quarks size = " << quarks.size() << std::endl;
    if (quarks.size() < 2) continue;
    // Fill vectors for b-quarks, leading and subleading quarks coming from tops      
    GenTops_BQuark.push_back(bquark);
    if(0) std::cout << "BQuark size = " << GenTops_BQuark.size() << std::endl;
    if (quarks.at(0).pt() > quarks.at(1).pt())
      {
        GenTops_LdgQuark.push_back(quarks.at(0));
        GenTops_SubldgQuark.push_back(quarks.at(1));
      }
    else
      {
        GenTops_LdgQuark.push_back(quarks.at(1));
        GenTops_SubldgQuark.push_back(quarks.at(0));
      }
    
  } // For-Loop over top quarks  
  
  
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
      h_HPus_dRCut_Pt -> Fill(Htb_HPlus_p4.pt());
      //cFatTau.increment();
    }
  
  bool doMatching = (GenTops_BQuark.size() == GenTops.size());
  const double twoSigmaDpt = 0.39;
  const double dRcut = 0.4;
  vector <double> dRminB;
  std::vector<math::XYZTLorentzVector> MGen_LdgJet, MGen_SubldgJet, MGen_Bjet;
  if (doMatching) {
    for (size_t i=0; i<GenTops.size(); i++)
      {
	math::XYZTLorentzVector mcMatchedT_BJet;
	math::XYZTLorentzVector BQuark;
	BQuark = GenTops_BQuark.at(i).p4();
	double dRmin   = 9999.9;
	double dRminIn = 0.4;
	// Do matching
       	for(size_t j=0; j < selJets_p4.size(); j++)
	  {
	    double dR = ROOT::Math::VectorUtil::DeltaR(selJets_p4.at(j), BQuark);
	    if (dR > dRminIn) continue;
	    double dPtOverPt = std::abs((selJets_p4.at(j).pt() - BQuark.pt())/BQuark.pt());
	    if (dPtOverPt > twoSigmaDpt) continue;
	    dRmin = dR;
	    dRminIn = dR;
	    mcMatchedT_BJet = selJets_p4.at(j);  
	  }
	dRminB.push_back(dRmin);
	MGen_Bjet.push_back(mcMatchedT_BJet);
      }
  }

  
  //h_GenuineTop_Bjet_Pt           -> Fill(mcMatchedT_BJet.pt());
  if (doMatching) {
    for (size_t i=0; i<GenTops.size(); i++)
      {
	math::XYZTLorentzVector LdgQuark;
	LdgQuark = GenTops_LdgQuark.at(i).p4();
	math::XYZTLorentzVector SubldgQuark;
	SubldgQuark = GenTops_SubldgQuark.at(i).p4();
	math::XYZTLorentzVector BJet;
	BJet = MGen_Bjet.at(i);
	math::XYZTLorentzVector mcMatchedT_LdgJet, mcMatchedT_SubldgJet;
	double dR1min, dR2min, dPtOverPt1min, dPtOverPt2min;
	dR1min = dR2min = dPtOverPt1min = dPtOverPt2min = 99999.9;
	
	for(size_t j=0; j < selJets_p4.size(); j++)
	  {
	    // Skip the jets that are matched with bquarks
	    double dRn = 9999.9;
	    for (size_t k=0; k<GenTops.size(); k++)
	      {
		if (dRminB.at(k) < dRcut)
		  {
		    dRn = ROOT::Math::VectorUtil::DeltaR(selJets_p4.at(j) , BJet);
		  }
	      }
	    if( dRn <= 0.8) continue;
	    
	    // Find dR for the two jets in top-decay dijet
	    
	    double dR1 = ROOT::Math::VectorUtil::DeltaR(selJets_p4.at(j), LdgQuark);
	    double dR2 = ROOT::Math::VectorUtil::DeltaR(selJets_p4.at(j), SubldgQuark);
	    
	    // Require both jets to be within dR <= dRcut 
	    // Calculate dPtOverPt for each jet in top-decay dijet   
	    
	    double dPtOverPt1 = std::abs((selJets_p4.at(j).pt() - LdgQuark.pt())/LdgQuark.pt());
	    double dPtOverPt2 = std::abs((selJets_p4.at(j).pt() - SubldgQuark.pt())/SubldgQuark.pt());
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
			mcMatchedT_LdgJet = selJets_p4.at(j);
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
			mcMatchedT_SubldgJet = selJets_p4.at(j);
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
			mcMatchedT_SubldgJet = selJets_p4.at(j);
		      }
		  }
		
		// Is Jet2 closer in eta-phi AND has smaller pT difference? 
		
		else if (dR1 < dR1min)                //at least two matched jets                                                         
		  {
		    if  (dPtOverPt1 < twoSigmaDpt)
		      {
			dR1min  = dR1;
			dPtOverPt1min = dPtOverPt1;
			mcMatchedT_LdgJet = selJets_p4.at(j);
		      }
		  }
	      }
	  }//for(auto Bjet: selectedBJets)
      
	// Check if TOP is genuine                                                                                                       
	// if (doMatching)      
	
	bool isGenuine = false;
	isGenuine = (dR1min<= dRcut && dR2min <= dRcut && dRminB.at(i) <= dRcut);
	bool twoJMatched = isGenuine;  //at least two matched jets two top decay products (used for plots)
	if (!twoJMatched) twoJMatched = max(dR1min, dR2min) <= dRcut || max(dR1min, dRminB.at(i)) <= dRcut || max(dR2min, dRminB.at(i)) <= dRcut;
	bool isGenuineW = false;
	isGenuineW = (dR1min<= dRcut && dR2min <= dRcut);
	
	if (isGenuineW)
	  {
	    math::XYZTLorentzVector WldgJet_p4;
	    math::XYZTLorentzVector WsubldgJet_p4;
	    math::XYZTLorentzVector Wdijet_p4;
	    WldgJet_p4  = mcMatchedT_LdgJet;
	    WsubldgJet_p4 = mcMatchedT_SubldgJet;
	    Wdijet_p4 = WldgJet_p4 + WsubldgJet_p4;
	    double dR_Wj1j2 = ROOT::Math::VectorUtil::DeltaR(WldgJet_p4, WsubldgJet_p4);
	    h_GenuineWfromTop_Dijet_M      -> Fill(Wdijet_p4.M());
	    h_GenuineWfromTop_ldgJet_Pt    -> Fill(WldgJet_p4.pt());
	    h_GenuineWfromTop_subldgJet_Pt -> Fill(WsubldgJet_p4.pt());
	    h_GenuineWfromTop_Dijet_Pt     -> Fill(Wdijet_p4.pt());
	    h_GenuineWfromTop_Dijet_dR     -> Fill(dR_Wj1j2);
	  }


	h_TopQuarkPt_isGenuineTop -> Fill(isGenuine, bHt_TQuark_p4.pt());
	if (isGenuine)
	  {
	    if(0) std::cout << "=== bag3" << std::endl;
	    math::XYZTLorentzVector ldgJet_p4;
	    math::XYZTLorentzVector subldgJet_p4;
	    math::XYZTLorentzVector Tbjet_p4, dijet_p4, trijet_p4;
	    if(0) std::cout << "=== bag511" << std::endl;
	    ldgJet_p4    = mcMatchedT_LdgJet;
	    if(0) std::cout << "=== bag512" << std::endl;
	    subldgJet_p4 = mcMatchedT_SubldgJet;
	    if(0) std::cout << "=== bag513" << std::endl;
	    Tbjet_p4 = BJet;
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
	    cgenTop.increment();
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
      }//For-loop: All top-quarks
  }
  
  //chris
  //second reco top
  const double SigmaWMass = 9.72589;
  const double SigmaTMass = 17.9292;
  const double SigmaWMass2 = SigmaWMass*SigmaWMass;
  const double SigmaTMass2 = SigmaTMass*SigmaTMass;
  const double WMass = 80.379;
  const double TMass = 173.0;
  std::vector<math::XYZTLorentzVector> v_WJet, v_TopJet;
  for (size_t l=0; l<GenTops.size(); l++)
    {
      double X2min = 99999.9;
      math::XYZTLorentzVector WJet_p4, TopJet_p4;
      // find the best combination for top candidate
      for(size_t i=0; i < selJets_p4.size(); i++)
	{
	  math::XYZTLorentzVector jet1 = selJets_p4.at(i);
	  for(size_t j=0; j < selJets_p4.size(); j++)
	    {
	      if (j == i) continue;
	      math::XYZTLorentzVector jet2 = selJets_p4.at(j);
	      math::XYZTLorentzVector dijetn = jet1 + jet2;
	      double dMW  = dijetn.M() - WMass;
	      double dMW2 = dMW*dMW;
	      double dMW2SWM2 = dMW2/SigmaWMass2;
	      for(size_t k=0; k < selJets_p4.size(); k++)
		{
		  if(k == i || k == j) continue;
		  math::XYZTLorentzVector jet3    = selJets_p4.at(k);
		  math::XYZTLorentzVector trijetn = dijetn + jet3;
		  double dMT  = trijetn.M() - TMass;
		  double dMT2 = dMT*dMT;
		  double dMT2SWM2 = dMT2/SigmaTMass2;
		  double X2 = dMW2SWM2 + dMT2SWM2;
		  if (X2 > X2min) continue;
		  X2min = X2;
		  WJet_p4   = dijetn;
		  TopJet_p4 = trijetn; 
		}// for(size_t k=0; k < selJets_p4.size(); k++)
	    }// for(size_t j=0; j < selJets_p4.size(); j++)
	}// for(size_t i=0; i < selJets_p4.size(); i++)
      h_SecApro_TQuarkCan_M -> Fill(TopJet_p4.M());
      h_SecApro_WBosonCan_M -> Fill(WJet_p4.M());
      v_WJet.push_back(WJet_p4);
      v_TopJet.push_back(TopJet_p4);
    }//for (size_t l=0; l<GenTops.size(); l++)


  // top candidate dR with GenTop
  if (GenTops.size() == v_TopJet.size() && GenTops.size() == v_WJet.size())
    {
      for (size_t i=0; i<GenTops.size(); i++)
	{
	  double dRgt = ROOT::Math::VectorUtil::DeltaR(GenTops.at(i).p4(), v_TopJet.at(i));
	  h_SecApro_TQuark_Topjet_dR -> Fill(dRgt);
	  if (dRgt > 0.3) continue;
	  h_SecApro_TQuark_M -> Fill(v_TopJet.at(i).M());
	  h_SecApro_WBoson_M -> Fill(v_WJet.at(i).M());
	}
    }//if (GenTops.size() == v_TopJet.size())


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

vector<genParticle> RecGenHToHW::GetGenParticlesRec(const vector<genParticle> genParticles, const int pdgId)
{
  std::vector<genParticle> particles;

  // For-loop: All genParticles
  for (auto& p: genParticles){
    // Find last copy of a given particle
    if (!p.isLastCopy()) continue;
    
    // Consider only particles
    if (std::abs(p.pdgId()) != pdgId) continue;

    // Save this particle
    particles.push_back(p);
  }
  return particles;
}



const genParticle RecGenHToHW::GetLastCopy(const vector<genParticle> genParticles, const genParticle &p){

  int gen_pdgId = p.pdgId();

  for (size_t i=0; i<p.daughters().size(); i++){

    const genParticle genDau = genParticles[p.daughters().at(i)];
    int genDau_pdgId   = genDau.pdgId();

    if (gen_pdgId == genDau_pdgId)  return GetLastCopy(genParticles, genDau);
  }
  return p;
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

//chris
Bool_t RecGenHToHW::IshadronicTop(const vector<genParticle> genParticles, const vector<genParticle> genP)
{
  vector<genParticle> d1;
  bool ifoundW = false;
  for (auto& p: genP){
    for (size_t i=0; i<p.daughters().size(); i++)
      {
        int dau_index = p.daughters().at(i);
        genParticle dau = genParticles[dau_index];
	if (abs(dau.pdgId()) == 24)
          {
	    genParticle W = GetLastCopy(fEvent.genparticles().getGenParticles(), dau);
            d1.push_back(W);
            ifoundW = true;
	  }
      } // for (size_t i=0; i<p.daughters().size(); i++)
    if (!ifoundW) return false;
    bool hadrd2 = AreHadroDau(genParticles, d1);
    if (!hadrd2) return false;
  } //for (auto& p: genP)
  return true;
}



Bool_t RecGenHToHW::AreHadroDau(const vector<genParticle> genParticles, const vector<genParticle>  genP)
{
  MCTools mcTools(fEvent);
  for (auto& p: genP){
    for (size_t i=0; i<p.daughters().size(); i++)
      {
	int dau_index = p.daughters().at(i);
	genParticle dau = genParticles[dau_index];
	if (mcTools.IsLepton(dau.pdgId())) return false;
      }
  }
  return true;
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
