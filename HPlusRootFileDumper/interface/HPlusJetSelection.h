#ifndef HPLUSANALYSISHPLUSJETSELECTION_H
#define HPLUSANALYSISHPLUSJETSELECTION_H

#include "HiggsAnalysis/HPlusRootFileDumper/interface/HPlusAnalysisBase.h"
#include "HiggsAnalysis/HPlusRootFileDumper/interface/HPlusSelectionBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "TH1F.h"


namespace HPlusAnalysis {
  
  //*************************************************************************************************************
  //   author     : Alexandros Attikis                                                                          *
  //   date       : 07 June 2010                                                                                *
  //   email      : attikis@cern.ch                                                                             *
  //   Description: Class for jet selection:                                                                    *
  //                1) Loop over all jets,                                                                      *
  //                2) Apply jet energy corrections,                                                            *
  //                3) Clean jet collection with specific cuts on Et, |Eta| and EMFraction,                     *
  //                4) Sort jets in jet Et (descending order),                                                  *
  //                5) Save to edm ntuple: a) the jets surviving selection, b) the number of jets surviving     *
  //                selection, c) Validation histograms including Jet Et, Eta of: "leading-Jet",                *
  // "2nd-Jet", "3rd-Jet", "4th-Jet", d) Informative counters.                                                  *
  //                6) Return boolean decision if Event satisfies "jet-selection criteria",                     *
  //                   i.e. At least N jets are left after the steps 1) -> 4).                                  *
  //*************************************************************************************************************
  
  class HPlusJetSelection : public HPlusAnalysis::HPlusAnalysisBase, public HPlusAnalysis::HPlusSelectionBase {
  public:
  /// Default constructor
  HPlusJetSelection(const edm::ParameterSet& iConfig);
  /// Default destructor
  ~HPlusJetSelection();
  
  void beginJob();
  void endJob();
  bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup);
  std::vector<reco::CaloJet> sortCaloJets(std::vector<reco::CaloJet> CaloJetsToSort, const size_t caloJetSize);
  std::vector<reco::CaloJet> filterCaloJets(std::vector<reco::CaloJet> caloJetsSorted, const size_t caloJetSize);
  std::vector<reco::CaloJet> eraseVectorElement( std::vector<reco::CaloJet> myJetVector, reco::CaloJet test);

  private:
  /// Name of jet collection
  edm::InputTag fJetCollectionName;
  /// Cut values (their values taken from py cfg file)
  double fCutMinNJets;
  double fCutMinJetEt;
  double fCutMaxAbsJetEta;
  double fCutMaxEMFraction;
  // ****************************//

  // Counters
  int fCounterTest;
  int fCounterJetsPriorSelection;
  int fCounterError;
  int fCounterJetsPostSelection;
  int fCounterJetCollectionHandleEmpty;
  // Histograms
  TH1F* hLeadJetMaxEt; /// Maximum muon pT in event (-1 if no muon in event)
  
  };
  
}

#endif
