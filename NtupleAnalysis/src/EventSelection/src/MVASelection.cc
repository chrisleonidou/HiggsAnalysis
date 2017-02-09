#include "EventSelection/interface/MVASelection.h"

#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"


//TMVA Stuff

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"


MVASelection::Data::Data() 
: bPassedSelection(false) { }

MVASelection::Data::~Data(){ }

MVASelection::MVASelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& prefix, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, prefix+postfix),
  cPassedMVASelection(fEventCounter.addCounter("passed MVA selection "+prefix+" ("+postfix+")")),
  cSubAll(fEventCounter.addSubCounter("MVA selection "+prefix+" ("+postfix+")", "All events"))
{
  initialize(config, postfix);
}


MVASelection::MVASelection(const ParameterSet& config, const std::string& postfix)
: BaseSelection(),
  cPassedMVASelection(fEventCounter.addCounter("passed MVA selection ("+postfix+")")),
  cSubAll(fEventCounter.addSubCounter("MVA selection ("+postfix+")", "All events"))
{
  initialize(config, postfix);
  bookHistograms(new TDirectory());
}

MVASelection::~MVASelection() {
  delete hMVAValueAll;
  delete reader;
}

void MVASelection::initialize(const ParameterSet& config, const std::string& postfix) {
  reader = new TMVA::Reader( "!Color:Silent" );
/*
  Float_t R_bb, MET, TransMass, TransMass_jj, TransMass_muEt;
  Float_t tjDist, jjDist, mujDist, ejDist;
  Float_t SelectedTau_pt, SelectedTau_eta, SelectedTau_phi;
  Float_t Muon_pt, Muon_eta, Muon_phi;
  Float_t Electron_pt, Electron_eta, Electron_phi;
  Float_t Jet1_pt, Jet1_phi, Jet1_eta, Jet1_IncSecVer_Btag, Jet1_MVA_Btag;
  Float_t Jet2_pt, Jet2_phi, Jet2_eta, Jet2_IncSecVer_Btag, Jet2_MVA_Btag;
  Float_t Jet3_pt, Jet3_phi, Jet3_eta, Jet3_IncSecVer_Btag, Jet3_MVA_Btag;
  Float_t nTaus, nJets, nMuons, nElectrons;
*/
  reader->AddVariable("R_bb",&R_bb);
  reader->AddVariable("MET",&MET);
  reader->AddVariable("TransMass",&TransMass);
  reader->AddVariable("TransMass_jj",&TransMass_jj);
  reader->AddVariable("TransMass_muEt",&TransMass_muEt);

  reader->AddVariable("tj1Dist",&tj1Dist);
  reader->AddVariable("tj2Dist",&tj2Dist);
  reader->AddVariable("tj3Dist",&tj3Dist);

  reader->AddVariable("j1j2Dist",&j1j2Dist);
  reader->AddVariable("j1j3Dist",&j1j3Dist);
  reader->AddVariable("j2j3Dist",&j2j3Dist);

  reader->AddVariable("muj1Dist",&muj1Dist);
  reader->AddVariable("muj2Dist",&muj2Dist);
  reader->AddVariable("muj3Dist",&muj3Dist);

  reader->AddVariable("ej1Dist",&ej1Dist);
  reader->AddVariable("ej2Dist",&ej2Dist);
  reader->AddVariable("ej3Dist",&ej3Dist);


  reader->AddVariable("tetaj1eta",&tetaj1eta);
  reader->AddVariable("tetaj2eta",&tetaj2eta);
  reader->AddVariable("tetaj3eta",&tetaj3eta);

  reader->AddVariable("j1etaj2eta",&j1etaj2eta);
  reader->AddVariable("j1etaj3eta",&j1etaj3eta);
  reader->AddVariable("j2etaj3eta",&j2etaj3eta);


//  reader->AddVariable("nJets",&nJets);

  reader->AddVariable("SelectedTau_pt",&SelectedTau_pt);
  reader->AddVariable("SelectedTau_phi",&SelectedTau_phi);
  reader->AddVariable("SelectedTau_eta",&SelectedTau_eta);
  reader->AddVariable("nTaus",&nTaus);

  reader->AddVariable("Jet1_pt",&Jet1_pt);
  reader->AddVariable("Jet1_phi",&Jet1_phi);
  reader->AddVariable("Jet1_eta",&Jet1_eta);
  reader->AddVariable("Jet1_IncSecVer_Btag",&Jet1_IncSecVer_Btag);
  reader->AddVariable("Jet1_MVA_Btag",&Jet1_MVA_Btag);

  reader->AddVariable("Jet2_pt",&Jet2_pt);
  reader->AddVariable("Jet2_phi",&Jet2_phi);
  reader->AddVariable("Jet2_eta",&Jet2_eta);
  reader->AddVariable("Jet2_IncSecVer_Btag",&Jet2_IncSecVer_Btag);
  reader->AddVariable("Jet2_MVA_Btag",&Jet2_MVA_Btag);

  reader->AddVariable("Jet3_pt",&Jet3_pt);
  reader->AddVariable("Jet3_phi",&Jet3_phi);
  reader->AddVariable("Jet3_eta",&Jet3_eta);
  reader->AddVariable("Jet3_IncSecVer_Btag",&Jet3_IncSecVer_Btag);
  reader->AddVariable("Jet3_MVA_Btag",&Jet3_MVA_Btag);

  reader->AddVariable("Electron_pt",&Electron_pt);
  reader->AddVariable("Electron_phi",&Electron_phi);
  reader->AddVariable("Electron_eta",&Electron_eta);
  reader->AddVariable("nElectrons",&nElectrons);

  reader->AddVariable("Muon_pt",&Muon_pt);
  reader->AddVariable("Muon_phi",&Muon_phi);
  reader->AddVariable("Muon_eta",&Muon_eta);
  reader->AddVariable("nMuons",&nMuons);

//  reader->BookMVA( "BDTG method","MyClassification_BDTG.weights.xml");
  reader->BookMVA("DNN method","MyClassification_DNN.weights.xml");

}

void MVASelection::bookHistograms(TDirectory* dir) {
  TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kDebug, dir, "MVASelection_"+sPostfix);
  hMVAValueAll = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "MVAValueAll", "MVAValue, all; BDT", 40, -1.0, 1.0);
}

MVASelection::Data MVASelection::silentAnalyze(const Event& event, TMVA::Reader& reader) {
  ensureSilentAnalyzeAllowed(event.eventID());
  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(event);
  enableHistogramsAndCounters();
  return myData;
}

MVASelection::Data MVASelection::analyze(const Event& event, TMVA::Reader& reader) {
  ensureAnalyzeAllowed(event.eventID());
  MVASelection::Data data = privateAnalyze(event);
  if (fCommonPlotsIsEnabled())
    fCommonPlots->fillControlPlotsAtMVASelection(event, data);
  return data;
}

MVASelection::Data MVASelection::privateAnalyze(const Event& event) {
  cSubAll.increment();
  bool passedMVA = false;
  Data output;

  SelectedTau_pt=event.taus()[0].pt();
  SelectedTau_eta=event.taus()[0].eta();
  SelectedTau_phi=event.taus()[0].phi();
  nTaus=event.taus().size();

  nMuons=event.muons().size();
  if(nMuons>0){
    Muon_pt=event.muons()[0].pt();
    Muon_eta=event.muons()[0].eta();
    Muon_phi=event.muons()[0].phi();
  }else{
    Muon_pt=0;
    Muon_eta=0;
    Muon_phi=0;
  }

  nElectrons=event.electrons().size();
  if(nElectrons>0){
    Electron_pt=event.electrons()[0].pt();
    Electron_eta=event.electrons()[0].eta();
    Electron_phi=event.electrons()[0].phi();
  }else{
    Electron_pt=0;
    Electron_eta=0;
    Electron_phi=0;
  }
  Jet1_pt=event.jets()[0].pt();
  Jet1_eta=event.jets()[0].eta();
  Jet1_phi=event.jets()[0].phi();
  Jet1_IncSecVer_Btag=event.jets()[0].pfCombinedInclusiveSecondaryVertexV2BJetTags();
  Jet1_MVA_Btag=event.jets()[0].pfCombinedMVAV2BJetTags();

  Jet2_pt=event.jets()[1].pt();
  Jet2_eta=event.jets()[1].eta();
  Jet2_phi=event.jets()[1].phi();
  Jet2_IncSecVer_Btag=event.jets()[1].pfCombinedInclusiveSecondaryVertexV2BJetTags();
  Jet2_MVA_Btag=event.jets()[1].pfCombinedMVAV2BJetTags();

  Jet3_pt=event.jets()[2].pt();
  Jet3_eta=event.jets()[2].eta();
  Jet3_phi=event.jets()[2].phi();
  Jet3_IncSecVer_Btag=event.jets()[2].pfCombinedInclusiveSecondaryVertexV2BJetTags();
  Jet3_MVA_Btag=event.jets()[2].pfCombinedMVAV2BJetTags();

  nJets=event.jets().size();

  Float_t met_x,met_y;
  met_x=event.met().x();
  met_y=event.met().y();
  MET=sqrt(pow(met_x,2)+pow(met_y,2));

  TransMass=sqrt(2*SelectedTau_pt*MET*(1-cos(atan2(met_y,met_x)+SelectedTau_phi)));
  TransMass_jj=sqrt(2*Jet1_pt*Jet2_pt*(1-cos(Jet2_phi+Jet1_phi)));
  TransMass_muEt=sqrt(2*Muon_pt*MET*(1-cos(atan2(met_y,met_x)+Muon_phi)));

  tj1Dist=sqrt(pow(SelectedTau_phi-Jet1_phi,2)+pow(SelectedTau_eta-Jet1_eta,2));
  tj2Dist=sqrt(pow(SelectedTau_phi-Jet2_phi,2)+pow(SelectedTau_eta-Jet2_eta,2));
  tj3Dist=sqrt(pow(SelectedTau_phi-Jet3_phi,2)+pow(SelectedTau_eta-Jet3_eta,2));

  j1j2Dist=sqrt(pow(Jet2_phi-Jet1_phi,2)+pow(Jet2_eta-Jet1_eta,2));
  j1j3Dist=sqrt(pow(Jet3_phi-Jet1_phi,2)+pow(Jet3_eta-Jet1_eta,2));  
  j2j3Dist=sqrt(pow(Jet2_phi-Jet3_phi,2)+pow(Jet2_eta-Jet3_eta,2));

  tetaj1eta=SelectedTau_eta*Jet1_eta;
  tetaj1eta=SelectedTau_eta*Jet2_eta;
  tetaj1eta=SelectedTau_eta*Jet3_eta;

  j1etaj2eta=Jet1_eta*Jet2_eta;
  j1etaj3eta=Jet1_eta*Jet3_eta;
  j2etaj3eta=Jet2_eta*Jet3_eta;


  if(nMuons>0){
    muj1Dist=sqrt(pow(Muon_phi-Jet1_phi,2)+pow(Muon_eta-Jet1_eta,2));
    muj2Dist=sqrt(pow(Muon_phi-Jet2_phi,2)+pow(Muon_eta-Jet2_eta,2));
    muj3Dist=sqrt(pow(Muon_phi-Jet3_phi,2)+pow(Muon_eta-Jet3_eta,2));
  }else{
    muj1Dist=0;
    muj2Dist=0;
    muj3Dist=0;
  }
  if(nElectrons>0){
    ej1Dist=sqrt(pow(Electron_phi-Jet1_phi,2)+pow(Electron_eta-Jet1_eta,2));
    ej2Dist=sqrt(pow(Electron_phi-Jet2_phi,2)+pow(Electron_eta-Jet2_eta,2));
    ej3Dist=sqrt(pow(Electron_phi-Jet3_phi,2)+pow(Electron_eta-Jet3_eta,2));
  }else{
    ej1Dist=0;
    ej2Dist=0;
    ej3Dist=0;
  }

//  output.setValue(reader->EvaluateMVA("BDTG method"));
  output.setValue(reader->EvaluateMVA("DNN method"));
  hMVAValueAll->Fill(output.mvaValue());
  passedMVA=(output.mvaValue()>0.6);
  if(passedMVA){
    output.setTrue();
    cPassedMVASelection.increment();
  }
  return output;
}
