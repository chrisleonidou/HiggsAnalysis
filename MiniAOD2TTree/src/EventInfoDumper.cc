#include "HiggsAnalysis/MiniAOD2TTree/interface/EventInfoDumper.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"


EventInfoDumper::EventInfoDumper(edm::ConsumesCollector&& iConsumesCollector, const edm::ParameterSet& pset)
: puSummaryToken(iConsumesCollector.consumes<std::vector<PileupSummaryInfo>>(pset.getParameter<edm::InputTag>("PileupSummaryInfoSrc"))),
  lheToken(iConsumesCollector.consumes<LHEEventProduct>(pset.getUntrackedParameter<edm::InputTag>("LHESrc", edm::InputTag("")))),
  vertexToken(iConsumesCollector.consumes<edm::View<reco::Vertex>>(pset.getParameter<edm::InputTag>("OfflinePrimaryVertexSrc"))),
  secvertexToken(iConsumesCollector.consumes<edm::View<reco::VertexCompositePtrCandidate>>(pset.getParameter<edm::InputTag>("SlimmedSecondaryVerticesSrc")))
//  topPtToken(iConsumesCollector.consumes<double>(pset.getParameter<edm::InputTag>("TopPtProducer")))
{
    bookTopPt = false;
    if(pset.exists("TopPtProducer")){
        bookTopPt = true;
        topPtToken = iConsumesCollector.consumes<double>(pset.getParameter<edm::InputTag>("TopPtProducer"));
    }
    bookLumiScalers = pset.exists("LumiScalersSrc");
    if(bookLumiScalers) 
        instLumiToken = iConsumesCollector.consumes<std::vector<LumiScalers>>(pset.getParameter<edm::InputTag>("LumiScalersSrc"));
}

EventInfoDumper::~EventInfoDumper(){}

void EventInfoDumper::book(TTree* tree){
    tree->Branch("event",&event);
    tree->Branch("run",&run);     
    tree->Branch("lumi",&lumi);
    tree->Branch("prescale",&prescale);
    if(bookLumiScalers) tree->Branch("instLumi",&instLumi);
    tree->Branch("nPUvertices",&nPU);
    tree->Branch("NUP",&NUP);
    tree->Branch("nGoodOfflineVertices",&nGoodOfflinePV);
    tree->Branch("pvX",&pvX);
    tree->Branch("pvY",&pvY);
    tree->Branch("pvZ",&pvZ);
    tree->Branch("pvDistanceToNextVertex",&distanceToNextPV);
    tree->Branch("pvDistanceToClosestVertex",&distanceToClosestPV);
    tree->Branch("pvPtSumRatioToNext",&ptSumRatio);
    if(bookTopPt)
    tree->Branch("topPtWeight", &topPtWeight);

    // Secondary Vertex Collection
    tree->Branch("svPt", &svPt);
    tree->Branch("svEta", &svEta);
    tree->Branch("svPhi", &svPhi);
    tree->Branch("svMass", &svMass);
    tree->Branch("svNTks", &svNTks);
    tree->Branch("svChi2", &svChi2);
    tree->Branch("svNdof", &svNdof);
    tree->Branch("svDxy", &svDxy);
    tree->Branch("svDxyErr", &svDxyErr);
    tree->Branch("svD3d", &svD3d);
    tree->Branch("svD3dErr", &svD3dErr);
    tree->Branch("costhetasvpv", &costhetasvpv);

}

bool EventInfoDumper::fill(edm::Event& iEvent, const edm::EventSetup& iSetup){
    event = iEvent.id().event();
    run   = iEvent.run();
    lumi  = iEvent.luminosityBlock();
    instLumi = -1;
    prescale = 1.0;
    nPU = -1;
    NUP = -1;
    // Amount of PU
    edm::Handle<std::vector<PileupSummaryInfo> > hpileup;
    iEvent.getByToken(puSummaryToken, hpileup);
    if(hpileup.isValid()) { // protection for data
        for(std::vector<PileupSummaryInfo>::const_iterator iPV = hpileup->begin(); iPV != hpileup->end(); ++iPV) {
            if(iPV->getBunchCrossing() == 0) {
                nPU = iPV->getTrueNumInteractions();
                break;
            }
        }
    }

    if(bookLumiScalers){
      edm::Handle<LumiScalersCollection> hinstLumi;
      iEvent.getByToken(instLumiToken, hinstLumi);
      if(hinstLumi.isValid()){ // this is in AOD, not in MINIAOD. If needed, use two-files-solution
	instLumi = hinstLumi->begin()->instantLumi();
      }   
    }

    // number of jets for combining W+Jets/Z+jets inclusive with exclusive
    edm::Handle<LHEEventProduct> lheHandle;
    iEvent.getByToken(lheToken, lheHandle);
    if (lheHandle.isValid()) {
        // Store NUP = number of partons
        NUP = lheHandle->hepeup().NUP;
    }

    // PV
    nGoodOfflinePV = -1;
    pvX = 0;
    pvY = 0;
    pvZ = 0;
    ptSumRatio = -1.0;
    distanceToNextPV = -1.0;
    distanceToClosestPV = -1.0;
    edm::Handle<edm::View<reco::Vertex> > hoffvertex;
    if(iEvent.getByToken(vertexToken, hoffvertex)){
      nGoodOfflinePV = hoffvertex->size();
      // Multiply by 10 to get mm
      pvX = hoffvertex->at(0).x()*10.0;
      pvY = hoffvertex->at(0).y()*10.0;
      pvZ = hoffvertex->at(0).z()*10.0;
      if (nGoodOfflinePV > 1) {
        distanceToNextPV = std::fabs(hoffvertex->at(0).z() - hoffvertex->at(1).z());
        for (size_t i = 1; i < hoffvertex->size(); ++i) {
          float delta = std::fabs(hoffvertex->at(0).z() - hoffvertex->at(i).z());
          if (delta < distanceToClosestPV || distanceToClosestPV < 0.0) {
            distanceToClosestPV = delta;
          }
        }
        double ptSum0 = 0.0;
        for (std::vector<reco::TrackBaseRef>::const_iterator iter = hoffvertex->at(0).tracks_begin(); iter != hoffvertex->at(0).tracks_end(); iter++) {
          ptSum0 += hoffvertex->at(0).trackWeight(*iter) * (*iter)->pt()*(*iter)->pt();
        }
        double ptSum1 = 0.0;
        for (std::vector<reco::TrackBaseRef>::const_iterator iter = hoffvertex->at(1).tracks_begin(); iter != hoffvertex->at(1).tracks_end(); iter++) {
          ptSum1 += hoffvertex->at(1).trackWeight(*iter) * (*iter)->pt()*(*iter)->pt();
        }
        if (ptSum0 > 0.0) {
          ptSumRatio = ptSum1 / ptSum0;
        }
      }
      distanceToNextPV *= 10.0;
      distanceToClosestPV *= 10.0;
    }
    
    // Top pt
    topPtWeight = 1.0;
    edm::Handle<double> topPtHandle;
    if (bookTopPt && iEvent.getByToken(topPtToken, topPtHandle)) {
      topPtWeight = *(topPtHandle.product());
    }

    // Secondary Vertex collection
    reco::Vertex PV = hoffvertex->at(0); //getting PV from above.
    edm::Handle<edm::View<reco::VertexCompositePtrCandidate> > secvertex;
    if(iEvent.getByToken(secvertexToken, secvertex)){
      for(size_t isv = 0; isv<secvertex->size(); isv++ ){
        const reco::VertexCompositePtrCandidate &sv = (*secvertex)[isv];
        svPt         = sv.pt();
        svEta        = sv.eta();
        svPhi        = sv.phi();
        svMass       = sv.mass();
        svNTks       = sv.numberOfDaughters();
        svChi2       = sv.vertexChi2();
        svNdof       = sv.vertexNdof();
        svDxy        = vertexDxy(sv,PV).value();
        svDxyErr     = vertexDxy(sv,PV).error();
        svD3d        = vertexD3d(sv,PV).value();
        svD3dErr     = vertexD3d(sv,PV).error();
        costhetasvpv = vertexDdotP(sv, PV);
      }//eof: loop on sec vertex
    }//eof: sec vertex handle
    
    return filter();
}

bool EventInfoDumper::filter(){
    return true;
}

void EventInfoDumper::reset(){
}

Measurement1D EventInfoDumper::vertexD3d(const reco::VertexCompositePtrCandidate &svcand,
                                         const reco::Vertex &pv) const{
  VertexDistance3D dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

Measurement1D EventInfoDumper::vertexDxy(const reco::VertexCompositePtrCandidate &svcand,
                                         const reco::Vertex &pv) const{
  VertexDistanceXY dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

float EventInfoDumper::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const
                                   reco::Vertex &pv) const{
  reco::Candidate::Vector p = sv.momentum();
  reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
  return p.Unit().Dot(d.Unit());
}
