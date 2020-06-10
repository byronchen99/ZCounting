
// CMSSW includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "ZCounting/TnPPairTreeProducer/plugins/MuonTreeProducer.h"

#include <TLorentzVector.h>

//
// constructors and destructor
//
MuonTreeProducer::MuonTreeProducer(const edm::ParameterSet& iConfig):
    triggerResultsInputTag_(iConfig.getParameter<edm::InputTag>("TriggerResults")),
    fTrackName(iConfig.getUntrackedParameter<std::string>("edmTrackName", "generalTracks")),
    fMuonName(iConfig.getUntrackedParameter<std::string>("edmMuonName", "muons")),
    fMuonHLTNames (iConfig.getParameter<std::vector<std::string>>("MuonTriggerNames")),
    PtCutL_(iConfig.getParameter<double>("PtCutL")),
    EtaCutL_(iConfig.getParameter<double>("EtaCutL")),
    fPVName(iConfig.getUntrackedParameter<std::string>("edmPVName", "offlinePrimaryVertices")),
    VtxNTracksFitCut_(iConfig.getParameter<double>("VtxNTracksFitMin")),
    VtxNdofCut_(iConfig.getParameter<double>("VtxNdofMin")),
    VtxAbsZCut_(iConfig.getParameter<double>("VtxAbsZMax")),
    VtxRhoCut_(iConfig.getParameter<double>("VtxRhoMax"))
{
    triggers = new triggertool();
    triggers->setTriggerResultsToken(consumes<edm::TriggerResults>(triggerResultsInputTag_));
    triggers->setTriggerEventToken(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerEvent")));

    edm::LogVerbatim("MuonTreeProducer") << "getInput: set trigger names";
    for(unsigned int i = 0; i < fMuonHLTNames.size(); ++i) {
        triggers->addTriggerRecord(fMuonHLTNames.at(i));
    }

    // additional trigger to be checked
    triggers->addTriggerRecord("HLT_IsoMu27_v*");

    fPVName_token = consumes<std::vector<reco::Vertex>>(fPVName);
    fMuonName_token = consumes<std::vector<reco::Muon>>(fMuonName);
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::LogInfo("MuonTreeProducer")<<"analyze";

    this->clearVariables();

    triggers->readEvent(iEvent);

    if(!isMuonTrigger()){
        edm::LogVerbatim("MuonTreeProducer") << "MuonTreeProducer::analyze - event did not pass any muon trigger";
        return;
    }

    run_ = iEvent.id().run();
    ls_ = iEvent.id().luminosityBlock();
    eventNumber_ = iEvent.id().event();

    //-------------------------------
    //--- Primary vertices
    //-------------------------------
    edm::Handle<reco::VertexCollection> hVertexProduct;
    iEvent.getByToken(fPVName_token, hVertexProduct);
    if (!hVertexProduct.isValid()) {
        edm::LogWarning("MuonTreeProducer") << "MuonTreeProducer::analyze - no valid pv product found" << std::endl;
        return;
    }

    const reco::VertexCollection* pvCol = hVertexProduct.product();
    const reco::Vertex* pv = &(*pvCol->begin());

    for (auto const& itVtx : *hVertexProduct) {
        if (itVtx.isFake())
            continue;
        if (itVtx.tracksSize() < VtxNTracksFitCut_)
            continue;
        if (itVtx.ndof() < VtxNdofCut_)
            continue;
        if (fabs(itVtx.z()) > VtxAbsZCut_)
            continue;
        if (itVtx.position().Rho() > VtxRhoCut_)
            continue;
        if (nPV_ == 0) {
            pv = &itVtx;
        }
        nPV_++;
    }

    //-------------------------------
    //--- Muons and Tracks
    //-------------------------------
    edm::Handle<std::vector<reco::Muon>> hMuonProduct;
    iEvent.getByToken(fMuonName_token, hMuonProduct);
    if (!hMuonProduct.isValid()){
        edm::LogWarning("MuonTreeProducer") << "MuonTreeProducer::analyze - no valid muon product found" << std::endl;
        return;
    }

    // Tag loop
    for (auto const& itMu : *hMuonProduct) {

        this->clearMuonVariables();

        pt_ = itMu.muonBestTrack()->pt();
        eta_ = itMu.muonBestTrack()->eta();
        phi_ = itMu.muonBestTrack()->phi();
        q_ = itMu.muonBestTrack()->charge();

        if (pt_ < PtCutL_)
            continue;
        if (fabs(eta_) > EtaCutL_)
            continue;
        if (!passBaseline(itMu))
            continue;

        is_HLT_IsoMu27_ = triggers->passObj("HLT_IsoMu27_v*", eta_, phi_);
        is_HLT_Mu17_ = triggers->passObj("HLT_Mu17_v*", eta_, phi_);



        if (!isMuonTriggerObj(eta_, phi_))
            continue;

        pfIso_ = getPFIso(itMu);
        tkIso_ = getTkIso(itMu);
        dxy_ = getDxy(itMu, *pv);
        dz_ = getDz(itMu, *pv);

        tree_->Fill();

    }

}

// ------------ method called once each run just before starting event loop  ------------
void MuonTreeProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    std::cout<<"now at "<<iRun.id()<<std::endl;

    bool hltChanged(true);
    if (hltConfigProvider_.init(iRun, iSetup, triggerResultsInputTag_.process(), hltChanged)) {
        edm::LogInfo("MuonTreeProducer")<< "[TriggerObjMatchValueMapsProducer::beginRun] HLTConfigProvider initialized [processName() = \""
            << hltConfigProvider_.processName() << "\", tableName() = \"" << hltConfigProvider_.tableName()
            << "\", size() = " << hltConfigProvider_.size() << "]";
    } else {
        edm::LogError("MuonTreeProducer") << "Initialization of HLTConfigProvider failed for Run=" << iRun.id() << " (process=\""
        << triggerResultsInputTag_.process() << "\") -> plugin will not produce outputs for this Run";
        return;
    }

    triggers->initHLTObjects(hltConfigProvider_);

}


// ------------ method called once each job just before beginRun ------------
void
MuonTreeProducer::beginJob()
{
    if( !fs ){
        edm::LogError("MuonTreeProducer") << "TFile Service is not registered in cfg file";
        return;
    }

    tree_=(fs->make<TTree>("tree" ,"tree" ));

    tree_->Branch("nPV", &nPV_,"nPV_/i");
    tree_->Branch("run", &run_,"run_/i");
    tree_->Branch("ls", &ls_,"ls_/i");
    tree_->Branch("eventNumber", &eventNumber_, "eventNumber_/i");

    tree_->Branch("pt", &pt_,"pt_/f");
    tree_->Branch("eta", &eta_,"eta_/f");
    tree_->Branch("phi", &phi_,"phi_/f");
    tree_->Branch("q", &q_,"q_/f");
    tree_->Branch("pfIso", &pfIso_,"pfIso_/f");
    tree_->Branch("tkIso", &tkIso_,"tkIso_/f");
    tree_->Branch("dxy", &dxy_,"dxy_/f");
    tree_->Branch("dz", &dz_,"dz_/f");

    tree_->Branch("is_HLT_Mu17", &is_HLT_Mu17_,"is_HLT_Mu17_/b");
    tree_->Branch("is_HLT_IsoMu27", &is_HLT_IsoMu27_,"is_HLT_IsoMu27_/b");


}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonTreeProducer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//--------------------------------------------------------------------------------------------------
void MuonTreeProducer::clearVariables(){
    edm::LogInfo("MuonTreeProducer")<<"clearVariables()";

    nPV_ = 0;
    run_ = 0;
    ls_ = 0;
    eventNumber_ = 0;

}

//--------------------------------------------------------------------------------------------------
void MuonTreeProducer::clearMuonVariables(){
    edm::LogInfo("MuonTreeProducer")<<"clearTagVariables()";

    pt_ = 0.;
    eta_ = 0.;
    phi_ = 0.;
    q_ = 0;
    pfIso_ = 0.;
    tkIso_ = 0.;
    dxy_ = 0.;
    dz_ = 0.;
    is_HLT_Mu17_ = false;
    is_HLT_IsoMu27_ = false;
}


//--------------------------------------------------------------------------------------------------
bool MuonTreeProducer::isMuonTrigger() {
    if (triggers->pass(fMuonHLTNames))
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
bool MuonTreeProducer::isMuonTriggerObj(const double eta, const double phi) {
    if (triggers->passObj(fMuonHLTNames, eta, phi))
        return true;
    return false;
}
