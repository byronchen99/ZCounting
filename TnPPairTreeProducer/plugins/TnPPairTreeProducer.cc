
// CMSSW includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "ZCounting/TnPPairTreeProducer/plugins/TnPPairTreeProducer.h"

#include <TLorentzVector.h>

//
// constructors and destructor
//
TnPPairTreeProducer::TnPPairTreeProducer(const edm::ParameterSet& iConfig):
    triggerResultsInputTag_(iConfig.getParameter<edm::InputTag>("TriggerResults"))
{
    fTrackName = iConfig.getUntrackedParameter<std::string>("edmTrackName", "generalTracks");
    fMuonName = iConfig.getUntrackedParameter<std::string>("edmMuonName", "muons");
    fMuonHLTNames  = iConfig.getParameter<std::vector<std::string>>("MuonTriggerNames");
    fMuonHLTDRMAX = iConfig.getParameter<double>("MuonTriggerDRMAX");
    fPVName = iConfig.getUntrackedParameter<std::string>("edmPVName", "offlinePrimaryVertices");

    VtxNTracksFitCut_ = iConfig.getUntrackedParameter<double>("VtxNTracksFitMin");
    VtxNdofCut_ = iConfig.getUntrackedParameter<double>("VtxNdofMin");
    VtxAbsZCut_ = iConfig.getUntrackedParameter<double>("VtxAbsZMax");
    VtxRhoCut_ = iConfig.getUntrackedParameter<double>("VtxRhoMax");

    // Muon-specific Cuts
    IDTypestr_ = iConfig.getUntrackedParameter<std::string>("IDType");
    IsoTypestr_ = iConfig.getUntrackedParameter<std::string>("IsoType");
    IsoCut_ = iConfig.getUntrackedParameter<double>("IsoCut");
    DxyCut_ = iConfig.getUntrackedParameter<double>("DxyCut");
    DzCut_ = iConfig.getUntrackedParameter<double>("DzCut");

    if (IDTypestr_ == "Loose")
        IDType_ = LooseID;
    else if (IDTypestr_ == "Medium")
        IDType_ = MediumID;
    else if (IDTypestr_ == "Tight")
        IDType_ = TightID;
    else if (IDTypestr_ == "Custom")
        IDType_ = CustomID;
    else if (IDTypestr_ == "None")
        IDType_ = NoneID;
    else
        edm::LogError("TnPPairTreeProducer") << "invalid IDType (choose one: \"Custom\", \"Loose\", \"Medium\", \"Tight\", \"None\" ) ";

    if (IsoTypestr_ == "Tracker-based")
        IsoType_ = TrackerIso;
    else if (IsoTypestr_ == "PF-based")
        IsoType_ = PFIso;
    else if (IsoTypestr_ == "None")
        IsoType_ = NoneIso;
    else
        edm::LogError("TnPPairTreeProducer") << "invalid IsoType (choose one: \"Tracker-based\", \"PF-based\", \"None\" ) ";


    PtCutL1_ = iConfig.getUntrackedParameter<double>("PtCutL1");
    EtaCutL1_ = iConfig.getUntrackedParameter<double>("EtaCutL1");
    PtCutL2_ = iConfig.getUntrackedParameter<double>("PtCutL2");
    EtaCutL2_ = iConfig.getUntrackedParameter<double>("EtaCutL2");

    MassMin_ = iConfig.getUntrackedParameter<double>("MassMin");
    MassMax_ = iConfig.getUntrackedParameter<double>("MassMax");

    triggers = new triggertool();
    triggers->setTriggerResultsToken(consumes<edm::TriggerResults>(triggerResultsInputTag_));
    triggers->setTriggerEventToken(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerEvent")));
    triggers->setDRMAX(fMuonHLTDRMAX);

    edm::LogVerbatim("TnPPairTreeProducer") << "getInput: set trigger names";
    for(unsigned int i = 0; i < fMuonHLTNames.size(); ++i) {
        triggers->addTriggerRecord(fMuonHLTNames.at(i));
    }
    //triggers->addTriggerRecord("HLT_IsoMu27_v*");

    fPVName_token = consumes<std::vector<reco::Vertex>>(fPVName);
    fMuonName_token = consumes<std::vector<reco::Muon>>(fMuonName);
    fTrackName_token = consumes<std::vector<reco::Track>>(fTrackName);

}


TnPPairTreeProducer::~TnPPairTreeProducer()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TnPPairTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::LogVerbatim("TnPPairTreeProducer")<<"analyze event "<<iEvent.id().event();

    this->clearVariables();

    triggers->readEvent(iEvent);

    if(!isMuonTrigger()){
        edm::LogVerbatim("TnPPairTreeProducer") << "TnPPairTreeProducer::analyze - event did not pass any muon trigger";
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
        edm::LogWarning("TnPPairTreeProducer") << "TnPPairTreeProducer::analyze - no valid pv product found" << std::endl;
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
        edm::LogWarning("TnPPairTreeProducer") << "TnPPairTreeProducer::analyze - no valid muon product found" << std::endl;
        return;
    }

    edm::Handle<reco::TrackCollection> hTrackProduct;
    iEvent.getByToken(fTrackName_token, hTrackProduct);
    if (!hTrackProduct.isValid()){
        edm::LogWarning("TnPPairTreeProducer") << "TnPPairTreeProducer::analyze - no valid track product found" << std::endl;
        return;
    }

    TLorentzVector vTag(0., 0., 0., 0.);
    TLorentzVector vProbe(0., 0., 0., 0.);

    // Tag loop
    for (auto const& itMu1 : *hMuonProduct) {

        this->clearTagVariables();

        pt1_ = itMu1.muonBestTrack()->pt();
        eta1_ = itMu1.muonBestTrack()->eta();
        phi1_ = itMu1.muonBestTrack()->phi();
        q1_ = itMu1.muonBestTrack()->charge();

        if (pt1_ < PtCutL1_)
            continue;
        if (fabs(eta1_) > EtaCutL1_)
            continue;
        if (!(passMuonID(itMu1, *pv) && passMuonIso(itMu1)))
            continue;
        if (!isMuonTriggerObj(eta1_, phi1_))
            continue;

        pfIso1_ = getPFIso(itMu1);
        tkIso1_ = getTkIso(itMu1);
        dxy1_ = getDxy(itMu1, *pv);
        dz1_ = getDz(itMu1, *pv);
        // is1IsoMu27_ = triggers->passObj("HLT_IsoMu27_v*", eta1_, phi1_);

        vTag.SetPtEtaPhiM(pt1_, eta1_, phi1_, MUON_MASS);

        // Probe loop
        for (auto const& itMu2 : *hMuonProduct) {
            if (&itMu2 == &itMu1)
                continue;

            this->clearProbeVariables();

            pt2_ = itMu2.muonBestTrack()->pt();
            eta2_ = itMu2.muonBestTrack()->eta();
            phi2_ = itMu2.muonBestTrack()->phi();
            q2_ = itMu2.muonBestTrack()->charge();

            if (pt2_ < PtCutL2_)
                continue;
            if (fabs(eta2_) > EtaCutL2_)
                continue;

            vProbe.SetPtEtaPhiM(pt2_, eta2_, phi2_, MUON_MASS);
            TLorentzVector vDilep = vTag + vProbe;
            dilepMass_ = vDilep.M();
            if ((dilepMass_ < MassMin_) || (dilepMass_ > MassMax_))
                continue;
            if (passMuonID(itMu2, *pv) && passMuonIso(itMu2)) {
                if (isMuonTriggerObj(eta2_, phi2_)) {
                    // category 2HLT: both muons passing trigger requirements
                    if (&itMu1 > &itMu2)
                        continue;  // make sure we don't double count MuMu2HLT category
                    is2HLT_ = true;
                }
                else // category 1HLT: probe passing selection but not trigger
                    isSel_ = true;
            }
            else if (itMu2.isGlobalMuon())
                isGlo_ = true;
            else if (itMu2.isStandAloneMuon())
                isSta_ = true;
            else if (itMu2.innerTrack()->hitPattern().trackerLayersWithMeasurement() >= 6
                      && itMu2.innerTrack()->hitPattern().numberOfValidPixelHits() >= 1)
                isTrk_ = true;
            else
                continue;

            pfIso2_ = getPFIso(itMu2);
            tkIso2_ = getTkIso(itMu2);
            dxy2_ = getDxy(itMu2, *pv);
            dz2_ = getDz(itMu2, *pv);
            // is2IsoMu27_ = triggers->passObj("HLT_IsoMu27_v*", eta2_, phi2_);

            tree_->Fill();
        }

        // Probe loop over tracks, only for standalone efficiency calculation
        for (auto const& itTrk : *hTrackProduct) {
            // Check track is not a muon
            bool isMuon = false;
            for (auto const& itMu : *hMuonProduct) {
                if (itMu.innerTrack().isNonnull() && itMu.innerTrack().get() == &itTrk) {
                    isMuon = true;
                    break;
                }
            }
            if (isMuon)
                continue;

            this->clearProbeVariables();

            pt2_ = itTrk.pt();
            eta2_ = itTrk.eta();
            phi2_ = itTrk.phi();
            q2_ = itTrk.charge();
            // is2IsoMu27_ = triggers->passObj("HLT_IsoMu27_v*", eta2_, phi2_);

            // Probe selection:  kinematic cuts and opposite charge requirement
            if (pt2_ < PtCutL2_)
                continue;
            if (fabs(eta2_) > EtaCutL2_)
                continue;

            vProbe.SetPtEtaPhiM(pt2_, eta2_, phi2_, MUON_MASS);

            TLorentzVector vDilep = vTag + vProbe;
            dilepMass_ = vDilep.M();

            if ((dilepMass_ < MassMin_) || (dilepMass_ > MassMax_))
                continue;

            if (itTrk.hitPattern().trackerLayersWithMeasurement() >= 6 && itTrk.hitPattern().numberOfValidPixelHits() >= 1) {
                isTrk_ = true;
                tree_->Fill();
            }

        }    //End of probe loop over tracks

    }

}


// ------------ method called once each run just before starting event loop  ------------
void TnPPairTreeProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    edm::LogVerbatim("TnPPairTreeProducer") << "now at "<<iRun.id();

    bool hltChanged(true);
    if (hltConfigProvider_.init(iRun, iSetup, triggerResultsInputTag_.process(), hltChanged)) {
        edm::LogVerbatim("TnPPairTreeProducer")<<" [TriggerObjMatchValueMapsProducer::beginRun] HLTConfigProvider initialized [processName() = \""
            << hltConfigProvider_.processName() << "\", tableName() = \"" << hltConfigProvider_.tableName()
            << "\", size() = " << hltConfigProvider_.size() << "]";
    } else {
        edm::LogError("TnPPairTreeProducer") << "Initialization of HLTConfigProvider failed for Run=" << iRun.id() << " (process=\""
        << triggerResultsInputTag_.process() << "\") -> plugin will not produce outputs for this Run";
        return;
    }
    edm::LogVerbatim("TnPPairTreeProducer") << "hlt: "<<hltChanged;

    triggers->initHLTObjects(hltConfigProvider_);

}

// ------------ method called once each job just before beginRun ------------
void
TnPPairTreeProducer::beginJob()
{
    if( !fs ){
        edm::LogError("TnPPairTreeProducer") << "TFile Service is not registered in cfg file";
        return;
    }

    tree_=(fs->make<TTree>("tree" ,"tree" ));

    tree_->Branch("nPV", &nPV_,"nPV_/i");
    tree_->Branch("run", &run_,"run_/i");
    tree_->Branch("ls", &ls_,"ls_/i");
    tree_->Branch("eventNumber", &eventNumber_, "eventNumber_/i");

    tree_->Branch("pt1", &pt1_,"pt1_/f");
    tree_->Branch("eta1", &eta1_,"eta1_/f");
    tree_->Branch("phi1", &phi1_,"phi1_/f");
    tree_->Branch("q1", &q1_,"q1_/f");
    tree_->Branch("pfIso1", &pfIso1_,"pfIso1_/f");
    tree_->Branch("tkIso1", &tkIso1_,"tkIso1_/f");
    tree_->Branch("dxy1", &dxy1_,"dxy1_/f");
    tree_->Branch("dz1", &dz1_,"dz1_/f");
    //tree_->Branch("is1IsoMu27", &is1IsoMu27_,"is1IsoMu27_/b");

    tree_->Branch("pt2", &pt2_,"pt2_/f");
    tree_->Branch("eta2", &eta2_,"eta2_/f");
    tree_->Branch("phi2", &phi2_,"phi2_/f");
    tree_->Branch("q2", &q2_,"q2_/f");
    tree_->Branch("pfIso2", &pfIso2_,"pfIso2_/f");
    tree_->Branch("tkIso2", &tkIso2_,"tkIso2_/f");
    tree_->Branch("dxy2", &dxy2_,"dxy2_/f");
    tree_->Branch("dz2", &dz2_,"dz2_/f");
    //tree_->Branch("is2IsoMu27", &is2IsoMu27_,"is2IsoMu27_/b");

    tree_->Branch("is2HLT", &is2HLT_,"is2HLT_/b");
    tree_->Branch("isSel", &isSel_,"isSel_/b");
    tree_->Branch("isGlo", &isGlo_,"isGlo_/b");
    tree_->Branch("isSta", &isSta_,"isSta_/b");
    tree_->Branch("isTrk", &isTrk_,"isTrk_/b");

    tree_->Branch("dilepMass", &dilepMass_,"dilepMass_/f");

}

// ------------ method called once each job just after ending the event loop  ------------
void
TnPPairTreeProducer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TnPPairTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//--------------------------------------------------------------------------------------------------
void TnPPairTreeProducer::clearVariables(){
    edm::LogVerbatim("TnPPairTreeProducer")<<"clearVariables()";

    nPV_ = 0;
    run_ = 0;
    ls_ = 0;
    eventNumber_ = 0;

}

//--------------------------------------------------------------------------------------------------
void TnPPairTreeProducer::clearTagVariables(){
    edm::LogVerbatim("TnPPairTreeProducer")<<"clearTagVariables()";

    pt1_ = 0.;
    eta1_ = 0.;
    phi1_ = 0.;
    q1_ = 0;
    pfIso1_ = 0.;
    tkIso1_ = 0.;
    dxy1_ = 0.;
    dz1_ = 0.;
    // is1IsoMu27_ = false;
}

//--------------------------------------------------------------------------------------------------
void TnPPairTreeProducer::clearProbeVariables(){
    //edm::LogVerbatim("TnPPairTreeProducer")<<"clearProbeVariables()";

    pt2_ = 0.;
    eta2_ = 0.;
    phi2_ = 0.;
    q2_ = 0;
    pfIso2_ = 0.;
    tkIso2_ = 0.;
    dxy2_ = 0.;
    dz2_ = 0.;
    //is2IsoMu27_ = false;

    is2HLT_ = false;
    isSel_ = false;
    isGlo_ = false;
    isSta_ = false;
    isTrk_ = false;

    dilepMass_ = 0.;
}

//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::isMuonTrigger() {
    if (triggers->pass(fMuonHLTNames)){
        return true;
    }
    return false;
}

//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::isMuonTriggerObj(const double eta, const double phi) {
    if (triggers->passObj(fMuonHLTNames, eta, phi))
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::isCustomID(const reco::Muon& muon, const reco::Vertex& vtx){
    return muon.isGlobalMuon()
                && muon.isPFMuon()
                && muon.globalTrack()->normalizedChi2() < 10.
                && muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
                && muon.numberOfMatchedStations() > 1
                && (getDxy(muon, vtx) < DxyCut_ || DxyCut_ <=0.)
                && (getDz(muon, vtx) < DzCut_ || DzCut_ <=0.)
                && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
                && muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
}

//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::passMuonID(const reco::Muon& muon, const reco::Vertex& vtx) {
    //Muon ID selection, using internal function "DataFormats/MuonReco/src/MuonSelectors.cc

    if (IDType_ == LooseID && muon::isLooseMuon(muon))
        return true;
    else if (IDType_ == MediumID && muon::isMediumMuon(muon))
        return true;
    else if (IDType_ == TightID && muon::isTightMuon(muon, vtx))
        return true;
    else if (IDType_ == CustomID && isCustomID(muon, vtx))
        return true;
    else if (IDType_ == NoneID)
        return true;

    return false;
}


//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::passMuonIso(const reco::Muon& muon) {
    //Muon isolation selection, up-to-date with MUO POG recommendation

    if (IsoType_ == TrackerIso && getTkIso(muon) < IsoCut_)
        return true;
    else if (IsoType_ == PFIso && getPFIso(muon) < IsoCut_)
        return true;
    else if (IsoType_ == NoneIso)
        return true;

    return false;
}

//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::passBaseline(const reco::Muon& muon){
    return muon.isGlobalMuon()
                && muon.isPFMuon()
                && muon.globalTrack()->normalizedChi2() < 10.
                && muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
                && muon.numberOfMatchedStations() > 1
                && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
                && muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
}


//--------------------------------------------------------------------------------------------------
float TnPPairTreeProducer::getPFIso(const reco::Muon& muon) {
    return (muon.pfIsolationR04().sumChargedHadronPt +
                   std::max(0.,
                            muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt -
                                0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();
}

//--------------------------------------------------------------------------------------------------
float TnPPairTreeProducer::getTkIso(const reco::Muon& muon) {
    return muon.isolationR03().sumPt / muon.pt();
}

//--------------------------------------------------------------------------------------------------
float TnPPairTreeProducer::getDxy(const reco::Muon& muon, const reco::Vertex& vtx) {
    return fabs(muon.muonBestTrack()->dxy(vtx.position()));
}

//--------------------------------------------------------------------------------------------------
float TnPPairTreeProducer::getDz(const reco::Muon& muon, const reco::Vertex& vtx) {
    return fabs(muon.muonBestTrack()->dz(vtx.position()));
}
