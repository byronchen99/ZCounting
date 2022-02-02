
// CMSSW includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "ZCounting/TnPPairTreeProducer/plugins/TnPPairTreeProducer.h"
#include "ZCounting/TnPPairTreeProducer/interface/getFilename.h"

#include "ZCounting/ZUtils/interface/triggertool.h"
#include "ZCounting/ZUtils/interface/Helper.h"

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

    // Muon-specific Cuts
    IDTypestr_ = iConfig.getUntrackedParameter<std::string>("IDType");
    IsoTypestr_ = iConfig.getUntrackedParameter<std::string>("IsoType");
    IsoCut_ = iConfig.getUntrackedParameter<double>("IsoCut");
    DxyCut_ = iConfig.getUntrackedParameter<double>("DxyCut");
    DzCut_ = iConfig.getUntrackedParameter<double>("DzCut");
    emulateTrigger_ = iConfig.getUntrackedParameter<bool>("emulateTrigger");

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
        if(!isGoodPV(itVtx))
            continue;

        if (nPV_ == 0) {
            pv = &itVtx;
        }
        nPV_++;
    }
    
    dzPV_ = pv->z();
    drhoPV_ = pv->position().Rho();


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
        
        const bool passTrigger1 = isMuonTriggerObj(eta1_, phi1_) || (emulateTrigger_ && isMuonTriggerObjEmulated(eta1_, phi1_, eventNumber_));
            
        if (!passTrigger1)
            continue;


        pfIso1_ = getPFIso(itMu1);
        tkIso1_ = getTkIso(itMu1);
        dxy1_ = itMu1.muonBestTrack()->dxy();
        dz1_ = itMu1.muonBestTrack()->dz();
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
                const bool passTrigger2 = isMuonTriggerObj(eta2_, phi2_) || (emulateTrigger_ && isMuonTriggerObjEmulated(eta2_, phi2_, eventNumber_));

                if (passTrigger2) {
                    // category 2HLT: both muons passing trigger requirements
                    if (&itMu1 > &itMu2)
                        continue;  // make sure we don't double count MuMu2HLT category
                    is2HLT_ = true;
                }
                else // category 1HLT: probe passing selection but not trigger
                    isSel_ = true;
            }
            else if (itMu2.isGlobalMuon()
            ){      // has valid global track
                isGlo_ = true;
            }

            if (itMu2.isStandAloneMuon()
            ){      // has valid outer track
                isSta_ = true;
            }
            if (itMu2.innerTrack().isNonnull()
            ){      // has valid inner track
                nTrackerLayers2_ = itMu2.innerTrack()->hitPattern().trackerLayersWithMeasurement();
                nValidPixelHits2_ = itMu2.innerTrack()->hitPattern().numberOfValidPixelHits();
                trackAlgo2_ = itMu2.innerTrack()->algo();
                isTrk_ = true;
            }

            pfIso2_ = getPFIso(itMu2);
            tkIso2_ = getTkIso(itMu2);
            dxy2_ = itMu2.muonBestTrack()->dxy();
            dz2_ = itMu2.muonBestTrack()->dz();
            // is2IsoMu27_ = triggers->passObj("HLT_IsoMu27_v*", eta2_, phi2_);
            delR_ = vTag.DeltaR(vProbe);

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
            delR_ = vTag.DeltaR(vProbe);

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

            dxy2_ = itTrk.dxy(pv->position());
            dz2_ = itTrk.dz(pv->position());
            dxy2_ = itTrk.dxy();
            dz2_ = itTrk.dz();
            nTrackerLayers2_ = itTrk.hitPattern().trackerLayersWithMeasurement();
            nValidPixelHits2_ = itTrk.hitPattern().numberOfValidPixelHits();
            trackAlgo2_ = itTrk.algo();
            isTrk_ = true;

            tree_->Fill();

        }    //End of probe loop over tracks

    }

}

// ------------ method called once each new LS  ------------
void
TnPPairTreeProducer::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup)
{
    ls_ = iLumi.id().luminosityBlock();

    std::cout<<"TnPPairTreeProducer::beginLuminosityBlock --- now at LS "<< ls_ <<std::endl;

    if(emulateTrigger_){
        // find and open file with emulated HLT information
        const std::string fNameHLT = getFilename(run_, ls_);

        if(_fileHLTEmulation != 0)
            _fileHLTEmulation->Close();

        _fileHLTEmulation = TFile::Open(("root://xrootd-cms.infn.it//"+fNameHLT).c_str());        
        fs->cd();
    }
}

// ------------ method called once each run just before starting event loop  ------------
void TnPPairTreeProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    edm::LogVerbatim("TnPPairTreeProducer") << "now at "<<iRun.id();
    run_ = iRun.id().run();

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
    tree_->Branch("eventNumber", &eventNumber_, "eventNumber_/l");

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
    tree_->Branch("nTrackerLayers2", &nTrackerLayers2_,"nTrackerLayers2_/f");
    tree_->Branch("nValidPixelHits2", &nValidPixelHits2_,"nValidPixelHits2_/f");
    tree_->Branch("trackAlgo2", &trackAlgo2_,"trackAlgo2_/i");

    //tree_->Branch("is2IsoMu27", &is2IsoMu27_,"is2IsoMu27_/b");

    tree_->Branch("is2HLT", &is2HLT_,"is2HLT_/b");
    tree_->Branch("isSel", &isSel_,"isSel_/b");
    tree_->Branch("isGlo", &isGlo_,"isGlo_/b");
    tree_->Branch("isSta", &isSta_,"isSta_/b");
    tree_->Branch("isTrk", &isTrk_,"isTrk_/b");

    tree_->Branch("dilepMass", &dilepMass_,"dilepMass_/f");
    tree_->Branch("delR", &delR_,"delR_/f");
    tree_->Branch("drhoPV", &drhoPV_,"drhoPV_/f");
    tree_->Branch("dzPV", &dzPV_,"dzPV_/f");

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
    dxy2_ = 0.;
    dz2_ = 0.;
    nTrackerLayers2_ = 0.;
    nValidPixelHits2_ = 0.;
    trackAlgo2_ = 0;
    //is2IsoMu27_ = false;

    is2HLT_ = false;
    isSel_ = false;
    isGlo_ = false;
    isSta_ = false;
    isTrk_ = false;

    dilepMass_ = 0.;
    delR_ = 0.;
}

//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::isMuonTrigger() {
    return triggers->pass(fMuonHLTNames);
}

//--------------------------------------------------------------------------------------------------
bool TnPPairTreeProducer::isMuonTriggerObj(const double eta, const double phi) {
    return triggers->passObj(fMuonHLTNames, eta, phi);
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
    else if (IDType_ == CustomID && isCustomTightMuon(muon))
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
// For trigger emulation in the 2017H (low PU) dataset
// We emulated the HLT_IsoMu24_v11 in separated samples and need to get the trigger objects from these separate samples

bool TnPPairTreeProducer::isMuonTriggerObjEmulated(const double eta, const double phi, const long long unsigned eventNumber) {

    // filter tag for HLT_IsoMu24_v11:
    const edm::InputTag filterTag("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07", "", "HLTX");
    
    TTreeReader myReader("Events", _fileHLTEmulation);
    TTreeReaderValue<edm::EventAuxiliary> eventAuxiliary_(myReader, "EventAuxiliary");
    TTreeReaderValue<trigger::TriggerEvent> triggerEvent_(myReader, "triggerTriggerEvent_hltTriggerSummaryAOD__HLTX.obj");

    // std::cout<<"Look for event "<< eventNumber <<std::endl;

    while(myReader.Next()){
        // find event
        if(eventNumber != eventAuxiliary_->event())
            continue;

        // std::cout<<"Found event!"<<std::endl;
        // look for trigger objects
        if(triggerEvent_->filterIndex(filterTag) < triggerEvent_->sizeFilters()){
            const trigger::Keys& trigKeys = triggerEvent_->filterKeys(triggerEvent_->filterIndex(filterTag));
            const trigger::TriggerObjectCollection & trigObjColl(triggerEvent_->getObjects());
            //now loop of the trigger objects passing filter
            for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
                const trigger::TriggerObject& obj = trigObjColl[*keyIt];
                // std::cout<<"Trigger object(pt | eta) = "<<obj.pt()<< " | "<<obj.eta()<<std::endl;
                if (reco::deltaR(eta, phi, obj.eta(), obj.phi()) < fMuonHLTDRMAX){
                    return true;
                }
            }
        }
        return false;
    }
    edm::LogWarning("isMuonTriggerObjEmulated")<<"Event was not found!"<<std::endl;
    return false;
}
