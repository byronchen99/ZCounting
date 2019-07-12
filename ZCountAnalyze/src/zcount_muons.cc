#include "ZCounting/ZCountAnalyze/interface/zcount_muons.h"

#include <TLorentzVector.h>

#include <algorithm>

zcount_muons::zcount_muons():zcount_module(){}

zcount_muons::~zcount_muons(){}

void zcount_muons::getInput(const edm::ParameterSet& iConfig){
    edm::LogVerbatim("zcount_muons") << " getInput";

    PtCutTag_    = iConfig.getUntrackedParameter<double>("MuonTagPtCut");
    PtCutProbe_  = iConfig.getUntrackedParameter<double>("MuonProbePtCut");
    EtaCutTag_   = iConfig.getUntrackedParameter<double>("MuonTagEtaCut");
    EtaCutProbe_ = iConfig.getUntrackedParameter<double>("MuonProbeEtaCut");

    IDTypestr_   = iConfig.getUntrackedParameter<std::string>("MuonIDType");
    IsoTypestr_  = iConfig.getUntrackedParameter<std::string>("MuonIsoType");
    IsoCut_      = iConfig.getUntrackedParameter<double>("MuonIsoCut");

    if (IDTypestr_ == "Loose")
        IDType_ = LooseID;
    else if (IDTypestr_ == "Medium")
        IDType_ = MediumID;
    else if (IDTypestr_ == "Tight")
        IDType_ = TightID;
    else
        IDType_ = NoneID;

    if (IsoTypestr_ == "Tracker-based")
        IsoType_ = TrackerIso;
    else if (IsoTypestr_ == "PF-based")
        IsoType_ = PFIso;
    else
        IsoType_ = NoneIso;

    fMuonHLTNames       = iConfig.getParameter<std::vector<std::string>>("MuonTriggerNames");
    fMuonHLTObjectNames = iConfig.getParameter<std::vector<std::string>>("MuonTriggerObjectNames");

    if (fMuonHLTNames.size() != fMuonHLTObjectNames.size()) {
        edm::LogError("zcount_muons") << "List of MuonTriggerNames and MuonTriggerObjectNames has to be the same length";
    }

    edm::LogVerbatim("zcount_muons") << "getInput: set trigger names";
    for(unsigned int i = 0; i < fMuonHLTNames.size(); ++i) {
        triggerModule->addTriggerRecord(fMuonHLTNames.at(i), fMuonHLTObjectNames.at(i));
    }

    MassMin_ = iConfig.getUntrackedParameter<double>("MassMin");
    MassMax_ = iConfig.getUntrackedParameter<double>("MassMax");

    edm::LogVerbatim("zcount_muons") << "getInput: done";
}

void zcount_muons::initBranches(TTree* tree){

    addBranch(tree,"nMuonTags",  &nTag_,  "nTag_/i");
    addBranch(tree,"MuonTagPt",  &tagPt_, "tagPt_[nTag_]/f");
    addBranch(tree,"MuonTagEta", &tagEta_,"tagEta_[nTag_]/f");
    addBranch(tree,"MuonTagPhi", &tagPhi_,"tagPhi_[nTag_]/f");

    addBranch(tree,"nMuonProbes",        &nProbe_,    "nProbe_/i");
    addBranch(tree,"MuonProbePt",        &probePt_,   "probePt_[nProbe_]/f");
    addBranch(tree,"MuonProbeEta",       &probeEta_,  "probeEta_[nProbe_]/f");
    addBranch(tree,"MuonProbePhi",       &probePhi_,  "probePhi_[nProbe_]/f");
    addBranch(tree,"MuonProbeCategory",  &probeCat_,  "probeCat_[nProbe_]/i");
    addBranch(tree,"MuonDilepMass",      &dilepMass_, "dilepMass_[nProbe_]/f");
    addBranch(tree,"MuonTagIndex",       &tagIndex_,  "tagIndex_[nProbe_]/i");

}

bool zcount_muons::readEvent(const edm::Event& iEvent){
    edm::LogVerbatim("zcount_muons") << "zcount_muons: readEvent";

    nTag_ = 0;
    nProbe_ = 0;
    recoMuons.clear();

    if(!isMuonTrigger()){
        edm::LogVerbatim("zcount_muons") << "zcount muons: event did not pass any muon trigger";
        return false;
    }

    iEvent.getByToken(fTrackName_token, hTrackProduct);
    if (!hTrackProduct.isValid()){
        edm::LogWarning("zcount_muons") << "zcount_muons: no valid track product";
        return false;
    }

    iEvent.getByToken(fMuonName_token, hMuonProduct);
    if (!hMuonProduct.isValid()){
        edm::LogWarning("zcount_muons") << "zcount_muons: no valid muon product";
        return false;
    }

    TLorentzVector vTag(0., 0., 0., 0.);
    TLorentzVector vProbe(0., 0., 0., 0.);
    TLorentzVector vTrack(0., 0., 0., 0.);

    const reco::Vertex& pv = *(pvModule->getPV());

    edm::LogVerbatim("zcount_muons") << "zcount_muons: nMuons = "<<hMuonProduct->size();
    edm::LogVerbatim("zcount_muons") << "zcount_muons: nTracks = "<<hTrackProduct->size();


    // Tag loop
    for (auto const& itMu1 : *hMuonProduct) {
        if((unsigned int)nTag_ >= max_num_tags || (unsigned int)nProbe_ >= max_num_probes)
            break;

        float pt1 = itMu1.muonBestTrack()->pt();
        float eta1 = itMu1.muonBestTrack()->eta();
        float phi1 = itMu1.muonBestTrack()->phi();
        float q1 = itMu1.muonBestTrack()->charge();

        // Tag selection: kinematic cuts, lepton selection and trigger matching
        if (pt1 < PtCutTag_)
            continue;
        if (fabs(eta1) > EtaCutTag_)
            continue;
        if (!(passMuonID(itMu1, pv, IDType_) && passMuonIso(itMu1, IsoType_, IsoCut_)))
            continue;
        if (!isMuonTriggerObj(eta1, phi1))
            continue;

        vTag.SetPtEtaPhiM(pt1, eta1, phi1, MUON_MASS);
        recoMuons.push_back({vTag,cHLT,q1});

        edm::LogVerbatim("zcount_muons") << "zcount_muons: good tag muon found";

        bool goodTPPair = false; // if at least one good probe is found to the corresponding tag

        // Probe loop over muons
        for (auto const& itMu2 : *hMuonProduct) {
            if((unsigned int)nProbe_ >= max_num_probes)
                break;

            probeCat_[nProbe_] = cNone;
            if (&itMu2 == &itMu1)
                continue;

            float pt2 = itMu2.muonBestTrack()->pt();
            float eta2 = itMu2.muonBestTrack()->eta();
            float phi2 = itMu2.muonBestTrack()->phi();
            float q2 = itMu2.muonBestTrack()->charge();

            edm::LogVerbatim("zcount_muons") << "zcount_muons: probe pt = "<<pt2;
            edm::LogVerbatim("zcount_muons") << "zcount_muons: probe eta = "<<eta2;
            edm::LogVerbatim("zcount_muons") << "zcount_muons: probe c1/c2 = "<<q1<<"/"<<q2;


            // Probe selection: kinematic cuts and opposite charge requirement
            if (pt2 < PtCutProbe_)
                continue;
            if (fabs(eta2) > EtaCutProbe_)
                continue;
            if (q1 == q2)
                continue;

            vProbe.SetPtEtaPhiM(pt2, eta2, phi2, MUON_MASS);

            // Mass window
            TLorentzVector vDilep = vTag + vProbe;
            float dilepMass = vDilep.M();
            edm::LogVerbatim("zcount_muons") << "zcount_muons: dilepMass = "<<dilepMass;
            if ((dilepMass < MassMin_) || (dilepMass > MassMax_))
                continue;


            // Determine event category for efficiency calculation
            if (passMuonID(itMu2, pv, IDType_) && passMuonIso(itMu2, IsoType_, IsoCut_)) {
                if (isMuonTriggerObj(eta2, phi2)) {

                    // category (2)HLT: both muons passing trigger requirements
                    if (&itMu1 > &itMu2)
                        continue;  // make sure we don't double count MuMu2HLT category

                    probeCat_[nProbe_] = cHLT;
                    recoMuons.push_back({vProbe,cHLT,q2});

                    edm::LogVerbatim("zcount_muons") << "zcount_muons: HLT probe";

                }
                else {
                    // category Sel: probe passing selection but not trigger
                    probeCat_[nProbe_] = cSel;
                    recoMuons.push_back({vProbe,cSel,q2});
                    edm::LogVerbatim("zcount_muons") << "zcount_muons: Sel probe";

                }

            }
            else if (itMu2.isGlobalMuon()) {
                // category Glo: probe is a Global muon but failing selection
                probeCat_[nProbe_] = cGlo;
                recoMuons.push_back({vProbe,cGlo,q2});
                edm::LogVerbatim("zcount_muons") << "zcount_muons: Glo probe";

            }
            else if (itMu2.isStandAloneMuon()) {
                // category Sta: probe is a Standalone muon
                probeCat_[nProbe_] = cSta;
                recoMuons.push_back({vProbe,cSta,q2});
                edm::LogVerbatim("zcount_muons") << "zcount_muons: Sta probe";

            }
            else if (itMu2.innerTrack()->hitPattern().trackerLayersWithMeasurement() >= 6 &&
                 itMu2.innerTrack()->hitPattern().numberOfValidPixelHits() >= 1) {
                // cateogry Trk: probe is a tracker track
                probeCat_[nProbe_] = cTrk;
                recoMuons.push_back({vProbe,cTrk,q2});
                edm::LogVerbatim("zcount_muons") << "zcount_muons: Trk probe";
            }

            if(probeCat_[nProbe_] != cNone){
                probePt_[nProbe_]   = pt2;
                probeEta_[nProbe_]  = eta2;
                probePhi_[nProbe_]  = phi2;
                dilepMass_[nProbe_] = dilepMass;
                tagIndex_[nProbe_]  = nTag_;

                nProbe_++;

                goodTPPair = true;

            }

        } // End of probe loop over muons

        // Probe loop over tracks, only for standalone efficiency calculation
        for (auto const& itTrk : *hTrackProduct) {
            if((unsigned int)nProbe_ >= max_num_probes)
                break;

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

            float pt2 = itTrk.pt();
            float eta2 = itTrk.eta();
            float phi2 = itTrk.phi();
            float q2 = itTrk.charge();

            // Probe selection:  kinematic cuts and opposite charge requirement
            if (pt2 < PtCutProbe_)
                continue;
            if (fabs(eta2) > EtaCutProbe_)
                continue;
            if (q1 == q2)
                continue;

            vTrack.SetPtEtaPhiM(pt2, eta2, phi2, MUON_MASS);

            TLorentzVector vDilep = vTag + vTrack;
            float dilepMass = vDilep.M();
            if ((dilepMass < MassMin_) || (dilepMass > MassMax_))
                continue;


            if (itTrk.hitPattern().trackerLayersWithMeasurement() >= 6 && itTrk.hitPattern().numberOfValidPixelHits() >= 1) {
                probeCat_[nProbe_] = cTrk;
                edm::LogVerbatim("zcount_muons") << "zcount_muons: Trk probe";
                recoMuons.push_back({vTrack,cTrk,q2});

                probePt_[nProbe_]   = pt2;
                probeEta_[nProbe_]  = eta2;
                probePhi_[nProbe_]  = phi2;
                dilepMass_[nProbe_] = dilepMass;
                tagIndex_[nProbe_]  = nTag_;

                nProbe_++;

                goodTPPair = true;
            }

        } //End of probe loop over tracks

        if(goodTPPair){
            tagPt_[nTag_] = pt1;
            tagEta_[nTag_] = eta1;
            tagPhi_[nTag_] = phi1;

            nTag_++;
        }

    } //End of tag loop

    if(nProbe_ > 0)
        edm::LogVerbatim("zcount_muons") << "zcount_muons: good tag-and-probe muon pair found";
    else
        return false;

    return true;
}

void zcount_muons::fillBranches(){
    return;
}

//--------------------------------------------------------------------------------------------------
bool zcount_muons::isMuonTrigger() {
    if (triggerModule->pass(fMuonHLTNames))
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
bool zcount_muons::isMuonTriggerObj(const double eta, const double phi) {
    if (triggerModule->passObj(fMuonHLTNames, fMuonHLTObjectNames, eta, phi))
        return true;
    return false;
}
//--------------------------------------------------------------------------------------------------
// Muon ID selection, using internal function "DataFormats/MuonReco/src/MuonSelectors.cc
bool zcount_muons::passMuonID(
    const reco::Muon& muon,
    const reco::Vertex& vtx,
    const MuonIDTypes& idType) {

    if (idType == LooseID && muon::isLooseMuon(muon))
        return true;
    else if (idType == MediumID && muon::isMediumMuon(muon))
        return true;
    else if (idType == TightID && muon::isTightMuon(muon, vtx))
        return true;
    else if (idType == NoneID)
        return true;
    else
        return false;
}
//--------------------------------------------------------------------------------------------------
//Muon isolation selection, up-to-date with MUO POG recommendation (FIXME: CHECK PLS)
bool zcount_muons::passMuonIso(const reco::Muon& muon,
                            const MuonIsoTypes& isoType,
                            const float isoCut) {

    if (isoType == TrackerIso && muon.isolationR03().sumPt < isoCut)
        return true;
    else if (isoType == PFIso &&
           muon.pfIsolationR04().sumChargedHadronPt +
                   std::max(0.,
                            muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt -
                                0.5 * muon.pfIsolationR04().sumPUPt) < isoCut)
        return true;
    else if (isoType == NoneIso)
        return true;
    else
        return false;
}

