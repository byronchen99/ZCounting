#include "ZCounting/ZCountAnalyze/interface/zcount_genInfo.h"

#include <algorithm>

zcount_genInfo::zcount_genInfo():zcount_module(){}

zcount_genInfo::~zcount_genInfo(){}

void zcount_genInfo::getInput(const edm::ParameterSet& iConfig){
    edm::LogVerbatim("zcount_genInfo") << "zcount_genInfo: getInput";
}

void zcount_genInfo::initBranches(TTree* tree){
    addBranch(tree,"eventWeight", &eventWeight_,"eventWeight_/f");

    // gen particle info
    addBranch(tree,"ZPt",         &ZPt_,"ZPt_/f");
    addBranch(tree,"ZEta",        &ZEta_,"ZEta_/f");
    addBranch(tree,"ZPhi",        &ZPhi_,"ZPhi_/f");
    addBranch(tree,"ZM",          &ZM_,"ZM_/f");
    addBranch(tree,"ZDecayMode",  &ZDecayMode_,"ZDecayMode_/i");

    addBranch(tree,"ZLeptonPt",   &ZLeptonPt_,"ZLeptonPt_/f");
    addBranch(tree,"ZLeptonEta",  &ZLeptonEta_,"ZLeptonEta_/f");
    addBranch(tree,"ZLeptonPhi",  &ZLeptonPhi_,"ZLeptonPhi_/f");

    addBranch(tree,"ZAntiLeptonPt",  &ZAntiLeptonPt_, "ZAntiLeptonPt_/f");
    addBranch(tree,"ZAntiLeptonEta", &ZAntiLeptonEta_,"ZAntiLeptonEta_/f");
    addBranch(tree,"ZAntiLeptonPhi", &ZAntiLeptonPhi_,"ZAntiLeptonPhi_/f");

    addBranch(tree,"ZStableMass",    &ZStableMass_,"ZStableMass_/f");

    // reco info of reco-gen matched particles
    addBranch(tree,"ZLeptonRecoCat", &ZLeptonRecoCat_,"ZLeptonRecoCat_[2]/i");
    addBranch(tree,"ZLeptonRecoPt",  &ZLeptonRecoPt_,"ZLeptonRecoPt_/f");
    addBranch(tree,"ZLeptonRecoEta", &ZLeptonRecoEta_,"ZLeptonRecoEta_/f");
    addBranch(tree,"ZLeptonRecoPhi", &ZLeptonRecoPhi_,"ZLeptonRecoPhi_/f");

    addBranch(tree,"ZAntiLeptonRecoCat", &ZAntiLeptonRecoCat_,"ZAntiLeptonRecoCat_[2]/i");
    addBranch(tree,"ZAntiLeptonRecoPt",  &ZAntiLeptonRecoPt_, "ZAntiLeptonRecoPt_/f");
    addBranch(tree,"ZAntiLeptonRecoEta", &ZAntiLeptonRecoEta_,"ZAntiLeptonRecoEta_/f");
    addBranch(tree,"ZAntiLeptonRecoPhi", &ZAntiLeptonRecoPhi_,"ZAntiLeptonRecoPhi_/f");

    addBranch(tree,"ZMassReco",  &ZMassReco_, "ZMassReco_/f");

    // vertex information from gen particle
    addBranch(tree,"ZVtxX",     &ZVtxX_, "ZVtxX_/f");
    addBranch(tree,"ZVtxY",     &ZVtxY_, "ZVtxY_/f");
    addBranch(tree,"ZVtxZ",     &ZVtxZ_,    "ZVtxZ_/f");
    addBranch(tree,"ZVtxRho",   &ZVtxRho_,  "ZVtxRho_/f");
    addBranch(tree,"ZVtxNDoF",  &ZVtxNDoF_, "ZVtxNDoF_/f");

    addBranch(tree,"matchRecoVtx",  &matchRecoVtx_, "matchRecoVtx_/b");
    addBranch(tree,"matchRecoGoodVtx",  &matchRecoGoodVtx_, "matchRecoGoodVtx_/b");

}

bool zcount_genInfo::readEvent(const edm::Event& iEvent){

    // ------------------ Get products from tokens ------------------- //
    // get reco vertices for vertex matching
    edm::Handle<reco::VertexCollection> hVertexProduct;
    iEvent.getByToken(fPVName_token, hVertexProduct);
    if (!hVertexProduct.isValid()) {
        edm::LogWarning("zcount_PV") << "No valid primary vertex product found" ;
        return false;
    }

    // get gen info for event weight
    edm::Handle<GenEventInfoProduct> hGenInfoProduct;
    iEvent.getByToken(fGenInfoName_token, hGenInfoProduct);
    if (!hGenInfoProduct.isValid()){
        edm::LogWarning("zcount_genInfo") << "zcount_genInfo: no valid gen info product";
        return false;
    }

    // get gen z info
    edm::Handle<std::vector<GenZDecayProperties>> hGenZInfoProduct;
    iEvent.getByToken(fGenZInfoName_token, hGenZInfoProduct);
    if (hGenZInfoProduct.failedToGet()){
         edm::LogWarning("zcount_genInfo") << "zcount_genZInfo: no valid gen Z info product";
         return false;
    }

    // get reco tracks for matching
    edm::Handle<reco::TrackCollection> hTrackProduct;
    iEvent.getByToken(fTrackName_token, hTrackProduct);
    if (!hTrackProduct.isValid()){
        edm::LogWarning("zcount_muons") << "zcount_muons: no valid track product";
        return false;
    }

    // get reco muons for matching
    edm::Handle<reco::MuonCollection> hMuonProduct;
    iEvent.getByToken(fMuonName_token, hMuonProduct);
    if (!hMuonProduct.isValid()){
        edm::LogWarning("zcount_muons") << "zcount_muons: no valid muon product";
        return false;
    }

    // ------------------ extract gen particles from gen Z token ------------------- //
    if(hGenZInfoProduct->size() < 1){
        // try GenZLeptonDecay
        iEvent.getByToken(fGenZLepInfoName_token, hGenZInfoProduct);
        if (hGenZInfoProduct.failedToGet()){
             edm::LogWarning("zcount_genInfo") << "zcount_genZInfo: no valid gen Z info product";
             return false;
        }
        else if(hGenZInfoProduct->size() < 1){
            edm::LogWarning("zcount_genInfo") << "zcount_genZInfo: no gen Z Product found";
            return false;
        }
    }
    if(hGenZInfoProduct->size() > 1)
        edm::LogWarning("zcount_genInfo") << "zcount_genZInfo: more then one gen Z found";

    if(hGenZInfoProduct->at(0).z()){
        ZGen = hGenZInfoProduct->at(0).z();
    }
    else{
        edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: no gen Z found";
        ZGen = 0;
    }
    ZDecayMode_ = hGenZInfoProduct->at(0).decayMode();

    if(hGenZInfoProduct->at(0).stableLepton()){
        ZLepton = hGenZInfoProduct->at(0).stableLepton();
    }
    else{
        edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: no stable gen Z lepton found";
        ZLepton = 0;
    }
    if(hGenZInfoProduct->at(0).stableAntiLepton()){
        ZAntiLepton = hGenZInfoProduct->at(0).stableAntiLepton();
    }
    else{
        edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: no stable gen Z antilepton found";
        ZAntiLepton = 0;
    }

    // ------------------ match gen leptons from gen Z to reco muons or reco tracks ------------------- //
    // 1. minimize delta R between gen muons and reco muons
    // 2. charge has to match
    // 3. relative pt difference has to be below threshold
    // 4. if no muon match is found, loop over track collection with same criterias 1-3
    double minDelRLep  = MUON_MAX_DELTAR_MATCH;
    double minDelRALep = MUON_MAX_DELTAR_MATCH;
    // first look if gen lepton can be matched with a reco muon
    ZLeptonRecoMuon = 0;
    ZAntiLeptonRecoMuon = 0;
    ZLeptonRecoTrack = 0;
    ZAntiLeptonRecoTrack = 0;
    for(auto const& itMu : *hMuonProduct) {
        TLorentzVector lvMuon(0., 0., 0., 0.);
        lvMuon.SetPtEtaPhiM(itMu.muonBestTrack()->pt(), itMu.muonBestTrack()->eta(), itMu.muonBestTrack()->phi(), MUON_MASS);

        if(ZLepton && itMu.muonBestTrack()->charge() == -1 && deltaR(ZLepton->polarP4(),lvMuon) < minDelRLep
                && std::abs(lvMuon.Pt() - ZLepton->pt())/ZLepton->pt() < MUON_MAX_DRELPT_MATCH ){
            ZLeptonRecoMuon = &itMu;
            minDelRLep = deltaR(lvMuon, ZLepton->polarP4());
        }
        else if(ZAntiLepton && itMu.muonBestTrack()->charge() == 1 && deltaR(lvMuon, ZAntiLepton->polarP4()) < minDelRALep
                && std::abs(lvMuon.Pt() - ZAntiLepton->pt())/ZAntiLepton->pt() < MUON_MAX_DRELPT_MATCH ){
            ZAntiLeptonRecoMuon = &itMu;
            minDelRALep = deltaR(lvMuon, ZAntiLepton->polarP4());
        }
    }
    // if gen lepton not matched with a reco muon, look if gen lepton can be matched with a reco track
    if(!ZLeptonRecoMuon || !ZAntiLeptonRecoMuon){
        for(auto const& itTrk : *hTrackProduct) {
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
            TLorentzVector lvTrk(0., 0., 0., 0.);
            lvTrk.SetPtEtaPhiM(itTrk.pt(), itTrk.eta(), itTrk.phi(), MUON_MASS);
            if(!ZLeptonRecoMuon && ZLepton && itTrk.charge() == -1 && deltaR(lvTrk, ZLepton->polarP4()) < minDelRLep
                    && std::abs(lvTrk.Pt() - ZLepton->pt())/ZLepton->pt() < MUON_MAX_DRELPT_MATCH ){
                ZLeptonRecoTrack = &itTrk;
                minDelRLep = deltaR(lvTrk, ZLepton->polarP4());
            }
            else if(!ZAntiLeptonRecoMuon && ZAntiLepton && itTrk.charge() == 1 && deltaR(lvTrk, ZAntiLepton->polarP4()) < minDelRALep
                    && std::abs(lvTrk.Pt() - ZAntiLepton->pt())/ZAntiLepton->pt() < MUON_MAX_DRELPT_MATCH ){
                ZAntiLeptonRecoTrack = &itTrk;
                minDelRALep = deltaR(lvTrk, ZAntiLepton->polarP4());
            }
        }
    }

    // ------------------ match gen vertex from gen Z to reco vertex from reco vertex collection ------------------- //
    const reco::VertexCollection* pvCol = hVertexProduct.product();
    ZRecoVtx = &(*pvCol->begin());      // reco vertex for Z in list of all vertices
    ZRecoGoodVtx = &(*pvCol->begin());  // reco vertex for Z in list of all good vertices
    if(ZVtxX_ != 0. || ZVtxY_ != 0. || ZVtxZ_ != 0){
        double minXYZ = 10.;
        for (auto const& itVtx : *hVertexProduct) {
            double distXYZ = std::sqrt( std::pow(itVtx.x() - ZVtxX_, 2) + std::pow(itVtx.y() - ZVtxY_, 2) + std::pow(itVtx.z() - ZVtxZ_, 2) );
            if(distXYZ < minXYZ){
                minXYZ = distXYZ;
                ZRecoVtx = &itVtx;
            }

            if(distXYZ < minXYZ && pvModule->isGoodPV(itVtx)){
                minXYZ = distXYZ;
                ZRecoGoodVtx = &itVtx;
            }
        }
    }

    eventWeight_ = hGenInfoProduct->weight();

    return true;
}


void zcount_genInfo::fillBranches(){
    edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: fill Branches";

    // ------------------ check if best reco vertex is z reco vertex ------------------- //
    const reco::Vertex *pv = pvModule->getPV();
    if(ZRecoVtx == pv)
        matchRecoVtx_ = true;
    else
        matchRecoVtx_ = false;

    if(ZRecoGoodVtx == pv)
        matchRecoGoodVtx_ = true;
    else
        matchRecoGoodVtx_ = false;

    // ------------------ categorization for the reco particle that is matched to the gen (anti) lepton ------------------- //
    // 1. using the ZRecoGoodVtx
    if(ZLeptonRecoMuon){
        ZLeptonRecoCat_[0] = muonModule->getMuonCategory(*ZLeptonRecoMuon, *ZRecoGoodVtx);
    }
    else if(ZLeptonRecoTrack && muonModule->isValidTrack(*ZLeptonRecoTrack)){
        ZLeptonRecoCat_[0] = cTrk;
    }
    else{
        ZLeptonRecoCat_[0] = cNone;
    }

    if(ZAntiLeptonRecoMuon){
        ZAntiLeptonRecoCat_[0] = muonModule->getMuonCategory(*ZAntiLeptonRecoMuon, *ZRecoGoodVtx);
    }
    else if(ZAntiLeptonRecoTrack && muonModule->isValidTrack(*ZAntiLeptonRecoTrack)){
        ZAntiLeptonRecoCat_[0] = cTrk;
    }
    else{
        ZAntiLeptonRecoCat_[0] = cNone;
    }
    // 2. using the first good PV
    if(ZLeptonRecoMuon){
        ZLeptonRecoCat_[1] = muonModule->getMuonCategory(*ZLeptonRecoMuon, *pv);
    }
    else if(ZLeptonRecoTrack && muonModule->isValidTrack(*ZLeptonRecoTrack)){
        ZLeptonRecoCat_[1] = cTrk;
    }
    else{
        ZLeptonRecoCat_[1] = cNone;
    }

    if(ZAntiLeptonRecoMuon){
        ZAntiLeptonRecoCat_[1] = muonModule->getMuonCategory(*ZAntiLeptonRecoMuon, *pv);
    }
    else if(ZAntiLeptonRecoTrack && muonModule->isValidTrack(*ZAntiLeptonRecoTrack)){
        ZAntiLeptonRecoCat_[1] = cTrk;
    }
    else{
        ZAntiLeptonRecoCat_[1] = cNone;
    }

    // ------------------ set output ------------------- //
    if(ZGen){
        ZPt_      = ZGen->pt();
        ZEta_     = ZGen->eta();
        ZPhi_     = ZGen->phi();
        ZM_       = ZGen->mass();

        ZVtxX_    = ZGen->vx();
        ZVtxY_    = ZGen->vy();
        ZVtxZ_    = ZGen->vz();
        ZVtxNDoF_ = ZGen->vertexNdof();
        ZVtxRho_  = ZGen->vertex().Rho();
    }
    else{
        ZPt_      = 0;
        ZEta_     = 0;
        ZPhi_     = 0;
        ZM_       = 0;

        ZVtxX_    = 0;
        ZVtxY_    = 0;
        ZVtxZ_    = 0;
        ZVtxNDoF_ = 0;
        ZVtxRho_  = 0;
    }

    if(ZLepton){
        ZLeptonPt_  = ZLepton->pt();
        ZLeptonEta_ = ZLepton->eta();
        ZLeptonPhi_ = ZLepton->phi();
    }
    else{
        ZLeptonPt_  = 0;
        ZLeptonEta_ = 0;
        ZLeptonPhi_ = 0;
    }

    if(ZAntiLepton){
        ZAntiLeptonPt_ = ZAntiLepton->pt();
        ZAntiLeptonEta_ = ZAntiLepton->eta();
        ZAntiLeptonPhi_ = ZAntiLepton->phi();
    }
    else{
        ZAntiLeptonPt_ = 0;
        ZAntiLeptonEta_ = 0;
        ZAntiLeptonPhi_ = 0;
    }

    if(ZLepton && ZAntiLepton){
        ZStableMass_ = (ZLepton->polarP4() + ZAntiLepton->polarP4()).M();
    }
    else{
        ZStableMass_ = 0;
    }

    TLorentzVector lvMuon(0., 0., 0., 0.);
    if(ZLeptonRecoMuon){
        lvMuon.SetPtEtaPhiM(ZLeptonRecoMuon->muonBestTrack()->pt(), ZLeptonRecoMuon->muonBestTrack()->eta(),
                              ZLeptonRecoMuon->muonBestTrack()->phi(), MUON_MASS);
    }
    else if(ZLeptonRecoTrack){
        lvMuon.SetPtEtaPhiM(ZLeptonRecoTrack->pt(), ZLeptonRecoTrack->eta(),
                              ZLeptonRecoTrack->phi(), MUON_MASS);
    }

    TLorentzVector lvAntiMuon(0., 0., 0., 0.);
    if(ZAntiLeptonRecoMuon){
        lvAntiMuon.SetPtEtaPhiM(ZAntiLeptonRecoMuon->muonBestTrack()->pt(), ZAntiLeptonRecoMuon->muonBestTrack()->eta(),
                              ZAntiLeptonRecoMuon->muonBestTrack()->phi(), MUON_MASS);
    }
    else if(ZAntiLeptonRecoTrack){
        lvAntiMuon.SetPtEtaPhiM(ZAntiLeptonRecoTrack->pt(), ZAntiLeptonRecoTrack->eta(),
                              ZAntiLeptonRecoTrack->phi(), MUON_MASS);
    }

    ZLeptonRecoPt_  = lvMuon.Pt();
    ZLeptonRecoEta_ = lvMuon.Eta();
    ZLeptonRecoPhi_ = lvMuon.Phi();

    ZAntiLeptonRecoPt_  = lvMuon.Pt();
    ZAntiLeptonRecoEta_ = lvMuon.Eta();
    ZAntiLeptonRecoPhi_ = lvMuon.Phi();

    ZMassReco_ = (lvAntiMuon + lvMuon).M();

    edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: end fill Branches";

    return;
}
