#include "ZCounting/ZCountAnalyze/interface/zcount_genInfo.h"

#include <algorithm>

zcount_genInfo::zcount_genInfo():zcount_module(){}

zcount_genInfo::~zcount_genInfo(){}

void zcount_genInfo::getInput(const edm::ParameterSet& iConfig){
    edm::LogVerbatim("zcount_genInfo") << "zcount_genInfo: getInput" ;

}

void zcount_genInfo::initBranches(TTree* tree){
    addBranch(tree,"eventWeight", &eventWeight_,"eventWeight_/f");
    addBranch(tree,"ZPt", &ZPt_,"ZPt_/f");
    addBranch(tree,"ZEta", &ZEta_,"ZEta_/f");
    addBranch(tree,"ZPhi", &ZPhi_,"ZPhi_/f");
    addBranch(tree,"ZM", &ZM_,"ZM_/f");
    addBranch(tree,"ZStableMass", &ZStableMass_,"ZStableMass_/f");
    addBranch(tree,"ZDecayMode",  &ZDecayMode_,"ZDecayMode_/i");

    addBranch(tree,"ZLeptonPt",  &ZLeptonPt_,"ZLeptonPt_/f");
    addBranch(tree,"ZLeptonEta", &ZLeptonEta_,"ZLeptonEta_/f");
    addBranch(tree,"ZLeptonPhi", &ZLeptonPhi_,"ZLeptonPhi_/f");
    addBranch(tree,"ZLeptonRecoCat", &ZLeptonRecoCat_,"ZLeptonRecoCat_/i");
    addBranch(tree,"ZLeptonRecoPt",  &ZLeptonRecoPt_,"ZLeptonRecoPt_/f");
    addBranch(tree,"ZLeptonRecoEta", &ZLeptonRecoEta_,"ZLeptonRecoEta_/f");
    addBranch(tree,"ZLeptonRecoPhi", &ZLeptonRecoPhi_,"ZLeptonRecoPhi_/f");
    addBranch(tree,"ZLeptonRecoDelR", &ZLeptonRecoDelR_,"ZLeptonRecoDelR_/f");

    addBranch(tree,"ZAntiLeptonPt",  &ZAntiLeptonPt_, "ZAntiLeptonPt_/f");
    addBranch(tree,"ZAntiLeptonEta", &ZAntiLeptonEta_,"ZAntiLeptonEta_/f");
    addBranch(tree,"ZAntiLeptonPhi", &ZAntiLeptonPhi_,"ZAntiLeptonPhi_/f");
    addBranch(tree,"ZAntiLeptonRecoCat", &ZAntiLeptonRecoCat_,"ZAntiLeptonRecoCat_/i");
    addBranch(tree,"ZAntiLeptonRecoPt",  &ZAntiLeptonRecoPt_, "ZAntiLeptonRecoPt_/f");
    addBranch(tree,"ZAntiLeptonRecoEta", &ZAntiLeptonRecoEta_,"ZAntiLeptonRecoEta_/f");
    addBranch(tree,"ZAntiLeptonRecoPhi", &ZAntiLeptonRecoPhi_,"ZAntiLeptonRecoPhi_/f");
    addBranch(tree,"ZAntiLeptonRecoDelR", &ZAntiLeptonRecoDelR_,"ZAntiLeptonRecoDelR_/f");

    addBranch(tree,"ZMassReco",  &ZMassReco_, "ZMassReco_/f");


}

bool zcount_genInfo::readEvent(const edm::Event& iEvent){

    edm::Handle<GenEventInfoProduct> hGenInfoProduct;
    iEvent.getByToken(fGenInfoName_token, hGenInfoProduct);
    if (!hGenInfoProduct.isValid()){
        edm::LogWarning("zcount_genInfo") << "zcount_genInfo: no valid gen info product";
        return false;
    }

    edm::Handle<std::vector<GenZDecayProperties>> hGenZInfoProduct;
    iEvent.getByToken(fGenZInfoName_token, hGenZInfoProduct);
    if (hGenZInfoProduct.failedToGet()){
         edm::LogWarning("zcount_genInfo") << "zcount_genZInfo: no valid gen Z info product";
         return false;
    }

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
        ZPt_        = hGenZInfoProduct->at(0).z()->pt();
        ZEta_       = hGenZInfoProduct->at(0).z()->eta();
        ZPhi_       = hGenZInfoProduct->at(0).z()->phi();
        ZM_         = hGenZInfoProduct->at(0).z()->mass();
    }
    else{
        edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: no gen Z found";
        ZPt_        = 0;
        ZEta_       = 0;
        ZPhi_       = 0;
        ZM_         = 0;
    }

    ZDecayMode_ = hGenZInfoProduct->at(0).decayMode();

    LV lep = nullP4_;
    if(hGenZInfoProduct->at(0).stableLepton())
        lep = hGenZInfoProduct->at(0).stableLepton()->polarP4();
    ZLepton_ = TLorentzVector(lep.X(),lep.Y(),lep.Z(),lep.T());

    lep = nullP4_;
    if(hGenZInfoProduct->at(0).stableAntiLepton())
        lep = hGenZInfoProduct->at(0).stableAntiLepton()->polarP4();
    ZAntiLepton_ = TLorentzVector(lep.X(),lep.Y(),lep.Z(),lep.T());

    if(hGenZInfoProduct->at(0).stableLepton() && hGenZInfoProduct->at(0).stableAntiLepton()){
        ZStableMass_ = (ZAntiLepton_ + ZLepton_).M();
    }
    else{
        ZStableMass_ = 0;
    }


    eventWeight_ = hGenInfoProduct->weight();

    return true;
}


void zcount_genInfo::fillBranches(){
    edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: fill Branches";

    const std::vector<recoMuon> *recoMuons = muonModule->getRecoMuons();

    ZLeptonRecoCat_ = 0;
    TLorentzVector ZLeptonReco_;

    ZAntiLeptonRecoCat_ = 0;
    TLorentzVector ZAntiLeptonReco_;

    double minDelRLep  = 0.4;
    double minDelRALep = 0.4;

    // delta R matching from reco particles to gen particles
    for (auto const& itMu : *recoMuons) {
        if(ZLepton_.Pt() && itMu.charge < 0. && itMu.lv.DeltaR(ZLepton_) < minDelRLep){
            minDelRLep = itMu.lv.DeltaR(ZLepton_);
            ZLeptonRecoCat_ = itMu.category;
            ZLeptonReco_ = itMu.lv;

        }
        else if(ZAntiLepton_.Pt() && itMu.charge > 0. && itMu.lv.DeltaR(ZAntiLepton_) < minDelRALep){
            minDelRALep = itMu.lv.DeltaR(ZAntiLepton_);
            ZAntiLeptonRecoCat_ = itMu.category;
            ZAntiLeptonReco_ = itMu.lv;
        }
    }

    ZLeptonPt_ = ZLepton_.Pt();
    ZLeptonEta_ = ZLepton_.Eta();
    ZLeptonPhi_ = ZLepton_.Phi();

    ZAntiLeptonPt_ = ZAntiLepton_.Pt();
    ZAntiLeptonEta_ = ZAntiLepton_.Eta();
    ZAntiLeptonPhi_ = ZAntiLepton_.Phi();

    ZLeptonRecoPt_  = ZLeptonReco_.Pt();
    ZLeptonRecoEta_ = ZLeptonReco_.Eta();
    ZLeptonRecoPhi_ = ZLeptonReco_.Phi();
    ZLeptonRecoDelR_ = minDelRLep;

    ZAntiLeptonRecoPt_  = ZAntiLeptonReco_.Pt();
    ZAntiLeptonRecoEta_ = ZAntiLeptonReco_.Eta();
    ZAntiLeptonRecoPhi_ = ZAntiLeptonReco_.Phi();
    ZAntiLeptonRecoDelR_ = minDelRALep;

    ZMassReco_ = (ZLeptonReco_ + ZAntiLeptonReco_).M();

    edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: end fill Branches";

    return;
}


