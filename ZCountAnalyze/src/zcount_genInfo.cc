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
    addBranch(tree,"ZDecayMode", &ZDecayMode_,"ZDecayMode_/i");

    addBranch(tree,"ZLeptonPt", &ZLeptonPt_,"ZLeptonPt_/f");
    addBranch(tree,"ZLeptonEta", &ZLeptonEta_,"ZLeptonEta_/f");
    addBranch(tree,"ZLeptonPhi", &ZLeptonPhi_,"ZLeptonPhi_/f");

    addBranch(tree,"ZAntiLeptonPt", &ZAntiLeptonPt_,"ZAntiLeptonPt_/f");
    addBranch(tree,"ZAntiLeptonEta", &ZAntiLeptonEta_,"ZAntiLeptonEta_/f");
    addBranch(tree,"ZAntiLeptonPhi", &ZAntiLeptonPhi_,"ZAntiLeptonPhi_/f");

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

    if(hGenZInfoProduct->at(0).stableLepton()){
        ZLeptonPt_        = hGenZInfoProduct->at(0).stableLepton()->pt();
        ZLeptonEta_       = hGenZInfoProduct->at(0).stableLepton()->eta();
        ZLeptonPhi_       = hGenZInfoProduct->at(0).stableLepton()->phi();
    }
    else{
        edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: no gen stable lepton found";
        ZLeptonPt_  = 0;
        ZLeptonEta_ = 0;
        ZLeptonPhi_ = 0;
    }
    if(hGenZInfoProduct->at(0).stableAntiLepton()){
        ZAntiLeptonPt_        = hGenZInfoProduct->at(0).stableAntiLepton()->pt();
        ZAntiLeptonEta_       = hGenZInfoProduct->at(0).stableAntiLepton()->eta();
        ZAntiLeptonPhi_       = hGenZInfoProduct->at(0).stableAntiLepton()->phi();
    }
    else{
        edm::LogVerbatim("zcount_genInfo") << "zcount_genZInfo: no gen stable antilepton found";
        ZAntiLeptonPt_  = 0;
        ZAntiLeptonEta_ = 0;
        ZAntiLeptonPhi_ = 0;
    }

    if(hGenZInfoProduct->at(0).stableLepton() && hGenZInfoProduct->at(0).stableAntiLepton()){
        ZStableMass_ = (hGenZInfoProduct->at(0).stableLepton()->polarP4() + hGenZInfoProduct->at(0).stableAntiLepton()->polarP4()).mass();
    }
    else{
        ZStableMass_ = 0;
    }


    eventWeight_ = hGenInfoProduct->weight();

    return true;
}


void zcount_genInfo::fillBranches(){

    return;
}

