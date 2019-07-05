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
    else if(hGenZInfoProduct->size() != 1){
        edm::LogWarning("zcount_genInfo") << "zcount_genZInfo: no or more then one gen Z found";
        return false;
    }
    else{
        ZPt_        = hGenZInfoProduct->at(0).z()->pt();
        ZEta_       = hGenZInfoProduct->at(0).z()->eta();
        ZPhi_       = hGenZInfoProduct->at(0).z()->phi();
        ZM_         = hGenZInfoProduct->at(0).z()->mass();
        ZDecayMode_ = hGenZInfoProduct->at(0).decayMode();

        ZLeptonPt_        = hGenZInfoProduct->at(0).stableLepton()->pt();
        ZLeptonEta_       = hGenZInfoProduct->at(0).stableLepton()->eta();
        ZLeptonPhi_       = hGenZInfoProduct->at(0).stableLepton()->phi();

        ZAntiLeptonPt_        = hGenZInfoProduct->at(0).stableAntiLepton()->pt();
        ZAntiLeptonEta_       = hGenZInfoProduct->at(0).stableAntiLepton()->eta();
        ZAntiLeptonPhi_       = hGenZInfoProduct->at(0).stableAntiLepton()->phi();

        ZStableMass_ = (hGenZInfoProduct->at(0).stableLepton()->polarP4() + hGenZInfoProduct->at(0).stableAntiLepton()->polarP4()).mass();

    }
    eventWeight_ = hGenInfoProduct->weight();

    return true;
}


void zcount_genInfo::fillBranches(){

    return;
}

