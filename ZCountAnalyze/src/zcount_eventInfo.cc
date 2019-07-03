#include "ZCounting/ZCountAnalyze/interface/zcount_eventInfo.h"

#include <algorithm>

zcount_eventInfo::zcount_eventInfo():zcount_module(){}

zcount_eventInfo::~zcount_eventInfo(){}

void zcount_eventInfo::getInput(const edm::ParameterSet& iConfig){

}

void zcount_eventInfo::initBranches(TTree* tree){
    addBranch(tree,"runNumber", &runNumber_,"runNumber_/i");
    addBranch(tree,"lumiSection", &lumiSection_,"lumiSection_/i");
    addBranch(tree,"eventNumber", &eventNumber_,"eventNumber_/i");
}

bool zcount_eventInfo::readEvent(const edm::Event& iEvent){

    runNumber_   = iEvent.id().run();
    lumiSection_ = iEvent.id().luminosityBlock();
    eventNumber_ = iEvent.id().event();

    return true;
}


void zcount_eventInfo::fillBranches(){
    return;
}

