#include "ZCounting/ZCountAnalyze/interface/zcount_PV.h"

#include <algorithm>

zcount_PV::zcount_PV():zcount_module(){}

zcount_PV::~zcount_PV(){}

void zcount_PV::getInput(const edm::ParameterSet& iConfig){
    VtxNTracksFitCut_ = iConfig.getUntrackedParameter<double>("VtxNTracksFitMin");
    VtxNdofCut_       = iConfig.getUntrackedParameter<double>("VtxNdofMin");
    VtxAbsZCut_       = iConfig.getUntrackedParameter<double>("VtxAbsZMax");
    VtxRhoCut_        = iConfig.getUntrackedParameter<double>("VtxRhoMax");
}

void zcount_PV::initBranches(TTree* tree){
    addBranch(tree,"nPV", &nPV_,"nPV_/f");
}

bool zcount_PV::readEvent(const edm::Event& iEvent){
    iEvent.getByToken(fPVName_token, hVertexProduct);

    if (!hVertexProduct.isValid()) {
        edm::LogWarning("ZCounting") << "No valid primary vertex product found" << std::endl;
        return false;
    }

    //const reco::VertexCollection* pvCol = hVertexProduct.product();
    //const reco::Vertex* pv = &(*pvCol->begin());
    nPV_ = 0;
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
        //if (nPV_ == 0) 
        //    pv = &itVtx;
        
        nPV_++;
    }
    if (nPV_ == 0)
        return false;

    edm::LogInfo("zcount_PV") << "good primary vertex ";

    return true;
}


void zcount_PV::fillBranches(){
    return;
}

