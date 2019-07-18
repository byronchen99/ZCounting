#include "ZCounting/ZCountAnalyze/interface/zcount_PV.h"

#include <algorithm>

zcount_PV::zcount_PV():zcount_module(){}

zcount_PV::~zcount_PV(){}

void zcount_PV::getInput(const edm::ParameterSet& iConfig){
    edm::LogVerbatim("zcount_PV") << " getInput" ;

    VtxNTracksFitCut_ = iConfig.getUntrackedParameter<double>("VtxNTracksFitMin");
    VtxNdofCut_       = iConfig.getUntrackedParameter<double>("VtxNdofMin");
    VtxAbsZCut_       = iConfig.getUntrackedParameter<double>("VtxAbsZMax");
    VtxRhoCut_        = iConfig.getUntrackedParameter<double>("VtxRhoMax");
}

void zcount_PV::initBranches(TTree* tree){
    addBranch(tree,"nPV", &nPV_,"nPV_/f");

    addBranch(tree,"VtxX",     &VtxX_,    "VtxX_/f");
    addBranch(tree,"VtxY",     &VtxY_,    "VtxY_/f");
    addBranch(tree,"VtxZ",     &VtxZ_,    "VtxZ_/f");
    addBranch(tree,"VtxRho",   &VtxRho_,  "VtxRho_/f");
    addBranch(tree,"VtxNDoF",  &VtxNDoF_, "VtxNDoF_/f");
}

bool zcount_PV::readEvent(const edm::Event& iEvent){
    edm::LogVerbatim("zcount_PV") << "zcount_PV: readEvent" ;
    iEvent.getByToken(fPVName_token, hVertexProduct);

    if (!hVertexProduct.isValid()) {
        edm::LogWarning("zcount_PV") << "No valid primary vertex product found" ;
        return false;
    }

    const reco::VertexCollection* pvCol = hVertexProduct.product();
    pv = &(*pvCol->begin());
    nPV_ = 0;
    for (auto const& itVtx : *hVertexProduct) {

        if(!isGoodPV(itVtx))
            continue;

        if (nPV_ == 0)
            pv = &itVtx;

        nPV_++;
    }
    if (nPV_ == 0)
        return false;

    edm::LogVerbatim("zcount_PV") << "zcount_PV: good primary vertex ";

    return true;
}


void zcount_PV::fillBranches(){

    VtxX_    = pv->x();
    VtxY_    = pv->y();
    VtxZ_    = pv->z();
    VtxRho_  = pv->position().Rho();
    VtxNDoF_ = pv->ndof();

    return;
}

//--------------------------------------------------------------------------------------------------
bool zcount_PV::isGoodPV(const reco::Vertex &vtx){
    if (vtx.isFake())
        return false;
    if (vtx.tracksSize() < VtxNTracksFitCut_)
        return false;
    if (vtx.ndof() < VtxNdofCut_)
        return false;
    if (fabs(vtx.z()) > VtxAbsZCut_)
        return false;
    if (vtx.position().Rho() > VtxRhoCut_)
        return false;

    return true;
}

