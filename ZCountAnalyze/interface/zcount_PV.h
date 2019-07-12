/**
 * * Primary Vertex selection and content handling
 *    
 * Created on: 29 June 2019
 *   Author: David Walter
 *  
 */

#ifndef ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_PV_H_
#define ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_PV_H_

#include "ZCounting/ZCountAnalyze/interface/zcount_module.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class zcount_PV: public zcount_module{
public:
    zcount_PV();
    ~zcount_PV();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    bool readEvent(const edm::Event& iEvent);

    void fillBranches();

    void setPVToken(edm::EDGetTokenT<reco::VertexCollection> pvToken) {
        fPVName_token = pvToken;
    }

    const reco::Vertex* getPV()const{
        return pv;
    }

private:

    const reco::Vertex* pv;

    edm::EDGetTokenT<reco::VertexCollection> fPVName_token;

    edm::Handle<reco::VertexCollection> hVertexProduct;

    // input 
    double VtxNTracksFitCut_;
    double VtxNdofCut_;
    double VtxAbsZCut_;
    double VtxRhoCut_;
    
    // output
    float nPV_;
};

#endif

