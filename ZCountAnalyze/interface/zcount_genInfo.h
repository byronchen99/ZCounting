/**
 * *
 *    
 * Created on: 04 July 2019
 *   Author: David Walter
 *  
 */

#ifndef ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_GENINFO_H_
#define ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_GENINFO_H_

#include "ZCounting/ZCountAnalyze/interface/zcount_module.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

class zcount_genInfo: public zcount_module{
public:
    zcount_genInfo();
    ~zcount_genInfo();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    bool readEvent(const edm::Event& iEvent);

    void fillBranches();

    void setGenInfoToken(edm::EDGetTokenT<GenEventInfoProduct> genInfoToken) {
        fGenInfoName_token = genInfoToken;
    }
    void setGenZInfoToken(edm::EDGetTokenT<std::vector<GenZDecayProperties>> genZInfoToken) {
        fGenZInfoName_token = genZInfoToken;
    }


private:

    edm::EDGetTokenT<GenEventInfoProduct> fGenInfoName_token;
    edm::EDGetTokenT<std::vector<GenZDecayProperties>> fGenZInfoName_token;

    // input
    
    // output
    float eventWeight_;

    float ZPt_;
    float ZEta_;
    float ZPhi_;
    float ZM_;
    float ZStableMass_;
    int ZDecayMode_;

    float ZLeptonPt_;
    float ZLeptonEta_;
    float ZLeptonPhi_;

    float ZAntiLeptonPt_;
    float ZAntiLeptonEta_;
    float ZAntiLeptonPhi_;

};

#endif

