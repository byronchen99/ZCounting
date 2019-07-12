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
#include "ZCounting/ZCountAnalyze/interface/zcount_muons.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

#include <TLorentzVector.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;


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
    void setGenZLepInfoToken(edm::EDGetTokenT<std::vector<GenZDecayProperties>> genZLepInfoToken) {
        fGenZLepInfoName_token = genZLepInfoToken;
    }

    void setMuonModule(zcount_muons &module){
        muonModule = &module;
    }


private:

    edm::EDGetTokenT<GenEventInfoProduct> fGenInfoName_token;
    edm::EDGetTokenT<std::vector<GenZDecayProperties>> fGenZInfoName_token;
    edm::EDGetTokenT<std::vector<GenZDecayProperties>> fGenZLepInfoName_token;

    zcount_muons *muonModule;

    TLorentzVector ZLepton_;
    TLorentzVector ZAntiLepton_;
    const LV nullP4_{0,0,0,0};
    // input
    
    // output
    float eventWeight_;

    // Gen particles
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

    // Associated Reco particles
    int ZLeptonRecoCat_;
    float ZLeptonRecoPt_;
    float ZLeptonRecoEta_;
    float ZLeptonRecoPhi_;
    float ZLeptonRecoDelR_;

    int ZAntiLeptonRecoCat_;
    float ZAntiLeptonRecoPt_;
    float ZAntiLeptonRecoEta_;
    float ZAntiLeptonRecoPhi_;
    float ZAntiLeptonRecoDelR_;


    float ZMassReco_;

};

#endif

