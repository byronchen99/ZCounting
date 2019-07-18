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
#include "ZCounting/ZCountAnalyze/interface/LorentzVector.h"

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

    void setGenInfoToken(edm::EDGetTokenT<GenEventInfoProduct> genInfoToken){
        fGenInfoName_token = genInfoToken;
    }
    void setGenZInfoToken(edm::EDGetTokenT<std::vector<GenZDecayProperties>> genZInfoToken){
        fGenZInfoName_token = genZInfoToken;
    }
    void setGenZLepInfoToken(edm::EDGetTokenT<std::vector<GenZDecayProperties>> genZLepInfoToken){
        fGenZLepInfoName_token = genZLepInfoToken;
    }
    void setPVToken(edm::EDGetTokenT<reco::VertexCollection> pvToken){
        fPVName_token = pvToken;
    }
    void setMuonToken(edm::EDGetTokenT<reco::MuonCollection> muonToken){
        fMuonName_token = muonToken;
    }
    void setTrackToken(edm::EDGetTokenT<reco::TrackCollection> trackToken){
        fTrackName_token = trackToken;
    }

    void setPVModule(zcount_PV &module){
        pvModule = &module;
    }
    void setMuonModule(zcount_muons &module){
        muonModule = &module;
    }

private:

    edm::EDGetTokenT<GenEventInfoProduct>              fGenInfoName_token;
    edm::EDGetTokenT<std::vector<GenZDecayProperties>> fGenZInfoName_token;
    edm::EDGetTokenT<std::vector<GenZDecayProperties>> fGenZLepInfoName_token;
    edm::EDGetTokenT<reco::VertexCollection>           fPVName_token;
    edm::EDGetTokenT<reco::MuonCollection>             fMuonName_token;
    edm::EDGetTokenT<reco::TrackCollection>            fTrackName_token;

    zcount_PV    *pvModule;
    zcount_muons *muonModule;

    const reco::GenParticle *ZGen;
    const reco::GenParticle *ZLepton;
    const reco::GenParticle *ZAntiLepton;

    const reco::Muon  *ZLeptonRecoMuon;
    const reco::Muon  *ZAntiLeptonRecoMuon;
    const reco::Track *ZLeptonRecoTrack;
    const reco::Track *ZAntiLeptonRecoTrack;

    const reco::Vertex* ZRecoVtx;
    const reco::Vertex* ZRecoGoodVtx;

    //----- input

    //----- output
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
    int ZLeptonRecoCat_[2];
    float ZLeptonRecoPt_;
    float ZLeptonRecoEta_;
    float ZLeptonRecoPhi_;

    int ZAntiLeptonRecoCat_[2];
    float ZAntiLeptonRecoPt_;
    float ZAntiLeptonRecoEta_;
    float ZAntiLeptonRecoPhi_;

    float ZMassReco_;

    // primary vertex from Z
    float ZVtxX_;
    float ZVtxY_;
    float ZVtxZ_;
    float ZVtxRho_;
    float ZVtxNDoF_;

    bool matchRecoVtx_;     // if reco vertex is matched with closest to Z vertex of reco vertex collection
    bool matchRecoGoodVtx_; // if reco vertex is matched with closest to Z vertex of good reco vertex collection

};

#endif

