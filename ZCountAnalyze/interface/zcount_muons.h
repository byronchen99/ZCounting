/**
 * *
 *    
 * Created on: 29 June 2019
 *   Author: David Walter
 *
 *  MuonProbeCategory {0,5}
 *    5: the probe passes the HLT trigger
 *    4: the probe passes the selection, but not the HLT trigger
 *    3: the probe passes is a global muon, but fails the selection
 *    2: the probe is a standalone muon, this means the (possible) muon does not match to a tracker track
 *    1: the probe is a tracker muon or a tracker track, this means that the (possible) muon does not match to a standalone muon
 *    0: the probe does not pass the lowest criteria - it should not be stored with this category
 */

#ifndef ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_MUONS_H_
#define ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_MUONS_H_

#include "ZCounting/ZCountAnalyze/interface/zcount_module.h"
#include "ZCounting/ZCountAnalyze/interface/zcount_PV.h"
#include "ZCounting/ZCountAnalyze/interface/zcount_trigger.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <TLorentzVector.h>


class zcount_muons: public zcount_module{
public:
    zcount_muons();
    ~zcount_muons();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    bool readEvent(const edm::Event& iEvent);

    void fillBranches();

    void setMuonToken(edm::EDGetTokenT<reco::MuonCollection> muonToken){
        fMuonName_token = muonToken;
    }
    void setTrackToken(edm::EDGetTokenT<reco::TrackCollection> trackToken){
        fTrackName_token = trackToken;
    }

    void setPVModule(zcount_PV &module){
        pvModule = &module;
    }
    void setTriggerModule(zcount_trigger &module){
        triggerModule = &module;
    }


    // functions
    int getMuonCategory(const reco::Muon &muon, const reco::Vertex &vtx);
    bool isValidTrack(const reco::Track &trk);

    bool isMuonTrigger();
    bool isMuonTriggerObj(const double eta, const double phi);

    bool passMuonID(const reco::Muon& muon, const reco::Vertex& vtx);
    bool passMuonIso(const reco::Muon& muon);

private:
   
    edm::EDGetTokenT<reco::MuonCollection>  fMuonName_token;
    edm::EDGetTokenT<reco::TrackCollection> fTrackName_token;

    edm::Handle<reco::MuonCollection>  hMuonProduct;
    edm::Handle<reco::TrackCollection> hTrackProduct;

    zcount_PV      *pvModule;
    zcount_trigger *triggerModule;

    static constexpr std::size_t max_num_tags   = 5;
    static constexpr std::size_t max_num_probes = 20;

    // input
    double PtCutTag_;
    double PtCutProbe_;
    double EtaCutTag_;
    double EtaCutProbe_;

    std::string IDTypestr_;
    std::string IsoTypestr_;
    MuonIDTypes IDType_{NoneID};
    MuonIsoTypes IsoType_{NoneIso};
    double IsoCut_;

    std::vector<std::string> fMuonHLTNames;
    std::vector<std::string> fMuonHLTObjectNames;

    double MassMin_;
    double MassMax_;

    // output
    int   nTag_;
    float tagPt_[max_num_tags];
    float tagEta_[max_num_tags];
    float tagPhi_[max_num_tags];

    int   nProbe_;
    int   probeCat_[max_num_probes];    // probe muon category
    float probePt_[max_num_probes];
    float probeEta_[max_num_probes];
    float probePhi_[max_num_probes];
    float dilepMass_[max_num_probes];   // tag-and-probe mass
    int   tagIndex_[max_num_probes];    // index of the corresponding tag muon

};

#endif

