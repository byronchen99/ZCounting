#ifndef ZCOUNTING_TNPPAIRTREEPRODUCER_INTERFACE_MUONDEFINITIONS_H
#define ZCOUNTING_TNPPAIRTREEPRODUCER_INTERFACE_MUONDEFINITIONS_H

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

class MuonDefinitions{
public:
    MuonDefinitions(){};
    ~MuonDefinitions(){};

protected:
    bool passBaseline(const reco::Muon& muon);
    float getPFIso(const reco::Muon& muon);
    float getTkIso(const reco::Muon& muon);
    float getDxy(const reco::Muon& muon, const reco::Vertex& vtx);
    float getDz(const reco::Muon& muon, const reco::Vertex& vtx);

};

//--------------------------------------------------------------------------------------------------
bool  MuonDefinitions::passBaseline(const reco::Muon& muon){
    return muon.isGlobalMuon()
                && muon.isPFMuon()
                && muon.globalTrack()->normalizedChi2() < 10.
                && muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
                && muon.numberOfMatchedStations() > 1
                && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
                && muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
}


//--------------------------------------------------------------------------------------------------
float  MuonDefinitions::getPFIso(const reco::Muon& muon) {
    return (muon.pfIsolationR04().sumChargedHadronPt +
                   std::max(0.,
                            muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt -
                                0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();
}

//--------------------------------------------------------------------------------------------------
float  MuonDefinitions::getTkIso(const reco::Muon& muon) {
    return muon.isolationR03().sumPt / muon.pt();
}

//--------------------------------------------------------------------------------------------------
float  MuonDefinitions::getDxy(const reco::Muon& muon, const reco::Vertex& vtx) {
    return fabs(muon.muonBestTrack()->dxy(vtx.position()));
}

//--------------------------------------------------------------------------------------------------
float MuonDefinitions::getDz(const reco::Muon& muon, const reco::Vertex& vtx) {
    return fabs(muon.muonBestTrack()->dz(vtx.position()));
}

#endif
