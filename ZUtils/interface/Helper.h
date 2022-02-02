#ifndef HELPER_H_
#define HELPER_H_

#include "DataFormats/Common/interface/AssociativeIterator.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

int doMatch(const edm::Handle<std::vector<reco::GenParticle>>& GenParticles, 
    float pt, float eta, float phi, int pdgId, 
    float &pt_gen, float &eta_gen, float &phi_gen, int &pdgId_gen, 
    int &isPromptFinalState_gen, int &isDirectPromptTauDecayProductFinalState_gen);

bool isTau(const reco::GenParticle* lepton);
const reco::GenParticle* tauDaughter(const reco::GenParticle* tau);
    
float getPFIso(const reco::Muon& muon);
float getTkIso(const reco::Muon& muon);
bool isCustomTightMuon(const reco::Muon& muon);

bool isGoodPV(const reco::Vertex &vtx);
bool isValidTrack(const reco::Track &trk);    

#endif
