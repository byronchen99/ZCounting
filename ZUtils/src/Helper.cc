#include "ZCounting/ZUtils/interface/Helper.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TMath.h"

//--------------------------------------------------------------------------------------------------
int doMatch(const edm::Handle<std::vector<reco::GenParticle>>& GenParticles, 
    float pt, float eta, float phi, int pdgId, 
    float &pt_gen, float &eta_gen, float &phi_gen, int &pdgId_gen, 
    int &isPromptFinalState_gen, int &isDirectPromptTauDecayProductFinalState_gen)
{
	int match_value = 0;
    // 0 = no MC match found
    // 1 = prompt lepton match found
    // 2 = also match charge 
    // 3 = GammaConv match found

	reco::GenParticleCollection genParticlesCollection = *GenParticles;
	reco::GenParticleCollection::const_iterator genParticleSrc;

	float drmin = 0.3;
	float minDPtRel = 0.5;

	int ipart = -1;
	for(genParticleSrc = genParticlesCollection.begin(); genParticleSrc != genParticlesCollection.end(); genParticleSrc++)
	{
		ipart++;

		reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

		float ptGen = mcp->pt();
		float etaGen = mcp->eta();
		float phiGen = mcp->phi();
		int idGen = mcp->pdgId();
		int statusGen = mcp->status();
		int isPromptFinalState = mcp->isPromptFinalState(); //final state (status 1) particle satisfying isPrompt() = not from hadron, muon or tau decay
		int isDirectPromptTauDecayProductFinalState = mcp->isDirectPromptTauDecayProductFinalState(); //final state (status 1) particle satisfying isDirectPromptTauDecayProduct() = direct decay product from a tau decay (ie no intermediate hadron), where the tau did not come from a hadron decay

		if(statusGen != 1) {continue;} //For ele and muons, ask particle to be stable (status=1)

		if((fabs(pt - ptGen) / ptGen) > minDPtRel) {continue;} //Pt requirement

		float dr = reco::deltaR(eta, phi, etaGen, phiGen);
		if(dr > drmin) {continue;} //Look for gen particle with smallest dR

		//-- Current best match found !
		if(abs(pdgId) == 11 && abs(idGen) == 22) //Gamma conversion matching
		{
			drmin = dr;
			match_value = 3;

			pt_gen = ptGen;
			eta_gen = etaGen;
			phi_gen = phiGen;
			pdgId_gen = idGen;
		}
		else if(abs(pdgId) == abs(idGen)) //Lepton matching
		{
			drmin = dr;
			match_value = 1;
			if(pdgId == idGen) {match_value = 2;} //Charge matched

			pt_gen = ptGen;
			eta_gen = etaGen;
			phi_gen = phiGen;
			pdgId_gen = idGen;
			isPromptFinalState_gen = isPromptFinalState;
			isDirectPromptTauDecayProductFinalState_gen = isDirectPromptTauDecayProductFinalState;
		}
	}

	return match_value;
}

//--------------------------------------------------------------------------------------------------
bool isTau(const reco::GenParticle* lepton)
{
    return std::abs(lepton->pdgId()) == 15;
}

//--------------------------------------------------------------------------------------------------
const reco::GenParticle* tauDaughter(const reco::GenParticle* tau)
{
    for(size_t iDaughter = 0; iDaughter < tau->numberOfDaughters(); ++iDaughter){
        const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(tau->daughter(iDaughter));
        if(std::abs(daughter->pdgId())==11 || std::abs(daughter->pdgId())==13) return daughter;
        else if(isTau(daughter)) return tauDaughter(daughter);
    }
    return tau;
}

//--------------------------------------------------------------------------------------------------
bool isGoodPV(const reco::Vertex &vtx){
    if (vtx.isFake())
        return false;
    if (vtx.tracksSize() < 0.)
        return false;
    if (vtx.ndof() < 4.)
        return false;
    if (fabs(vtx.z()) > 24.)
        return false;
    if (vtx.position().Rho() > 2.)
        return false;

    return true;
}

//--------------------------------------------------------------------------------------------------
bool isValidTrack(const reco::Track &trk){
    if(trk.hitPattern().trackerLayersWithMeasurement() >= 6 && trk.hitPattern().numberOfValidPixelHits() >= 1)
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
bool isCustomTightMuon(const reco::Muon& muon){
    return muon.isGlobalMuon()
        && muon.isPFMuon()
        && muon.globalTrack()->normalizedChi2() < 10.
        && muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
        && muon.numberOfMatchedStations() > 1
        && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
        && muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
}

//--------------------------------------------------------------------------------------------------
float getPFIso(const reco::Muon& muon) {
    return (muon.pfIsolationR04().sumChargedHadronPt +
                   std::max(0.,
                            muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt -
                                0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();
}

//--------------------------------------------------------------------------------------------------
float getTkIso(const reco::Muon& muon) {
    return muon.isolationR03().sumPt / muon.pt();
}
