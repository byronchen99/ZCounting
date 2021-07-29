#include "TMath.h"

#ifndef HELPER_H_
#define HELPER_H_

float GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
   float DeltaPhi = TMath::Abs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}


int doMatch(const edm::Handle<std::vector<reco::GenParticle>>& GenParticles, float pt, float eta, float phi, int pdgId, float &pt_gen, float &eta_gen, float &phi_gen, int &pdgId_gen, int &isPromptFinalState_gen, int &isDirectPromptTauDecayProductFinalState_gen)
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

		float dr = GetDeltaR(eta, phi, etaGen, phiGen);
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

#endif
