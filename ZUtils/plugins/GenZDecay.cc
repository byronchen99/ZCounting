// -*- C++ -*-
//
// Package:    GenZDecay
// Class:      GenZDecay
//
/**\class GenZDecay GenZDecay.cc ZCounting/ZUtils/plugins/GenZDecay.cc

 Description: Producer for information of generator-level Z bosons and their decay products

 Implementation:
     The decayMode of the Z is assigned as follows:
     + For each particle x of the decay Z->xxbar the PDG ID of x is assigned, except for the tau,
       i.e. these are values like 1,2,3,4,5,11,12,13,14,16
     + For decays to taus, the further decay is identified for leptonic decays, for hadronic decays it is set to 0,
       i.e. the number is 15aabb, where aa is for the tau (11, 13, 00) and bb is for the antitau (11, 13, 00)
*/
//
// Original Author:  Johannes Hauk,01a/O2.115,3139,
//         Created:  Mon Oct  7 16:09:18 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

//
// class declaration
//

class GenZDecay : public edm::one::EDProducer<>
{
public:
    explicit GenZDecay(const edm::ParameterSet&);
    ~GenZDecay() = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void endRun(edm::Run&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

    /// Find Z's in list of input particles
    std::vector<const reco::GenParticle*> findZs(const reco::GenParticleCollection&)const;

    /// Find status3 (Matrix-Element) daughters of Zs
    std::pair<const reco::GenParticle*, const reco::GenParticle*> findMeDaughters(const reco::GenParticle* z)const;

    /// Find status1 (stable) lepton daughters, in case of leptonic decays
    std::pair<const reco::GenParticle*, const reco::GenParticle*> findStableLeptons(std::pair<const reco::GenParticle*, const reco::GenParticle*>& meDaughters)const;

    /// Find stable (anti)lepton from a ME (anti)lepton (for (anti)taus only if it exists)
    const reco::GenParticle* findStableLepton(const reco::GenParticle* particle)const;

    /// Assign IDs to identify Z decay modes
    int decayId();

    // ----------member data ---------------------------

    /// input token for the genParticle source
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
GenZDecay::GenZDecay(const edm::ParameterSet& iConfig):
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("src")))
{
    produces<std::vector<GenZDecayProperties> >();
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenZDecay::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::unique_ptr<std::vector<GenZDecayProperties> > v_genZDecayProperties(new std::vector<GenZDecayProperties>());

    // Get source collection
    edm::Handle<reco::GenParticleCollection> src;
    iEvent.getByToken(genParticlesToken_, src);

    // Find Z's in list of input particles
    std::vector<const reco::GenParticle*> v_z = this->findZs(*src);

    // Loop over Z's and get decay particles
    for(std::vector<const reco::GenParticle*>::const_iterator i_z = v_z.begin(); i_z != v_z.end(); ++i_z){

        // Find status3 (matrix element) daughters of Z
        std::pair<const reco::GenParticle*, const reco::GenParticle*> meDaughters = this->findMeDaughters(*i_z);

        // Find status1 (stable) daughters in case of e, mu, tau --> electrons or muons
        std::pair<const reco::GenParticle*, const reco::GenParticle*> stableLeptons = this->findStableLeptons(meDaughters);

        // Assign the ID for identifying the type of Z decay
        int decayId(0);
        const int meParticlePdgId = meDaughters.first->pdgId();
        //const int meAntiParticlePdgId = meDaughters.second->pdgId();
        const int stableLeptonPdgId = stableLeptons.first ? stableLeptons.first->pdgId() : 0;
        const int stableAntiLeptonPdgId = stableLeptons.second ? stableLeptons.second->pdgId() : 0;
        if(meParticlePdgId != 15) decayId = meParticlePdgId;
        else decayId = 150000 + 100*stableLeptonPdgId - stableAntiLeptonPdgId;

        // Put all information into the container
        GenZDecayProperties genZDecayProperties(*i_z,
                                                meDaughters.first, meDaughters.second,
                                                stableLeptons.first, stableLeptons.second,
                                                decayId
                                               );
        //genZDecayProperties.print();
        v_genZDecayProperties->push_back(genZDecayProperties);
    }

    iEvent.put(std::move(v_genZDecayProperties));

}



std::vector<const reco::GenParticle*> GenZDecay::findZs(const reco::GenParticleCollection& v_particle)const
{
    std::vector<const reco::GenParticle*> v_target;

    for(reco::GenParticleCollection::const_iterator i_particle=v_particle.begin(); i_particle!=v_particle.end(); ++i_particle){
        if(i_particle->pdgId()==23 && (i_particle->status()==3 || i_particle->status()==62)){
            v_target.push_back(&(*i_particle));
        }
    }

    return v_target;
}



std::pair<const reco::GenParticle*, const reco::GenParticle*> GenZDecay::findMeDaughters(const reco::GenParticle* z)const
{
    std::pair<const reco::GenParticle*, const reco::GenParticle*> result;

    // Call Z the mother, check for its daughters
    const reco::GenParticle* mother = z;
    const reco::GenParticle* daughterParticle(0);
    const reco::GenParticle* daughterAntiParticle(0);
    for(size_t iDaughter = 0; iDaughter < mother->numberOfDaughters(); ++iDaughter){
        const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(mother->daughter(iDaughter));
        if(daughter->pdgId()==23) continue;
        if(daughter->pdgId()>0) daughterParticle = daughter;
        else daughterAntiParticle = daughter;
    }
    result = std::make_pair(daughterParticle, daughterAntiParticle);

    return result;
}



std::pair<const reco::GenParticle*, const reco::GenParticle*> GenZDecay::findStableLeptons(std::pair<const reco::GenParticle*, const reco::GenParticle*>& meDaughters)const
{
    std::pair<const reco::GenParticle*, const reco::GenParticle*> result;

    const reco::GenParticle* lepton = this->findStableLepton(meDaughters.first);
    const reco::GenParticle* antiLepton = this->findStableLepton(meDaughters.second);
    result = std::make_pair(lepton, antiLepton);

    return result;
}



const reco::GenParticle* GenZDecay::findStableLepton(const reco::GenParticle* particle)const
{
    // Do only for charged leptons
    while(std::abs(particle->pdgId())==11 || std::abs(particle->pdgId())==13 || std::abs(particle->pdgId())==15){
        // Is it the final stable particle
        if(particle->status() == 1){
            return particle;
        }

        bool foundLeptonDaughter(false);
        for(size_t i = 0; i < particle->numberOfDaughters(); ++i){
            const reco::GenParticle* daugther = dynamic_cast<const reco::GenParticle*>(particle->daughter(i));
            const int absPdgId = std::abs(daugther->pdgId());
            if(absPdgId==11 || absPdgId==13 || absPdgId==15){
                particle = daugther;
                foundLeptonDaughter = true;
                break;
            }
        }
        if(!foundLeptonDaughter) return 0;
    }

    return 0;
}



// ------------ method called once each job just before starting event loop  ------------
void
GenZDecay::beginJob()
{
}



// ------------ method called once each job just after ending the event loop  ------------
void
GenZDecay::endJob()
{
}



// ------------ method called when starting to processes a run  ------------
void
GenZDecay::beginRun(edm::Run&, edm::EventSetup const&)
{
}



// ------------ method called when ending the processing of a run  ------------
void
GenZDecay::endRun(edm::Run&, edm::EventSetup const&)
{
}



// ------------ method called when starting to processes a luminosity block  ------------
void
GenZDecay::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}



// ------------ method called when ending the processing of a luminosity block  ------------
void
GenZDecay::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenZDecay::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(GenZDecay);
