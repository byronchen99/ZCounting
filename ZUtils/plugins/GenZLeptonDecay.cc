// -*- C++ -*-
//
// Package:    TopAnalysis/HiggsUtils
// Class:      GenZLeptonDecay
//
/**\class GenZLeptonDecay GenZLeptonDecay.cc TopAnalysis/HiggsUtils/plugins/GenZLeptonDecay.cc

 Description: Info about generated leptons in Z/gamma* decays using new generator access via GenStatusFlags

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk
//         Created:  Mon, 17 Aug 2015 09:45:34 GMT
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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"



//
// class declaration
//

class GenZLeptonDecay : public edm::one::EDProducer<>
{
public:
    explicit GenZLeptonDecay(const edm::ParameterSet&);
    ~GenZLeptonDecay() = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

    /// input tag for the genParticle source
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

GenZLeptonDecay::GenZLeptonDecay(const edm::ParameterSet& iConfig):
genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("src")))
{
    produces<std::vector<GenZDecayProperties> >();
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenZLeptonDecay::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::unique_ptr<std::vector<GenZDecayProperties> > v_genZDecayProperties(new std::vector<GenZDecayProperties>());

    // Get source collection
    edm::Handle<reco::GenParticleCollection> v_particle;
    iEvent.getByToken(genParticlesToken_, v_particle);

    // Find stable leptons, from prompt or tau decay
    std::vector<const reco::GenParticle*> v_stableLepton;
    for(reco::GenParticleCollection::const_iterator i_particle = v_particle->begin(); i_particle != v_particle->end(); ++i_particle){
        if(i_particle->status() != 1) continue;
        const int absPdgId = std::abs(i_particle->pdgId());
        if(absPdgId<11 || absPdgId>16 || absPdgId%2 == 0) continue;
        if(!i_particle->isPromptFinalState() && !i_particle->isDirectPromptTauDecayProductFinalState()) continue;
        v_stableLepton.push_back(&(*i_particle));
    }

    // Find ME mother leptons
    std::vector<const reco::GenParticle*> v_meDaughter;
    for(std::vector<const reco::GenParticle*>::const_iterator i_particle = v_stableLepton.begin(); i_particle != v_stableLepton.end();){
        const reco::GenParticle* particle = *i_particle;
        while(particle->mother()->pdgId() == particle->pdgId()){
            particle = dynamic_cast<const reco::GenParticle*>(particle->mother());
        }
        if((*i_particle)->isDirectPromptTauDecayProductFinalState()){
            while(std::abs(particle->mother()->pdgId()) == 15){
                particle = dynamic_cast<const reco::GenParticle*>(particle->mother());
            }
        }
        if(std::abs(particle->mother()->pdgId())==24) {
            i_particle = v_stableLepton.erase(i_particle);
        } else {
            v_meDaughter.push_back(&(*particle));
            ++i_particle;
        }
    }

    // Set decay ID and access final output info
    if(v_stableLepton.size() != v_meDaughter.size()){
        edm::LogError("GenZLeptonDecay")<<"Could not find for each final state lepton the initial one";
    }
    else{
        int decayId(0);
        const reco::GenParticle* stableLepton(0);
        const reco::GenParticle* stableAntiLepton(0);
        const reco::GenParticle* meDaughterParticle(0);
        const reco::GenParticle* meDaughterAntiParticle(0);


        if(v_stableLepton.size() == 2){
            if(v_stableLepton.at(0)->pdgId() > 0){
                stableLepton = v_stableLepton.at(0);
                stableAntiLepton = v_stableLepton.at(1);
                meDaughterParticle = v_meDaughter.at(0);
                meDaughterAntiParticle = v_meDaughter.at(1);
            }
            else{
                stableLepton = v_stableLepton.at(1);
                stableAntiLepton = v_stableLepton.at(0);
                meDaughterParticle = v_meDaughter.at(1);
                meDaughterAntiParticle = v_meDaughter.at(0);
            }
            const int meParticlePdgId = meDaughterParticle->pdgId();
            if(meParticlePdgId != 15) {
                decayId = meParticlePdgId;
            }
            else {
                decayId = 150000 + 100*(stableLepton->pdgId()) - stableAntiLepton->pdgId();
            }

        }
        else if(v_stableLepton.size() == 1){
            decayId = 150000;
            if(v_stableLepton.at(0)->pdgId() > 0){
                stableLepton = v_stableLepton.at(0);
                meDaughterParticle = v_meDaughter.at(0);
                decayId += 100*stableLepton->pdgId();
            }
            else{
                stableAntiLepton = v_stableLepton.at(0);
                meDaughterAntiParticle = v_meDaughter.at(0);
                decayId -= stableAntiLepton->pdgId();
            }
        }
        else{
            decayId = 150000;
        }

        if(stableLepton && stableLepton->pdgId()!=11 && stableLepton->pdgId()!=13) stableLepton = 0;
        if(stableAntiLepton && stableAntiLepton->pdgId()!=-11 && stableAntiLepton->pdgId()!=-13) stableAntiLepton = 0;

        // Put all information into the container
        GenZDecayProperties genZDecayProperties(0,
                                                meDaughterParticle, meDaughterAntiParticle,
                                                stableLepton, stableAntiLepton,
                                                decayId
                                               );
        // genZDecayProperties.print();
        v_genZDecayProperties->push_back(genZDecayProperties);
    }

    iEvent.put(std::move(v_genZDecayProperties));
}



// ------------ method called once each job just before starting event loop  ------------
void
GenZLeptonDecay::beginJob()
{}



// ------------ method called once each job just after ending the event loop  ------------
void
GenZLeptonDecay::endJob()
{}



// ------------ method called when starting to processes a run  ------------
/*
void
GenZLeptonDecay::beginRun(edm::Run const&, edm::EventSetup const&)
{}
*/



// ------------ method called when ending the processing of a run  ------------
/*
void
GenZLeptonDecay::endRun(edm::Run const&, edm::EventSetup const&)
{}
*/



// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenZLeptonDecay::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}
*/



// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenZLeptonDecay::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}
*/



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenZLeptonDecay::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(GenZLeptonDecay);
