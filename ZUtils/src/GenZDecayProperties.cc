#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

#include <iostream>
#include <string>


GenZDecayProperties::GenZDecayProperties()
{
z_ = 0;
meDaughterParticle_ = 0;
stableAntiLepton_ = 0;
meDaughterAntiParticle_ = 0;
stableLepton_ = 0;
decayMode_ = -999;
}


GenZDecayProperties::GenZDecayProperties(const reco::GenParticle* _z,
                               const reco::GenParticle* _meDaughterParticle, const reco::GenParticle* _meDaughterAntiParticle,
                               const reco::GenParticle* _stableLepton, const reco::GenParticle* _stableAntiLepton,
                               const int _decayMode):
z_(_z),
meDaughterParticle_(_meDaughterParticle),
meDaughterAntiParticle_(_meDaughterAntiParticle),
stableLepton_(_stableLepton),
stableAntiLepton_(_stableAntiLepton),
decayMode_(_decayMode)
{}

void print_particle(std::string name, const reco::GenParticle* particle) {
    if(particle) {
        std::cout << name << ": " << particle->pdgId() << ", " << particle->status()
                  << ", " << particle->polarP4().M()
                  << ", " << particle->polarP4().Pt()
                  << ", " << particle->polarP4().Eta()
                  << ", " << particle->polarP4().Phi() << "\n";
    } else {
        std::cout << name << " not found!\n";
    }
}

void
GenZDecayProperties::print()const
{
    // std::cout<<"\n"
    //    <<"--------------------------------------\n"
    //    <<"- Dump GenZDecay Content             -\n"
    //    <<"--------------------------------------\n";
    // std::cout<<"     Particle (PDG ID, Status, mass, pt, eta, phi)\n";
    print_particle("Z", z_);
    print_particle("ME daughter particle", meDaughterParticle_);
    print_particle("ME daughter antiparticle", meDaughterAntiParticle_);
    print_particle("Stable lepton", stableLepton_);
    print_particle("Stable antilepton", stableAntiLepton_);
}

const reco::GenParticle* GenZDecayProperties::z()const{return z_;}
const reco::GenParticle* GenZDecayProperties::meDaughterParticle()const{return meDaughterParticle_;}
const reco::GenParticle* GenZDecayProperties::meDaughterAntiParticle()const{return meDaughterAntiParticle_;}
const reco::GenParticle* GenZDecayProperties::stableLepton()const{return stableLepton_;}
const reco::GenParticle* GenZDecayProperties::stableAntiLepton()const{return stableAntiLepton_;}
int GenZDecayProperties::decayMode()const{return decayMode_;}
