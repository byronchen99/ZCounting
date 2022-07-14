#ifndef GenZDecayProperties_h
#define GenZDecayProperties_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class GenZDecayProperties{

public:

    /// Empty constructor
    GenZDecayProperties();
    /// Default constructor
    GenZDecayProperties(const reco::GenParticle*,
        const reco::GenParticle*, const reco::GenParticle*,
        const reco::GenParticle*, const reco::GenParticle*,
        const int
    );
    /// Default destructor
    virtual ~GenZDecayProperties(){};

    const reco::GenParticle* z()const;
    const reco::GenParticle* meDaughterParticle()const;
    const reco::GenParticle* meDaughterAntiParticle()const;
    const reco::GenParticle* stableLepton()const;
    const reco::GenParticle* stableAntiLepton()const;

    int decayMode()const;

    /// Print content of the Z decay
    void print()const;



private:

    /// Pointer to the ME-level Z
    const reco::GenParticle* z_;
    /// Pointer to the ME-level decay particle of the Z
    /// (always existent independent of decay mode)
    const reco::GenParticle* meDaughterParticle_;
    /// Pointer to the ME-level decay anti-particle of the Z
    /// (always existent independent of decay mode)
    const reco::GenParticle* meDaughterAntiParticle_;
    /// Pointer to the final stable charged lepton from Z->ee or Z->mumu or Z->tautau
    /// (in case it exists, i.e. the decay mode has leptons in the final state)
    const reco::GenParticle* stableLepton_;
    /// Pointer to the final stable charged anti-lepton from Z->ee or Z->mumu or Z->tautau
    /// (in case it exists, i.e. the decay mode has leptons in the final state)
    const reco::GenParticle* stableAntiLepton_;

    /// Number specifying the decay mode (desciption see above)
    int decayMode_;
};


#endif
