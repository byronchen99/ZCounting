import FWCore.ParameterSet.Config as cms

genZLeptonDecay = cms.EDProducer('GenZLeptonDecay',
    ## input particle collection of type edm::View<reco::GenParticle>
    src = cms.InputTag("genParticles")
)
