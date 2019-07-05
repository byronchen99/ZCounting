import FWCore.ParameterSet.Config as cms

genZDecay = cms.EDProducer('GenZDecay',
    ## input particle collection of type edm::View<reco::GenParticle>
    src = cms.InputTag("genParticles")
)