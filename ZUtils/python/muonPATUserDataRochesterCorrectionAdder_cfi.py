import FWCore.ParameterSet.Config as cms

muonPATUserDataRochesterCorrectionAdder = cms.EDProducer('MuonPATUserDataRochesterCorrectionAdder',
    ## input particle collection of type edm::View<reco::GenParticle>
    src = cms.InputTag("slimmedMuons"),
    path = cms.string("ZCounting/ZUtils/data/MuonRochCorr.txt"),
    debug = cms.bool(False),
    applyEnergyCorrections = cms.bool(False)
)
