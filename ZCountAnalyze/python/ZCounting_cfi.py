import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCounting',
    edmPVName = cms.InputTag("offlinePrimaryVertices"),

    #Primary vertex
    VtxNTracksFitMin = cms.untracked.double(0.),
    VtxNdofMin       = cms.untracked.double(4.),
    VtxAbsZMax       = cms.untracked.double(24.),
    VtxRhoMax        = cms.untracked.double(2.)
)
