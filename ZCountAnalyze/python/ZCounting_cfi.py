import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCounting',
    # Primary vertex - source
    edmPVName = cms.InputTag("offlineSlimmedPrimaryVertices"),
    # Primary vertex - settings
    VtxNTracksFitMin = cms.untracked.double(0.),
    VtxNdofMin       = cms.untracked.double(4.),
    VtxAbsZMax       = cms.untracked.double(24.),
    VtxRhoMax        = cms.untracked.double(2.),

    # Trigger - sources
    trigger_bits = cms.InputTag("TriggerResults","","HLT"),
    trigger_objects = cms.InputTag("slimmedPatTrigger"),

    # Muons - source
    pat_muons = cms.InputTag("slimmedMuons","","PAT"),

    # Muon - triggers
    muon_trigger_patterns       = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"),

    # For MC only
    genEventInfo         = cms.InputTag('generator'),
    genZLeptonCollection = cms.InputTag('genZLeptonDecay'),

)
