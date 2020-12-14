import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCounting',
    pileupSummaryInfoCollection=cms.InputTag("slimmedAddPileupInfo"),

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
    muon_trigger_patterns=cms.vstring("HLT_L1SingleMu18_v*",
                                      "HLT_L1SingleMu25_v*",
                                      "HLT_IsoMu24_v*",
                                      "HLT_IsoMu27_v*",
                                      "HLT_IsoMu30_v*"),
    muon_trigger_DRMAX=cms.untracked.double(1.0),

    # For MC only
    hasGenZ = cms.untracked.bool(False),
    hasGenTt = cms.untracked.bool(False),

    genEventInfo         = cms.InputTag('generator'),
    genZLeptonCollection = cms.InputTag(''),
    genTtCollection = cms.InputTag(''),

)
