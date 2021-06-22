import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCounting',
    pileupSummaryInfoCollection=cms.InputTag("slimmedAddPileupInfo"),

    era = cms.string("None"),
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
    # 2^{x}
    muon_trigger_patterns=cms.vstring(
                                # x =
        "HLT_L1SingleMu18_v*",  # 0
        "HLT_L1SingleMu25_v*",  # 1
        "HLT_IsoMu24_v*",       # 2
        "HLT_IsoTkMu24_v*",     # 3
        "HLT_IsoMu27_v*",       # 4
        "HLT_IsoTkMu27_v*",     # 5
        "HLT_IsoMu30_v*",       # 6
        "HLT_IsoTkMu30_v*",     # 7
        "HLT_Mu17_v*",          # 8
        "HLT_TkMu17_v*",        # 9
        "HLT_Mu19_v*",          # 10
        "HLT_TkMu19_v*",        # 11
        "HLT_Mu20_v*",          # 12
        "HLT_TkMu20_v*",        # 13
        "HLT_Mu27_v*",          # 14
        "HLT_TkMu27_v*",        # 15
        "HLT_IsoMu22_eta2p1",   # 16
        "HLT_IsoTkMu22_eta2p1", # 17
        "HLT_IsoMu24_eta2p1",   # 18
        "HLT_IsoTkMu24_eta2p1", # 19
        "HLT_Mu50_v*",          # 20
        "HLT_TkMu50_v*"         # 21
    ),
    muon_trigger_DRMAX=cms.untracked.double(0.1),

    # For MC only
    hasGenZ = cms.untracked.bool(False),
    hasGenTt = cms.untracked.bool(False),

    genEventInfo         = cms.InputTag('generator'),
    genZLeptonCollection = cms.InputTag(''),
    genTtCollection = cms.InputTag(''),

)
