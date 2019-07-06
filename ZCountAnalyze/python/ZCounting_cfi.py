import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCounting',

    # Primary vertex - source
    edmPVName = cms.InputTag("offlinePrimaryVertices"),
    # Primary vertex - settings
    VtxNTracksFitMin = cms.untracked.double(0.),
    VtxNdofMin       = cms.untracked.double(4.),
    VtxAbsZMax       = cms.untracked.double(24.),
    VtxRhoMax        = cms.untracked.double(2.),

    # Trigger - sources
    TriggerEvent   = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
    TriggerResults = cms.InputTag('TriggerResults','','HLT'),

    # Muons - source
    edmMuonName    = cms.InputTag('muons'),
    edmTrackName   = cms.InputTag('generalTracks'),
    # Muons - triggers
    MuonTriggerNames       = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"),
    MuonTriggerObjectNames = cms.vstring("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                         "hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"),
    # Muons - settings
    MuonTagPtCut    = cms.untracked.double(27.0),
    MuonProbePtCut  = cms.untracked.double(27.0),
    MuonTagEtaCut   = cms.untracked.double(2.4),
    MuonProbeEtaCut = cms.untracked.double(2.4),
    MuonIDType      = cms.untracked.string("Tight"),
    MuonIsoType     = cms.untracked.string("NULL"),
    MuonIsoCut      = cms.untracked.double(0.),

    # General information
    MassMin = cms.untracked.double(66.0),
    MassMax = cms.untracked.double(116.0),

    isData       = cms.untracked.bool(True),
    selectEvents = cms.untracked.bool(True),    #If 'False', store all events that are processed

    # For MC only
    genEventInfo         = cms.InputTag('generator'),
    genZCollection       = cms.InputTag('genZDecay'),
    genZLeptonCollection = cms.InputTag('genZLeptonDecay'),

)
