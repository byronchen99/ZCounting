import FWCore.ParameterSet.Config as cms

tnpPairTreeProducer = cms.EDAnalyzer('TnPPairTreeProducer',
    # Trigger - sources
    TriggerEvent   = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
    TriggerResults = cms.InputTag('TriggerResults', '', 'HLT'),

    # PV source
    edmPVName=cms.untracked.string('offlinePrimaryVertices'),

    # PV selection
    VtxNTracksFitMin=cms.untracked.double(0.),
    VtxNdofMin=cms.untracked.double(4.),
    VtxAbsZMax=cms.untracked.double(24.),
    VtxRhoMax=cms.untracked.double(2.),

    # Muon and Track sources
    edmMuonName  = cms.untracked.string('muons'),
    edmTrackName = cms.untracked.string('generalTracks'),

    # Muons - triggers
    MuonTriggerNames=cms.vstring("HLT_IsoMu27_v*"),
    MuonTriggerObjectNames=cms.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", ),

    # Muon selection
    IDType   = cms.untracked.string("Custom"),  # Custom, Tight, Medium, Loose
    DxyCut=cms.untracked.double(-1.),
    DzCut=cms.untracked.double(-1.),
    IsoType  = cms.untracked.string("None"),  # Tracker-based, PF-based
    IsoCut   = cms.untracked.double(0.),

    PtCutL1  = cms.untracked.double(30.0),
    EtaCutL1 = cms.untracked.double(2.4),
    PtCutL2=cms.untracked.double(30.0),
    EtaCutL2=cms.untracked.double(2.4),

    MassMin  = cms.untracked.double(50.0),
    MassMax  = cms.untracked.double(150.0),

)
