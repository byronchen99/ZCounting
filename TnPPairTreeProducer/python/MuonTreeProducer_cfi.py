import FWCore.ParameterSet.Config as cms

muonTreeProducer = cms.EDAnalyzer('MuonTreeProducer',
    # Trigger - sources
    TriggerEvent   = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
    TriggerResults = cms.InputTag('TriggerResults', '', 'HLT'),

    # PV source
    edmPVName=cms.untracked.string('offlinePrimaryVertices'),

    # PV selection
    VtxNTracksFitMin=cms.double(0.),
    VtxNdofMin=cms.double(4.),
    VtxAbsZMax=cms.double(24.),
    VtxRhoMax=cms.double(2.),

    # Muon and Track sources
    edmMuonName  = cms.untracked.string('muons'),
    edmTrackName = cms.untracked.string('generalTracks'),

    # Muons - triggers
    MuonTriggerNames=cms.vstring("HLT_Mu17_v*"),

    # Muon selection
    PtCutL  = cms.double(27.0),
    EtaCutL = cms.double(2.4)
)
