import FWCore.ParameterSet.Config as cms

beforeSelectionInfo = cms.EDAnalyzer('BeforeSelectionInfo',
  src = cms.InputTag('slimmedAddPileupInfo'),
)
