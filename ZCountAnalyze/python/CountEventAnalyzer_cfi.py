import FWCore.ParameterSet.Config as cms

countEvents = cms.EDAnalyzer("CountEventAnalyzer",
    genEventInfo = cms.InputTag("generator")
)
