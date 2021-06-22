import FWCore.ParameterSet.Config as cms

pileupMCTemplateMaker = cms.EDAnalyzer('PileupMCTemplateMaker',
  src = cms.InputTag('slimmedAddPileupInfo'),
)
