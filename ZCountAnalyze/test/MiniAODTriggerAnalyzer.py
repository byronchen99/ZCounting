import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/FEF559B1-95DF-E811-9528-D4AE52900EF9.root'
        #'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/0289A3FC-E8BE-604A-9BFA-33FD036DFBE0.root'
        # data - SingleMuon - 2018
        #'/store/data/Run2018B/SingleMuon/MINIAOD/17Sep2018-v1/60000/FF47BB90-FC1A-CC44-A635-2B8B8C64AA39.root'
        # data - SingleMuon - 2017H LowPU
        #'/store/data/Run2017H/SingleMuon/MINIAOD/17Nov2017-v2/90000/FA9FA831-8B34-E811-BA1D-008CFAC93CFC.root'
        # DYJets - Fall17 Flat PU
        '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP1_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU0to75TuneCP1_12Apr2018_94X_mc2017_realistic_v14-v2/90000/F8ACEB3C-072C-E911-B67D-002590A82B8E.root'
    )
)

process.demo = cms.EDAnalyzer("MiniAODTriggerAnalyzer",
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("slimmedPatTrigger"),
)

process.p = cms.Path(process.demo)
