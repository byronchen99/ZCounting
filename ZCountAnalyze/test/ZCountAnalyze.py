import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('isData', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch off generator Info")
options.register('selectEvents', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch off cuts ")
options.register('max', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "specify number of events")


process = cms.Process("zcounting")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('ZCounting')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.max))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #data
        #'root://xrootd-cms.infn.it//store/data/Run2016B/SingleMuon/AOD/PromptReco-v2/000/275/376/00000/FE9B738C-A839-E611-A098-02163E01472F.root'
        #MC
        #'root://xrootd-cms.infn.it///store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/FED3B50D-8AB1-E611-B230-FA163E71DC21.root'
        #'root://xrootd-cms.infn.it///store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/FE3F6B2D-6CB1-E611-B673-001E67457E7C.root'
        #local
        'file:/eos/home-d/dwalter/data/ZCounting/FE3F6B2D-6CB1-E611-B673-001E67457E7C.root'
    )
)


outFileName = options.outputFile + '_' + str(options.job) +  '.root'
print ('Using output file ' + outFileName)
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(outFileName))

# if a good Z candidate is found
from ZCounting.ZUtils.GenZDecay_cfi import genZDecay
process.genZDecay = genZDecay.clone(
    src = "genParticles",
)
# if no good Z candidate is found. E.g. if gamma is found instead
from ZCounting.ZUtils.GenZLeptonDecay_cfi import genZLeptonDecay
process.genZLeptonDecay = genZLeptonDecay.clone(
    src = "genParticles",
)

process.load("ZCounting.ZCountAnalyze.ZCounting_cfi")
process.zcounting.isData = cms.untracked.bool(options.isData)
process.zcounting.selectEvents = cms.untracked.bool(options.selectEvents)
process.zcounting.genZCollection = cms.InputTag("genZDecay")
process.zcounting.genZLeptonCollection = cms.InputTag("genZLeptonDecay")

process.p = cms.Path(process.genZDecay + process.genZLeptonDecay + process.zcounting)
