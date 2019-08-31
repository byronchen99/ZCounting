import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nEvents', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "specify number of events")


process = cms.Process("zcounting")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('ZCounting')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
        '/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/0289A3FC-E8BE-604A-9BFA-33FD036DFBE0.root'
    )
)


outFileName = options.outputFile + '_' + str(options.job) + '.root'
print('Using output file ' + outFileName)
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFileName))

# if no good Z candidate is found. E.g. if gamma is found instead
from ZCounting.ZUtils.GenZLeptonDecay_cfi import genZLeptonDecay
process.genZLeptonDecay = genZLeptonDecay.clone(
    src="prunedGenParticles",
)

process.load("ZCounting.ZCountAnalyze.ZCounting_cfi")
process.zcounting.genZLeptonCollection = cms.InputTag("genZLeptonDecay")

process.p = cms.Path(process.genZLeptonDecay + process.zcounting)
