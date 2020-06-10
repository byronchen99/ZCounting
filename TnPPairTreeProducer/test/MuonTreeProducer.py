import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')

options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "output File (w/o .root)")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nMax', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "specify number of events")
options.parseArguments()

process = cms.Process("zcounting")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
#process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('ZCounting')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nMax))

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
     # local
     'file:/nfs/dust/cms/user/dwalter/data/Lumi/AOD/FE9FE5BF-03BF-834D-8720-1E917906FDD0.root'
     # Single Mu 2017 Low PU -
     #'/store/data/Run2017H/SingleMuon/AOD/17Nov2017-v2/90000/FEE0C793-6A34-E811-BC0A-1866DA7F95AE.root'
     # Single Mu 2017 B UL-Reco
     #'/store/data/Run2017B/SingleMuon/AOD/09Aug2019_UL2017-v1/50010/6053BF5B-683C-8C4B-8EF3-B630DC182CAA.root'
     #'/store/data/Run2017B/SingleMuon/AOD/09Aug2019_UL2017-v1/50010/6053BF5B-683C-8C4B-8EF3-B630DC182CAA.root'
    )
)


outFileName = options.outputFile + '_' + str(options.job) + '.root'
print('Using output file ' + outFileName)
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFileName))

process.load("ZCounting.TnPPairTreeProducer.MuonTreeProducer_cfi")


process.p = cms.Path(process.muonTreeProducer)
