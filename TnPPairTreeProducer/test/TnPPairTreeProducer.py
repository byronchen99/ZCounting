import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')

options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "output File (w/o .root)")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nMax', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "specify number of events")
options.register('era', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "year + era")
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
     #'file:/nfs/dust/cms/user/dwalter/data/Lumi/AOD_files/Run2017B/SingleMuon/AOD/17Nov2017-v1/C6DAEE21-37D8-E711-8AC0-02163E0145B8.root'
     # remote test
     '/store/data/Run2017F/SingleMuon/AOD/09Aug2019_UL2017-v1/270004/90E36450-548B-BF41-BA70-438F46B99A7D.root'
     # Single Mu 2017 Low PU -
     #'/store/data/Run2017H/SingleMuon/AOD/17Nov2017-v2/90000/FEE0C793-6A34-E811-BC0A-1866DA7F95AE.root'
     # Single Mu 2017 B UL-Reco
     #'/store/data/Run2017B/SingleMuon/AOD/09Aug2019_UL2017-v1/50010/6053BF5B-683C-8C4B-8EF3-B630DC182CAA.root'
    )
)


outFileName = options.outputFile + '_' + str(options.job) + '.root'
print('Using output file ' + outFileName)
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFileName))

print("producing for era "+options.era)
process.load("ZCounting.TnPPairTreeProducer.TnPPairTreeProducer_cfi")
if options.era == '2017':
    print("set 2017 configuration")
    process.tnpPairTreeProducer.MuonTriggerNames = cms.vstring("HLT_IsoMu27_v*")
    process.tnpPairTreeProducer.PtCutL1 = cms.untracked.double(30.0)
    process.tnpPairTreeProducer.PtCutL2 = cms.untracked.double(30.0)

elif options.era == '2017H':
    print("set 2017 Low PU configuration")
    process.tnpPairTreeProducer.MuonTriggerNames = cms.vstring("HLT_HIMu17_v*")
    process.tnpPairTreeProducer.PtCutL1 = cms.untracked.double(27.0)
    process.tnpPairTreeProducer.PtCutL2 = cms.untracked.double(27.0)

process.p = cms.Path(process.tnpPairTreeProducer)
