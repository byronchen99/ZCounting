import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('verbosity', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "verbosity level ()")
options.register('maxEvents', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "specify number of events")
options.register('samplename', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "specify sample name (dy, tt) ")
options.register('year', 2017, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "specify year (2016, 2017, 2018) ")
options.parseArguments()


if options.samplename == '':
    raise RuntimeError('ZCountAnalyze.py: cannot run without specifying a samplename')
elif options.samplename not in ('tt', 'dy', 'w', 'qcd', 'zz', 'wz', 'ww'):
    raise RuntimeError('ZCountAnalyze.py: unknown samplename '+options.samplename)

if options.year not in (2016, 2017, 2018):
    raise RuntimeError('ZCountAnalyze.py: unknown year '+options.year)

print("processing {0} events".format(options.maxEvents))

process = cms.Process("zcounting")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

if options.verbosity == 0:
    process.MessageLogger.cerr.threshold = 'WARNING'
elif options.verbosity == 1:
    process.MessageLogger.cerr.threshold = 'INFO'
elif options.verbosity == 2:
    process.MessageLogger.cerr.threshold = 'DEBUG'
    process.MessageLogger.cerr.debugModules = cms.untracked.vstring(
        'zcounting',
        )
else:
    process.MessageLogger.cerr.threshold = 'ERROR'

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('ZCounting')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
        # '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/00000/CE864285-6D1C-E911-B09D-34E6D7BDDECE.root'
        # '/store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0075FEFF-D341-E811-AD4B-001E673D35A9.root'
        # '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU0to75TuneCUETP8M1_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/128C7A83-4F17-E711-A485-FA163E23B354.root'
        '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU28to62HcalNZSRAW_94X_mcRun2_asymptotic_v3-v1/250000/1A702941-9A6C-E911-A2BD-AC1F6BAB6860.root'
    )
)

outFileName = options.outputFile + '_' + str(options.job) + '.root'
print('Using output file ' + outFileName)
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFileName))

process.load("ZCounting.ZCountAnalyze.ZCounting_cfi")

tasks = cms.Task()

if options.year == 2016:
    process.zcounting.muon_trigger_patterns=cms.vstring(
        "HLT_IsoMu24_v*",
        "HLT_IsoTkMu24_v*")

# if no good Z candidate is found. E.g. if gamma is found instead
if options.samplename == 'dy':
    from ZCounting.ZUtils.GenZLeptonDecay_cfi import genZLeptonDecay
    process.genZLeptonDecay = genZLeptonDecay.clone(src="prunedGenParticles",)
    tasks.add(process.genZLeptonDecay)

    process.zcounting.hasGenZ = True
    process.zcounting.genZLeptonCollection = cms.InputTag("genZLeptonDecay")

elif options.samplename == 'tt':
    process.load('TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff')
    process.initSubset.src = cms.InputTag("prunedGenParticles")
    process.decaySubset.src = cms.InputTag("prunedGenParticles")
    process.decaySubset.runMode = 'Run2'
    process.decaySubset.fillMode = 'kStable' # Top before Decay, after Radiation
    tasks.add(process.makeGenEvtTask)

    process.zcounting.hasGenTt = True
    process.zcounting.genTtCollection = cms.InputTag("genEvt")


process.p = cms.Path(process.zcounting, tasks)
