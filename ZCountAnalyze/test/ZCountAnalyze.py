import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('outputFile','output',
    VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,
    "output File (w/o .root)"
    )
options.register('job', 0,
    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
    "job number"
    )
options.register('verbosity', 1,
    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
    "verbosity level ()"
    )
options.register('maxEvents', -1,
    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
    "specify number of events"
    )
options.register('samplename', '',
    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,
    "specify sample name (dy, tt) "
    )
options.register('year', '2017',
    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,
    "specify era (2016preVFP, 2016postVFP, 2017, 2018) "
    )
options.parseArguments()


if options.samplename == '':
    raise RuntimeError('ZCountAnalyze.py: cannot run without specifying a samplename')
elif options.samplename not in ('tt', 'dy', 'w', 'qcd', 'zz', 'wz', 'ww'):
    raise RuntimeError('ZCountAnalyze.py: unknown samplename '+options.samplename)

if options.year not in ("2016preVFP", "2016postVFP", "2017", "2017H", "2018"):
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
        # '/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU0to75TuneCP5_102X_upgrade2018_realistic_v15-v1/00000/01C7EF7D-6860-7F41-A815-5CCC6856AB28.root'
        # '/store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0075FEFF-D341-E811-AD4B-001E673D35A9.root'
        # '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU0to75TuneCUETP8M1_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/128C7A83-4F17-E711-A485-FA163E23B354.root'
        # '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU28to62HcalNZSRAW_94X_mcRun2_asymptotic_v3-v1/250000/1A702941-9A6C-E911-A2BD-AC1F6BAB6860.root'

        # UL 2016
        # '/store/mc/RunIISummer20UL16MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU0to75_106X_mcRun2_asymptotic_v13-v2/230000/0ACCBC6C-6308-2D42-8D55-18F3522FE776.root'
        # 'root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL16MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU0to75_106X_mcRun2_asymptotic_v13-v2/230000/10E0D710-F696-7248-9A7C-80B0F127C6B5.root'
        # UL 2016 APV
        # '/store/mc/RunIISummer20UL16MiniAODAPV/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/FlatPU0to75_106X_mcRun2_asymptotic_preVFP_v8-v2/00000/0040C830-7251-BA46-BE12-86EE8CFD0907.root'

        # UL 2017
        '/store/mc/RunIISummer20UL17MiniAODv2/WW_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/280000/8E5C41E3-0E1B-B54E-A18C-09C2244E6B5F.root'
        # UL 2018
        # '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/120000/02FD7B88-3EDC-C64D-8E18-F9A8F9E7E7DF.root'
        # '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/230000/0271F29A-6EA5-054E-B733-EB5EC7A80F2C.root'
        # '/store/mc/RunIISummer20UL18MiniAODv2/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E17A5D77-2DC3-F24A-A13D-C65191D2BDCC.root'
        # '/store/mc/RunIISummer20UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/260000/0071F930-6376-7A48-89F1-74E189BD3BFC.root'
    )
)

if options.year == "2016preVFP":
    globalTag = '106X_mcRun2_asymptotic_preVPF_v8'
    l1PrefireECAL = "UL2016preVFP"
    l1PrefireMuon = "2016preVFP"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2016aUL.txt'
elif options.year == "2016postVFP":
    globalTag = '106X_mcRun2_asymptotic_v13'
    l1PrefireECAL = "UL2016postVFP"
    l1PrefireMuon = "2016BG"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2016bUL.txt'
elif options.year == "2017":
    globalTag = '106X_mc2017_realistic_v6'
    l1PrefireECAL = "UL2017BtoF"
    l1PrefireMuon = "20172018"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2017UL.txt'
elif options.year == "2018":
    globalTag = '106X_upgrade2018_realistic_v11_L1v1'
    l1PrefireECAL = "None"
    l1PrefireMuon = "20172018"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2018UL.txt'

# ## Geometry and Detector Conditions
# process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
# process.load("Configuration.StandardSequences.MagneticField_cff")
# process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
# process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
# process.load("Configuration.Geometry.GeometryRecoDB_cff")

print("Using global tag: "+globalTag)
process.GlobalTag.globaltag = globalTag

outFileName = options.outputFile + '_' + str(options.job) + '.root'
print('Using output file ' + outFileName)
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFileName))

tasks = cms.Task()

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

# main analyzer
process.load("ZCounting.ZCountAnalyze.ZCounting_cfi")
# events before selection
process.load("ZCounting.ZCountAnalyze.CountEventAnalyzer_cfi")
# pileup MC template maker
process.load("ZCounting.ZCountAnalyze.PileupMCTemplateMaker_cfi")

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer

process.prefiringweight = l1PrefiringWeightProducer.clone(
TheJets = cms.InputTag("updatedPatJetsUpdatedJEC"),
# L1Maps = cms.string("../../../ZCounting/ZCountAnalyze/data/All2017Gand2017HPrefiringMaps.root"),
DataEraECAL = cms.string(l1PrefireECAL),
DataEraMuon = cms.string(l1PrefireMuon),
UseJetEMPt = cms.bool(False),
PrefiringRateSystematicUnctyECAL = cms.double(0.2),
PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)

if options.year == "2016postVFP":
    process.prefiringweight2016H = l1PrefiringWeightProducer.clone(
    DataEraECAL = cms.string("None"),
    DataEraMuon = cms.string("2016H"),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
    )
    process.prefireSequence = cms.Sequence(process.prefiringweight + process.prefiringweight2016H)
else:
    process.prefireSequence = cms.Sequence(process.prefiringweight)

### apply rochester corrections
from ZCounting.ZUtils.muonPATUserDataRochesterCorrectionAdder_cfi import muonPATUserDataRochesterCorrectionAdder
process.selectedMuonsWithEnCorrInfo = muonPATUserDataRochesterCorrectionAdder.clone(
    src = 'slimmedMuons',
    path = roccoFile,
    applyEnergyCorrections = False,
    debug = False,
)
tasks.add(process.selectedMuonsWithEnCorrInfo)

process.zcounting.era = options.year
process.zcounting.pat_muons = cms.InputTag('selectedMuonsWithEnCorrInfo')

# if no good Z candidate is found. E.g. if gamma is found instead
if options.samplename in ('dy', 'zz', 'wz'):
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


process.p = cms.Path(
    process.countEvents *
    process.pileupMCTemplateMaker *
    process.jecSequence *
    process.prefireSequence *
    process.zcounting,
    tasks)
