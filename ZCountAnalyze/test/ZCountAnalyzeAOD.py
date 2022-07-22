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
options.register('samplename', 'dy',
    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,
    "specify sample name ('dy', 'tt', 'w', 'qcd', 'zz', 'wz', 'ww', 't', 'ttz', 'ttw', 'met', 'smu', 'se') "
    )
options.register('year', '2017',
    VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,
    "specify era (2016preVFP, 2016postVFP, 2017, 2017H, 2018, 2022) "
    )
options.parseArguments()


if options.samplename == '':
    raise RuntimeError('ZCountAnalyze.py: cannot run without specifying a samplename')
elif options.samplename not in ('dy', 'tt', 'w', 'qcd', 'zz', 'wz', 'ww', 't', 'ttz', 'ttw', 'met', 'smu', 'se'):
    raise RuntimeError('ZCountAnalyze.py: unknown samplename '+options.samplename)

if options.year not in ("2016preVFP", "2016postVFP", "2017", "2017H", "2018", "2022"):
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
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
        # Data
        # Single Mu 2017 D UL-Reco
        # 'file:/pnfs/desy.de/cms/tier2/store/data/Run2017E/SingleMuon/AOD/09Aug2019_UL2017-v1/260000/0005DF00-5EE0-C84E-9241-08FA62D9EFF7.root'
        # "file:/pnfs/desy.de/cms/tier2/store/data/Run2017H/SingleMuon/AOD/09Aug2019_UL2017_LowPU-v1/10000/02FD55DE-EA02-FD4C-86D3-0E93C4E9E280.root"
        #'file:/eos/cms/store/data/Run2018D/SingleMuon/AOD/12Nov2019_UL2018-v6/230000/1217FEF3-168A-7C46-96EA-F5D8B6C6BC4A.root'
        # SIMULATION: AODSIM
        #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer20UL17RECO/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/AODSIM/106X_mc2017_realistic_v6-v1/00000/0C0FD46D-6DAE-6F45-B9F4-129A6ABC69C4.root'
        # "/store/mc/RunIISummer20UL17RECO/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/106X_mc2017_realistic_v6-v1/130000/2678F234-2D42-AB4F-89E4-80AFDE9EFB82.root"
        # "/store/mc/RunIISummer20UL17RECO/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/AODSIM/106X_mc2017_realistic_v6-v1/00000/0C0FD46D-6DAE-6F45-B9F4-129A6ABC69C4.root"
        # '/store/mc/Run3Winter22DR/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/AODSIM/FlatPU0to70_122X_mcRun3_2021_realistic_v9_ext1-v2/2520000/002de84f-b906-4509-8cb6-872ee0cf4029.root'
        'file:/eos/home-d/dwalter/Lumi/ZCounting_Standalone/CMSSW_12_4_3/src/ZCounting/ZCountAnalyze/002de84f-b906-4509-8cb6-872ee0cf4029.root'
    )
)

genWeights = cms.vint32(1005, 1009, 1004, 1007, 1002, 1003),
pdfWeights = cms.vint32(1252, 1354)

if options.year == "2016preVFP":
    globalTag = '106X_mcRun2_asymptotic_preVFP_v11'
    l1PrefireECAL = "UL2016preVFP"
    l1PrefireMuon = "2016preVFP"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2016aUL.txt'
elif options.year == "2016postVFP":
    globalTag = '106X_mcRun2_asymptotic_v17'
    l1PrefireECAL = "UL2016postVFP"
    l1PrefireMuon = "2016BG"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2016bUL.txt'
elif options.year in ("2017", "2017H"):
    globalTag = '106X_mc2017_realistic_v6'
    l1PrefireECAL = "UL2017BtoF"
    l1PrefireMuon = "20172018"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2017UL.txt'
elif options.year == "2018":
    globalTag = '106X_upgrade2018_realistic_v11_L1v1'
    l1PrefireECAL = "None"
    l1PrefireMuon = "20172018"
    roccoFile = 'ZCounting/ZUtils/data/RoccoR2018UL.txt'
elif options.year == "2022":
    globalTag = '112X_mcRun3_2021_realistic_v16'
    genWeights = cms.vint32(1021, 1041, 1006, 1011, 1016, 1031)
    pdfWeights = cms.vint32(1252, 1354)
    l1PrefireECAL = "None"
    l1PrefireMuon = "None"
    roccoFile = "None"

# TODO: prefiring is currently not working in AOD... need to setup

if options.samplename in ('smu','met'):
# set global tag (INFO: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions?redirectedfrom=CMS.SWGuideFrontierConditions)
    if options.year in ('2016', '2017', '2018'):
        globalTag = '106X_dataRun2_v32'
    elif options.year in ('2022', ):
        globalTag = '123X_dataRun3_Prompt_v12'
    else:
        globalTag = '124X_dataRun3_Prompt_v4'

# ## Geometry and Detector Conditions
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

print("Using global tag: "+globalTag)
process.GlobalTag.globaltag = globalTag

outFileName = options.outputFile + '_' + str(options.job) + '.root'
print('Using output file ' + outFileName)
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFileName))

tasks = cms.Task()

# main analyzer
process.load("ZCounting.ZCountAnalyze.ZCountingAOD_cfi")
# events before selection
process.load("ZCounting.ZCountAnalyze.CountEventAnalyzer_cfi")
# pileup MC template maker
process.load("ZCounting.ZCountAnalyze.PileupMCTemplateMaker_cfi")
process.pileupMCTemplateMaker.src = cms.InputTag('addPileupInfo')

process.zcounting.era = options.year
process.zcounting.roccorFile = cms.string(roccoFile)

# if no good Z candidate is found. E.g. if gamma is found instead
if options.samplename in ('dy', 'zz', 'wz', 'ttz'):
    from ZCounting.ZUtils.GenZLeptonDecay_cfi import genZLeptonDecay
    process.genZLeptonDecay = genZLeptonDecay.clone(src="genParticles",)
    tasks.add(process.genZLeptonDecay)
    process.zcounting.hasGenZ = True
    process.zcounting.genZLeptonCollection = cms.InputTag("genZLeptonDecay")
    
    process.zcounting.genWeights = genWeights
    process.zcounting.pdfWeights = pdfWeights
    process.zcounting.printLHE = True
    
    ### particle level definition of gen events
    from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
    process.genParticles2HepMC = genParticles2HepMC.clone()
    tasks.add(process.genParticles2HepMC)

    from GeneratorInterface.RivetInterface.particleLevel_cfi import particleLevel
    process.particleLevel = particleLevel.clone()
    tasks.add(process.particleLevel)
    
elif options.samplename in ('smu','se','met'):
    # only store necessary information for data
    process.zcounting.isData = True
    if options.samplename in ('smu', 'se'):
        process.zcounting.met_trigger_patterns = cms.vstring()
        if options.samplename in ('smu',):
            process.zcounting.electron_trigger_patterns = cms.vstring()
            process.zcounting.store_electrons = False
        elif options.samplename in ('se',):
            process.zcounting.muon_trigger_patterns = cms.vstring()
            process.zcounting.store_muons = False

if options.year == '2017H':
    print("set 2017 Low PU configuration")
    # trigger emulation of HLT_IsoMu24
    # process.load("HLTrigger.Configuration.HLT_User_cff")
    process.zcounting.muon_trigger_patterns = cms.vstring("HLT_IsoMu24_v*") #"HLT_HIMu17_v*"
    process.zcounting.emulateTrigger = cms.untracked.bool(True)

process.p = cms.Path(
    process.countEvents *
    process.pileupMCTemplateMaker *
    process.zcounting,
    tasks)
