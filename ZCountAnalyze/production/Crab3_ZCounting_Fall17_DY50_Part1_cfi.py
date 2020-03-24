from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = '200318_ZCountAnalyze_V10_DY50_Fall17_Part1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/nfs/dust/cms/user/dwalter/CMSSW_10_2_20_UL/src/ZCounting/ZCountAnalyze/test/ZCountAnalyze.py'
config.JobType.pyCfgParams = ['samplename=dy']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ['']

config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
config.Data.publication = False
config.Data.outputDatasetTag = '200318_ZCountAnalyze_V10'

config.Site.storageSite = 'T2_DE_DESY'
#config.Site.blacklist = ['T3_US_Baylor','T2_US_Vanderbilt','T2_US_UCSD','T2_US_Florida' ]
