from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'TnPPairTrees_UL_2017B_V13_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/nfs/dust/cms/user/dwalter/CMSSW_10_6_13/src/ZCounting/TnPPairTreeProducer/test/TnPPairTreeProducer.py'
config.JobType.pyCfgParams = ['era=2017']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ['']

config.Data.inputDataset = '/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
# config.Data.runRange = '297020-299329'
config.Data.lumiMask = '/nfs/dust/cms/user/dwalter/ZCounting/CMSSW_10_6_26/src/ZCounting/TnPPairTreeProducer/production/res/Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt'
config.Data.publication = False
config.Data.outputDatasetTag = 'TnPPairTrees_V13_UL2017B'

config.Site.storageSite = 'T2_DE_DESY'
#config.Site.blacklist = ['T3_US_Baylor','T2_US_Vanderbilt','T2_US_UCSD','T2_US_Florida' ]
config.Site.whitelist = ['T2_IT_Pisa', ]
