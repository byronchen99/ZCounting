from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'TnPPairTrees_UL_2016H_V13_v4'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/nfs/dust/cms/user/dwalter/CMSSW_10_6_13/src/ZCounting/TnPPairTreeProducer/test/TnPPairTreeProducer.py'
config.JobType.pyCfgParams = ['era=2016']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ['']

config.Data.inputDataset = '/SingleMuon/Run2016H-21Feb2020_UL2016-v1/AOD'
# config.Data.inputDataset = '/SingleMuon/Run2016H-21Feb2020_UL2016_WMass-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
# config.Data.runRange = '280919-284044'
config.Data.lumiMask = '/nfs/dust/cms/user/dwalter/ZCounting/CMSSW_10_6_26/src/ZCounting/TnPPairTreeProducer/production/res/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt'
config.Data.publication = False
config.Data.outputDatasetTag = 'TnPPairTrees_V13_UL2016H'

config.Site.storageSite = 'T2_DE_DESY'
config.Site.blacklist = ['T1_IT_CNAF',]
