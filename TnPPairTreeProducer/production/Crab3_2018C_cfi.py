from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'TnPPairTrees_UL_2018C_V13'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/nfs/dust/cms/user/dwalter/CMSSW_10_6_13/src/ZCounting/TnPPairTreeProducer/test/TnPPairTreeProducer.py'
config.JobType.pyCfgParams = ['era=2018']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ['']

config.Data.inputDataset = '/SingleMuon/Run2018C-12Nov2019_UL2018-v3/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
# config.Data.runRange = '319337-320065'
config.Data.lumiMask = '/nfs/dust/cms/user/dwalter/ZCounting/CMSSW_10_6_26/src/ZCounting/TnPPairTreeProducer/production/res/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys_mod.txt'
config.Data.publication = False
config.Data.outputDatasetTag = 'TnPPairTrees_V13_UL2018C'

config.Site.storageSite = 'T2_DE_DESY'
#config.Site.blacklist = ['T3_US_Baylor','T2_US_Vanderbilt','T2_US_UCSD','T2_US_Florida' ]
