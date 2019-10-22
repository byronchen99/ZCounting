from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TnPPairTrees_2017D_V02'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/eos/home-d/dwalter/CMSSW_10_6_4/src/ZCounting/TnPPairTreeProducer/test/TnPPairTreeProducer.py'
config.JobType.pyCfgParams = ['era=2017']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ['']

config.Data.inputDataset = '/SingleMuon/Run2017D-09Aug2019_UL2017-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.runRange = '302030-303434'
config.Data.outLFNDirBase = '/store/user/%s/TnPPairTrees_2017D_V02/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'TnPPairTrees_2017D_V02'

config.Site.storageSite = 'T2_DE_DESY'
#config.Site.blacklist = ['T3_US_Baylor','T2_US_Vanderbilt','T2_US_UCSD','T2_US_Florida' ]
