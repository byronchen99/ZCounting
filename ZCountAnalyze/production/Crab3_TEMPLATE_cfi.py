# This is a template file for submitting crab3 jobs
# edit fields with <...>
# to setup crab, do
#     $ cmsenv
#     $ source /cvmfs/cms.cern.ch/crab3/crab.sh
#     $ voms-proxy-init --voms cms
# submit jobs with
#     $ crab submit -c <this file>

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = '<name>'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pyCfgParams = ['isData=<False/True>','selectEvents=<False/True>']
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '<total path to project cfg file, e.g. /eos/home-d/<user>/CMSSW_9_4_13/src/ZCounting/ZCountAnalyze/test/ZCountAnalyze.py'
#config.JobType.inputFiles = ['<specify, if resource files are needed>']

config.Data.inputDataset = '<source file from https://cmsweb.cern.ch/das has to be /AODSIM (MC) /AOD (data)>'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/<project name>/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = '<output name>'

config.Site.storageSite = 'T2_DE_DESY'
#config.Site.blacklist = ['T3_US_Baylor','T2_US_Vanderbilt','T2_US_UCSD','T2_US_Florida' ]
