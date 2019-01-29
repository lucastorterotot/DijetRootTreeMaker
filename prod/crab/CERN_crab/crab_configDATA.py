from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = 'TMP-TEST-2017-data'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/l/ltortero/JEC-task/CMSSW_9_4_10/src/CMSDIJET/DijetRootTreeMaker/prod/Config_RunonMiniaod.py'
config.JobType.pyCfgParams = ['isRunonData=True']

config.Data.inputDataset = '/SinglePhoton/Run2017B-17Nov2017-v1/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

config.Site.storageSite = 'T3_FR_IPNL'

