from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = 'TMP-TEST-2017-mc'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/l/ltortero/JEC-task/CMSSW_9_4_10/src/CMSDIJET/DijetRootTreeMaker/prod/Config_RunonMiniaod.py'

config.Data.inputDataset = '/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6-v1/MINIAODSIM'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 500000
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

config.Site.storageSite = 'T3_FR_IPNL'

