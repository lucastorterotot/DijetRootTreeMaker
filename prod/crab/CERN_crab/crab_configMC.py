from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = '2018-MC'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/l/ltortero/JEC-task/CMSSW_10_2_5/src/CMSDIJET/DijetRootTreeMaker/prod/Config_RunonMiniaod.py'
config.JobType.priority = 6

config.Data.inputDataset = '/GJet_Pt-15To6000_TuneCP5-Flat_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.splitting = 'Automatic' #  'EventAwareLumiBased'
#config.Data.unitsPerJob = 500000
config.Data.publication = False
config.Data.outLFNDirBase =  '/store/user/%s/' % (getUsernameFromSiteDB())

config.Site.storageSite = 'T3_FR_IPNL'

