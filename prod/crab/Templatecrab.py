from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'name of the crab job'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'absolute path to your config file'
config.JobType.pyCfgParams = ['runOnData=True']

config.Data.inputDataset = '/SinglePhoton/Run2016B-03Feb2017_ver2-v2/MINIAOD'
config.Data.lumiMask = 'JSON file of your run era '
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = 'Your output directory at your storage site '

config.Site.storageSite = 'your storage site'
