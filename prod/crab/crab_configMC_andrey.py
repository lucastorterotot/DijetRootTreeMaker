from CRABClient.UserUtilities import config

config = config()

config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../MultijetJEC_cfg.py'
config.JobType.pyCfgParams = ['runOnData=False', 'saveGenJets=True']

config.Data.inputDataset = ''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 500000
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/aapopov/JetMET/160910/'

config.Site.storageSite = 'T3_FR_IPNL'
