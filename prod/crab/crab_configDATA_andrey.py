from CRABClient.UserUtilities import config

config = config()

config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../MultijetJEC_cfg.py'
config.JobType.pyCfgParams = ['runOnData=True']

config.Data.inputDataset = ''
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/aapopov/JetMET/160910/'

config.Site.storageSite = 'T3_FR_IPNL'
