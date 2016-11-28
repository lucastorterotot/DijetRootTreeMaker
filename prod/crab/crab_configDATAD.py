from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'GJet-ReReco-Run2016D_v2'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../localtest_flat-dataBCD-cfg_miniAOD.py'
config.JobType.pyCfgParams = ['runOnData=True']

config.Data.inputDataset = '/SinglePhoton/Run2016D-23Sep2016-v1/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/hlattaud/GammaJet/GJet-ReReco-Run2016D'

config.Site.storageSite = 'T3_FR_IPNL'
config.Site.whitelist = ['T3_FR_IPNL']
