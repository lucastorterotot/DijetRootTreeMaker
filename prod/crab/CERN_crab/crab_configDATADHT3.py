from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'GJet-Feb03-Run2016H_ver2_V1_74X_T3'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/h/hlattaud/private/production_GJet/CMSSW_8_0_26_patch1/src/CMSDIJET/DijetRootTreeMaker/prod/localtest_flat-dataH-cfg_miniAOD.py'
config.JobType.pyCfgParams = ['runOnData=True']

config.Data.inputDataset = '/SinglePhoton/Run2016H-03Feb2017_ver2-v1/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/hlattaud/GammaJet/GJet-ReReco-Run2016H'

config.Site.storageSite = 'T3_FR_IPNL'
