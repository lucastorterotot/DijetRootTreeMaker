from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'GJet-Feb03-Run2016Flate_scale_74X'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/h/hlattaud/private/production_GJet/CMSSW_8_0_26_patch1/src/CMSDIJET/DijetRootTreeMaker/prod/localtest_flat-dataG-cfg_miniAOD.py'
config.JobType.pyCfgParams = ['runOnData=True']

config.Data.inputDataset = '/SinglePhoton/Run2016F-03Feb2017-v1/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/work/h/hlattaud/private/production_GJet/CMSSW_8_0_26_patch1/src/CMSDIJET/DijetRootTreeMaker/prod/crab/CERN_crab/lateFJSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/hlattaud/GammaJet/GJet-ReReco-Run2016F'

config.Site.storageSite = 'T3_FR_IPNL'
