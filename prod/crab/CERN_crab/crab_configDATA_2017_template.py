from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'GJet-promptreco-Run2017A'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/h/hlattaud/private/production_GJet/CMSSW_9_2_0_patch4/src/CMSDIJET/DijetRootTreeMaker/prod/localtest_flat-dataA2017-cfg_miniAOD.py'
config.JobType.pyCfgParams = ['runOnData=True']

config.Data.inputDataset = '/SinglePhoton/Run2016B-03Feb2017_ver2-v2/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-302343_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/hlattaud/GammaJet/Gjet_promptreco_2017_runA'

config.Site.storageSite = 'T3_FR_IPNL'

