from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'GJet-Feb03-Run2016B'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/h/hlattaud/private/CMSSW_8_0_26_patch1/src/CMSDIJET/DijetRootTreeMaker/prod/localtest_flat-dataBCD-cfg_miniAOD.py'
config.JobType.pyCfgParams = ['runOnData=True']

config.Data.inputDataset = '/SinglePhoton/Run2016B-03Feb2017_ver2-v2/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_jetmet/hlattaud/GammaJet/GJet-ReReco-Run2016B'

config.Site.storageSite = 'T2_CH_CERN'

