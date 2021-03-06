from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'NEWGJet-RunIISummer16MiniAODv2_V1val'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/user/h/hlattaud/private/tutorial/CMSSW_8_0_24_patch1/src/CMSDIJET/DijetRootTreeMaker/prod/localtest_flat-MC-cfg_miniAOD.py'

config.Data.inputDataset = '/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_20M/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 500000
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_jetmet/hlattaud/GammaJet/GJet-RunIISummer16MiniAODv2'

config.Site.storageSite = 'T2_CH_CERN'
