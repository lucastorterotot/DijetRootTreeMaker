from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'GJet-RunIISpring16MiniAODv1_v2_endcapPhoton_20M'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/h/hlattaud/private/production_GJet/CMSSW_8_0_29/src/CMSDIJET/DijetRootTreeMaker/prod/Config_RunonMiniaod.py'

config.Data.inputDataset = '/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_20M/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 500000
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/hlattaud/GammaJet/'

config.Site.storageSite = 'T3_FR_IPNL'

