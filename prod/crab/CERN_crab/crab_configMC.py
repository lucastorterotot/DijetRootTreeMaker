from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'GJet-RunIISpring16MiniAODv1_v2'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/user/h/hlattaud/private/tutorial/CMSSW_8_0_11/src/CMSDIJET/DijetRootTreeMaker/prod/localtest_flat-MC-cfg_miniAOD.py'

config.Data.inputDataset = '/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 500000
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_jetmet/hlattaud/GammaJet/GJet-RunIISpring16MiniAODv1'

config.Site.storageSite = 'T2_CH_CERN'

