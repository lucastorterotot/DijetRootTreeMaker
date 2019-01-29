from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = 'TMP-TEST-2018-data'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/l/ltortero/JEC-task/CMSSW_10_2_5/src/CMSDIJET/DijetRootTreeMaker/prod/Config_RunonMiniaod.py'
config.JobType.pyCfgParams = ['isRunonData=True']

config.Data.inputDataset = '/EGamma/Run2018B-26Sep2018-v1/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Era/Prompt/Cert_316998-319312_13TeV_PromptReco_Collisions18_JSON_eraB.txt' # TODO : update if available
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

config.Site.storageSite = 'T3_FR_IPNL'

