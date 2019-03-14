from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = '2018-DATA'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/l/ltortero/JEC-task/CMSSW_10_2_5/src/CMSDIJET/DijetRootTreeMaker/prod/Config_RunonMiniaod.py'
config.JobType.pyCfgParams = ['isRunonData=True']
config.JobType.priority = 5

config.Data.inputDataset = '/EGamma/Run2018B-26Sep2018-v1/MINIAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

config.Site.storageSite = 'T3_FR_IPNL'
#config.Site.whitelist = ['T2_CH_CERN', 'T3_FR_IPNL']

