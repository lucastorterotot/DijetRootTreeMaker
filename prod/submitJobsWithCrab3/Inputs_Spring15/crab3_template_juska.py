from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'WORKINGAREA'
config.General.requestName = 'WORKINGDIR'
config.section_('JobType')
config.JobType.psetName = 'CMSSWCFG'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['OUTFILENAME']
config.section_('Data')
config.Data.inputDataset = 'INPUTDATASET'
config.Data.unitsPerJob = FILESPERJOB #without '' since it must be an int
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/juska/Spring15_AK4cors_v2'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ['T2_UA_KIPT','T2_US_Caltech']
