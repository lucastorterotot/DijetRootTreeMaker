import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')


process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
qgDatabaseVersion = 'v1' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

from CondCore.DBCommon.CondDBSetup_cfi import *
QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDBSetup,
      toGet = cms.VPSet(),
      connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
)

for type in ['AK4PFchs','AK4PFchs_antib']:
  QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
    record = cms.string('QGLikelihoodRcd'),
    tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
    label  = cms.untracked.string('QGL_'+type)
  )))





## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v4'
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' #80X_mcRun2_asymptotic_2016_miniAODv2

#--------------------- Report and output ---------------------------
# Note: in grid runs this parameter is not used.
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(50))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string('mylocaltest_Run2016B_10.root'),
                                 #fileName=cms.string(THISROOTFILE),
                                 closeFileFast = cms.untracked.bool(True)
                                 )

## --- suppress long output ---> wantSummary = cms.untracked.bool(False) 

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)

############## output  edm format ###############
process.out = cms.OutputModule('PoolOutputModule',                                                                                                                  
                               fileName = cms.untracked.string('jettoolbox.root'),                                                                              
                               outputCommands = cms.untracked.vstring([                                                                 
                                                                      'keep *_slimmedJets_*_*',                                                                  
                                                                      'keep *_slimmedJetsAK8_*_*',                                                                  
                                                                       ])                                                                                           
                               )


# Added 'vertexRef().isNonnull() &&' check for 80X data compatibility. Juska
process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string(' fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.slimmedGenJetsAK8 = ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)





#-------------------------------------------------------
# Gen Particles Pruner
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.prunedGenParticlesDijet = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
    "drop  *  ", # by default
    "keep ( status = 3 || (status>=21 && status<=29) )", # keep hard process particles
    )
)



#------------- Recluster Gen Jets to access the constituents -------
#already in toolbox, just add keep statements

process.out.outputCommands.append("keep *_slimmedGenJets_*_*")
process.out.outputCommands.append("keep *_slimmedGenJetsAK8_*_*")

##-------------------- Define the source  ----------------------------



# Note: for grid running it does not matter what's here, as input data is
# handled separately there.

process.source = cms.Source("PoolSource",
    #  fileNames = cms.untracked.vstring("/store/data/Run2016H/SinglePhoton/MINIAOD/PromptReco-v3/000/284/036/00000/726BEBFC-619F-E611-862A-02163E0121A2.root")
    fileNames = cms.untracked.vstring("/store/data/Run2016H/SinglePhoton/MINIAOD/03Feb2017_ver2-v1/100000/0027C019-EFEA-E611-8E79-7845C4FC35E1.root","file:/afs/cern.ch/work/h/hlattaud/private/production_GJet/CMSSW_8_0_26_patch1/src/CMSDIJET/DijetRootTreeMaker/pickevents.root")
   # fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/h/hlattaud/private/production_GJet/CMSSW_8_0_26_patch1/src/CMSDIJET/DijetRootTreeMaker/pickevents_oldreco.root")
    #fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/h/hlattaud/private/CMSSW_9_1_0/src/pickevents.root")
    
)
#process.source.eventsToProcess = cms.untracked.VEventRange("281613:11018807")#,"283884:939706570","283884:870499187","283885:16020018","274316:389398083")

#---keep un reg photon slimmed and before gx fix -------------

#process.slimmedPhotons_noreg = slimmedPhotonsclone()
#process.slimmedPhotonsBeforeGSFix_noreg = slimmedPhotonsBeforeGSFix.clone()
#------------------------------------------------------------
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
for idmod in my_id_modules:
         setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection) 


process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)

from EgammaAnalysis.ElectronTools.calibrationTablesRun2 import files

process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
#process.load('EgammaAnalysis.ElectronTools.calibratedPatbeforeGXPhotonsRun2_cfi')


#process.regressionApplication

#process.calibratedPatPhotons
process.calibratedPatPhotons.isMC = cms.bool(False)# this is 74X
#process.calibratedPatPhotons.correctionFile = cms.string(files["Moriond2017_JEC"])

#process.calibratedPatPhotonsbeforeGS.correctionFile = cms.string(files["Moriond2017_JEC"])
#process.calibratedPatPhotonsbeforeGS.isMC = cms.bool(False)

process.calibratedPatPhotons80X.isMC = cms.bool(False)


#process.calibratedPatbeforegxPhotons 
#process.calibratedPatbeforegxPhotons.isMC = cms.bool(False)


#process.EGMRegression = cms.Path(process.regressionApplication)

#process.EGMSmearerPhotons   = cms.Path(process.calibratedPatPhotons)



process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
    src = cms.InputTag('calibratedPatPhotons80X'), # cms.InputTag('slimmedphoton74X'),# this is 74X regression 
    cut = cms.string('pt>5 && abs(eta)')
)
srcViD = "selectedPhotons"#"slimmedPhotons"
process.egmPhotonIDs.physicsObjectSrc = cms.InputTag(srcViD)
process.egmPhotonIsolation.srcToIsolate = cms.InputTag(srcViD)
process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)
process.photonRegressionValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)
process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)


##-------------------- User analyzer  --------------------------------


##-------------------- MET --------------------------------



#---------------NewMet ReminiAOD recipe ---------------------




# EG cleaned

	## Following lines are for default MET for Type1 corrections.
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

   # If you only want to re-correct for JEC and get the proper uncertainties for the default MET
#runMetCorAndUncFromMiniAOD(process,
#                          isData=True ,
#                         )

   # Now you are creating the e/g corrected MET on top of the bad muon corrected MET (on re-miniaod)
#from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
#corMETFromMuonAndEG(process,
#                    pfCandCollection="", #not needed                                                                                                                                                                                                                                                                                                                      
#                    electronCollection="slimmedElectronsBeforeGSFix",
#                    photonCollection="slimmedPhotonsBeforeGSFix",
#                    corElectronCollection="slimmedElectrons",
#                    corPhotonCollection="slimmedPhotons", #slimmedPhotons
#                    allMETEGCorrected=True,
#                    muCorrection=False,
#                    eGCorrection=True,
#                    runOnMiniAOD=True,
#                    postfix="MuEGClean"
#                    )
#process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
#process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
#process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
#process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
#del process.slimmedMETsMuEGClean.caloMET
 
     # If you are running in the scheduled mode:
#process.egcorrMET = cms.Sequence(
#        process.cleanedPhotonsMuEGClean+process.cleanedCorPhotonsMuEGClean+
#        process.matchedPhotonsMuEGClean + process.matchedElectronsMuEGClean +
#        process.corMETPhotonMuEGClean+process.corMETElectronMuEGClean+
 #       process.patPFMetT1MuEGClean+process.patPFMetRawMuEGClean+
 #       process.patPFMetT1SmearMuEGClean+process.patPFMetT1TxyMuEGClean+
#        process.patPFMetTxyMuEGClean+process.patPFMetT1JetEnUpMuEGClean+
#        process.patPFMetT1JetResUpMuEGClean+process.patPFMetT1SmearJetResUpMuEGClean+
#        process.patPFMetT1ElectronEnUpMuEGClean+process.patPFMetT1PhotonEnUpMuEGClean+
#        process.patPFMetT1MuonEnUpMuEGClean+process.patPFMetT1TauEnUpMuEGClean+
#        process.patPFMetT1UnclusteredEnUpMuEGClean+process.patPFMetT1JetEnDownMuEGClean+
#        process.patPFMetT1JetResDownMuEGClean+process.patPFMetT1SmearJetResDownMuEGClean+
#        process.patPFMetT1ElectronEnDownMuEGClean+process.patPFMetT1PhotonEnDownMuEGClean+
#        process.patPFMetT1MuonEnDownMuEGClean+process.patPFMetT1TauEnDownMuEGClean+
#        process.patPFMetT1UnclusteredEnDownMuEGClean+process.slimmedMETsMuEGClean)









##-------------Add quarkGluon tagging---------------------------



process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag("slimmedJets")       # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion









# Residue from AOD and RECO running
calo_collection=''
cluster_collection=''
pfcalo_collection=''
   

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',

  # There's no avoiding this in Consumes era
  isData          = cms.bool(True),
  isreMiniAOD     = cms.bool(False),
  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETsUncorrected'),
  metforggen       = cms.InputTag('slimmedMETsUncorrected'),
  metEGcleaned     = cms.InputTag('slimmedMETsMuEGClean',processName = "jetToolbox"),   
  metpuppi              = cms.InputTag('slimmedMETsPuppi'),
  PFCands = cms.InputTag('packedPFCandidates'),
  
 # QGT              = cms.InputTag('QGTagger'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  
  ## PHOTONS ########################################
  ptMinPhoton               = cms.double(10),
  Photon                    = cms.InputTag('selectedPhotons'),
  Photonsmeared             = cms.InputTag('calibratedPatPhotons80X'),
 # Photonsmeared_nofix       = cms.InputTag('slimmedPhotonsBeforeGSFix',processName=cms.InputTag.skipCurrentProcess()),
  GenPhoton                 = cms.InputTag('slimmedGenPhotons'),
  full5x5SigmaIEtaIEtaMap   = cms.InputTag('photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta'),
  phoChargedIsolation       = cms.InputTag('photonIDValueMapProducer:phoChargedIsolation'),
  phoNeutralHadronIsolation = cms.InputTag('photonIDValueMapProducer:phoNeutralHadronIsolation'),
  phoPhotonIsolation        = cms.InputTag('photonIDValueMapProducer:phoPhotonIsolation'),
  #PhotonUncorr              = cms.InputTag('slimmedPhotonsBeforeGSFix'),
 # PhotonUncorr              = cms.InputTag('calibratedPatbeforegxPhotons'),
  eb               = cms.InputTag('reducedEgamma:reducedEBRecHits'),
  ee               = cms.InputTag('reducedEgamma:reducedEERecHits'),
  
  ## MC ########################################
  pu                        = cms.untracked.InputTag('slimmedAddPileupInfo'), # Updated from untracked to 80X by Juska
  ptHat                     = cms.untracked.InputTag('generator'), # Why do these need to be 'untracked' anyway?
  genParticles              = cms.InputTag('prunedGenParticlesDijet'),
  genJetsAK4                = cms.InputTag('slimmedGenJets'), 
  genJetsAK8                = cms.InputTag('slimmedGenJetsAK8'),  
 
 
   ## electrons ######################################## 

  Electrons                 = cms.InputTag('slimmedElectrons'),
  Electronssmeared          = cms.InputTag('slimmedElectrons'),
  ## muons ########################################

  Muons                     = cms.InputTag('slimmedMuons'),
  
  ## trigger ###################################
  #triggerAlias     = cms.vstring('Fat','PFHT650','PFNoPUHT650','HT750','HT550'),
  ##### For 0T data  #####
  #triggerAlias     = cms.vstring('L1Jet68','L1Jet36','L1Jet16','L1EG20','L1EG5'),
  ##### For JetHT PD ##### 
  triggerAlias     = cms.vstring('HLTPhoton30','HLTPhoton50','HLTPhoton75','HLTPhoton90','HLTPhoton120','HLTPhoton165'),                                
  triggerSelection = cms.vstring(

     
     ###for SinglePhotons
     'HLT_Photon30_R9Id90_HE10_IsoM_v*',
     'HLT_Photon50_R9Id90_HE10_IsoM_v*',
     'HLT_Photon75_R9Id90_HE10_IsoM_v*',
     'HLT_Photon90_R9Id90_HE10_IsoM_v*',
     'HLT_Photon120_R9Id90_HE10_IsoM_v*',
     'HLT_Photon165_R9Id90_HE10_IsoM_v*',
     
  ),
  
  prescalesTag          = cms.InputTag("patTrigger"),
  triggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
  triggerConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMask         = cms.bool(False),
    l1techIgnorePrescales = cms.bool(False),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    throw                 = cms.bool(False)
  ),

  triggerObjects = cms.InputTag('selectedPatTrigger'),
  filters = cms.vstring(
        'hltEG30R9Id90HE10IsoMHcalIsoFilter','hltEG50R9Id90HE10IsoMHcalIsoFilter', 'hltEG75R9Id90HE10IsoMHcalIsoFilter','hltEG90R9Id90HE10IsoMHcalIsoFilter','hltEG120R9Id90HE10IsoMTrackIsoFilter','hltEG165R9Id90HE10IsoMTrackIsoFilter'),

  ## JECs ################
  redoJECs  = cms.bool(True),

  
  L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFPuppi.txt'),
  L3corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFPuppi.txt'),
  ResCorrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFPuppi.txt'),
  L1corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFchs.txt'),
  L3corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFchs.txt'),
  ResCorrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt'),
  L1corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFPuppi.txt'),
  L3corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFPuppi.txt'),
  L1corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt'),
  L3corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt'),
  L1RCcorr_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1RC_AK4PFchs.txt')
)


# ------------------ path --------------------------

process.p = cms.Path()  
process.p +=                      process.chs
process.p +=                      process.regressionApplication 
process.p +=                      process.calibratedPatPhotons*process.calibratedPatPhotons80X#*process.calibratedPatPhotonsbeforeGS
process.p +=                      process.selectedPhotons
process.p +=                      process.egmPhotonIDSequence
#process.p +=                      process.fullPatMetSequence  # If you are re-correctign the default MET
#process.p +=                      process.egcorrMET
process.p +=                      process.dijets
