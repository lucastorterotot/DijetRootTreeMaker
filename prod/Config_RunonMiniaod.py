import random
import re
import string


import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500



#Parser to the command line option.
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

options.setDefault('maxEvents', 1001)


options.register(
    'globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
    'Global tag to be used'
)

options.register(
    'isRunonData', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    'Indicates whether the job processes data or simulation'
)

options.register(
    'isLegacy', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    'In case of data, distinguishes Reminiaod and legacy rereco. Ignored for simulation'
)

options.register(
    'isEndcaps', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    'specify wether or not you run on endcaps photons'
)

options.register(
    'UseReclusteredJets', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    'specify wether or not you run with reclustered Jets from reco'
)

options.register(
    'JECData', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
    'Which JEC do we use for data'
)

options.register(
    'JECMC', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
    'Which JEC do we use for MC '
)

options.parseArguments()

runOnData            = options.isRunonData
runOnLegacy          = options.isLegacy
runOnEndcaps         = options.isEndcaps
runOnReclusteredJets = options.UseReclusteredJets

print(runOnData)
print(runOnLegacy)
if not runOnEndcaps:
    print(runOnEndcaps)

if runOnData:
   print('You are running on datas')
else:
   print('You are running on MC')

# Provide a default global tag if user has not given any.  Chosen as
# according to recommendations for JEC [1].
# [1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC?rev=125
if len(options.globalTag) == 0:
    if runOnData:
        if runOnLegacy:
            options.globalTag = '80X_dataRun2_2016LegacyRepro_v4'
        else:
            options.globalTag = '80X_dataRun2_2016SeptRepro_v7'
    else:
        options.globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
    
    print 'WARNING: No global tag provided. Will use the default one: {}.'.format(
        options.globalTag
)




if len(options.JECData) == 0:
    if runOnLegacy:
            options.JECData = 'Summer16_07Aug2017BCD_V1_DATA'
    else:
            options.JECData = 'Summer16_03Feb2017BCD_V3'
    
    
    print 'WARNING: No data JEC provided. Will use the default one: {}.'.format(
        options.JECData
)

if len(options.JECMC) == 0:
    if runOnLegacy:
            options.JECMC = 'Summer16_07Aug2017_V1_MC'
    else:
            options.JECMC = 'Summer16_03Feb2017_V1_MC'
    
    
    print 'WARNING: No MC JEC provided. Will use the default one: {}.'.format(
        options.JECMC
)




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

## This local test conf is updated to match the grid production version
## 13 June 2016 by Juska.



## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag)
#--------------------- Report and output ---------------------------
# Note: in grid runs this parameter is not used.
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))




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


if runOnData:
        if runOnLegacy: 
           process.source = cms.Source("PoolSource",
           fileNames = cms.untracked.vstring("/store/data/Run2016C/SinglePhoton/MINIAOD/07Aug17-v1/50000/0825D2D7-C89E-E711-A061-008CFAC94258.root")
           )
        else:
             process.source = cms.Source("PoolSource",
             fileNames = cms.untracked.vstring("/store/data/Run2016H/SinglePhoton/MINIAOD/03Feb2017_ver2-v1/100000/0027C019-EFEA-E611-8E79-7845C4FC35E1.root")
             )
else:
       process.source = cms.Source("PoolSource",
           fileNames = cms.untracked.vstring("/store/mc/RunIISummer16MiniAODv2/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_20M/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/08AC6DB1-19B7-E611-A8F8-001E67E71E20.root")
           )
#process.source.eventsToProcess = cms.untracked.VEventRange("281613:11018807")#,"283884:939706570","283884:870499187","283885:16020018","274316:389398083") 


# Add map producer for photon identification 
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


if not runOnData:
     ##Photon Energy smearer------------------------------------------
      process.load('Configuration.StandardSequences.Services_cff')
      process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",

                                                       calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                                                 engineName = cms.untracked.string('TRandom3'),
                                                                                           ),
                                                       
                                                       calibratedPatPhotons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                                                 engineName = cms.untracked.string('TRandom3'),
                                                                                           ),
                                                       ) 



if not runOnLegacy:
      process.load('EgammaAnalysis.ElectronTools.calibratedPatbeforeGXPhotonsRun2_cfi')
      if runOnData: 
            process.calibratedPatPhotons80X.isMC = cms.bool(False)
      else:
            process.calibratedPatPhotons80X.isMC = cms.bool(True)       

if runOnData: 
      process.calibratedPatPhotons.isMC = cms.bool(False)
else: 
      process.calibratedPatPhotons.isMC = cms.bool(True)
      


if runOnLegacy:
        
        process.calibratedPatPhotons.correctionFile = cms.string(files["Legacy2016_17AUg17"])
        process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
              src = cms.InputTag('calibratedPatPhotons'),
              cut = cms.string('pt>5 && abs(eta)'))
else:
        process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
              src = cms.InputTag('calibratedPatPhotons80X'),
              cut = cms.string('pt>5 && abs(eta)'))
              
srcViD = "selectedPhotons"
process.egmPhotonIDs.physicsObjectSrc = cms.InputTag(srcViD)
process.egmPhotonIsolation.srcToIsolate = cms.InputTag(srcViD)
process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)
process.photonRegressionValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)
process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)









#---------------Met ReminiAOD recipe ---------------------

if not runOnLegacy:
         
         # EG cleaned

	## Following lines are for default MET for Type1 corrections.
          from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
          runMetCorAndUncFromMiniAOD(process,
                                    isData=True ,
                                   )
          from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
          corMETFromMuonAndEG(process,
                             pfCandCollection="", #not needed                                                                                                                                                                                                                                                                                                                      
                              electronCollection="slimmedElectronsBeforeGSFix",
                              photonCollection="slimmedPhotonsBeforeGSFix",
                              corElectronCollection="slimmedElectrons",
                              corPhotonCollection="slimmedPhotons", 
                              allMETEGCorrected=True,
                              muCorrection=False,
                              eGCorrection=True,
                              runOnMiniAOD=True,
                              postfix="MuEGClean"
                              )
          process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
          process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
          process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
          process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
          del process.slimmedMETsMuEGClean.caloMET
           
     # If you are running in the scheduled mode:
          process.egcorrMET = cms.Sequence(
                  process.cleanedPhotonsMuEGClean+process.cleanedCorPhotonsMuEGClean+
                  process.matchedPhotonsMuEGClean + process.matchedElectronsMuEGClean +
                  process.corMETPhotonMuEGClean+process.corMETElectronMuEGClean+
                  process.patPFMetT1MuEGClean+process.patPFMetRawMuEGClean+
                  process.patPFMetT1SmearMuEGClean+process.patPFMetT1TxyMuEGClean+
                  process.patPFMetTxyMuEGClean+process.patPFMetT1JetEnUpMuEGClean+
                  process.patPFMetT1JetResUpMuEGClean+process.patPFMetT1SmearJetResUpMuEGClean+
                  process.patPFMetT1ElectronEnUpMuEGClean+process.patPFMetT1PhotonEnUpMuEGClean+
                  process.patPFMetT1MuonEnUpMuEGClean+process.patPFMetT1TauEnUpMuEGClean+
                  process.patPFMetT1UnclusteredEnUpMuEGClean+process.patPFMetT1JetEnDownMuEGClean+
                  process.patPFMetT1JetResDownMuEGClean+process.patPFMetT1SmearJetResDownMuEGClean+
                  process.patPFMetT1ElectronEnDownMuEGClean+process.patPFMetT1PhotonEnDownMuEGClean+
                  process.patPFMetT1MuonEnDownMuEGClean+process.patPFMetT1TauEnDownMuEGClean+
                  process.patPFMetT1UnclusteredEnDownMuEGClean+process.slimmedMETsMuEGClean)


if not runOnLegacy and not runOnData:
         
         from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

   # If you only want to re-correct for JEC and get the proper uncertainties for the default MET
         runMetCorAndUncFromMiniAOD(process,
                          isData= False,
                         )

   # Now you are creating the bad muon corrected MET
         process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
         process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
         process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)

         from PhysicsTools.PatUtils.tools.muonRecoMitigation import muonRecoMitigation

         muonRecoMitigation(
                                process = process,
                       pfCandCollection = "packedPFCandidates", #input PF Candidate Collection
                                runOnMiniAOD = True, #To determine if you are running on AOD or MiniAOD
                                selection="", #You can use a custom selection for your bad muons. Leave empty if you would like to use the bad muon recipe definition.
                                muonCollection="", #The muon collection name where your custom selection will be applied to. Leave empty if you would like to use the bad muon recipe definition.
                                cleanCollName="cleanMuonsPFCandidates", #output pf candidate collection ame
                                cleaningScheme="computeAllApplyClone", #Options are: "all", "computeAllApplyBad","computeAllApplyClone". Decides which (or both) bad muon collections to be used for MET cleaning coming from the bad muon recipe.
                                postfix="" #Use if you would like to add a post fix to your muon / pf collections
                                )

         runMetCorAndUncFromMiniAOD(process,
                                    isData=False,
                                    pfCandColl="cleanMuonsPFCandidates",
                                    recoMetFromPFCs=True,
                                    postfix="MuClean"
                                    )


         process.mucorMET = cms.Sequence(                     
                 process.badGlobalMuonTaggerMAOD *
                 process.cloneGlobalMuonTaggerMAOD *
               #process.badMuons * # If you are using cleaning mode "all", uncomment this line
                 process.cleanMuonsPFCandidates *
                 process.fullPatMetSequenceMuClean
                 )




##-------------Add quarkGluon tagging---------------------------



process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag("slimmedJets")       # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion



# Recluster the ak4 jets or not

#if runOnReclusteredJets:
   





# Residue from AOD and RECO running
calo_collection=''
cluster_collection=''
pfcalo_collection=''
   

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',

  # There's no avoiding this in Consumes era
  isData          = cms.bool(True) if runOnData else cms.bool(False),
  isreMiniAOD     = cms.bool(False) if runOnLegacy else cms.bool(True),
  Endcaps_photon  = cms.bool(True) if runOnEndcaps else cms.bool(False),
  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  metforggen       = cms.InputTag('slimmedMETs'),
  metEGcleaned     = cms.InputTag('slimmedMETs') if runOnLegacy or not runOnData else cms.InputTag('slimmedMETsMuEGClean',processName = "jetToolbox"),   
  metpuppi              = cms.InputTag('slimmedMETsPuppi'),
  PFCands = cms.InputTag('packedPFCandidates'),
  
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  
  ## PHOTONS ########################################
  ptMinPhoton               = cms.double(10),
  Photon                    = cms.InputTag('selectedPhotons'),
  Photonsmeared             = cms.InputTag('calibratedPatPhotons') if runOnLegacy else cms.InputTag('calibratedPatPhotons80X'),
  GenPhoton                 = cms.InputTag('slimmedGenPhotons'),
  full5x5SigmaIEtaIEtaMap   = cms.InputTag('photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta'),
  phoChargedIsolation       = cms.InputTag('photonIDValueMapProducer:phoChargedIsolation'),
  phoNeutralHadronIsolation = cms.InputTag('photonIDValueMapProducer:phoNeutralHadronIsolation'),
  phoPhotonIsolation        = cms.InputTag('photonIDValueMapProducer:phoPhotonIsolation'),
  PhotonUncorr              = cms.InputTag('calibratedPatbeforegxPhotons') if runOnLegacy else cms.InputTag('selectedPhotons'),
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

  
  L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L2Relative_AK4PFPuppi.txt'),
  L3corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L3Absolute_AK4PFPuppi.txt'),
  ResCorrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L2L3Residual_AK4PFPuppi.txt'),
  L1corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L2Relative_AK8PFchs.txt'),
  L3corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L3Absolute_AK8PFchs.txt'),
  ResCorrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L3Absolute_AK4PFchs.txt'),
  L1corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L2Relative_AK4PFPuppi.txt'),
  L3corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L3Absolute_AK4PFPuppi.txt'),
  L1corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L2Relative_AK8PFchs.txt'),
  L3corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECMC+'/'+options.JECMC+'_L3Absolute_AK8PFchs.txt'),
  L1RCcorr_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/'+options.JECData+'/'+options.JECData+'_L1RC_AK4PFchs.txt')
)


# ------------------ path --------------------------
process.p = cms.Path()  
process.p +=                      process.chs
if not runOnData:
       process.p +=                     process.prunedGenParticlesDijet
       process.p +=                     process.slimmedGenJetsAK8
if not runOnData and not runOnLegacy:
       process.p +=                     process.mucorMET
       process.p +=                     process.fullPatMetSequence

process.p +=                      process.regressionApplication
if runOnLegacy: 
   process.p +=                      process.calibratedPatPhotons
else:
   process.p +=                      process.calibratedPatPhotons*process.calibratedPatPhotons80X*process.calibratedPatPhotonsbeforeGS
process.p +=                      process.selectedPhotons
process.p +=                      process.egmPhotonIDSequence
if not runOnLegacy and runOnData:
       process.p +=                      process.fullPatMetSequence  # If you are re-correctign the default MET
       process.p +=                      process.egcorrMET
process.p +=                      process.dijets
