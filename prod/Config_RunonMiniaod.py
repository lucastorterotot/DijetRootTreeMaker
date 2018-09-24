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
runOnEndcaps         = options.isEndcaps
runOnReclusteredJets = options.UseReclusteredJets

print(runOnData)

if runOnData:
   print('You are running on datas')
else:
   print('You are running on MC')

# Provide a default global tag if user has not given any.  Chosen as
# according to recommendations for JEC [1].
# [1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC?rev=125
if len(options.globalTag) == 0:
    if runOnData:
            options.globalTag = '94X_dataRun2_ReReco_EOY17_v2'
    
    print 'WARNING: No global tag provided. Will use the default one: {}.'.format(
        options.globalTag
)




if len(options.JECData) == 0:
    
    options.JECData = 'Summer16_07Aug2017BCD_V1_DATA'
   
    
    
    print 'WARNING: No data JEC provided. Will use the default one: {}.'.format(
        options.JECData
)

if len(options.JECMC) == 0:
    options.JECMC = 'Summer16_07Aug2017_V1_MC'
    
    
    
    print 'WARNING: No MC JEC provided. Will use the default one: {}.'.format(
        options.JECData
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

#adding the Quark gluon discriminant to ak4 jets
for type in ['AK4PFchs','AK4PFchs_antib']:
  QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
    record = cms.string('QGLikelihoodRcd'),
    tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
    label  = cms.untracked.string('QGL_'+type)
  )))






## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = options.globalTag

#--------------------- Report and output ---------------------------
# Note: in grid runs this parameter is not used.
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))



process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string('mylocaltest_Run2017A_10.root'),
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


# apply the charge hadron substraction on Pf candidate
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
     process.source = cms.Source("PoolSource",
         
         fileNames = cms.untracked.vstring("/store/data/Run2017F/SinglePhoton/MINIAOD/17Nov2017-v1/50000/004425E2-EBDE-E711-8F80-008CFAF5543A.root"),
    
         inputCommands =  cms.untracked.vstring("keep *","drop  *_isolatedTracks__*")
    
     )
else: 
     process.source = cms.Source("PoolSource",
         
         fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAOD/GJet_Pt-15To6000_TuneCP5-Flat_13TeV_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/40000/2C41D5C8-5DDA-E711-876C-1866DA87C2CD.root"),
    
         inputCommands =  cms.untracked.vstring("keep *","drop  *_isolatedTracks__*")
    
     ) 



# Photon Energy scale/smearing for 94X (Feb 2018) have to be updated 

useMiniAOD=True # set for false for AOD

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if useMiniAOD:
    switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
    switchOnVIDPhotonIdProducer(process,DataFormat.MiniAOD)
else:
    switchOnVIDElectronIdProducer(process,DataFormat.AOD)
    switchOnVIDPhotonIdProducer(process,DataFormat.AOD)


# define which IDs we want to produce and add them to VID 
ele_id_modules =  [ 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
                  ]
pho_id_modules =  [ 'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_TrueVtx_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff', 
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff', 
                    'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff'
                  ]
for idmod in ele_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in pho_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


if useMiniAOD:
    process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
    process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
    process.calibratedPatElectrons.electrons = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess())
    process.calibratedPatPhotons.photons = cms.InputTag('slimmedPhotons',processName=cms.InputTag.skipCurrentProcess())
    
    process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
    src = cms.InputTag('calibratedPatPhotons'), 
    cut = cms.string('pt>5 && abs(eta)')
)
    srcViD = "selectedPhotons"
    process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons')
    process.egmPhotonIDs.physicsObjectSrc = cms.InputTag(srcViD)
    process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('calibratedPatElectrons') 
    process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)
    process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag(srcViD)
    process.egmPhotonIsolation.srcToIsolate = cms.InputTag(srcViD)
    if hasattr(process,'heepIDVarValueMaps'):
        process.heepIDVarValueMaps.elesMiniAOD = cms.InputTag('calibratedPatElectrons')

    from RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff import egamma_modifications
    from RecoEgamma.EgammaTools.egammaObjectModifications_tools import makeVIDBitsModifier,makeVIDinPATIDsModifier,makeEnergyScaleAndSmearingSysModifier                                     
    egamma_modifications.append(makeVIDBitsModifier(process,"egmGsfElectronIDs","egmPhotonIDs"))
    egamma_modifications.append(makeVIDinPATIDsModifier(process,"egmGsfElectronIDs","egmPhotonIDs"))
    egamma_modifications.append(makeEnergyScaleAndSmearingSysModifier("calibratedPatElectrons",srcViD))

    #add the HEEP trk isol to the slimmed electron
    egamma_modifications[0].electron_config.heepV70TrkPtIso = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso")
    for pset in egamma_modifications:
        pset.overrideExistingValues = cms.bool(True)
        if hasattr(pset,"electron_config"): pset.electron_config.electronSrc = cms.InputTag("calibratedPatElectrons")
        if hasattr(pset,"photon_config"): pset.photon_config.photonSrc = cms.InputTag(srcViD)

    process.slimmedElectrons = cms.EDProducer("ModifiedElectronProducer",
                                              src=cms.InputTag("calibratedPatElectrons"),
                                              modifierConfig = cms.PSet(
            modifications = egamma_modifications
            )
                                              )
    process.slimmedPhotons = cms.EDProducer("ModifiedPhotonProducer",
                                              src=cms.InputTag(srcViD),
                                              modifierConfig = cms.PSet(
            modifications = cms.VPSet(egamma_modifications)
            )
                                            )

else:
    process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
    process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
    #we will rename the collections to the standard names so VID and other things pick it up

    process.gedGsfElectrons = process.calibratedElectrons.clone(electrons=cms.InputTag("gedGsfElectrons",processName=cms.InputTag.skipCurrentProcess()))
    process.gedPhotons = process.calibratedPhotons.clone(photons=cms.InputTag("gedPhotons",processName=cms.InputTag.skipCurrentProcess())) 

if useMiniAOD:
    process.egammaScaleSmearTask = cms.Task(process.calibratedPatElectrons,process.slimmedElectrons,
                                            process.selectedPhotons,process.slimmedPhotons
                                            )
else:
    process.egammaScaleSmearTask = cms.Task(process.gedGsfElectrons,
                                            process.gedPhotons
                                            )
                                          

process.egammaScaleSmearSeq = cms.Sequence( process.egammaScaleSmearTask)
process.egammaScaleSmearAndVIDSeq = cms.Sequence(process.egammaScaleSmearSeq*
    process.egmGsfElectronIDSequence*
    process.egmPhotonIDSequence)



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
  isData = cms.bool(True) if runOnData else cms.bool(False),
  isreMiniAOD     = cms.bool(False),
  Endcaps_photon = cms.bool(True) if runOnEndcaps else cms.bool(False),
  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  metforggen       = cms.InputTag('slimmedMETs'),
  metEGcleaned     = cms.InputTag('slimmedMETs'),   
  metpuppi              = cms.InputTag('slimmedMETsPuppi'),
  PFCands = cms.InputTag('packedPFCandidates'),
  
 # QGT              = cms.InputTag('QGTagger'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  
  ## PHOTONS ########################################
  ptMinPhoton               = cms.double(10),
  Photon                    = cms.InputTag('selectedPhotons'),
  Photonsmeared             = cms.InputTag('selectedPhotons'),
  GenPhoton                 = cms.InputTag('slimmedGenPhotons'),
  full5x5SigmaIEtaIEtaMap   = cms.InputTag('photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta'),
  phoChargedIsolation       = cms.InputTag('photonIDValueMapProducer:phoChargedIsolation'),
  phoNeutralHadronIsolation = cms.InputTag('photonIDValueMapProducer:phoNeutralHadronIsolation'),
  phoPhotonIsolation        = cms.InputTag('photonIDValueMapProducer:phoPhotonIsolation'),
  PhotonUncorr              = cms.InputTag('selectedPhotons'),
  eb               = cms.InputTag('reducedEgamma:reducedEBRecHits'),
  ee               = cms.InputTag('reducedEgamma:reducedEERecHits'),
  
  ## MC ########################################
  pu                        = cms.untracked.InputTag('slimmedAddPileupInfo'), 
  ptHat                     = cms.untracked.InputTag('generator'), 
  genParticles              = cms.InputTag('prunedGenParticlesDijet'),
  genJetsAK4                = cms.InputTag('slimmedGenJets'), 
  genJetsAK8                = cms.InputTag('slimmedGenJetsAK8'),  
 
 
   ## electrons ######################################## 

  Electrons                 = cms.InputTag('slimmedElectrons'),
  Electronssmeared          = cms.InputTag('slimmedElectrons'),
  ## muons ########################################

  Muons                     = cms.InputTag('slimmedMuons'),
  
  ## trigger ###################################
 
  ##### For single Photon ##### 
  triggerAlias     = cms.vstring('HLTPhoton30','HLTPhoton50','HLTPhoton75','HLTPhoton90','HLTPhoton120','HLTPhoton165','HLTPhoton200'),                         
  triggerSelection = cms.vstring(

     
     ###for SinglePhotons
     'HLT_Photon33_v*',
     'HLT_Photon50_R9Id90_HE10_IsoM_v*',
     'HLT_Photon75_R9Id90_HE10_IsoM_v*',
     'HLT_Photon90_R9Id90_HE10_IsoM_v*',
     'HLT_Photon120_R9Id90_HE10_IsoM_v*',
     'HLT_Photon165_R9Id90_HE10_IsoM_v*',
     'HLT_Photon200_v*'
     
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
  triggerObjects = cms.InputTag('slimmedPatTrigger'),
  filters = cms.vstring(
        'hltEG33L1EG26HEFilter','hltEG50R9Id90HE10IsoMTrackIsoFilter', 'hltEG75R9Id90HE10IsoMTrackIsoFilter','hltEG90R9Id90HE10IsoMTrackIsoFilter','hltEG120R9Id90HE10IsoMTrackIsoFilter','hltEG165R9Id90HE10IsoMTrackIsoFilter','hltEG200HEFilter'),#

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
process.p +=                      process.calibratedPatPhotons
process.p +=                      process.selectedPhotons
process.p +=                      process.egammaScaleSmearAndVIDSeq
process.p +=                      process.QGTagger
process.p +=                      process.dijets
