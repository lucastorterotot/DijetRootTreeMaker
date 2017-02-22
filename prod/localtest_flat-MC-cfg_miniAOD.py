import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
qgDatabaseVersion = 'v1' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion


###################################### Run on AOD instead of MiniAOD? ########
runOnAOD=False
###################################### Run on RECO instead of MiniAOD? ########
runOnRECO=False
if runOnRECO: runOnAOD=True



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
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'



#--------------------- Report and output ---------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string('mylocaltest_10.root'),
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



# ----------------------- Jet Tool Box  -----------------
# ----- giulia test: do not recluster ak4 and ca8 jets to save time --------


process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.slimmedGenJetsAK8 = ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)


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



process.source = cms.Source("PoolSource",
    # 2016B data "file:/afs/cern.ch/user/j/juska/eos/cms/store/data/Run2016B/JetHT/MINIAOD/PromptReco-v1/000/272/771/00000/B4A77EBA-DB15-E611-A15E-02163E013590.root")
	fileNames = cms.untracked.vstring("/store/mc/RunIISpring16MiniAODv1/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/02ED80EA-5012-E611-BC71-842B2B76670F.root") # mAODv2 
	# mAODv1 "file:/afs/cern.ch/user/j/juska/eos/cms/store/mc/RunIISpring16MiniAODv1/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/4EC1A37D-840D-E611-957D-0025905C543A.root")
    )




#-----------------------------------------------------------
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
process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.calibratedPatPhotons.isMC=cms.bool(True)
process.calibratedPatElectrons.isMC=cms.bool(True)
#process.calibratedPatPhotons 
#process.calibratedPatElectrons 
##-------------Add quarkGluon tagging---------------------------



process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag("slimmedJets")       # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

#------------------------------------------------------------
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

calo_collection=''
cluster_collection=''
pfcalo_collection=''
   
   
   
#-------------------------To RUN ON ReminiAOD ------------------------   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',

  # There's no avoiding this in Consumes era
  isData          = cms.bool(False),

  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  metforggen       = cms.InputTag('slimmedMETs'),  
  metpuppi         = cms.InputTag('slimmedMETsPuppi'),
  eb               = cms.InputTag('reducedEBRecHits'),
  ee               = cms.InputTag('reducedEERecHits'),
  metTypeI         = cms.InputTag('slimmedMETs'),
  PFCands          = cms.InputTag('packedPFCandidates'),

  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  
   ## PHOTONS ########################################
  ptMinPhoton               = cms.double(10),
  Photon                    = cms.InputTag('slimmedPhotons'),
  Photonsmeared             = cms.InputTag('calibratedPatPhotons'),
  
 # GenPhoton                 = cms.InputTag('slimmedGenPhotons'),
  full5x5SigmaIEtaIEtaMap   = cms.InputTag('photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta'),
  phoChargedIsolation       = cms.InputTag('photonIDValueMapProducer:phoChargedIsolation'),
  phoNeutralHadronIsolation = cms.InputTag('photonIDValueMapProducer:phoNeutralHadronIsolation'),
  phoPhotonIsolation        = cms.InputTag('photonIDValueMapProducer:phoPhotonIsolation'),
  
  ## MC ########################################
  pu               = cms.untracked.InputTag('slimmedAddPileupInfo'), # Updated from untracked to 80X by Juska
  ptHat            = cms.untracked.InputTag('generator'), # Why do these need to be 'untracked' anyway?
  genParticles     = cms.InputTag('prunedGenParticlesDijet'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
  genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),


  ## trigger ###################################
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

  ## electrons ######################################## 

  Electrons                 = cms.InputTag('slimmedElectrons'),
  Electronssmeared          = cms.InputTag('calibratedPatElectrons'),
  ## muons ########################################

  Muons                     = cms.InputTag('slimmedMuons'),

  ## JECs ################
  redoJECs  = cms.bool(True),

  ## Version Summer15_25nsV3
   L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFPuppi.txt'),
  L3corrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFPuppi.txt'),
  ResCorrAK4PUPPI_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFPuppi.txt'),
  L1corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt'),
  L3corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt'),
  ResCorrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt'),
  L1corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFPuppi.txt'),
  L2corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFPuppi.txt'),
  L3corrAK4PUPPI_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFPuppi.txt'),
  L1corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt'),
  L3corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt'),
  L1RCcorr_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1RC_AK4PFchs.txt')


)


# ------------------ path --------------------------

process.p = cms.Path()

if runOnRECO:
   process.p += process.pfClusterRefsForJets_step
                                                        
process.p +=                     process.prunedGenParticlesDijet
process.p +=                     process.chs 
process.p +=                     process.slimmedGenJetsAK8                      
process.p +=                     process.dijets 
