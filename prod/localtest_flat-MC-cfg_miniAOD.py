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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))

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
	fileNames = cms.untracked.vstring("/store/mc/RunIISpring16MiniAODv2/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/1CE70D6E-334C-E611-8BB4-002590E50AF2.root") # mAODv2 
	# mAODv1 "file:/afs/cern.ch/user/j/juska/eos/cms/store/mc/RunIISpring16MiniAODv1/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/4EC1A37D-840D-E611-957D-0025905C543A.root")
    )



##-------------------- User analyzer  --------------------------------
##-------------------- MET --------------------------------

## MET CHS (not available as slimmedMET collection)
## copied from https://github.com/cms-jet/JMEValidator/blob/CMSSW_7_6_X/python/FrameworkConfiguration.py
def clean_met_(met):
     del met.t01Variation
     del met.t1Uncertainties
     del met.t1SmearedVarsAndUncs
     del met.tXYUncForRaw
     del met.tXYUncForT1
     del met.tXYUncForT01
     del met.tXYUncForT1Smear
     del met.tXYUncForT01Smear

from PhysicsTools.PatAlgos.tools.metTools import addMETCollection

## Raw PF METs
process.load('RecoMET.METProducers.PFMET_cfi')

process.pfMet.src = cms.InputTag('packedPFCandidates')
addMETCollection(process, labelName='patPFMet', metSource='pfMet') # RAW MET
process.patPFMet.addGenMET = False

process.pfMetCHS = process.pfMet.clone()
process.pfMetCHS.src = cms.InputTag("chs")
process.pfMetCHS.alias = cms.string('pfMetCHS')
addMETCollection(process, labelName='patPFMetCHS', metSource='pfMetCHS')
# RAW CHS MET
process.patPFMetCHS.addGenMET = False


## Slimmed METs
from PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi import slimmedMETs
#### CaloMET is not available in MiniAOD
del slimmedMETs.caloMET

### CHS
process.slimmedMETsCHS = slimmedMETs.clone()
if hasattr(process, "patPFMetCHS"):
     # Create MET from Type 1 PF collection
     process.patPFMetCHS.addGenMET = False
     process.slimmedMETsCHS.src = cms.InputTag("patPFMetCHS")
     process.slimmedMETsCHS.rawUncertainties = cms.InputTag("patPFMetCHS") # only central value
else:
     # Create MET from RAW PF collection
     process.patPFMetCHS.addGenMET = False
     process.slimmedMETsCHS.src = cms.InputTag("patPFMetCHS")
     del process.slimmedMETsCHS.rawUncertainties # not available

clean_met_(process.slimmedMETsCHS)
addMETCollection(process, labelName="slMETsCHS", metSource="slimmedMETsCHS")
process.slMETsCHS.addGenMET = False
#-----------------------------------------------------------

process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                       
                                                       calibratedPatPhotons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                                                 engineName = cms.untracked.string('TRandom3'),
                                                                                           ),
                                                       )
process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
process.calibratedPatPhotons 

##-------------Add quarkGluon tagging---------------------------



process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag("slimmedJets")       # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

#------------------------------------------------------------
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

calo_collection=''
cluster_collection=''
pfcalo_collection=''
   

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',

  # There's no avoiding this in Consumes era
  isData          = cms.bool(False),

  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slMETsCHS'),
  metTypeI         = cms.InputTag('slMETsCHS'),

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
    throw                 = cms.bool(False)
  ),



  ## JECs ################
  redoJECs  = cms.bool(True),

  ## Version Summer15_25nsV3
  L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2Relative_AK8PFchs.txt'),
  L3corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L3Absolute_AK8PFchs.txt'),
  ResCorrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK4PFchs.txt'),
  L1corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK8PFchs.txt'),
  L3corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK8PFchs.txt'),
  L1RCcorr_DATA = cms.FileInPath('CMSDIJET/responsecomputing/DijetRootTreeAnalyzer/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L1RC_AK4PFchs.txt')

)


# ------------------ path --------------------------

process.p = cms.Path()

if runOnRECO:
   process.p += process.pfClusterRefsForJets_step
                                                        
process.p +=                     process.prunedGenParticlesDijet
process.p +=                     process.chs 
process.p +=                     process.slimmedGenJetsAK8                      
process.p +=                     process.dijets 
