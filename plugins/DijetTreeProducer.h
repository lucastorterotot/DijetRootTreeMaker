#ifndef DijetTreeProducer_h
#define DijetTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "TTree.h"
#include "TH1F.h"
#include "TParameter.h"

// For JECs
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

// For 80X updates (for consumes system)
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// EGM smearer package
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

enum JetAlgorithm {
  AK4,
  AK8,
  PUPPI
};

struct JetInfos {
  JetAlgorithm algo;
  edm::EDGetTokenT<pat::JetCollection> inputTag;

};


class DijetTreeProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit DijetTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual bool isValidPhotonLoose(const pat::PhotonRef& photonRef,const edm::Event& event, double generatorWeight);
    virtual bool isValidPhotonMedium(const pat::PhotonRef& photonRef,const edm::Event& event, double generatorWeight);
    virtual bool isValidPhotonTight(const pat::PhotonRef& photonRef,const edm::Event& event, double generatorWeight); 
    virtual void endJob();
    virtual ~DijetTreeProducer();
    

  private: 
    
    void initialize();
    // For JECs
    bool redoJECs_;
    edm::FileInPath L1corrAK4_DATA_, L2corrAK4_DATA_, L3corrAK4_DATA_, ResCorrAK4_DATA_,L1corrPUPPI_DATA_, L2corrPUPPI_DATA_, L3corrPUPPI_DATA_, ResCorrPUPPI_DATA_, L1corrAK8_DATA_, L2corrAK8_DATA_, L3corrAK8_DATA_, ResCorrAK8_DATA_,L1RCcorr_DATA_;
    edm::FileInPath L1corrAK4_MC_, L2corrAK4_MC_, L3corrAK4_MC_,L1corrPUPPI_MC_, L2corrPUPPI_MC_, L3corrPUPPI_MC_, L1corrAK8_MC_, L2corrAK8_MC_, L3corrAK8_MC_;
    JetCorrectorParameters *L1ParAK4_DATA;
    JetCorrectorParameters *L2ParAK4_DATA;
    JetCorrectorParameters *L3ParAK4_DATA;
    JetCorrectorParameters *L2L3ResAK4_DATA;
    FactorizedJetCorrector *JetCorrectorAK4_DATA;
    JetCorrectorParameters *L1ParPUPPI_DATA;
    JetCorrectorParameters *L2ParPUPPI_DATA;
    JetCorrectorParameters *L3ParPUPPI_DATA;
    JetCorrectorParameters *L2L3ResPUPPI_DATA;
    FactorizedJetCorrector *JetCorrectorPUPPI_DATA;
    JetCorrectorParameters *L1JetParForTypeI;
    FactorizedJetCorrector *jetCorrectorForTypeI;    
    JetCorrectorParameters *L1ParAK4_MC;
    JetCorrectorParameters *L2ParAK4_MC;
    JetCorrectorParameters *L3ParAK4_MC;
    FactorizedJetCorrector *JetCorrectorAK4_MC;
    JetCorrectorParameters *L1ParPUPPI_MC;
    JetCorrectorParameters *L2ParPUPPI_MC;
    JetCorrectorParameters *L3ParPUPPI_MC;
    FactorizedJetCorrector *JetCorrectorPUPPI_MC;
    JetCorrectorParameters *L1ParAK8_DATA;
    JetCorrectorParameters *L2ParAK8_DATA;
    JetCorrectorParameters *L3ParAK8_DATA;
    JetCorrectorParameters *L2L3ResAK8_DATA;
    FactorizedJetCorrector *JetCorrectorAK8_DATA;
    JetCorrectorParameters *L1ParAK8_MC;
    JetCorrectorParameters *L2ParAK8_MC;
    JetCorrectorParameters *L3ParAK8_MC;
    FactorizedJetCorrector *JetCorrectorAK8_MC;
    //---- configurable parameters --------   
    double ptMinAK4_,ptMinAK8_;
    double ptMinPhoton_;
    bool mVerbose = false;
    bool isData_;
    bool goodPVtx_;
    
    
    // Migrate to consumes-system for running in 80X

    edm::EDGetTokenT<pat::PhotonCollection> srcPhoton_;
    edm::EDGetTokenT<pat::PhotonCollection> srcPhotonsmeared_;
    edm::EDGetTokenT<pat::ElectronCollection> srcElectron_;
    edm::EDGetTokenT<pat::ElectronCollection> srcElectronsmeared_;
    edm::EDGetTokenT<pat::MuonCollection> srcMuon_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> srcPfCands_;
    
    
    edm::EDGetTokenT<edm::ValueMap<float>> full5x5SigmaIEtaIEtaMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> phoChargedIsolationToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> phoNeutralHadronIsolationToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> phoPhotonIsolationToken_;
    
    edm::InputTag barrelRecHitCollection_;
    edm::InputTag endcapRecHitCollection_;
    
    edm::EDGetTokenT<EcalRecHitCollection> srcebrechit_;
    edm::EDGetTokenT<EcalRecHitCollection> srceerechit_;

    
    edm::EDGetTokenT<pat::JetCollection> srcJetsAK4_;
    edm::EDGetTokenT<edm::View<pat::Jet>> srcJetsAK4View_;
    edm::EDGetTokenT<edm::ValueMap<float>> qgToken ;
    edm::EDGetTokenT<pat::JetCollection> srcJetsAK8_;
    edm::EDGetTokenT<pat::JetCollection> srcJetsPUPPI_;
    edm::EDGetTokenT<edm::View<pat::Jet>> srcJetsAK4puppiView_;
    

    edm::EDGetTokenT<double> srcRho_;
    edm::EDGetTokenT<std::vector<pat::MET> > srcMET_      ;
    edm::EDGetTokenT<std::vector<pat::MET> > srcMETforgen_;
    edm::EDGetTokenT<std::vector<pat::MET> > srcMETpuppi_ ;

    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<reco::VertexCollection> srcVrtx_;

    edm::EDGetTokenT<reco::GenJetCollection> srcGenJetsAK4_;
    edm::EDGetTokenT<reco::GenJetCollection> srcGenJetsAK8_;
    edm::EDGetTokenT<reco::GenParticleCollection> srcPrunedGenParticles_;
    
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > srcPU_;
    edm::EDGetTokenT<GenEventInfoProduct> srcGenInfo_;

    edm::Service<TFileService> fs_;
    TTree *outTree_;

    //---- TRIGGER -------------------------
    edm::EDGetTokenT<edm::TriggerResults> srcTriggerResults_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> srcTriggerPrescale_;
    std::vector<bool>  *triggerResultfromToken_;
    std::vector<double> *triggerPrescale_; 
    std::vector<std::string> *triggerName_;
    
    
    triggerExpression::Data triggerCache_;
    std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
    std::vector<std::string> vtriggerAlias_,vtriggerSelection_;
    TH1F *triggerPassHisto_,*triggerNamesHisto_,*puHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int   run_,evt_,nVtx_,lumi_,BXnumber_;
    int   nJetsAK4_,nJetsPUPPI_, nJetsAK8_, nGenJetsAK4_, nGenJetsAK8_;
    int   nPhotons_, nPhotonsLoose_,nPhotonsMedium_,nPhotonsTight_, nGenphotons_;
    int   nMuonsLoose_;
    TParameter<int>*  Nevent_tot;

    float rho_,metEnergy_,metPt_,metPhi_,metEta_,metEnergypuppi_,metPtpuppi_,metPhipuppi_,metEtapuppi_,metSig_,metcorrected_;
    float metEnergyGen_,metPtGen_,metPhiGen_,metEtaGen_,metEnergypuppiGen_,metPtpuppiGen_,metPhipuppiGen_,metEtapuppiGen_;
    float htAK4_;
    float htAK8_;
    std::vector<bool>  *triggerResult_;
    
    //---- photon variables --------------
    std::vector<float> *ptphoton_,*etaphoton_,*phiphoton_,*energyphoton_,*full5x5SigmaIEtaIEtaMapTokenphoton_,*phoChargedIsolationTokenphoton_,*phoNeutralHadronIsolationTokenphoton_,*phoPhotonIsolationTokenphoton_,*hadTowOverEm_;
    std::vector<float> *ptsmearedphoton_,*etasmearedphoton_,*phismearedphoton_,*energysmearedphoton_,*full5x5SigmaIEtaIEtaMapTokensmearedphoton_,*phoChargedIsolationTokensmearedphoton_,*phoNeutralHadronIsolationTokensmearedphoton_,*phosmearedphotonIsolationTokensmearedphoton_;
    std::vector<bool>  *isPhotonLoose_,*isPhotonMedium_,*isPhotonTight_,  *HaspixelSeed_ , *electronconvVeto_;
    std::vector<float> *ptphotonSC_,*etaphotonSC_,*phiphotonSC_,*energyphotonSC_;
    std::vector<float> *ptGenphoton_,*etaGenphoton_,*phiGenphoton_,*energyGenphoton_;
    std::vector<double> *Ecorrbump_;
    
    
    //  std::vector<float> *PFpx, *PFpy, *PFpz, 
    
    
    std::vector<float> *elecPt_, *elecEta_, *elecPhi_, *elecEnergy_, *elecID_, *elecISO_;
    std::vector<float> *elecPtsmeared_, *elecEtasmeared_, *elecPhismeared_, *elecEnergysmeared_, *elecIDsmeared_, *elecISOsmeared_;    
    std::vector<float> *muPt_, *muEta_, *muPhi_, *muEnergy_;
    
    //------------variable for met correction--------
    FactorizedJetCorrector *jetCorrector;
  //  FactorizedJetCorrector *jetCorrectorForTypeI;
    
    //---- jet and genJet variables --------------
    std::vector<float> *ptAK4_,*jecAK4_,*etaAK4_,*phiAK4_,*massAK4_,*energyAK4_,*areaAK4_,*csvAK4_,*qgdAK4_,*chfAK4_,*nhfAK4_,*phfAK4_,*elfAK4_,*mufAK4_,*nemfAK4_,*cemfAK4_;
    std::vector<float> *ptAK4raw_,*etaAK4raw_,*phiAK4raw_,*massAK4raw_,*energyAK4raw_;
    
    std::vector<int> *idLAK4_,*idTAK4_, *chHadMultAK4_, *chMultAK4_, *neHadMultAK4_, *neMultAK4_, *phoMultAK4_,*pdgIDGenAK4_;
    std::vector<float> *hf_hfAK4_, *hf_emfAK4_, *hofAK4_;
    std::vector<float> *ptGenAK4_,*etaGenAK4_,*phiGenAK4_,*massGenAK4_,*energyGenAK4_;
    
    std::vector<float> *ptPUPPI_,*jecPUPPI_,*etaPUPPI_,*phiPUPPI_,*massPUPPI_,*energyPUPPI_,*areaPUPPI_,*csvPUPPI_,*qgdPUPPI_,*chfPUPPI_,*nhfPUPPI_,*phfPUPPI_,*elfPUPPI_,*mufPUPPI_,*nemfPUPPI_,*cemfPUPPI_;
    std::vector<float> *ptPUPPIraw_,*etaPUPPIraw_,*phiPUPPIraw_,*massPUPPIraw_,*energyPUPPIraw_;
    
    std::vector<int> *idLPUPPI_,*idTPUPPI_, *chHadMultPUPPI_, *chMultPUPPI_, *neHadMultPUPPI_, *neMultPUPPI_, *phoMultPUPPI_,*pdgIDGenPUPPI_;
    std::vector<float> *hf_hfPUPPI_, *hf_emfPUPPI_, *hofPUPPI_;
    std::vector<float> *ptGenPUPPI_,*etaGenPUPPI_,*phiGenPUPPI_,*massGenPUPPI_,*energyGenPUPPI_;
    
    std::vector<float>*ptAK8_,*jecAK8_,*etaAK8_,*phiAK8_,*massAK8_,*energyAK8_,*areaAK8_,*csvAK8_,*qgdAK8_,*chfAK8_,*nhfAK8_,*phfAK8_,*elfAK8_,*mufAK8_,*nemfAK8_,*cemfAK8_, *massPrunedAK8_, *massSoftDropAK8_, *dR_AK8_,*tau1AK8_,*tau2AK8_, *tau3AK8_ ;
    std::vector<int> *idLAK8_,*idTAK8_, *chHadMultAK8_, *chMultAK8_, *neHadMultAK8_, *neMultAK8_, *phoMultAK8_;
    std::vector<float> *hf_hfAK8_, *hf_emfAK8_, *hofAK8_;
    std::vector<float> *ptGenAK8_,*etaGenAK8_,*phiGenAK8_,*massGenAK8_,*energyGenAK8_;

    //---- MC variables ---------------
    std::vector<float> *npu_; 
    std::vector<int> *Number_interactions;
    std::vector <int> *OriginBX;
   
    //-----gen particles from hard scattering ------
  
    float ptHat_; 
    int processID_;
    float weight_;
    std::vector<float> *gen_eta, *gen_phi, *gen_p, *gen_px, *gen_py, *gen_pz, *gen_pt, *gen_energy,  *gen_vx, *gen_vy, *gen_vz;
    std::vector<int> *gen_numDaught, *gen_status, *gen_index, *gen_motherIndex, *gen_pdgId;  

};

#endif
