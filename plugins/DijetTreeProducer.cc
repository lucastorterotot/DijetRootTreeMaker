#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include <map>
#include <unordered_map>
#include "TFile.h"
#include "TH1D.h"
#include "TParameter.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <boost/regex.hpp>
#include "CMSDIJET/DijetRootTreeMaker/plugins/DijetTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
// EGM smearer package
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"




#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"

using namespace std;
using namespace reco;
using namespace pat;
using namespace edm;


DijetTreeProducer::DijetTreeProducer(edm::ParameterSet const&cfg):srcJetsAK4View_(consumes<edm::View<pat::Jet>>(cfg.getParameter<edm::InputTag>("jetsAK4")))
{
  
  // Migrate to Consumes-system. Skip Calo-stuff

  //ADD photons 
  srcPhoton_         = (consumes<pat::PhotonCollection>(cfg.getParameter<InputTag>("Photon")));
  srcPhotonsmeared_  = (consumes<pat::PhotonCollection>(cfg.getParameter<InputTag>("Photonsmeared")));
  
  srcPhotonUncorr_   = (consumes<pat::PhotonCollection>(cfg.getParameter<InputTag>("Photon")));
  ptMinPhoton_       = cfg.getParameter<double>                    ("ptMinPhoton");
  // full5x5SigmaIEtaIEtaMapToken_    =   (consumes <edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap")));
 /* phoChargedIsolationToken_        =   (consumes <edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("phoChargedIsolation")));
  phoNeutralHadronIsolationToken_  =   (consumes <edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("phoNeutralHadronIsolation")));
  phoPhotonIsolationToken_         =   (consumes <edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("phoPhotonIsolation")));*/

  barrelRecHitCollection_ = cfg.getParameter<InputTag>("eb");
  srcebrechit_            = consumes<EcalRecHitCollection>(barrelRecHitCollection_);
  
  endcapRecHitCollection_ = cfg.getParameter<InputTag>("ee");
  srceerechit_            = consumes<EcalRecHitCollection>(endcapRecHitCollection_);
 
  isData_ = cfg.getParameter<bool>("isData");
  isReminiAOD_ = cfg.getParameter<bool>("isreMiniAOD");
  RunEndcapPhoton_ = cfg.getParameter<bool>("Endcaps_photon");
  srcJetsAK4_   = (consumes<pat::JetCollection>(cfg.getParameter<InputTag>    ("jetsAK4")));
//  qgToken       = (consumes<edm::ValueMap<float> >(edm::InputTag              ("QGTagger", "qgLikelihood")));
  srcJetsPUPPI_ = (consumes<pat::JetCollection>(cfg.getParameter<InputTag>    ("jetsPUPPI")));
  
  srcRho_             = (consumes<double>(cfg.getParameter<edm::InputTag>                          ("rho")));
  srcMET_             = (consumes<pat::METCollection >(cfg.getParameter<edm::InputTag>             ("met")));
  srcMETforgen_       = (consumes<pat::METCollection >(cfg.getParameter<edm::InputTag>             ("metforggen")));
  srcMETpuppi_        = (consumes<pat::METCollection >(cfg.getParameter<edm::InputTag>             ("metpuppi")));
  if (isData_ && isReminiAOD_){
  srcMETEGcleaned_    = (consumes<pat::METCollection >(cfg.getParameter<edm::InputTag>             ("metEGcleaned")));
  srcPhotonUncorr_    = (consumes<pat::PhotonCollection>(cfg.getParameter<InputTag>                ("PhotonUncorr")));
  }else {
  
  srcMETEGcleaned_    = (consumes<pat::METCollection >(cfg.getParameter<edm::InputTag>             ("met")));
  }
  srcVrtx_             = (consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>         ("vtx")));
  srcJetsAK4puppiView_ =(consumes<edm::View<pat::Jet>>(cfg.getParameter<edm::InputTag>             ("jetsPUPPI")));
  
  ptMinAK4_           = cfg.getParameter<double>                                                   ("ptMinAK4");
  
  srcElectron_         = (consumes<pat::ElectronCollection>(cfg.getParameter<InputTag>             ("Electrons")));
  srcElectronsmeared_  = (consumes<pat::ElectronCollection>(cfg.getParameter<InputTag>             ("Electronssmeared")));
  srcMuon_             = (consumes<pat::MuonCollection>(cfg.getParameter<InputTag>                 ("Muons")));
  
  srcPU_              = consumes<std::vector<PileupSummaryInfo> >(cfg.getUntrackedParameter<edm::InputTag>    ("pu"));
  srcPfCands_         = consumes<pat::PackedCandidateCollection>(cfg.getParameter<InputTag>    ("PFCands"));

  if (!isData_){
    

    

     srcGenJetsAK4_      = (consumes<GenJetCollection>(cfg.getParameter<edm::InputTag>("genJetsAK4")));
     srcPrunedGenParticles_ = (consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles")));
     srcGenInfo_           = consumes<GenEventInfoProduct>(cfg.getUntrackedParameter<edm::InputTag>  ("ptHat"));
     }

  srcTriggerResults_  = (consumes<edm::TriggerResults >(cfg.getParameter<edm::InputTag>             ("triggerResultsTag")));
  srcTriggerPrescale_ = (consumes<pat::PackedTriggerPrescales >(cfg.getParameter<edm::InputTag>             ("prescalesTag")));
  triggerCache_       = triggerExpression::Data(cfg.getParameterSet("triggerConfiguration"),consumesCollector());
  vtriggerAlias_      = cfg.getParameter<std::vector<std::string> > ("triggerAlias");
  vtriggerSelection_  = cfg.getParameter<std::vector<std::string> > ("triggerSelection");
  triggerObjectsToken = consumes<edm::View<pat::TriggerObjectStandAlone>>(cfg.getParameter<edm::InputTag>("triggerObjects"));
  filters_name        = cfg.getParameter<std::vector<std::string> >("filters");
  

  if (vtriggerAlias_.size() != vtriggerSelection_.size()) {
    cout<<"ERROR: the number of trigger aliases does not match the number of trigger names !!!"<<endl;
    return;
  }
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    vtriggerSelector_.push_back(triggerExpression::parse(vtriggerSelection_[i]));
  }
  
  // For JECs
  redoJECs_ = cfg.getParameter<bool>("redoJECs");
  // AK4 DATA
  L1corrAK4_DATA_ = cfg.getParameter<edm::FileInPath>("L1corrAK4_DATA");
  L2corrAK4_DATA_ = cfg.getParameter<edm::FileInPath>("L2corrAK4_DATA");
  L3corrAK4_DATA_ = cfg.getParameter<edm::FileInPath>("L3corrAK4_DATA");
  ResCorrAK4_DATA_ = cfg.getParameter<edm::FileInPath>("ResCorrAK4_DATA");
  
  // AK4 MC 
  L1corrAK4_MC_ = cfg.getParameter<edm::FileInPath>("L1corrAK4_MC");
  L2corrAK4_MC_ = cfg.getParameter<edm::FileInPath>("L2corrAK4_MC");
  L3corrAK4_MC_ = cfg.getParameter<edm::FileInPath>("L3corrAK4_MC");
  
  L1RCcorr_DATA_ = cfg.getParameter<edm::FileInPath>("L1RCcorr_DATA");
  
  // AK4PUPPI DATA
  L1corrPUPPI_DATA_ = cfg.getParameter<edm::FileInPath>("L1corrAK4PUPPI_DATA");
  L2corrPUPPI_DATA_ = cfg.getParameter<edm::FileInPath>("L2corrAK4PUPPI_DATA");
  L3corrPUPPI_DATA_ = cfg.getParameter<edm::FileInPath>("L3corrAK4PUPPI_DATA");
  ResCorrPUPPI_DATA_ = cfg.getParameter<edm::FileInPath>("ResCorrAK4PUPPI_DATA");
  
  // AK4PUPPI MC 
  L1corrPUPPI_MC_ = cfg.getParameter<edm::FileInPath>("L1corrAK4PUPPI_MC");
  L2corrPUPPI_MC_ = cfg.getParameter<edm::FileInPath>("L2corrAK4PUPPI_MC");
  L3corrPUPPI_MC_ = cfg.getParameter<edm::FileInPath>("L3corrAK4PUPPI_MC");
  
  
  if(redoJECs_)
  {
    // AK4
    L1ParAK4_DATA = new JetCorrectorParameters(L1corrAK4_DATA_.fullPath());
    L2ParAK4_DATA = new JetCorrectorParameters(L2corrAK4_DATA_.fullPath());
    L3ParAK4_DATA = new JetCorrectorParameters(L3corrAK4_DATA_.fullPath());
    L2L3ResAK4_DATA = new JetCorrectorParameters(ResCorrAK4_DATA_.fullPath());
    L1ParAK4_MC = new JetCorrectorParameters(L1corrAK4_MC_.fullPath());
    L2ParAK4_MC = new JetCorrectorParameters(L2corrAK4_MC_.fullPath());
    L3ParAK4_MC = new JetCorrectorParameters(L3corrAK4_MC_.fullPath());
    L1JetParForTypeI = new JetCorrectorParameters(L1RCcorr_DATA_.fullPath());

    std::vector<JetCorrectorParameters> vParAK4_DATA;
    std::vector<JetCorrectorParameters> vParAK4_MC;
    std::vector<JetCorrectorParameters> vParTypeI;
    vParAK4_DATA.push_back(*L1ParAK4_DATA);
    vParAK4_DATA.push_back(*L2ParAK4_DATA);
    vParAK4_DATA.push_back(*L3ParAK4_DATA);
    vParAK4_DATA.push_back(*L2L3ResAK4_DATA);
    
    vParTypeI.push_back(*L1JetParForTypeI);
    
    vParAK4_MC.push_back(*L1ParAK4_MC);
    vParAK4_MC.push_back(*L2ParAK4_MC);
    vParAK4_MC.push_back(*L3ParAK4_MC);

    JetCorrectorAK4_DATA = new FactorizedJetCorrector(vParAK4_DATA);
    JetCorrectorAK4_MC = new FactorizedJetCorrector(vParAK4_MC);
    jetCorrectorForTypeI = new FactorizedJetCorrector(vParTypeI);
    
    // AK4PUPPI
    L1ParPUPPI_DATA = new JetCorrectorParameters(L1corrPUPPI_DATA_.fullPath());
    L2ParPUPPI_DATA = new JetCorrectorParameters(L2corrPUPPI_DATA_.fullPath());
    L3ParPUPPI_DATA = new JetCorrectorParameters(L3corrPUPPI_DATA_.fullPath());
    L2L3ResPUPPI_DATA = new JetCorrectorParameters(ResCorrPUPPI_DATA_.fullPath());
    L1ParPUPPI_MC = new JetCorrectorParameters(L1corrPUPPI_MC_.fullPath());
    L2ParPUPPI_MC = new JetCorrectorParameters(L2corrPUPPI_MC_.fullPath());
    L3ParPUPPI_MC = new JetCorrectorParameters(L3corrPUPPI_MC_.fullPath());
    L1JetParForTypeI = new JetCorrectorParameters(L1RCcorr_DATA_.fullPath());

    std::vector<JetCorrectorParameters> vParPUPPI_DATA;
    std::vector<JetCorrectorParameters> vParPUPPI_MC;
    vParPUPPI_DATA.push_back(*L1ParPUPPI_DATA);
    vParPUPPI_DATA.push_back(*L2ParPUPPI_DATA);
    vParPUPPI_DATA.push_back(*L3ParPUPPI_DATA);
    vParPUPPI_DATA.push_back(*L2L3ResPUPPI_DATA);
    
    
    vParPUPPI_MC.push_back(*L1ParPUPPI_MC);
    vParPUPPI_MC.push_back(*L2ParPUPPI_MC);
    vParPUPPI_MC.push_back(*L3ParPUPPI_MC);

    JetCorrectorPUPPI_DATA = new FactorizedJetCorrector(vParPUPPI_DATA);
    JetCorrectorPUPPI_MC = new FactorizedJetCorrector(vParPUPPI_MC);

  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enum class IsolationType {
  CHARGED_HADRONS,
    NEUTRAL_HADRONS,
    PHOTONS
    };
    

float getEffectiveArea(float eta, IsolationType type) {
  eta = fabs(eta);
  switch (type) {
  case IsolationType::CHARGED_HADRONS:
    if (eta < 1.0)
      return 0.0112;
    else if (eta < 1.479)
      return 0.0108;
    else if (eta < 2.0)
      return 0.0106;
    else if (eta < 2.2)
      return 0.01002;
    else if (eta < 2.3)
      return 0.0098;
    else if (eta < 2.4)
      return 0.0089;
    else
      return 0.0087;
    break;
    
  case IsolationType::NEUTRAL_HADRONS:
    if (eta < 1.0)
      return 0.0668;
    else if (eta < 1.479)
      return 0.1054;
    else if (eta < 2.0)
      return 0.0786;
    else if (eta < 2.2)
      return 0.0233;
    else if (eta < 2.3)
      return 0.0078;
    else if (eta < 2.4)
      return 0.0028;
    else
      return 0.0137;
    break;

    //Official  
  case IsolationType::PHOTONS:
    if (eta < 1.0)
      return 0.1113;
    else if (eta < 1.479)
      return 0.0953;
    else if (eta < 2.0)
      return 0.0619;
    else if (eta < 2.2)
      return 0.0837;
    else if (eta < 2.3)
      return 0.1070;
    else if (eta < 2.4)
      return 0.1212;
    else
      return 0.1466;
    break;
  }
  
  return -1;
}

double getCorrectedPFIsolation(double isolation, double rho, float eta, IsolationType type) 
{
  float effectiveArea = getEffectiveArea(eta, type); 
  return std::max(isolation - rho*effectiveArea, 0.);
}


//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  
  
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    triggerNamesHisto_->Fill(vtriggerSelection_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  
  
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("BXnumber"             ,&BXnumber_          ,"BXnumber_/I");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("metEnergy"            ,&metEnergy_         ,"metEnergy_/F");
  outTree_->Branch("metPt"                ,&metPt_             ,"metPt_/F");
  outTree_->Branch("metEta"               ,&metEta_            ,"metEta_/F");
  outTree_->Branch("metPhi"               ,&metPhi_            ,"metPhi_/F");
  outTree_->Branch("metEnergyPUPPI"            ,&metEnergypuppi_         ,"metEnergypuppi_/F");
  outTree_->Branch("metPtPUPPI"                ,&metPtpuppi_             ,"metPtpuppi_/F");
  outTree_->Branch("metEtaPUPPI"               ,&metEtapuppi_            ,"metEtapuppi_/F");
  outTree_->Branch("metPhiPUPPI"               ,&metPhipuppi_            ,"metPhipuppi_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  outTree_->Branch("metTypeI"             ,&metcorrected_      ,"metcorrected_/F");
  outTree_->Branch("goodPVtx"             ,&goodPVtx_      ,"goodPVtx_/B");
  outTree_->Branch("NgoodPV"             ,&NgoodPV_      ,"NgoodPV_/I");
  
  
  
  outTree_->Branch("deltaNHfootprintX"               ,&deltaNHfootprintX_            ,"deltaNHfootprintX_/F");
  outTree_->Branch("deltaNHfootprintY"               ,&deltaNHfootprintY_            ,"deltaNHfootprintY_/F");

  
  outTree_->Branch("metEnergyGen"            ,&metEnergyGen_         ,"metEnergyGen_/F");
  outTree_->Branch("metPtGen"                ,&metPtGen_             ,"metPtGen_/F");
  outTree_->Branch("metEtaGen"               ,&metEtaGen_            ,"metEtaGen_/F");
  outTree_->Branch("metPhiGen"               ,&metPhiGen_            ,"metPhiGen_/F");
  outTree_->Branch("metEnergyPUPPIGen"            ,&metEnergypuppiGen_         ,"metEnergypuppiGen_/F");
  outTree_->Branch("metPtPUPPIGen"                ,&metPtpuppiGen_             ,"metPtpuppiGen_/F");
  outTree_->Branch("metEtaPUPPIGen"               ,&metEtapuppiGen_            ,"metEtapuppiGen_/F");
  outTree_->Branch("metPhiPUPPIGen"               ,&metPhipuppiGen_            ,"metPhipuppiGen_/F");
  
  
  
  
  gen_eta          = new std::vector<float>;
  gen_phi          = new std::vector<float>;
  gen_p            = new std::vector<float>;
  gen_px           = new std::vector<float>;
  gen_py           = new std::vector<float>;
  gen_pz           = new std::vector<float>;
  gen_pt           = new std::vector<float>;
  gen_energy	   = new std::vector<float>; 
  gen_pdgId	   = new std::vector<int>; 
  gen_vx	   = new std::vector<float>; 
  gen_vy	   = new std::vector<float>; 
  gen_vz	   = new std::vector<float>;	 
  gen_numDaught    = new std::vector<int>;      
  gen_status	   = new std::vector<int>; 
  gen_index	   = new std::vector<int>; 
  gen_motherIndex  = new std::vector<int>; 

  outTree_->Branch("gen_eta"		,"vector<float>" , &gen_eta      	);
  outTree_->Branch("gen_phi"		,"vector<float>" , &gen_phi 		);
  outTree_->Branch("gen_p"		,"vector<float>",   &gen_p 		);
  outTree_->Branch("gen_px"		,"vector<float>",  &gen_px 		);
  outTree_->Branch("gen_py"		,"vector<float>",  &gen_py 		);
  outTree_->Branch("gen_pz"		,"vector<float>",  &gen_pz 		);
  outTree_->Branch("gen_pt"		,"vector<float>",  &gen_pt 		);
  outTree_->Branch("gen_energy"	    	,"vector<float>",  &gen_energy		);
  outTree_->Branch("gen_pdgId"	    	,"vector<int>",   &gen_pdgId  		);
  outTree_->Branch("gen_vx"	    	,"vector<float>", &gen_vx		);
  outTree_->Branch("gen_vy"	    	,"vector<float>", &gen_vy		);
  outTree_->Branch("gen_vz"	    	,"vector<float>", &gen_vz		);
  outTree_->Branch("gen_numDaught"  	,"vector<int>",    &gen_numDaught       );
  outTree_->Branch("gen_status"	    	,"vector<int>",  &gen_status     	);
  outTree_->Branch("gen_index"	    	,"vector<int>",  &gen_index      	);
  outTree_->Branch("gen_motherIndex"	,"vector<int>", &gen_motherIndex 	);
  outTree_->Branch("nJetsAK4"           ,&nJetsAK4_          ,"nJetsAK4_/I"		);
  outTree_->Branch("nJetsPUPPI"           ,&nJetsPUPPI_          ,"nJetsPUPPI_/I"		);
  outTree_->Branch("htAK4"              ,&htAK4_             ,"htAK4_/F"		);
  
  outTree_->Branch("nPhoton"                  ,&nPhotons_                ,"nPhotons_/I");
  outTree_->Branch("nPhotonLoose"             ,&nPhotonsLoose_           ,"nPhotonsLoose_/I");
  outTree_->Branch("nPhotonMedium"            ,&nPhotonsMedium_          ,"nPhotonsMedium_/I");
  outTree_->Branch("nPhotonTight"             ,&nPhotonsTight_           ,"nPhotonsTight_/I");
  //------------------------------------------------------------- --------
  ptphoton_            = new std::vector<float>;
  ptsmearedphoton_     = new std::vector<float>;
  ptphotonSC_          = new std::vector<float>;  
  etaphoton_           = new std::vector<float>;
  etasmearedphoton_    = new std::vector<float>;
  etaphotonSC_         = new std::vector<float>;  
  phiphoton_           = new std::vector<float>;
  phismearedphoton_    = new std::vector<float>;
  phiphotonSC_          = new std::vector<float>;
  energyphoton_        = new std::vector<float>;
  energysmearedphoton_ = new std::vector<float>;
  energyphotonSC_       = new std::vector<float>; 
  HaspixelSeed_      = new std::vector<bool>;
  electronconvVeto_  = new std::vector<bool>;
  hadTowOverEm_      = new std::vector<float>;  
  isPhotonLoose_       = new std::vector<bool>;
  isPhotonMedium_      = new std::vector<bool>;
  isPhotonTight_       = new std::vector<bool>;
  Ecorrbump_           = new std::vector<double>;
  
  isMatch30_        = new std::vector<bool>;
  isMatch50_        = new std::vector<bool>;
  isMatch75_        = new std::vector<bool>;
  isMatch90_        = new std::vector<bool>;
  isMatch120_       = new std::vector<bool>;
  isMatch165_       = new std::vector<bool>;
  isGenMatch_       = new std::vector<bool>;
  
  outTree_->Branch("isMatch30"                                 ,"vector<bool>"   ,&isMatch30_);
  outTree_->Branch("isMatch50"                                 ,"vector<bool>"   ,&isMatch50_);
  outTree_->Branch("isMatch75"                                 ,"vector<bool>"   ,&isMatch75_);
  outTree_->Branch("isMatch90"                                 ,"vector<bool>"   ,&isMatch90_);
  outTree_->Branch("isMatch120"                                ,"vector<bool>"   ,&isMatch120_);
  outTree_->Branch("isMatch165"                                ,"vector<bool>"   ,&isMatch165_);
  outTree_->Branch("isGenMatch"                                ,"vector<bool>"   ,&isGenMatch_);
      
 
  full5x5SigmaIEtaIEtaMapTokenphoton_   = new std::vector<float>;
  phoChargedIsolationTokenphoton_       = new std::vector<float>;
  phoNeutralHadronIsolationTokenphoton_ = new std::vector<float>;
  phoPhotonIsolationTokenphoton_        = new std::vector<float>;
  
  
  outTree_->Branch("PhotonLoosePt"         ,"vector<float>"   ,&ptphoton_);
  outTree_->Branch("PhotonsmearPt"         ,"vector<float>"   ,&ptsmearedphoton_);
  outTree_->Branch("PhotonSCPt"           ,"vector<float>"   ,&ptphotonSC_);  
  outTree_->Branch("PhotonLooseEta"        ,"vector<float>"   ,&etaphoton_);
  outTree_->Branch("PhotonsmearEta"        ,"vector<float>"   ,&etasmearedphoton_);
  outTree_->Branch("PhotonSCEta"           ,"vector<float>"   ,&etaphotonSC_);    
  outTree_->Branch("PhotonLoosePhi"        ,"vector<float>"   ,&phiphoton_);
  outTree_->Branch("PhotonsmearPhi"        ,"vector<float>"   ,&phismearedphoton_);
  outTree_->Branch("PhotonSCPhi"           ,"vector<float>"   ,&phiphotonSC_); 
  outTree_->Branch("PhotonLooseEnergy"     ,"vector<float>"   ,&energyphoton_);
  outTree_->Branch("PhotonsmearEnergy"     ,"vector<float>"   ,&energysmearedphoton_);
  outTree_->Branch("PhotonSCEnergy"        ,"vector<float>"   ,&energyphotonSC_);
  outTree_->Branch("PhotonEcorrBump"        ,"vector<double>" ,&Ecorrbump_); 
  outTree_->Branch("Photonfull5x5SigmaIEtaIEtaMapToken"           ,"vector<float>"   ,&full5x5SigmaIEtaIEtaMapTokenphoton_);
  outTree_->Branch("PhotonphoChargedIsolationToken"               ,"vector<float>"   ,&phoChargedIsolationTokenphoton_);
  outTree_->Branch("PhotonphoNeutralHadronIsolationToken"         ,"vector<float>"   ,&phoNeutralHadronIsolationTokenphoton_);
  outTree_->Branch("PhotonphoPhotonIsolationToken"                ,"vector<float>"   ,&phoPhotonIsolationTokenphoton_);
  outTree_->Branch("HaspixelSeed"                                 ,"vector<bool>"   ,&HaspixelSeed_);
  outTree_->Branch("ElectronVeto"                                 ,"vector<bool>"   ,&electronconvVeto_);
  outTree_->Branch("hadTowOverEm"                                 ,"vector<float>"   ,&hadTowOverEm_);
  outTree_->Branch("isPhotonLoose"                                 ,"vector<bool>"   ,&isPhotonLoose_);
  outTree_->Branch("isPhotonMedium"                                ,"vector<bool>"   ,&isPhotonMedium_);
  outTree_->Branch("isPhotonTight"                                 ,"vector<bool>"   ,&isPhotonTight_);
  
  
  
  
  //--------------Photon scale ---------------
  
 photon_scaleUNC_gainup_= new std::vector<float>;
 photon_scaleUNC_gaindown_= new std::vector<float>; 
 photon_scaleUNC_systup_= new std::vector<float>; 
 photon_scaleUNC_systdown_= new std::vector<float>; 
 photon_scaleUNC_statup_= new std::vector<float>; 
 photon_scaleUNC_statdown_= new std::vector<float>; 
 photon_scaleUNC_ETup_= new std::vector<float>; 
 photon_scaleUNC_ETdown_= new std::vector<float>; 
 photon_smearUNC_phiup_= new std::vector<float>; 
 photon_smearUNC_rhodown_= new std::vector<float>; 
 photon_smearUNC_rhoup_= new std::vector<float>; 
 photon_smear_central_= new std::vector<float>; 
 photon_scale_central_= new std::vector<float>; 
     
  outTree_->Branch("Photon_scaleUNC_gainup"         ,"vector<float>"   ,&photon_scaleUNC_gainup_);
  outTree_->Branch("Photon_scaleUNC_gaindown"         ,"vector<float>"   ,&photon_scaleUNC_gaindown_);
  outTree_->Branch("Photon_scaleUNC_systup"         ,"vector<float>"   ,&photon_scaleUNC_systup_);
  outTree_->Branch("Photon_scaleUNC_systdown"         ,"vector<float>"   ,&photon_scaleUNC_systdown_);
  outTree_->Branch("Photon_scaleUNC_statup"         ,"vector<float>"   ,&photon_scaleUNC_statup_);
  outTree_->Branch("Photon_scaleUNC_statdown"         ,"vector<float>"   ,&photon_scaleUNC_statdown_);
  outTree_->Branch("Photon_scaleUNC_ETup"         ,"vector<float>"   ,&photon_scaleUNC_ETup_);
  outTree_->Branch("Photon_scaleUNC_ETdown"         ,"vector<float>"   ,&photon_scaleUNC_ETdown_);
  outTree_->Branch("Photon_smearUNC_phiup"         ,"vector<float>"   ,&photon_smearUNC_phiup_);
  outTree_->Branch("Photon_smearUNC_rhodown"         ,"vector<float>"   ,&photon_smearUNC_rhodown_);
  outTree_->Branch("Photon_smearUNC_rhoup"         ,"vector<float>"   ,&photon_smearUNC_rhoup_);
  
  outTree_->Branch("photon_scale_central"         ,"vector<float>"   ,&photon_scale_central_);
  outTree_->Branch("photon_smear_central"         ,"vector<float>"   ,&photon_smear_central_);
  
  //---------------Electrons----------------------------------
  
  
  elecPt_        = new std::vector<float>;
  elecEta_       = new std::vector<float>;
  elecPhi_      = new std::vector<float>;
  elecEnergy_    = new std::vector<float>;
  elecID_        = new std::vector<float>;
  elecISO_        = new std::vector<float>;
  
  elecPtsmeared_        = new std::vector<float>;
  elecEtasmeared_       = new std::vector<float>;
  elecPhismeared_      = new std::vector<float>;
  elecEnergysmeared_    = new std::vector<float>;
  elecIDsmeared_        = new std::vector<float>;
  elecISOsmeared_        = new std::vector<float>;
  
  
  outTree_->Branch("electronPt"         ,"vector<float>"   ,&elecPt_);
  outTree_->Branch("electronEta"        ,"vector<float>"   ,&elecEta_);
  outTree_->Branch("electronPhi"        ,"vector<float>"   ,&elecPhi_);
  outTree_->Branch("electronEnergy"     ,"vector<float>"   ,&elecEnergy_);
  outTree_->Branch("electronID"         ,"vector<float>"   ,&elecID_);
  outTree_->Branch("electronISO"        ,"vector<float>"   ,&elecISO_);
  
  outTree_->Branch("electronPtsmeared"         ,"vector<float>"   ,&elecPtsmeared_);
  outTree_->Branch("electronEtasmeared"        ,"vector<float>"   ,&elecEtasmeared_);
  outTree_->Branch("electronPhismeared"        ,"vector<float>"   ,&elecPhismeared_);
  outTree_->Branch("electronEnergysmeared"     ,"vector<float>"   ,&elecEnergysmeared_);
  outTree_->Branch("electronIDsmeared"         ,"vector<float>"   ,&elecIDsmeared_);
  outTree_->Branch("electronISOsmeared"        ,"vector<float>"   ,&elecISOsmeared_);
  
  
  
  //-------------Muons-------------------------------
  
  
  muPt_       = new std::vector<float>;
  muEta_      = new std::vector<float>;
  muPhi_      = new std::vector<float>;
  muEnergy_   = new std::vector<float>;
  
  outTree_->Branch("muonPt"         ,"vector<float>"   ,&muPt_);
  outTree_->Branch("muonEta"        ,"vector<float>"   ,&muEta_);
  outTree_->Branch("muonPhi"        ,"vector<float>"   ,&muPhi_);
  outTree_->Branch("muonEnergy"     ,"vector<float>"   ,&muEnergy_);
  outTree_->Branch("nMuonsLoose"             ,&nMuonsLoose_          ,"nMuonsLoose_/I");
  //------------------------------------------------------------------
  ptAK4_             = new std::vector<float>;
  jecAK4_            = new std::vector<float>;
  etaAK4_            = new std::vector<float>;
  phiAK4_            = new std::vector<float>;
  massAK4_           = new std::vector<float>;
  energyAK4_         = new std::vector<float>;  
  ptAK4raw_             = new std::vector<float>;
  etaAK4raw_            = new std::vector<float>;
  phiAK4raw_            = new std::vector<float>;
  massAK4raw_           = new std::vector<float>;
  energyAK4raw_         = new std::vector<float>;  
  areaAK4_           = new std::vector<float>;
  csvAK4_            = new std::vector<float>;
  qgdAK4_            = new std::vector<float>; 
  chfAK4_            = new std::vector<float>;
  nhfAK4_            = new std::vector<float>;
  phfAK4_            = new std::vector<float>;
  mufAK4_            = new std::vector<float>;
  elfAK4_            = new std::vector<float>;
  nemfAK4_           = new std::vector<float>;
  cemfAK4_           = new std::vector<float>;
  // Hadronic forward hadrons
  hf_hfAK4_          = new std::vector<float>;
  // Hadronic forward electromagnetic fraction
  hf_emfAK4_         = new std::vector<float>;
  hofAK4_            = new std::vector<float>;
  idLAK4_            = new std::vector<int>;
  idTAK4_            = new std::vector<int>;
  chHadMultAK4_     = new std::vector<int>;
  chMultAK4_         = new std::vector<int>;
  neHadMultAK4_      = new std::vector<int>;
  neMultAK4_         = new std::vector<int>;
  phoMultAK4_        = new std::vector<int>;
  hadronflavour_ = new std::vector<int>;  


  outTree_->Branch("jetPtAK4"                ,"vector<float>"     ,&ptAK4_);
  outTree_->Branch("jetJecAK4"               ,"vector<float>"     ,&jecAK4_);
  outTree_->Branch("jetEtaAK4"               ,"vector<float>"     ,&etaAK4_);
  outTree_->Branch("jetPhiAK4"               ,"vector<float>"     ,&phiAK4_);
  outTree_->Branch("jetMassAK4"              ,"vector<float>"     ,&massAK4_);
  outTree_->Branch("jetEnergyAK4"            ,"vector<float>"     ,&energyAK4_);  
  outTree_->Branch("jetPtAK4RC"                ,"vector<float>"     ,&ptAK4raw_);
  outTree_->Branch("jetEtaAK4RC"               ,"vector<float>"     ,&etaAK4raw_);
  outTree_->Branch("jetPhiAK4RC"               ,"vector<float>"     ,&phiAK4raw_);
  outTree_->Branch("jetMassAK4RC"              ,"vector<float>"     ,&massAK4raw_);
  outTree_->Branch("jetEnergyAK4RC"            ,"vector<float>"     ,&energyAK4raw_);  
  outTree_->Branch("jetAreaAK4"              ,"vector<float>"     ,&areaAK4_);
  outTree_->Branch("jetCSVAK4"               ,"vector<float>"     ,&csvAK4_);  
  outTree_->Branch("jetQGDAK4"               ,"vector<float>"     ,&qgdAK4_);  
  outTree_->Branch("jetChfAK4"               ,"vector<float>"     ,&chfAK4_);
  outTree_->Branch("jetNhfAK4"               ,"vector<float>"     ,&nhfAK4_);
  outTree_->Branch("jetPhfAK4"               ,"vector<float>"     ,&phfAK4_);
  outTree_->Branch("jetMufAK4"               ,"vector<float>"     ,&mufAK4_);
  outTree_->Branch("jetElfAK4"               ,"vector<float>"     ,&elfAK4_);
  outTree_->Branch("jetNemfAK4"              ,"vector<float>"     ,&nemfAK4_);
  outTree_->Branch("jetCemfAK4"              ,"vector<float>"     ,&cemfAK4_);
  outTree_->Branch("jetHf_hfAK4"             ,"vector<float>"     ,&hf_hfAK4_);
  outTree_->Branch("jetHf_emfAK4"            ,"vector<float>"    ,&hf_emfAK4_);
  outTree_->Branch("jetHofAK4"               ,"vector<float>"    ,&hofAK4_);
  outTree_->Branch("idLAK4"                  ,"vector<int>"      ,&idLAK4_);   
  outTree_->Branch("idTAK4"                  ,"vector<int>"      ,&idTAK4_);   
  outTree_->Branch("chHadMultAK4"           ,"vector<int>"      ,&chHadMultAK4_);   
  outTree_->Branch("chMultAK4"              ,"vector<int>"      ,&chMultAK4_);   
  outTree_->Branch("neHadMultAK4"           ,"vector<int>"      ,&neHadMultAK4_);   
  outTree_->Branch("neMultAK4"              ,"vector<int>"      ,&neMultAK4_);   
  outTree_->Branch("phoMultAK4"             ,"vector<int>"      ,&phoMultAK4_);   
  outTree_->Branch("hadronflavour"          ,"vector<int>"      ,&hadronflavour_);   
  
  //-------Jet PUPPI-----------------
  ptPUPPI_             = new std::vector<float>;
  jecPUPPI_            = new std::vector<float>;
  etaPUPPI_            = new std::vector<float>;
  phiPUPPI_            = new std::vector<float>;
  massPUPPI_           = new std::vector<float>;
  energyPUPPI_         = new std::vector<float>;  
  ptPUPPIraw_             = new std::vector<float>;
  etaPUPPIraw_            = new std::vector<float>;
  phiPUPPIraw_            = new std::vector<float>;
  massPUPPIraw_           = new std::vector<float>;
  energyPUPPIraw_         = new std::vector<float>;  
  areaPUPPI_           = new std::vector<float>;
  csvPUPPI_            = new std::vector<float>;
  qgdPUPPI_            = new std::vector<float>; 
  chfPUPPI_            = new std::vector<float>;
  nhfPUPPI_            = new std::vector<float>;
  phfPUPPI_            = new std::vector<float>;
  mufPUPPI_            = new std::vector<float>;
  elfPUPPI_            = new std::vector<float>;
  nemfPUPPI_           = new std::vector<float>;
  cemfPUPPI_           = new std::vector<float>;
  // Hadronic forward hadrons
  hf_hfPUPPI_          = new std::vector<float>;
  // Hadronic forward electromagnetic fraction
  hf_emfPUPPI_         = new std::vector<float>;
  hofPUPPI_            = new std::vector<float>;
  idLPUPPI_            = new std::vector<int>;
  idTPUPPI_            = new std::vector<int>;
  chHadMultPUPPI_     = new std::vector<int>;
  chMultPUPPI_         = new std::vector<int>;
  neHadMultPUPPI_      = new std::vector<int>;
  neMultPUPPI_         = new std::vector<int>;
  phoMultPUPPI_        = new std::vector<int>;
  
  outTree_->Branch("jetPtPUPPI"                ,"vector<float>"     ,&ptPUPPI_);  
  outTree_->Branch("jetJecPUPPI"               ,"vector<float>"     ,&jecPUPPI_);
  outTree_->Branch("jetEtaPUPPI"               ,"vector<float>"     ,&etaPUPPI_);
  outTree_->Branch("jetPhiPUPPI"               ,"vector<float>"     ,&phiPUPPI_);
  outTree_->Branch("jetMassPUPPI"              ,"vector<float>"     ,&massPUPPI_);
  outTree_->Branch("jetEnergyPUPPI"            ,"vector<float>"     ,&energyPUPPI_);  
  outTree_->Branch("jetPtPUPPIRC"                ,"vector<float>"     ,&ptPUPPIraw_);
  outTree_->Branch("jetEtaPUPPIRC"               ,"vector<float>"     ,&etaPUPPIraw_);
  outTree_->Branch("jetPhiPUPPIRC"               ,"vector<float>"     ,&phiPUPPIraw_);
  outTree_->Branch("jetMassPUPPIRC"              ,"vector<float>"     ,&massPUPPIraw_);
  outTree_->Branch("jetEnergyPUPPIRC"            ,"vector<float>"     ,&energyPUPPIraw_);  
  outTree_->Branch("jetAreaPUPPI"              ,"vector<float>"     ,&areaPUPPI_);
  outTree_->Branch("jetCSVPUPPI"               ,"vector<float>"     ,&csvPUPPI_);  
  outTree_->Branch("jetQGDPUPPI"               ,"vector<float>"     ,&qgdPUPPI_);  
  outTree_->Branch("jetChfPUPPI"               ,"vector<float>"     ,&chfPUPPI_);
  outTree_->Branch("jetNhfPUPPI"               ,"vector<float>"     ,&nhfPUPPI_);
  outTree_->Branch("jetPhfPUPPI"               ,"vector<float>"     ,&phfPUPPI_);
  outTree_->Branch("jetMufPUPPI"               ,"vector<float>"     ,&mufPUPPI_);
  outTree_->Branch("jetElfPUPPI"               ,"vector<float>"     ,&elfPUPPI_);
  outTree_->Branch("jetNemfPUPPI"              ,"vector<float>"     ,&nemfPUPPI_);
  outTree_->Branch("jetCemfPUPPI"              ,"vector<float>"     ,&cemfPUPPI_);
  outTree_->Branch("jetHf_hfPUPPI"             ,"vector<float>"     ,&hf_hfPUPPI_);
  outTree_->Branch("jetHf_emfPUPPI"            ,"vector<float>"    ,&hf_emfPUPPI_);
  outTree_->Branch("jetHofPUPPI"               ,"vector<float>"    ,&hofPUPPI_);
  outTree_->Branch("idLPUPPI"                  ,"vector<int>"      ,&idLPUPPI_);   
  outTree_->Branch("idTPUPPI"                  ,"vector<int>"      ,&idTPUPPI_);   
  outTree_->Branch("chHadMultPUPPI"          ,"vector<int>"      ,&chHadMultPUPPI_);   
  outTree_->Branch("chMultPUPPI"              ,"vector<int>"      ,&chMultPUPPI_);   
  outTree_->Branch("neHadMultPUPPI"           ,"vector<int>"      ,&neHadMultPUPPI_);   
  outTree_->Branch("neMultPUPPI"              ,"vector<int>"      ,&neMultPUPPI_);   
  outTree_->Branch("phoMultPUPPI"             ,"vector<int>"      ,&phoMultPUPPI_); 
  
  
  
  //----------end PUPPI----------------

  triggerResult_  = new std::vector<bool>;
  triggerPrescale_ = new std::vector<int>;
  triggerName_     = new std::vector<std::string>;
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);
  outTree_->Branch("triggerPrescale","vector<int>",&triggerPrescale_);
  outTree_->Branch("triggerName","vector<string>",&triggerName_);



  //------------------- MC ---------------------------------
  npu_                = new std::vector<float>;  
  Number_interactions = new std::vector<int>;
  OriginBX            = new std::vector<int>;
 
  outTree_->Branch("npu"                  ,"vector<float>"     , &npu_ );
  outTree_->Branch("PileupInteractions"   ,"vector<int>"       , &Number_interactions );
  outTree_->Branch("PileupOriginBX"       ,"vector<int>"       , &OriginBX );
  outTree_->Branch("ptHat"                ,&ptHat_             ,"ptHat_/F");
  outTree_->Branch("processID"            ,&processID_         ,"processID_/I");
  outTree_->Branch("weight"               ,&weight_            ,"weight_/F");

  outTree_->Branch("nGenJetsAK4"             ,&nGenJetsAK4_          ,"nGenJetsAK4_/I");
  outTree_->Branch("nGenPhoton"             ,&nGenphotons_          ,"nGenPhotons_/I");
  
  
  
  //------------------GenPhotons----------------------------------
  ptGenphoton_            = new std::vector<float>;
  etaGenphoton_           = new std::vector<float>;
  phiGenphoton_           = new std::vector<float>;
  energyGenphoton_        = new std::vector<float>;
  outTree_->Branch("photonPtGen"                 ,"vector<float>"     ,&ptGenphoton_);
  outTree_->Branch("photonEtaGen"                ,"vector<float>"     ,&etaGenphoton_);
  outTree_->Branch("photonPhiGen"                ,"vector<float>"     ,&phiGenphoton_);
  outTree_->Branch("photonEnergyGen"             ,"vector<float>"     ,&energyGenphoton_);
  
  
  

  //------------------GenJets-------------------------------------
  ptGenAK4_             = new std::vector<float>;
  etaGenAK4_            = new std::vector<float>;
  phiGenAK4_            = new std::vector<float>;
  massGenAK4_           = new std::vector<float>;
  pdgIDGenAK4_          = new std::vector<int>;
  energyGenAK4_         = new std::vector<float>;

  outTree_->Branch("jetPtGenAK4"                ,"vector<float>"     ,&ptGenAK4_);
  outTree_->Branch("jetEtaGenAK4"               ,"vector<float>"     ,&etaGenAK4_);
  outTree_->Branch("jetPhiGenAK4"               ,"vector<float>"     ,&phiGenAK4_);
  outTree_->Branch("jetMassGenAK4"              ,"vector<float>"     ,&massGenAK4_);
  outTree_->Branch("jetEnergyGenAK4"            ,"vector<float>"     ,&energyGenAK4_);
  outTree_->Branch("jetpdgIDGenAK4"             ,"vector<int>"     ,&pdgIDGenAK4_);
  
  ptGenPUPPI_             = new std::vector<float>;
  etaGenPUPPI_            = new std::vector<float>;
  phiGenPUPPI_            = new std::vector<float>;
  massGenPUPPI_           = new std::vector<float>;
  pdgIDGenPUPPI_          = new std::vector<int>;
  energyGenPUPPI_         = new std::vector<float>;

  outTree_->Branch("jetPtGenPUPPI"                ,"vector<float>"     ,&ptGenPUPPI_);
  outTree_->Branch("jetEtaGenPUPPI"               ,"vector<float>"     ,&etaGenPUPPI_);
  outTree_->Branch("jetPhiGenPUPPI"               ,"vector<float>"     ,&phiGenPUPPI_);
  outTree_->Branch("jetMassGenPUPPI"              ,"vector<float>"     ,&massGenPUPPI_);
  outTree_->Branch("jetEnergyGenPUPPI"            ,"vector<float>"     ,&energyGenPUPPI_);
  outTree_->Branch("jetpdgIDGenPUPPI"             ,"vector<int>"     ,&pdgIDGenPUPPI_);
  
  
     

}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::endJob() 
{  
  delete triggerResult_;
  delete triggerPrescale_;
  delete triggerName_;


  delete gen_eta	;
  delete gen_phi	;
  delete gen_p		;
  delete gen_px	;
  delete gen_py	;
  delete gen_pz	;
  delete gen_pt	;
  delete gen_energy    ;
  delete gen_pdgId	;
  delete gen_vx	;
  delete gen_vy	;
  delete gen_vz	;
  delete gen_numDaught	;
  delete gen_status	;
  delete gen_index   	;
  delete gen_motherIndex;

  delete ptAK4_;
  delete jecAK4_;
  delete etaAK4_;
  delete phiAK4_;
  delete massAK4_;
  delete energyAK4_;
  
  delete ptAK4raw_;
  delete etaAK4raw_;
  delete phiAK4raw_;
  delete massAK4raw_;
  delete energyAK4raw_;
  
  delete areaAK4_;
  delete csvAK4_;
  delete qgdAK4_;
  delete chfAK4_;
  delete nhfAK4_;
  delete phfAK4_;
  delete mufAK4_;
  delete elfAK4_;
  delete nemfAK4_;
  delete cemfAK4_;
  delete hf_hfAK4_;
  delete hf_emfAK4_;
  delete hofAK4_;
  delete idLAK4_;
  delete idTAK4_;
  delete chHadMultAK4_ ;
  delete chMultAK4_    ;
  delete neHadMultAK4_ ;
  delete neMultAK4_    ;
  delete phoMultAK4_   ;
  delete hadronflavour_;

 delete ptPUPPI_;
 delete jecPUPPI_;
  delete etaPUPPI_;
  delete phiPUPPI_;
  delete massPUPPI_;
  delete energyPUPPI_;
  
  delete ptPUPPIraw_;
  delete etaPUPPIraw_;
  delete phiPUPPIraw_;
  delete massPUPPIraw_;
  delete energyPUPPIraw_;
  
  delete areaPUPPI_;
  delete csvPUPPI_;
  delete qgdPUPPI_;
  delete chfPUPPI_;
  delete nhfPUPPI_;
  delete phfPUPPI_;
  delete mufPUPPI_;
  delete elfPUPPI_;
  delete nemfPUPPI_;
  delete cemfPUPPI_;
  delete hf_hfPUPPI_;
  delete hf_emfPUPPI_;
  delete hofPUPPI_;
  delete idLPUPPI_;
  delete idTPUPPI_;
  delete chHadMultPUPPI_ ;
  delete chMultPUPPI_    ;
  delete neHadMultPUPPI_ ;
  delete neMultPUPPI_    ;
  delete phoMultPUPPI_   ;

  delete ptphoton_         ;
  delete ptsmearedphoton_  ;
  delete ptphotonSC_       ;

  delete etaphoton_        ;
  delete etasmearedphoton_ ;
  delete etaphotonSC_      ;

  delete phiphoton_        ;
  delete phismearedphoton_ ;
  delete phiphotonSC_      ;

  delete energyphoton_     ;
  delete energysmearedphoton_ ;
  delete energyphotonSC_    ;

  delete ptGenphoton_      ;
  delete etaGenphoton_     ;
  delete phiGenphoton_     ;
  delete energyGenphoton_  ;
  delete full5x5SigmaIEtaIEtaMapTokenphoton_    ;
  delete phoChargedIsolationTokenphoton_        ;
  delete phoNeutralHadronIsolationTokenphoton_  ;
  delete phoPhotonIsolationTokenphoton_         ;
  delete isPhotonLoose_                         ;
  delete isPhotonMedium_                        ;
  delete isPhotonTight_                         ;
  delete HaspixelSeed_                          ;
  delete hadTowOverEm_                          ;
  delete pdgIDGenAK4_                           ;
  delete pdgIDGenPUPPI_                         ;
  delete electronconvVeto_                      ; 
  
 delete elecPt_ ;       
 delete elecEta_  ;     
 delete elecPhi_;     
 delete elecEnergy_ ;
 delete elecID_;
 delete elecISO_;   
  
 delete elecPtsmeared_ ;       
 delete elecEtasmeared_  ;     
 delete elecPhismeared_  ;    
 delete elecEnergysmeared_ ;
 delete elecIDsmeared_;
 delete elecISOsmeared_;  
  
  
 delete muPt_   ;    
 delete muEta_  ;    
 delete muPhi_  ;    
 delete muEnergy_;
 delete Ecorrbump_;
 
 delete isMatch30_     ;
 delete isMatch50_     ;    
 delete isMatch75_     ;    
 delete isMatch90_     ;    
 delete isMatch120_    ;    
 delete isMatch165_    ;
 delete isGenMatch_    ;
 
 delete photon_scaleUNC_gainup_;
 delete photon_scaleUNC_gaindown_; 
 delete photon_scaleUNC_systup_; 
 delete photon_scaleUNC_systdown_; 
 delete photon_scaleUNC_statup_; 
 delete photon_scaleUNC_statdown_; 
 delete photon_scaleUNC_ETup_; 
 delete photon_scaleUNC_ETdown_; 
 delete photon_smearUNC_phiup_; 
 delete photon_smearUNC_rhodown_; 
 delete photon_smearUNC_rhoup_ ; 
 delete photon_smear_central_ ; 
 delete photon_scale_central_;
 
 
// delete ptphotonnofix_;
  
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }

}
//////////////////////////////////////////////////////////////////////////////////////////




void DijetTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  double generatorWeight = 1.;
  initialize();
 
  
  Handle<pat::JetCollection> jetsAK4;
  iEvent.getByToken(srcJetsAK4_,jetsAK4);
  pat::JetCollection jetrawRC = *jetsAK4;
  
  Handle<pat::JetCollection> jetsPUPPI;
  iEvent.getByToken(srcJetsPUPPI_,jetsPUPPI);
  pat::JetCollection jetrawRCpuppi = *jetsPUPPI;
  
  edm::Handle<edm::View<pat::Jet>> jetsview;     
  iEvent.getByToken(srcJetsAK4View_,jetsview);
  

  edm::Handle<edm::View<pat::Jet>> jetsviewpuppi;     
  iEvent.getByToken(srcJetsAK4puppiView_,jetsviewpuppi);
  
  
  
 /* edm::Handle<edm::ValueMap<float>> qgHandle;
  iEvent.getByToken(qgToken, qgHandle);*/


  Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(srcPhoton_,photons);
  
  Handle<pat::PhotonCollection> photonsUncorr;
  iEvent.getByToken(srcPhotonUncorr_,photonsUncorr);
  
  
  Handle<pat::PhotonCollection> photonssmear;
  iEvent.getByToken(srcPhotonsmeared_,photonssmear);
  
 
  
//   //------------------ Genphoton ----------------------------------- 

  Handle<reco::GenJetCollection> handle_genJetsAK4;
  if (!iEvent.isRealData())
    iEvent.getByToken(srcGenJetsAK4_,handle_genJetsAK4);

  Handle<double>  rho;
  iEvent.getByToken(srcRho_,rho);
  
  

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(srcVrtx_,recVtxs);
  if(!recVtxs.isValid() || recVtxs->size() == 0 ||  recVtxs->front().isFake()) goodPVtx_ = false ;
  
  NgoodPV_ = 0;
  for(reco::VertexCollection::const_iterator ivec = recVtxs->begin(); ivec != recVtxs->end() ; ++ivec ){
  
  if(!ivec->isFake() && ivec->ndof() > 4 && ivec->z() <= 24 && ivec->position().rho() <= 2)NgoodPV_ ++;
  //pv_rho[i] <= 2 && pv_z[i] <= 24 && pv_ndof[i] > 4
  }
  
  
  
  const reco::Vertex& primaryVertex = recVtxs->at(0);
  
  edm::Handle<EcalRecHitCollection> _ebrechits;
  iEvent.getByToken(srcebrechit_,_ebrechits);
  

  
  edm::Handle< EcalRecHitCollection> _eerechits;
  iEvent.getByToken(srceerechit_,_eerechits);
  
  
  //-------------- Event Info -----------------------------------
  rho_    = *rho;


  nVtx_   = recVtxs->size();
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();
  BXnumber_ = iEvent.bunchCrossing();

  //---------- pu -----------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (!iEvent.isRealData()) {
    iEvent.getByToken(srcPU_,PupInfo);
    

    if(PupInfo.isValid()) {
      for( std::vector<PileupSummaryInfo>::const_iterator it = PupInfo->begin(); it != PupInfo->end(); ++it ) {
	npu_ -> push_back ( it -> getTrueNumInteractions() );
	Number_interactions -> push_back ( it->getPU_NumInteractions() ); 
	OriginBX -> push_back ( it -> getBunchCrossing());                
	
      }
    }
    else {
      cout << "an edm::LogError call for PileUpError used to be here, but that does not work anymore -Juska" << endl;
    }
  }// if MC
  
  //-------------- Gen Event Info -----------------------------------
  if (!iEvent.isRealData()) {

    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(srcGenInfo_,genEvtInfo);
    
    if( !genEvtInfo.isValid() ) {
      edm::LogInfo("GenEvtInfo") << "ERROR: genEvtInfo not valid! " << genEvtInfo;
    }

    if( genEvtInfo.isValid() ) {
      edm::LogInfo("GenEvtInfo") << "Successfully obtained " << genEvtInfo;
      ptHat_ = (genEvtInfo->hasBinningValues() ? genEvtInfo->binningValues()[0] : -999.);
      processID_ = genEvtInfo->signalProcessID();
      weight_ = genEvtInfo->weight();            
    }
      generatorWeight = genEvtInfo->weight();
      if (generatorWeight == 0.) {
	generatorWeight = 1.;
      }

    //------------------ Gen particles hard scattering -------------------
    
    // to be saved only for partons that start the jet -> from genJets take the costituents -> 
    //see hypernews https://hypernews.cern.ch/HyperNews/CMS/get/csa14/49/2.html
    //and https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Advanced_topics_re_clustering_ev 
    
   
      
    
    edm::Handle<reco::GenParticleCollection> prunedGenParticles;
    if (!iEvent.isRealData())
      iEvent.getByToken(srcPrunedGenParticles_, prunedGenParticles);
    

    
    if( prunedGenParticles.isValid() ) {
            
      for( reco::GenParticleCollection::const_iterator it = prunedGenParticles->begin(); it != prunedGenParticles->end(); ++it ) {
       
    	//save only particles from hard scattering 
	//already done from the pruner

    	int idx = std::distance(prunedGenParticles->begin(),it);

        // fill in all the vectors
        gen_eta		->push_back( it->eta() );
        gen_phi		->push_back( it->phi() );
        gen_p		->push_back( it->p() );
        gen_px		->push_back( it->px() );
        gen_py		->push_back( it->py() );
        gen_pz		->push_back( it->pz() );
        gen_pt		->push_back( it->pt() );
        gen_energy	->push_back( it->energy() );
        gen_pdgId	->push_back( it->pdgId() );
        gen_vx		->push_back( it->vx() );
        gen_vy		->push_back( it->vy() );
        gen_vz		->push_back( it->vz() );
        gen_numDaught	->push_back( it->numberOfDaughters() );
        gen_status	->push_back( it->status() );
    	gen_index   	->push_back( idx );
  
    	int midx = -1;

	for( reco::GenParticleCollection::const_iterator mit = prunedGenParticles->begin(); mit != prunedGenParticles->end(); ++mit ) {
	
    	  if( it->mother()==&(*mit) ) {
    	    midx = std::distance(prunedGenParticles->begin(),mit);
    	    break;
    	  }
    	}
    	gen_motherIndex->push_back( midx ); 

      }//loop over genParticles

    }
    
  }// if MC
  
  //-------------- Trigger Info -----------------------------------
  Handle<edm::TriggerResults> Trigger_result;
  iEvent.getByToken(srcTriggerResults_,Trigger_result);
  
  const edm::TriggerResults Static_trig = (*Trigger_result.product());
  Handle<pat::PackedTriggerPrescales> Trigger_prescale;
  iEvent.getByToken(srcTriggerPrescale_,Trigger_prescale);
  
  size_t sizetrigger = Trigger_result->size();
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*Trigger_result);
  int prescalefactor = 1 ;
  std::string triggerName ;
  std::string ValidTriggerregex;
  std::vector<std::string> ValidTrigger ;

  for(unsigned itrig=0;itrig<vtriggerSelector_.size();itrig++){
         ValidTriggerregex = vtriggerSelection_[itrig];
         if(ValidTriggerregex.size() == 32) ValidTriggerregex.replace(31,2,".*");
         if(ValidTriggerregex.size() == 33) ValidTriggerregex.replace(32,1,".*");
         ValidTrigger.push_back(ValidTriggerregex);
        }
  
  
  triggerPassHisto_->Fill("totalEvents",1);
  
  if (triggerCache_.setEvent(iEvent,iSetup)) {
    for(unsigned itrig=0;itrig<vtriggerSelector_.size();itrig++) {
      bool result(false);
      if (vtriggerSelector_[itrig]) {
        if (triggerCache_.configurationUpdated()) {
          vtriggerSelector_[itrig]->init(triggerCache_);
        }

        result = (*(vtriggerSelector_[itrig]))(triggerCache_);

        for(unsigned int i = 0 ; i<sizetrigger; i++){
        
         triggerName = triggerNames.triggerName(i);         
         std::vector<boost::regex> validTriggers = { boost::regex(ValidTrigger.at(itrig), boost::regex_constants::icase) };
                  
         for (boost::regex& validTrigger: validTriggers) {
	  if (boost::regex_match(triggerName, validTrigger)) {	 
           unsigned int index = triggerNames.triggerIndex(triggerName);
           prescalefactor = Trigger_prescale->getPrescaleForIndex(index);	   
	  }
         } 
        }
      }
      if (result) {
        triggerPassHisto_->Fill(vtriggerAlias_[itrig].c_str(),1);
      }
      triggerResult_->push_back(result);
      triggerPrescale_->push_back(prescalefactor);
      triggerName_->push_back(vtriggerSelection_[itrig].c_str());
    }
  }
  
  //------- store the object that launch the last filter in the HLT photon
  
   edm::Handle<edm::View<pat::TriggerObjectStandAlone>> triggerObjects;
   iEvent.getByToken(triggerObjectsToken, triggerObjects);
   
   std::vector<TLorentzVector> candhlt30;
   std::vector<TLorentzVector> candhlt50;
   std::vector<TLorentzVector> candhlt75;
   std::vector<TLorentzVector> candhlt90;
   std::vector<TLorentzVector> candhlt120;
   std::vector<TLorentzVector> candhlt165;
   
   for (auto /*const*/ &obj: *triggerObjects)
    {
        TLorentzVector cand;
        cand.SetPtEtaPhiE(obj.pt(),obj.eta(),obj.phi(),obj.energy());
        pat::TriggerObjectStandAlone To_unpack= obj ;        
        To_unpack.unpackFilterLabels(iEvent,Static_trig);
        if(To_unpack.hasFilterLabel(filters_name.at(0)))
        {
        	candhlt30.push_back(cand);
        }
        
        if(To_unpack.hasFilterLabel(filters_name.at(1)))
        {
        	candhlt50.push_back(cand);
        }
        if(To_unpack.hasFilterLabel(filters_name.at(2)))
        {
        	candhlt75.push_back(cand);
        }
        if(To_unpack.hasFilterLabel(filters_name.at(3)))
        {
        	candhlt90.push_back(cand);
        }
        if(To_unpack.hasFilterLabel(filters_name.at(4)))
        {
        	candhlt120.push_back(cand);
        }
        if(To_unpack.hasFilterLabel(filters_name.at(5)))
        {
        	candhlt165.push_back(cand);
        }
     
}
  
  
      
  //----------------------electron-----------------------
  
    Handle<pat::ElectronCollection> electron;
    iEvent.getByToken(srcElectron_,electron);
    
    Handle<pat::ElectronCollection> electronsmeared;
    iEvent.getByToken(srcElectronsmeared_,electronsmeared);
  
   pat::ElectronCollection::const_iterator ielectron = electron->begin();
   pat::ElectronCollection::const_iterator ielectronsmeared = electronsmeared->begin();
  int ie = 0;
  for(; ielectron != electron->end() ;++ielectron,++ielectronsmeared, ie++){
  
   if (ie >= 30)
   break;
   bool elecID = fabs(primaryVertex.z() - ielectron->vertex().z()) < 1.;
    elecID     &= ielectron->et() > 30.;
    elecID     &= fabs(ielectron->eta()) < 2.5 && (ielectron->superCluster()->eta() > 1.4442 && ielectron->superCluster()->eta() < 1.5660);
    elecID     &= ielectron->dB() < 0.02;
    elecID     &= ((int) ielectron->electronID("cutBasedElectronID-Fall17-94X-V2-loose") & 0x1);
    
    elecID_           ->push_back(elecID);
    elecIDsmeared_    ->push_back(elecID);
    
    float iso = (ielectron->dr03TkSumPt() + ielectron->dr03EcalRecHitSumEt() + ielectron->dr03HcalTowerSumEt()) / ielectron->et();
    
    elecISO_   ->push_back(iso);
    elecISOsmeared_   ->push_back(iso);
    
    elecPt_      ->push_back(ielectron->pt());
    elecEta_     ->push_back(ielectron->eta());
    elecPhi_     ->push_back(ielectron->phi());
    elecEnergy_  ->push_back(ielectron->energy());
    
    elecPtsmeared_      ->push_back(ielectronsmeared->pt());
    elecEtasmeared_     ->push_back(ielectronsmeared->eta());
    elecPhismeared_     ->push_back(ielectronsmeared->phi());
    elecEnergysmeared_  ->push_back(ielectronsmeared->energy());
    
  }
  
  //---------------------------muons
  
   Handle<pat::MuonCollection> muon;
   iEvent.getByToken(srcMuon_,muon);
   pat::MuonCollection::const_iterator imuon= muon->begin();
   nMuonsLoose_ = 0;
   for(; imuon != muon->end() ;++imuon){
  
     if(imuon->isLooseMuon())
     {
       nMuonsLoose_++;
       muPt_            ->push_back(imuon->pt());
       muEta_           ->push_back(imuon->eta());
       muPhi_           ->push_back(imuon->phi());
       muEnergy_        ->push_back(imuon->energy());
     }
  
  
    }
  //------------------------photons---------------

      uint32_t index = 0;
      pat::PhotonCollection::const_iterator iphoton = photons->begin();
      pat::PhotonCollection::const_iterator iphotonsmear= photonssmear->begin();
      pat::PhotonCollection::const_iterator iphotonsuncorr= photonsUncorr->begin();
      pat::Photon pho;
      nPhotons_ = 0;
      nPhotonsLoose_ = 0;
      nPhotonsMedium_ = 0;
      nPhotonsTight_ = 0;
      nGenphotons_=0;
      
    //  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
     // iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
     // edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
      //iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    //  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
     // iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    //  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    //  iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);     
      double rhod = *rho;


     std::vector<pat::Photon> PhotonT_vec;
     std::vector<pat::Photon> PhotonTOoB_vec ;
     pat::Photon PhotonT ;
     pat::Photon PhotonTOoB ;
      for(; iphoton != photons->end(); ++iphoton,++iphotonsmear, index++ ) 
      { 
           nPhotons_++;
           pho =*iphoton;
           PhotonT_vec.push_back((*iphoton));
          if(!RunEndcapPhoton_){ 
	   if (fabs(iphoton->eta()) <= 1.3) 
	   {               
	       pat::PhotonRef PhotonReftmp(photons, index);	      
	        if (isValidPhotonLoose(PhotonReftmp, iEvent, generatorWeight)) 
		{			
		  
		         
		  nPhotonsLoose_++;                           
                  ptphoton_             ->push_back( iphoton->pt()         );
		  ptphotonSC_           ->push_back( iphoton->superCluster()->rawEnergy()/ cosh(iphoton->superCluster()->eta()));
                  phiphoton_            ->push_back( iphoton->phi()        );
		  phiphotonSC_          ->push_back( iphoton->superCluster()->phi());


                  etaphoton_            ->push_back( iphoton->eta()        );
		  etaphotonSC_          ->push_back( iphoton->superCluster()->eta());

                  energyphoton_         ->push_back( iphoton->energy() );
		  energyphotonSC_       ->push_back( iphoton->superCluster()->rawEnergy() );
		  full5x5SigmaIEtaIEtaMapTokenphoton_   ->push_back(iphoton->full5x5_sigmaIetaIeta());
		  phoChargedIsolationTokenphoton_       ->push_back(iphoton-> userFloat("phoChargedIsolation"));
		  phoNeutralHadronIsolationTokenphoton_ ->push_back(getCorrectedPFIsolation(iphoton-> userFloat("phoNeutralHadronIsolation"), rhod, PhotonReftmp->eta(), IsolationType::NEUTRAL_HADRONS));
		  phoPhotonIsolationTokenphoton_        ->push_back(getCorrectedPFIsolation(iphoton-> userFloat("phoChargedIsolation"), rhod, PhotonReftmp->eta(), IsolationType::PHOTONS));
		  
		  isPhotonLoose_        ->push_back( true) ;
		  HaspixelSeed_         ->push_back(iphoton->hasPixelSeed());
		  hadTowOverEm_         ->push_back(iphoton->hadronicOverEm());
		  electronconvVeto_     ->push_back(iphoton->passElectronVeto());		  				  
		  double Ecorr=1;
		  Ecorrbump_            ->push_back( Ecorr                 );
		  phismearedphoton_     ->push_back( iphotonsmear->phi()   );		  
		  ptsmearedphoton_      ->push_back( iphotonsmear->pt()    );		  
		  etasmearedphoton_     ->push_back( iphotonsmear->eta()   );		  
		  energysmearedphoton_  ->push_back( iphotonsmear->energy());
		  
		  photon_scaleUNC_gainup_->push_back(iphoton-> userFloat("energyScaleGainUp") );
                  photon_scaleUNC_gaindown_->push_back(iphoton-> userFloat("energyScaleGainDown") );
                  photon_scaleUNC_systup_->push_back(iphoton-> userFloat("energyScaleSystUp") );
                  photon_scaleUNC_systdown_->push_back(iphoton-> userFloat("energyScaleSystDown") ); 
                  photon_scaleUNC_statup_->push_back(iphoton-> userFloat("energyScaleStatUp") ); 
                  photon_scaleUNC_statdown_->push_back(iphoton-> userFloat("energyScaleStatDown") ); 
                  photon_scaleUNC_ETup_->push_back(iphoton-> userFloat("energyScaleEtUp") ); 
                  photon_scaleUNC_ETdown_->push_back(iphoton-> userFloat("energyScaleEtDown") ); 
                  photon_scale_central_->push_back(iphoton-> userFloat("energyScaleValue") );
		  
		  
		  		  
		  if (!iEvent.isRealData() ) {
		  
		  
		        photon_smearUNC_phiup_->push_back(iphoton-> userFloat("energySigmaPhiUp") ); 
                        photon_smearUNC_rhodown_->push_back(iphoton-> userFloat("energySigmaRhoDown") );
                        photon_smearUNC_rhoup_ ->push_back(iphoton-> userFloat("energySigmaRhoUp") ); 
                        photon_smear_central_ ->push_back(iphoton->userFloat("energySigmaValue") );
                        
		  	if(iphoton->genPhoton()){
		  		ptGenphoton_     ->push_back( iphoton->genPhoton()->pt() );
                  		phiGenphoton_    ->push_back( iphoton->genPhoton()->phi() );
                  		etaGenphoton_    ->push_back( iphoton->genPhoton()->eta() );
                  		energyGenphoton_ ->push_back( iphoton->genPhoton()->energy());
                  		isGenMatch_      ->push_back(true);
	          		nGenphotons_++;
	          	}else{
	          		ptGenphoton_     ->push_back( -999 );
                  		phiGenphoton_    ->push_back( -999 );
                  		etaGenphoton_    ->push_back( -999 );
                  		energyGenphoton_ ->push_back( -999 );
                  		isGenMatch_      ->push_back(false);
	          	}
		  }
		  
                }
		
		if (isValidPhotonLoose(PhotonReftmp, iEvent, generatorWeight) && !isValidPhotonMedium(PhotonReftmp, iEvent, generatorWeight) ) 
		{ 
		  isPhotonMedium_     ->push_back(0);
		}
		
	       if (isValidPhotonLoose(PhotonReftmp, iEvent, generatorWeight) && !isValidPhotonTight(PhotonReftmp, iEvent, generatorWeight) ) 
		{ 
		  isPhotonTight_     ->push_back(0);
		}
		
		
		if (isValidPhotonMedium(PhotonReftmp, iEvent, generatorWeight) ) 
		{
		    nPhotonsMedium_++;
		    isPhotonMedium_     ->push_back(true);
		    
		    
		 }
		if (isValidPhotonTight(PhotonReftmp, iEvent, generatorWeight)) 
		{
		     nPhotonsTight_++;
		     PhotonT = (*iphoton);
		     PhotonTOoB = (*iphotonsmear);
		     PhotonT_vec.push_back((*iphoton));
                     PhotonTOoB_vec.push_back((*iphotonsmear)); 
		     isPhotonTight_   ->push_back(true); 
		 }             
          }}else{
             if (fabs(iphoton->eta()) >= 1.305 && fabs(iphoton->eta()) <= 2.5) 
	   {               
	       pat::PhotonRef PhotonReftmp(photons, index);	      
	        if (isValidEndcapPhotonLoose(PhotonReftmp, iEvent, generatorWeight)) 
		{			
		  
		         
		  nPhotonsLoose_++;                           
                  ptphoton_             ->push_back( iphoton->pt()         );
		  ptphotonSC_           ->push_back( iphoton->superCluster()->rawEnergy()/ cosh(iphoton->superCluster()->eta()));
                  phiphoton_            ->push_back( iphoton->phi()        );
		  phiphotonSC_          ->push_back( iphoton->superCluster()->eta());


                  etaphoton_            ->push_back( iphoton->eta()        );
		  etaphotonSC_          ->push_back( iphoton->superCluster()->phi());

                  energyphoton_         ->push_back( iphoton->energy() );
		  energyphotonSC_       ->push_back( iphoton->superCluster()->rawEnergy() );

		  full5x5SigmaIEtaIEtaMapTokenphoton_   ->push_back(iphoton->full5x5_sigmaIetaIeta());
		  phoChargedIsolationTokenphoton_       ->push_back(iphoton-> userFloat("phoChargedIsolation"));
		  phoNeutralHadronIsolationTokenphoton_ ->push_back(getCorrectedPFIsolation(iphoton-> userFloat("phoNeutralHadronIsolation"), rhod, PhotonReftmp->eta(), IsolationType::NEUTRAL_HADRONS));
		  phoPhotonIsolationTokenphoton_        ->push_back(getCorrectedPFIsolation(iphoton-> userFloat("phoChargedIsolation"), rhod, PhotonReftmp->eta(), IsolationType::PHOTONS));
		  
		  isPhotonLoose_        ->push_back( true) ;
		  HaspixelSeed_         ->push_back(iphoton->hasPixelSeed());
		  hadTowOverEm_         ->push_back(iphoton->hadronicOverEm());
		  electronconvVeto_     ->push_back(iphoton->passElectronVeto());		  				  
		  double Ecorr=1;
		  Ecorrbump_            ->push_back( Ecorr                 );
		  phismearedphoton_     ->push_back( iphotonsmear->phi()   );		  
		  ptsmearedphoton_      ->push_back( iphotonsmear->pt()    );		  
		  etasmearedphoton_     ->push_back( iphotonsmear->eta()   );		  
		  energysmearedphoton_  ->push_back( iphotonsmear->energy());		  
		  if (!iEvent.isRealData() ) {
		  	if(iphoton->genPhoton()){
		  		ptGenphoton_     ->push_back( iphoton->genPhoton()->pt() );
                  		phiGenphoton_    ->push_back( iphoton->genPhoton()->phi() );
                  		etaGenphoton_    ->push_back( iphoton->genPhoton()->eta() );
                  		energyGenphoton_ ->push_back( iphoton->genPhoton()->energy());
                  		isGenMatch_      ->push_back(true);
	          		nGenphotons_++;
	          	}else{
	          		ptGenphoton_     ->push_back( -999 );
                  		phiGenphoton_    ->push_back( -999 );
                  		etaGenphoton_    ->push_back( -999 );
                  		energyGenphoton_ ->push_back( -999 );
                  		isGenMatch_      ->push_back(false);
	          	}
		  }
		  
                }
		
		if (isValidEndcapPhotonLoose(PhotonReftmp, iEvent, generatorWeight) && !isValidEndcapPhotonMedium(PhotonReftmp, iEvent, generatorWeight) ) 
		{ 
		  isPhotonMedium_     ->push_back(0);
		}
		
	       if (isValidEndcapPhotonLoose(PhotonReftmp, iEvent, generatorWeight) && !isValidEndcapPhotonTight(PhotonReftmp, iEvent, generatorWeight) ) 
		{ 
		  isPhotonTight_     ->push_back(0);
		}
		
		
		if (isValidEndcapPhotonMedium(PhotonReftmp, iEvent, generatorWeight) ) 
		{
		    nPhotonsMedium_++;
		    isPhotonMedium_     ->push_back(true);
		    
		    
		 }

		if (isValidEndcapPhotonTight(PhotonReftmp, iEvent, generatorWeight)) 
		{
		     nPhotonsTight_++;
		     PhotonT = (*iphoton);
		     PhotonTOoB = (*iphotonsmear);
		     PhotonT_vec.push_back((*iphoton));
                     PhotonTOoB_vec.push_back((*iphotonsmear)); 
		     isPhotonTight_   ->push_back(true); 
		 }
          
          }}
      }//end loop over photon collection
      
      
      edm::Handle<pat::PackedCandidateCollection> pfs;
      iEvent.getByToken(srcPfCands_, pfs);
      
      
      Handle<vector<pat::MET> > met;
      iEvent.getByToken(srcMET_,met);
      
      Handle<vector<pat::MET> > rawmet;
      iEvent.getByToken(srcMET_,rawmet);
  
      Handle<vector<pat::MET> > metforgenchs;
      iEvent.getByToken(srcMETforgen_,metforgenchs);
      
      Handle<vector<pat::MET> > metEGcleaned;
      iEvent.getByToken(srcMETEGcleaned_, metEGcleaned);
      
  
     Handle<vector<pat::MET> > metpuppi;
     iEvent.getByToken(srcMETpuppi_,metpuppi);
     
     
     pat::MET Met = (*met)[0];
     
     pat::MET Metraw = (*met)[0];
     Metraw.setP4( Met.uncorP4() );
     
     
     pat::MET rawMet = (*rawmet)[0];
     rawMet.setP4( Met.uncorP4() );
     
     pat::MET rawMet74 = (*rawmet)[0];
     rawMet74.setP4( Met.uncorP4() );
     
     
     pat::MET EGcleanedMet = (*metEGcleaned)[0];
     EGcleanedMet.setP4( EGcleanedMet.uncorP4() );
     pat::MET EGcleanedMetRAW = (*metEGcleaned)[0];
     EGcleanedMetRAW.setP4(EGcleanedMet.uncorP4());
      
      if(nPhotonsTight_ == 1){
      double deltar = 0 ;
      int nlim = 0 ;
      pat::Photon Photonuncorr ;
           for(;iphotonsuncorr != photonsUncorr->end(); ++iphotonsuncorr){
            
            
            deltar = std::hypot((PhotonT.eta()-iphotonsuncorr->eta()),(PhotonT.phi()-iphotonsuncorr->phi())   );
            if(deltar <= 0.1 && nlim < 1){
            nlim ++ ;  
            Photonuncorr = (*iphotonsuncorr);
            }   
      
      
          }
      
      
      
      
       float FootprintMEx = 0;
       float FootprintMEy = 0;
       
       float FootprintfromNHMEx = 0;
       float FootprintfromNHMEy = 0;
       
       float FootprintMEx74 = 0;
       float FootprintMEy74 = 0;
       
  
  std::vector<reco::CandidatePtr> footprint;

  for (unsigned int i = 0, n = Photonuncorr.numberOfSourceCandidatePtrs(); i < n; ++i) {
    footprint.push_back(Photonuncorr.sourceCandidatePtr(i) );
  }
  
  std::vector<reco::CandidatePtr> footprint_rereco;

  for (unsigned int i = 0, n = PhotonTOoB.numberOfSourceCandidatePtrs(); i < n; ++i) {
    footprint_rereco.push_back(PhotonTOoB.sourceCandidatePtr(i) );
  }
  
  
  
  // now loop on pf candidates
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate &pf = (*pfs)[i];
    // pfcandidate-based footprint removal 
    
    if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end())  continue;
    
    
    if(std::hypot((PhotonT.eta()-pf.eta()),(PhotonT.phi()-pf.phi())  < 0.3  && pf.pdgId() == 130)){
        FootprintfromNHMEx += -1. * pf.px();
        FootprintfromNHMEy += -1. * pf.py();
        
    }// need to improve with removal of neutral hadron and photon id particle
    
    


    if(pf.fromPV() == 0) continue;

    FootprintMEx += -1.* pf.px();
    FootprintMEy += -1.* pf.py();
    
    FootprintMEx74 += -1.* pf.px();
    FootprintMEy74 += -1.* pf.py();
      
  }// loop over pfCand
  
  
  // Re-adding  photon but reco 
  if(footprint.size() > 0 ){
  FootprintMEx += -1.* PhotonT.px();//
  FootprintMEy += -1.* PhotonT.py();//
  FootprintMEx74 += -1.* PhotonTOoB.px();
  FootprintMEy74 += -1.* PhotonTOoB.py();
  
  }


  
  double FootprintMEPt = sqrt(FootprintMEx * FootprintMEx + FootprintMEy * FootprintMEy) ;
  double FootprintMEpt74 = sqrt(FootprintMEx74 * FootprintMEx74 + FootprintMEy74 * FootprintMEy74) ;
  
rawMet.setP4(reco::Candidate::LorentzVector(FootprintMEx, FootprintMEy, 0., FootprintMEPt));
rawMet74.setP4(reco::Candidate::LorentzVector(FootprintMEx74, FootprintMEy74, 0., FootprintMEpt74));
   
  metEnergy_    = rawMet.energy();
  metEta_       = rawMet.eta(); 
  metPhi_       = rawMet.phi();
  metPt_        = rawMet.pt();
  
  deltaNHfootprintX_ = FootprintfromNHMEx;
  deltaNHfootprintY_ = FootprintfromNHMEy;
  
  if(!iEvent.isRealData() &&(*metforgenchs)[0].genMET ())
  {
      metEnergyGen_    = rawMet.genMET ()->energy();
      metEtaGen_       = rawMet.genMET ()->eta();      
      metPhiGen_       = rawMet.genMET ()->phi();
      metPtGen_        = rawMet.genMET ()->pt();
  }
  
  pat::MET MetrawPuppi = (*metpuppi)[0];
  MetrawPuppi.setP4( MetrawPuppi.uncorP4() );
  
  metEnergypuppi_    = MetrawPuppi.energy();
  metEtapuppi_       = MetrawPuppi.eta();      
  metPhipuppi_       = MetrawPuppi.phi();
  metPtpuppi_        = MetrawPuppi.pt(); 
  if(!iEvent.isRealData() &&(*metpuppi)[0].genMET ())
  {
     metEnergypuppiGen_    = (*metpuppi)[0].genMET ()->energy();
     metEtapuppiGen_       = (*metpuppi)[0].genMET ()->eta();      
     metPhipuppiGen_       = (*metpuppi)[0].genMET ()->phi();
     metPtpuppiGen_        = (*metpuppi)[0].genMET ()->pt();
  }
  
  
  
  
  if ((*met)[0].sumEt() > 0) {
    metSig_ = (*met)[0].et()/(*met)[0].sumEt();
  }
      
      
    //---- match the photon tight to the trigger object;
    if(candhlt30.size() != 0 || candhlt50.size() != 0 || candhlt75.size() != 0 || candhlt90.size() != 0 || candhlt120.size() != 0 || candhlt165.size() != 0){
    for(size_t itrig = 0 ; itrig < candhlt30.size(); ++itrig ){
    
    if(std::hypot((PhotonT.eta()-candhlt30.at(itrig).Eta()),(PhotonT.phi()-candhlt30.at(itrig).Phi())) < 0.3 && candhlt30.at(itrig).Pt()/PhotonT.pt() > 0.5 && candhlt30.at(itrig).Pt()/PhotonT.pt() < 1.5) { isMatch30_ -> push_back(true);} else{isMatch30_ -> push_back(false); } 
      
     } 
     
     for(size_t itrig = 0 ; itrig < candhlt50.size(); ++itrig ){
    
    if(std::hypot((PhotonT.eta()-candhlt50.at(itrig).Eta()),(PhotonT.phi()-candhlt50.at(itrig).Phi())) < 0.3  && candhlt50.at(itrig).Pt()/PhotonT.pt() > 0.5 && candhlt50.at(itrig).Pt()/PhotonT.pt() < 1.5){ isMatch50_ -> push_back(true);} else{isMatch50_ -> push_back(false); } 
      
     }
     
     for(size_t itrig = 0 ; itrig < candhlt75.size(); ++itrig ){
    
    if(std::hypot((PhotonT.eta()-candhlt75.at(itrig).Eta()),(PhotonT.phi()-candhlt75.at(itrig).Phi())) < 0.3 && candhlt75.at(itrig).Pt()/PhotonT.pt() > 0.5 && candhlt75.at(itrig).Pt()/PhotonT.pt() < 1.5) { isMatch75_ -> push_back(true);} else{isMatch75_ -> push_back(false); } 
      
     }
     
     for(size_t itrig = 0 ; itrig < candhlt90.size(); ++itrig ){
    
    if(std::hypot((PhotonT.eta()-candhlt90.at(itrig).Eta()),(PhotonT.phi()-candhlt90.at(itrig).Phi())) < 0.3 && candhlt90.at(itrig).Pt()/PhotonT.pt() > 0.5 && candhlt90.at(itrig).Pt()/PhotonT.pt() < 1.5) { isMatch90_ -> push_back(true);} else{isMatch90_ -> push_back(false); } 
      
     }
     
     for(size_t itrig = 0 ; itrig < candhlt120.size(); ++itrig ){
    
    if(std::hypot((PhotonT.eta()-candhlt120.at(itrig).Eta()),(PhotonT.phi()-candhlt120.at(itrig).Phi())) < 0.3  && candhlt120.at(itrig).Pt()/PhotonT.pt() > 0.5 && candhlt120.at(itrig).Pt()/PhotonT.pt() < 1.5){ isMatch120_ -> push_back(true);} else{isMatch120_ -> push_back(false); } 
      
     }
     
     for(size_t itrig = 0 ; itrig < candhlt165.size(); ++itrig ){
    
    if(std::hypot((PhotonT.eta()-candhlt165.at(itrig).Eta()),(PhotonT.phi()-candhlt165.at(itrig).Phi())) < 0.3  && candhlt165.at(itrig).Pt()/PhotonT.pt() > 0.5 && candhlt165.at(itrig).Pt()/PhotonT.pt() < 1.5){ isMatch165_ -> push_back(true);} else{isMatch165_ -> push_back(false); } 
      
     }
     

    }else{
    
       isMatch30_ -> push_back(false);
       isMatch50_ -> push_back(false);
       isMatch75_ -> push_back(false);
       isMatch90_ -> push_back(false);
       isMatch120_ -> push_back(false);
       isMatch165_ -> push_back(false);
    
    
    }
    
   }   
    
  
  
  // AK4
  std::vector<double> jecFactorsAK4;
  std::vector<unsigned> sortedAK4JetIdx;

  uint32_t indexjet = 0;
  if(redoJECs_)
    {
      // sort AK4 jets by increasing pT
      std::multimap<double, unsigned> sortedAK4Jets;
     
    

      for(pat::JetCollection::const_iterator ijet = jetsAK4->begin();ijet != jetsAK4->end(); ++ijet)
	{
	  double correction = 1.;
	  JetCorrectorAK4_DATA->setJetEta(ijet->eta());
	  JetCorrectorAK4_DATA->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorAK4_DATA->setJetA(ijet->jetArea());
	  JetCorrectorAK4_DATA->setRho(rho_);
	  JetCorrectorAK4_MC->setJetEta(ijet->eta());
	  JetCorrectorAK4_MC->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorAK4_MC->setJetA(ijet->jetArea());
	  JetCorrectorAK4_MC->setRho(rho_);
	  if (iEvent.isRealData()) 
	    correction = JetCorrectorAK4_DATA->getCorrection();
	  else
	    correction = JetCorrectorAK4_MC->getCorrection();
	    
           
	  jecFactorsAK4.push_back(correction);
	  sortedAK4Jets.insert(std::make_pair(ijet->correctedJet(0).pt()*correction, ijet - jetsAK4->begin()));
	
	}
      // get jet indices in decreasing pT order

      for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedAK4Jets.rbegin(); it != sortedAK4Jets.rend(); ++it)
        sortedAK4JetIdx.push_back(it->second);

    }
  else
    {
      for(pat::JetCollection::const_iterator ijet = jetsAK4->begin();ijet != jetsAK4->end(); ++ijet)
	{
	  jecFactorsAK4.push_back(1./ijet->jecFactor(0));
	}
    }
    
  pat::Jet firstJet ;
  nJetsAK4_ = 0;
  float htAK4(0.0);
  int jetotresp =  0 ;
  vector<TLorentzVector> vP4AK4;
  for(std::vector<unsigned>::const_iterator i = sortedAK4JetIdx.begin(); i != sortedAK4JetIdx.end(); ++i, indexjet++) {

    pat::JetCollection::const_iterator ijet = (jetsAK4->begin() + *i);
    pat::JetCollection::const_iterator ijetrawRCstore = (jetrawRC.begin()+ *i);
    
    jetotresp ++;
    if (jetotresp == 1 ){
    firstJet = (*ijet);}
    
    edm::View<pat::Jet>::const_iterator ijetview = (jetsview->begin() + *i);
    double chf = ijet->chargedHadronEnergyFraction();
    
    double nhf = ijet->neutralHadronEnergyFraction(); 
    double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
    double muf = ijet->muonEnergyFraction();

    double hf_hf = ijet->HFHadronEnergyFraction();
    double hf_emf= ijet->HFEMEnergyFraction();
    double hof   = ijet->hoEnergyFraction();

    int chm    = ijet->chargedHadronMultiplicity();
      
    int chMult = ijet->chargedMultiplicity();
    int neMult = ijet->neutralMultiplicity();
    int npr    = chMult + neMult;

    int chHadMult = chm; 
    int neHadMult = ijet->neutralHadronMultiplicity();
    int phoMult = ijet->photonMultiplicity();
      
    double nemf = ijet->neutralEmEnergyFraction();
    double cemf = ijet->chargedEmEnergyFraction();
    int NumConst = npr;
    
    float eta  = ijet->eta(); 
    float pt   = ijet->correctedJet("Uncorrected").pt()*jecFactorsAK4.at(*i);   
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
    int idL = -999 ; 
    int idT = -999 ; 

    if(fabs(eta) > 3.0)
      {
        idL = ( nemf<0.90 && neMult>10 ) ;
        idT = ( nemf<0.90 && nhf>0.2 && neMult>10 ) ;
      }else{
      if(fabs(eta) > 2.7)
        {
          idL = ( nemf>0.01 && nhf<0.98 && neMult > 2 );
          idT = ( nemf>0.02 && nemf<0.99 && neMult>2);
        }else{
        if(fabs(eta)>2.6)
          {
            idL = ( nemf>0.01 && nhf<0.98 && neMult > 2 );
            idT = ( cemf<0.8 && chm>0 && nemf<0.99 && muf <0.8 && nhf < 0.9 ) ;
          }else{
          idL = ( nemf>0.01 && nhf<0.98 && neMult > 2 ) ;
          idT = ( cemf<0.8 && chm>0 && chf>0 && NumConst>1 && nemf<0.9 && muf <0.8 && nhf < 0.9 ) ;
        }
      }
    }
  
    edm::RefToBase<pat::Jet> jetReftmp(edm::Ref<edm::View<pat::Jet> >(jetsview, ijetview - jetsview->begin()));

       
      
    if (pt > ptMinAK4_) {
      htAK4 += pt;
      nJetsAK4_++;
      vP4AK4.push_back(TLorentzVector(ijet->correctedJet(0).px()*jecFactorsAK4.at(*i),ijet->correctedJet(0).py()*jecFactorsAK4.at(*i),ijet->correctedJet(0).pz()*jecFactorsAK4.at(*i),ijet->correctedJet(0).energy()*jecFactorsAK4.at(*i)));
      chfAK4_           ->push_back(chf);
      nhfAK4_           ->push_back(nhf);
      phfAK4_           ->push_back(phf);
      elfAK4_           ->push_back(elf);
      mufAK4_           ->push_back(muf);
      nemfAK4_          ->push_back(nemf);
      cemfAK4_          ->push_back(cemf);
      hf_hfAK4_         ->push_back(hf_hf);
      hf_emfAK4_        ->push_back(hf_emf);
      hofAK4_           ->push_back(hof);
      jecAK4_           ->push_back(jecFactorsAK4.at(*i));
      ptAK4_            ->push_back(pt);
      phiAK4_           ->push_back(ijet->phi());
      etaAK4_           ->push_back(ijet->eta());
      massAK4_          ->push_back(ijet->correctedJet(0).mass()*jecFactorsAK4.at(*i));
      energyAK4_        ->push_back(ijet->correctedJet(0).energy()*jecFactorsAK4.at(*i));
      
      
      if (!iEvent.isRealData() && ijet->genJet()) {
        ptGenAK4_            ->push_back(ijet->genJet()->pt());
 	phiGenAK4_           ->push_back(ijet->genJet()->phi());
 	etaGenAK4_           ->push_back(ijet->genJet()->eta());
 	massGenAK4_          ->push_back(ijet->genJet()->mass());
 	energyGenAK4_        ->push_back(ijet->genJet()->energy());
	if(ijet->genParton())
	{
 	  pdgIDGenAK4_         ->push_back(ijet->genParton()->pdgId());
	}else{pdgIDGenAK4_         ->push_back(-999.);}
      }else
      {
        
        ptGenAK4_            ->push_back(-999.);
 	phiGenAK4_           ->push_back(-999.);
 	etaGenAK4_           ->push_back(-999.);
 	massGenAK4_          ->push_back(-999.);
 	energyGenAK4_        ->push_back(-999.);
	pdgIDGenAK4_         ->push_back(-999.);
 	
      }
      
      
      ptAK4raw_            ->push_back(ijetrawRCstore->pt());
      phiAK4raw_           ->push_back(ijetrawRCstore->phi());
      etaAK4raw_           ->push_back(ijetrawRCstore->eta());
      massAK4raw_          ->push_back(ijetrawRCstore->mass());
      energyAK4raw_        ->push_back(ijetrawRCstore->energy());
      
      areaAK4_          ->push_back(ijet->jetArea());
      csvAK4_           ->push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      qgdAK4_           ->push_back(0./*(*qgHandle)[jetReftmp]*/);
      

      idLAK4_           ->push_back(idL);
      idTAK4_           ->push_back(idT);
      chHadMultAK4_     ->push_back(chHadMult);
      chMultAK4_        ->push_back(chMult);
      neHadMultAK4_     ->push_back(neHadMult);  
      neMultAK4_        ->push_back(neMult);
      phoMultAK4_       ->push_back(phoMult); 
      hadronflavour_    ->push_back(ijet->hadronFlavour());
    }

  }// jet loop  
  htAK4_     = htAK4;


// PUPPI


  std::vector<double> jecFactorsPUPPI;
  std::vector<unsigned> sortedPUPPIJetIdx;


  uint32_t indexjetpuppi = 0;
  if(redoJECs_)
    {
      // sort PUPPI jets by increasing pT
      std::multimap<double, unsigned> sortedPUPPIJets;
     
    

      for(pat::JetCollection::const_iterator ijet = jetsPUPPI->begin();ijet != jetsPUPPI->end(); ++ijet)
	{
	  double correction = 1.;
	  JetCorrectorPUPPI_DATA->setJetEta(ijet->eta());
	  JetCorrectorPUPPI_DATA->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorPUPPI_DATA->setJetA(ijet->jetArea());
	  JetCorrectorPUPPI_DATA->setRho(rho_);
	  JetCorrectorPUPPI_MC->setJetEta(ijet->eta());
	  JetCorrectorPUPPI_MC->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorPUPPI_MC->setJetA(ijet->jetArea());
	  JetCorrectorPUPPI_MC->setRho(rho_);
	  if (iEvent.isRealData()) 
	    correction = JetCorrectorPUPPI_DATA->getCorrection();
	  else
	    correction = JetCorrectorPUPPI_MC->getCorrection();
	    
           
	  jecFactorsPUPPI.push_back(correction);
	  sortedPUPPIJets.insert(std::make_pair(ijet->correctedJet(0).pt()*correction, ijet - jetsPUPPI->begin()));
	
	}
      // get jet indices in decreasing pT order

      for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedPUPPIJets.rbegin(); it != sortedPUPPIJets.rend(); ++it)
        sortedPUPPIJetIdx.push_back(it->second);

    }
  else
    {
      for(pat::JetCollection::const_iterator ijet = jetsPUPPI->begin();ijet != jetsPUPPI->end(); ++ijet)
	{
	  jecFactorsPUPPI.push_back(1./ijet->jecFactor(0));
	}
    }
    
  nJetsPUPPI_ = 0;
  vector<TLorentzVector> vP4PUPPI;
  for(std::vector<unsigned>::const_iterator i = sortedPUPPIJetIdx.begin(); i != sortedPUPPIJetIdx.end(); ++i, indexjetpuppi++) {

    pat::JetCollection::const_iterator ijet = (jetsPUPPI->begin() + *i);
    pat::JetCollection::const_iterator ijetrawRCstore = (jetrawRCpuppi.begin()+ *i);
    
    edm::View<pat::Jet>::const_iterator ijetview = (jetsviewpuppi->begin() + *i);
    double chf = ijet->chargedHadronEnergyFraction();
    double nhf = ijet->neutralHadronEnergyFraction(); 
    double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());

    double muf = ijet->muonEnergyFraction();

    double hf_hf = ijet->HFHadronEnergyFraction();
    double hf_emf= ijet->HFEMEnergyFraction();
    double hof   = ijet->hoEnergyFraction();

    int chm    = ijet->chargedHadronMultiplicity();
      
    int chMult = ijet->chargedMultiplicity();
    int neMult = ijet->neutralMultiplicity();
    int npr    = chMult + neMult;

    int chHadMult = chm; 
    int neHadMult = ijet->neutralHadronMultiplicity();
    int phoMult = ijet->photonMultiplicity();
      
    double nemf = ijet->neutralEmEnergyFraction();
    double cemf = ijet->chargedEmEnergyFraction();
    int NumConst = npr;

    float eta  = ijet->eta(); 
    float pt   = ijet->correctedJet(0).pt()*jecFactorsPUPPI.at(*i); 
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
    int idL = -999 ; 
    int idT = -999 ; 
    if(fabs(eta) < 3.0)
    {
      idL = ( nemf>0.01 && nhf<0.98 && neMult > 2);
      if(fabs(eta) <= 2.7){
      
         idL =(nhf<0.99 && nemf<0.99 && NumConst>1);
         if(fabs(eta) <= 2.4){
         
            idL = ( chf>0. && chMult>0 && cemf<0.99 && nhf<0.99 && nemf<0.99 && NumConst>1 );
         
         }
      
      }
    
    }else{
       idL = ( nemf<0.90 && neMult>10)   ;
       
       
       idT = (nhf<0.90 && nemf<0.90 && NumConst>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.90) || fabs(eta)>2.4)      ;
    }
   
    edm::RefToBase<pat::Jet> jetReftmp(edm::Ref<edm::View<pat::Jet> >(jetsviewpuppi, ijetview - jetsviewpuppi->begin()));

       
      
    if (pt > ptMinAK4_) {
      nJetsPUPPI_++;
      vP4PUPPI.push_back(TLorentzVector(ijet->correctedJet(0).px()*jecFactorsPUPPI.at(*i),ijet->correctedJet(0).py()*jecFactorsPUPPI.at(*i),ijet->correctedJet(0).pz()*jecFactorsPUPPI.at(*i),ijet->correctedJet(0).energy()*jecFactorsPUPPI.at(*i)));
      chfPUPPI_           ->push_back(chf);
      nhfPUPPI_           ->push_back(nhf);
      phfPUPPI_           ->push_back(phf);
      elfPUPPI_           ->push_back(elf);
      mufPUPPI_           ->push_back(muf);
      nemfPUPPI_          ->push_back(nemf);
      cemfPUPPI_          ->push_back(cemf);
      hf_hfPUPPI_         ->push_back(hf_hf);
      hf_emfPUPPI_        ->push_back(hf_emf);
      hofPUPPI_           ->push_back(hof);
      jecPUPPI_           ->push_back(jecFactorsPUPPI.at(*i));
      ptPUPPI_            ->push_back(pt);
      phiPUPPI_           ->push_back(ijet->phi());
      etaPUPPI_           ->push_back(ijet->eta());
      massPUPPI_          ->push_back(ijet->correctedJet(0).mass()*jecFactorsPUPPI.at(*i));
      energyPUPPI_        ->push_back(ijet->correctedJet(0).energy()*jecFactorsPUPPI.at(*i));

      
      if (!iEvent.isRealData() && ijet->genJet()) {
        ptGenPUPPI_            ->push_back(ijet->genJet()->pt());
 	phiGenPUPPI_           ->push_back(ijet->genJet()->phi());
 	etaGenPUPPI_           ->push_back(ijet->genJet()->eta());
 	massGenPUPPI_          ->push_back(ijet->genJet()->mass());
 	energyGenPUPPI_        ->push_back(ijet->genJet()->energy());
	if(ijet->genParton())
	{
 	  pdgIDGenPUPPI_         ->push_back(ijet->genParton()->pdgId());
	}else{pdgIDGenPUPPI_         ->push_back(-999.);}
      }else
      {
        
        ptGenPUPPI_            ->push_back(-999.);
 	phiGenPUPPI_           ->push_back(-999.);
 	etaGenPUPPI_           ->push_back(-999.);
 	massGenPUPPI_          ->push_back(-999.);
 	energyGenPUPPI_        ->push_back(-999.);
	pdgIDGenPUPPI_         ->push_back(-999.);
 	
      }

      
      ptPUPPIraw_            ->push_back(ijetrawRCstore->pt());
      phiPUPPIraw_           ->push_back(ijetrawRCstore->phi());
      etaPUPPIraw_           ->push_back(ijetrawRCstore->eta());
      massPUPPIraw_          ->push_back(ijetrawRCstore->mass());
      energyPUPPIraw_        ->push_back(ijetrawRCstore->energy());
      
      areaPUPPI_          ->push_back(ijet->jetArea());      
      idLPUPPI_           ->push_back(idL);
      idTPUPPI_           ->push_back(idT);
      chHadMultPUPPI_     ->push_back(chHadMult);
      chMultPUPPI_        ->push_back(chMult);
      neHadMultPUPPI_     ->push_back(neHadMult);  
      neMultPUPPI_        ->push_back(neMult);
      phoMultPUPPI_       ->push_back(phoMult); 
      
    }

  }	
      
  //---- Fill Tree --- 
  outTree_->Fill();     
  //------------------
  
  
}//end analyze for each event

//--------------------------Photon ID and selection--------------------------

// Barrel photon here, end cap further down
// see https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
bool DijetTreeProducer::isValidPhotonLoose(const pat::PhotonRef& photonRef, const edm::Event& event, double generatorWeight)
{
    bool isValid = true;
    //edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  //  event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
   // if (!isData_ && !photonRef->genPhoton())   
   // return false;
   // #1: H/E
    isValid &= photonRef->hadronicOverEm() < 0.04596;
    if (! isValid)
    return false;
    //#2: sigma ietaieta
    isValid &= photonRef->full5x5_sigmaIetaIeta() < 0.0106; //Official    
    if (! isValid)
    return false;
    
    edm::Handle<double> rhos;
    event.getByToken( srcRho_, rhos);
    double rho = *rhos;
    /*edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
    event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
    event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
    */
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoChargedIsolation"), rho, photonRef->eta(), IsolationType::CHARGED_HADRONS)< 1.694;
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoNeutralHadronIsolation"), rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (24.032 + 0.01512*photonRef->pt()+0.00002259*(photonRef->pt()*photonRef->pt() ) );
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoPhotonIsolation"), rho, photonRef->eta(), IsolationType::PHOTONS) < (2.876 + 0.004017*photonRef->pt());
    isValid &= photonRef->passElectronVeto();
    if (! isValid)
    return false;
    isValid &= photonRef->r9() >0.90;
    if (! isValid)
    return false;
    return isValid;
    
    
}

bool DijetTreeProducer::isValidPhotonMedium(const pat::PhotonRef& photonRef, const edm::Event& event, double generatorWeight)
{
    bool isValid = true;
   // edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
   // event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
   // if (!isData_ && !photonRef->genPhoton())   
   // return false;
   // #1: H/E
    isValid &= photonRef->hadronicOverEm() < 0.02197;
    if (! isValid)
    return false;
    //#2: sigma ietaieta
    isValid &= photonRef->full5x5_sigmaIetaIeta()< 0.01015; //Official    
    if (! isValid)
    return false;
    
    edm::Handle<double> rhos;
    event.getByToken( srcRho_, rhos);
    double rho = *rhos;
    /*edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
    event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
    event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);*/
    
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoChargedIsolation"), rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 1.141;
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoNeutralHadronIsolation"), rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (1.189 + 0.01512*photonRef->pt()+0.00002259*(photonRef->pt()*photonRef->pt() ) );
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoPhotonIsolation"), rho, photonRef->eta(), IsolationType::PHOTONS) < (2.08 + 0.004017*photonRef->pt());
    isValid &= photonRef->passElectronVeto();
    if (! isValid)
    return false;
    isValid &= photonRef->r9() >0.90;
    if (! isValid)
    return false;
    return isValid;
    
    
}
bool DijetTreeProducer::isValidPhotonTight(const pat::PhotonRef& photonRef, const edm::Event& event, double generatorWeight)
{
    bool isValid = true;
   // edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
    //event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
    //if (!isData_ && !photonRef->genPhoton())   
   // return false;
   // #1: H/E
    isValid &= photonRef->hadronicOverEm() < 0.02148;
    if (! isValid)
    return false;
    //#2: sigma ietaieta
    isValid &= photonRef->full5x5_sigmaIetaIeta() < 0.00996; //Official    
    if (! isValid)
    return false;
    
    edm::Handle<double> rhos;
    event.getByToken( srcRho_, rhos);
    double rho = *rhos;
  /*  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
    event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
    event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);*/
    
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoChargedIsolation"), rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 0.65;
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoNeutralHadronIsolation"), rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (0.317 + 0.01512*photonRef->pt()+0.00002259*(photonRef->pt()*photonRef->pt() ) );
    isValid &= getCorrectedPFIsolation(photonRef->userFloat("phoPhotonIsolation"), rho, photonRef->eta(), IsolationType::PHOTONS) < (2.044+0.004017*photonRef->pt());
    isValid &= photonRef->passElectronVeto();
    if (! isValid)
    return false;
    isValid &= photonRef->r9() >0.90;
    if (! isValid)
    return false;
    return isValid;
    
    
}


bool DijetTreeProducer::isValidEndcapPhotonLoose(const pat::PhotonRef& photonRef, const edm::Event& event, double generatorWeight)
{
    bool isValid = true;
    //edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
   // event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
   // if (!isData_ && !photonRef->genPhoton())   
   // return false;
   // #1: H/E
    isValid &= photonRef->hadronicOverEm() < 0.0590;
    if (! isValid)
    return false;
    //#2: sigma ietaieta
    isValid &= photonRef->full5x5_sigmaIetaIeta()< 0.0272; //Official    
    if (! isValid)
    return false;
    
    edm::Handle<double> rhos;
    event.getByToken( srcRho_, rhos);
    double rho = *rhos;
    edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
    event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
    event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
    
    isValid &= getCorrectedPFIsolation((*phoChargedIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS)< 2.089;
    isValid &= getCorrectedPFIsolation((*phoNeutralHadronIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (19.722 + 0.0117*photonRef->pt()+0.000023*(photonRef->pt()*photonRef->pt() ) );
    isValid &= getCorrectedPFIsolation((*phoPhotonIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) < (4.162 + 0.0037*photonRef->pt());
    isValid &= photonRef->passElectronVeto();
    if (! isValid)
    return false;
    isValid &= photonRef->r9() >0.90;
    if (! isValid)
    return false;
    return isValid;
    
    
}

bool DijetTreeProducer::isValidEndcapPhotonMedium(const pat::PhotonRef& photonRef, const edm::Event& event, double generatorWeight)
{
    bool isValid = true;
  //  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
   // event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
   // if (!isData_ && !photonRef->genPhoton())   
   // return false;
   // #1: H/E
    isValid &= photonRef->hadronicOverEm() < 0.0326;
    if (! isValid)
    return false;
    //#2: sigma ietaieta
    isValid &= photonRef->full5x5_sigmaIetaIeta() < 0.0272; //Official    
    if (! isValid)
    return false;
    
    edm::Handle<double> rhos;
    event.getByToken( srcRho_, rhos);
    double rho = *rhos;
    edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
    event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
    event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
    
    isValid &= getCorrectedPFIsolation((*phoChargedIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 1.051;
    isValid &= getCorrectedPFIsolation((*phoNeutralHadronIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (2.718 + 0.0117*photonRef->pt()+0.000023*(photonRef->pt()*photonRef->pt() ) );
    isValid &= getCorrectedPFIsolation((*phoPhotonIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) < (3.867 + 0.0037*photonRef->pt());
    isValid &= photonRef->passElectronVeto();
    if (! isValid)
    return false;
    isValid &= photonRef->r9() >0.90;
    if (! isValid)
    return false;
    return isValid;
    
    
}
bool DijetTreeProducer::isValidEndcapPhotonTight(const pat::PhotonRef& photonRef, const edm::Event& event, double generatorWeight)
{
    bool isValid = true;
    //edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
    //event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
    //if (!isData_ && !photonRef->genPhoton())   
   // return false;
   // #1: H/E
    isValid &= photonRef->hadronicOverEm() < 0.0321;
    if (! isValid)
    return false;
    //#2: sigma ietaieta
    isValid &= photonRef->full5x5_sigmaIetaIeta() < 0.0271; //Official    
    if (! isValid)
    return false;
    
    edm::Handle<double> rhos;
    event.getByToken( srcRho_, rhos);
    double rho = *rhos;
    edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
    event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
    event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
    
    isValid &= getCorrectedPFIsolation((*phoChargedIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 0.517;
    isValid &= getCorrectedPFIsolation((*phoNeutralHadronIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (2.716 + 0.0117*photonRef->pt()+0.000023*(photonRef->pt()*photonRef->pt() ) );
    isValid &= getCorrectedPFIsolation((*phoPhotonIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) < (3.032+0.0037*photonRef->pt());
    isValid &= photonRef->passElectronVeto();
    if (! isValid)
    return false;
    isValid &= photonRef->r9() >0.90;
    if (! isValid)
    return false;
    return isValid;
    
    
}

//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::initialize()
{
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  BXnumber_       = -999;
  rho_            = -999;
  metEnergy_      = -999;
  metEta_         = -999;
  metPt_          = -999;
  metPhi_         = -999;
  goodPVtx_       = true;
  NgoodPV_        = -999;
  
  metEnergypuppi_      = -999;
  metEtapuppi_         = -999;
  metPtpuppi_          = -999;
  metPhipuppi_         = -999;
  
  metEnergyGen_      = -999;
  metEtaGen_         = -999;
  metPtGen_          = -999;
  metPhiGen_         = -999;
  
  deltaNHfootprintX_ = -999;
  deltaNHfootprintY_ = -999;
  /*
  
  PFmetX_ = -999  ;
  PFmetY_ = -999  ;
  CHSmetX_= -999  ;
  CHSmetY_= -999  ;
  EGmetX_ = -999  ; 
  EGmetY_ = -999  ;
  */
  metEnergypuppiGen_      = -999;
  metEtapuppiGen_         = -999;
  metPtpuppiGen_          = -999;
  metPhipuppiGen_         = -999;
  
  metSig_         = -999;
  metcorrected_   = -999;
  nJetsAK4_          = -999;
  htAK4_             = -999;
  ptAK4_             ->clear();
  etaAK4_            ->clear();
  phiAK4_            ->clear();
  massAK4_           ->clear();
  energyAK4_         ->clear();
  
  ptAK4raw_             ->clear();
  etaAK4raw_            ->clear();
  phiAK4raw_            ->clear();
  massAK4raw_           ->clear();
  energyAK4raw_         ->clear();
  
  areaAK4_           ->clear();
  csvAK4_            ->clear();
  qgdAK4_            ->clear();
  chfAK4_            ->clear();
  nhfAK4_            ->clear();
  phfAK4_            ->clear();
  elfAK4_            ->clear();
  mufAK4_            ->clear();
  nemfAK4_           ->clear();
  cemfAK4_           ->clear();
  hf_hfAK4_             ->clear();
  hf_emfAK4_            ->clear();
  hofAK4_            ->clear();
  jecAK4_            ->clear();
  jecAK4_            ->clear();
  idLAK4_            ->clear();
  idTAK4_            ->clear();
  // Juska's fix
  chHadMultAK4_     ->clear();
  chMultAK4_        ->clear();
  neHadMultAK4_     ->clear();
  neMultAK4_        ->clear();
  phoMultAK4_        ->clear();
  hadronflavour_     ->clear();
 
  nJetsPUPPI_ = -999;
 
  ptPUPPI_             ->clear();
  etaPUPPI_            ->clear();
  phiPUPPI_            ->clear();
  massPUPPI_           ->clear();
  energyPUPPI_         ->clear();
  
  ptPUPPIraw_             ->clear();
  etaPUPPIraw_            ->clear();
  phiPUPPIraw_            ->clear();
  massPUPPIraw_           ->clear();
  energyPUPPIraw_         ->clear();
  
  areaPUPPI_           ->clear();
  csvPUPPI_            ->clear();
  qgdPUPPI_            ->clear();
  chfPUPPI_            ->clear();
  nhfPUPPI_            ->clear();
  phfPUPPI_            ->clear();
  elfPUPPI_            ->clear();
  mufPUPPI_            ->clear();
  nemfPUPPI_           ->clear();
  cemfPUPPI_           ->clear();
  hf_hfPUPPI_             ->clear();
  hf_emfPUPPI_            ->clear();
  hofPUPPI_            ->clear();
  jecPUPPI_            ->clear();
  jecPUPPI_            ->clear();
  idLPUPPI_            ->clear();
  idTPUPPI_            ->clear();
  // Juska's fix
  chHadMultPUPPI_     ->clear();
  chMultPUPPI_        ->clear();
  neHadMultPUPPI_     ->clear();
  neMultPUPPI_        ->clear();
  phoMultPUPPI_        ->clear();
  
  triggerResult_     ->clear();
  triggerPrescale_   ->clear();
  triggerName_       ->clear();

  //----- MC -------
  npu_ ->clear();
  Number_interactions ->clear();
  OriginBX            -> clear();
  
  ptHat_     = -999; 
  processID_ = -999; 
  weight_    = -999;

  nGenJetsAK4_ = -999;
  
  ptGenAK4_    ->clear();
  phiGenAK4_   ->clear();
  etaGenAK4_   ->clear();
  massGenAK4_  ->clear();
  energyGenAK4_->clear();
  pdgIDGenAK4_ ->clear();
  ptGenPUPPI_    ->clear();
  phiGenPUPPI_   ->clear();
  etaGenPUPPI_   ->clear();
  massGenPUPPI_  ->clear();
  energyGenPUPPI_->clear();
  pdgIDGenPUPPI_ ->clear();
  
  gen_eta		->clear();
  gen_phi		->clear();
  gen_p		        ->clear();
  gen_px		->clear();
  gen_py		->clear();
  gen_pz		->clear();
  gen_pt		->clear();
  gen_energy    	->clear();
  gen_pdgId	        ->clear();
  gen_vx		->clear();
  gen_vy		->clear();
  gen_vz		->clear();
  gen_numDaught	        ->clear();
  gen_status	        ->clear();
  gen_index   	        ->clear();
  gen_motherIndex       ->clear();  
  
  
  nPhotons_    =-999;
  nGenphotons_ =-999;
  nPhotonsLoose_  =-999;
  nPhotonsMedium_ =-999;
  nPhotonsTight_  =-999;
  nMuonsLoose_     =-999; 
  
  ptphoton_              ->clear();
  ptphotonSC_            ->clear();
  ptsmearedphoton_       ->clear();

  etaphoton_             ->clear();
  etaphotonSC_           ->clear();
  etasmearedphoton_      ->clear();

  phiphoton_             ->clear();
  phiphotonSC_           ->clear();
  phismearedphoton_      ->clear();

  energyphoton_          ->clear();
  energyphotonSC_        ->clear();
  energysmearedphoton_   ->clear();

  full5x5SigmaIEtaIEtaMapTokenphoton_      ->clear();
  phoChargedIsolationTokenphoton_          ->clear();
  phoNeutralHadronIsolationTokenphoton_    ->clear();
  phoPhotonIsolationTokenphoton_           ->clear();
  isPhotonLoose_                           ->clear();
  isPhotonMedium_                          ->clear();
  isPhotonTight_                           ->clear();
  HaspixelSeed_                            ->clear();
  hadTowOverEm_                            ->clear();
  electronconvVeto_                        ->clear();
  
  ptGenphoton_            ->clear();
  etaGenphoton_           ->clear();
  phiGenphoton_           ->clear();
  energyGenphoton_        ->clear();
  
  
  
   elecPt_          ->clear();       
  elecEta_          ->clear();     
  elecPhi_          ->clear();     
  elecEnergy_       ->clear();
  elecID_           ->clear();   
  elecISO_          ->clear();
  elecPtsmeared_           ->clear();       
  elecEtasmeared_          ->clear();     
  elecPhismeared_          ->clear();    
  elecEnergysmeared_       ->clear();
  elecIDsmeared_           ->clear();
  elecISOsmeared_          ->clear();
  
  muPt_        ->clear();    
  muEta_       ->clear();    
  muPhi_       ->clear();    
  muEnergy_    ->clear();
  
  isMatch30_    ->clear() ;
  isMatch50_    ->clear() ;    
  isMatch75_    ->clear() ;    
  isMatch90_    ->clear() ;    
  isMatch120_   ->clear() ;    
  isMatch165_   ->clear() ;
  isGenMatch_   ->clear() ;
      
  Ecorrbump_   ->clear();
  
  
  photon_scaleUNC_gainup_->clear();
  photon_scaleUNC_gaindown_->clear();
  photon_scaleUNC_systup_->clear();
  photon_scaleUNC_systdown_->clear(); 
  photon_scaleUNC_statup_->clear(); 
  photon_scaleUNC_statdown_->clear(); 
  photon_scaleUNC_ETup_->clear(); 
  photon_scaleUNC_ETdown_->clear(); 
  photon_smearUNC_phiup_->clear(); 
  photon_smearUNC_rhodown_->clear();
  photon_smearUNC_rhoup_ ->clear(); 
  photon_smear_central_ ->clear();
  photon_scale_central_->clear();
  
}
//////////////////////////////////////////////////////////////////////////////////////////
DijetTreeProducer::~DijetTreeProducer() 
{
}

DEFINE_FWK_MODULE(DijetTreeProducer);

