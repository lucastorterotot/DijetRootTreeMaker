# Instructions for producing Gamma + jet step 1 root file 

**Target** CMS Gamma + Jet  - JEC calculation

**Authors** Hugues Lattaud and Lucas Torterotot

**Last update** 20 Mar 20120

## Installation recipe for era 2016
More informations on [this twiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2)
### CMSSW release
You may want to install this in your work space, for example
```
cd /afs/cern.ch/work/${USER:0:1}/$USER/
mkdir -p JEC && cd JEC
```

If you run on 2016 legacy datasets, all the step needed to set up an environment with calibrated photon are sumarize in [this twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2016DataRecommendations#Running_the_legacy_scale_and_sme).
Following this twikis you should have a CMSSW environment ready for Egamma tuple production.
```
cmsrel CMSSW_8_0_31
cd CMSSW_8_0_31/src
cmsenv
git cms-init
```

### PhotonID and EGamma regression smearer etc
```
git cms-merge-topic ikrav:egm_id_80X_v3_photons
git cms-merge-topic cms-egamma:EGScaleAndSmearLeg_8026 #a backport of the 94X code for 80X
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data/ -rf  #removes the data directory so can replace with the cms-data repo for this package
git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
```

### SetUp the Gamma + Jet frame work and compile
```
cd $CMSSW_BASE/src
git clone -o lucas git@github.com:lucastorterotot/DijetRootTreeMaker.git -b JEC_JER_2016_master CMSDIJET/DijetRootTreeMaker/
scram b -j 8
```

You can now go back to [the general instructions](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md) to continue.