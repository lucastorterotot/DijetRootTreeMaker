# Instructions for producing Gamma + jet step 1 root file 

**Target** CMS Gamma + Jet  - JEC calculation

**Authors** Hugues Lattaud and Lucas Torterotot

**Last update** 20 Mar 2020

## Installation recipe for era 2017 UL
More informations on [this twiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2)
### CMSSW release
You may want to install this in your work space, for example
```
cd /afs/cern.ch/work/${USER:0:1}/$USER/
mkdir -p JEC && cd JEC
```

```
cmsrel CMSSW_10_6_3
cd CMSSW_10_6_3/src
cmsenv
git cms-init
git cms-merge-topic jainshilpi:UL2017SSV3
git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data -b UL2017SSV2
```

### SetUp the Gamma + Jet frame work and compile
```
cd $CMSSW_BASE/src
git clone -o lucas git@github.com:lucastorterotot/DijetRootTreeMaker.git -b JEC_JER_2017UL_CMSSW_10_6_3_master CMSDIJET/DijetRootTreeMaker/
scram b -j 8
```

You can now go back to [the general instructions](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md) to continue.