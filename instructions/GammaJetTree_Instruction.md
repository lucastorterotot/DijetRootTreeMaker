# Instructions for producing Gamma + jet step 1 root file 

**Target** CMS Gamma + Jet  - JEC calculation

**Authors** Hugues Lattaud and Lucas Torterotot

**Last update** 20 Mar 2020

## Installation recipe for era 2017
More informations on [this twiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2)
### CMSSW release
You may want to install this in your work space, for example
```
cd /afs/cern.ch/work/${USER:0:1}/$USER/
mkdir -p JEC && cd JEC
```

```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaID_949 #if you want the FallV2 IDs, otherwise skip
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940 #just adds in an extra file to have a setup function to make things easier
```

### SetUp the Gamma + Jet frame work and compile
```
cd $CMSSW_BASE/src
git clone -o lucas git@github.com:lucastorterotot/DijetRootTreeMaker.git -b JEC_JER_2017_CMSSW_9_4_10_master CMSDIJET/DijetRootTreeMaker/
scram b -j 8
```

You can now go back to [the general instructions](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md) to continue.