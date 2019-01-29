# Instructions for producing Gamma + jet step 1 root file 

**Target** CMS Gamma + Jet  - JEC calculation

**Authors** Hugues Lattaud and Lucas Torterotot

**Last update** 29 Jan 2019

## Installation recipe
More informations on [this twiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2)
### CMSSW release
You may want to install this in your work space, for example
```
cd /afs/cern.ch/work/${USER:0:1}/$USER/
mkdir -p JEC && cd JEC
```

Check if one has a good compiler version in use `echo $SCRAM_ARCH`.
`slc6_amd64_gcc491` seems to be incompatible with CMSSW 8.0.26.patch1, onemay update update with
```
export SCRAM_ARCH=slc6_amd64_gcc530
```

#### For 2016
If you run on 2016 legacy datasets, all the step needed to set up an environment with calibrated photon are sumarize in [this twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2016DataRecommendations#Running_the_legacy_scale_and_sme).
Following this twikis you should have a CMSSW environment ready for Egamma tuple production.
```
cmsrel CMSSW_8_0_31
cd CMSSW_8_0_31/src
cmsenv
git cms-init
git cms-merge-topic ikrav:egm_id_80X_v3_photons
git cms-merge-topic cms-egamma:EGScaleAndSmearLeg_8026 #a backport of the 94X code for 80X
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data/ -rf  #removes the data directory so can replace with the cms-data repo for this package
git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
scram b -j 8
```

#### For 2017
```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaID_949 #if you want the FallV2 IDs, otherwise skip
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940 #just adds in an extra file to have a setup function to make things easier
scram b -j 8
```

#### For 2018
```
cmsrel CMSSW_10_2_5
cd CMSSW_10_2_5/src/
cmsenv

git cms-init
git cms-merge-topic cms-egamma:EgammaID_1023 #if you want the V2 IDs, otherwise skip
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier

scram b -j 8
```

### SetUp the Gamma + Jet frame work
```
cd $CMSSW_BASE/src
git clone git@github.com:lucastorterotot/DijetRootTreeMaker.git CMSDIJET/DijetRootTreeMaker/
```
- To run on 2016 data, switch to `JEC_JER_2016_master` branch.
- To run on 2017 data, switch to `JEC_JER_2017_CMSSW_9_4_10_master` branch.
- To run on 2017 data, switch to `JEC_JER_2018_CMSSW_10_2_5_master` branch.
```
scram b -j 8
```

### Get the JEC and JER databases
```
cd $CMSSW_BASE/..
git clone git@github.com:cms-jet/JECDatabase.git JECDatabase
git clone git@github.com:cms-jet/JRDatabase.git JRDatabase
```
This allows to run local test if you past the JERC text file you want to use in the `data` directory.

## Test run
You'll find the CMSSW config file `Config_RunonMiniaod.py` in the `prod` repository.
    
You can now launch a local test:
```
cmsRun Config_RunonMiniaod.py globalTag=the_wanted_global_tag  isRunonData=True_or_False isEndcaps=True_or_False JECData=the_wanted_JEC JECMC=the_wanted_JEC maxEvents=numberofeventstoprocess
```
Successful run looks something like this:
```
Begin processing the 8th record. Run 273730, Event 3146673520, LumiSection 2623 at 13-Jun-2016 10:49:45.182 CEST
Begin processing the 9th record. Run 273730, Event 3146111224, LumiSection 2623 at 13-Jun-2016 10:49:45.187 CEST
Begin processing the 10th record. Run 273730, Event 3146571947, LumiSection 2623 at 13-Jun-2016 10:49:45.189 CEST
13-Jun-2016 10:49:45 CEST  Closed file file:/afs/cern.ch/user/j/juska/eos/cms/store/data/Run2016B/JetHT/MINIAOD/PromptReco-v2/000/273/730/00000/EA345ED4-B821-E611-BEA5-02163E0138E2.root
```

You can now look at your output rootfile. A Healthy rootfile look like a file with a single tree called events that store all the variables your are interested in.

For any developpment in the tuples production you can modify the plugin files in the plugin repository : `plugin/DijetTreeProducer.cc` and `plugin/DijetTreeProducer.h`.


If all went correctly you can now begin crab production.

## Production
Make sure to have your grid certificate installed on your lxplus account
```
cd prod/crab/CERN_crab/
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --voms cms --valid 168:00
```

To launch the job processing on crab3 use the `crabSubmit.py` with the required crab config file template and the json file that list the sample you want to process.
* to process data: 
```
python crabSubmit.py -t crab_configDATA.py Sample_data.json 
```
* to process MC: 
```
python crabSubmit.py -t crab_configMC.py Sample_NLO.json  
```

Your are now ready to go to the next step called [DijettreeAnalyser](https://github.com/lucastorterotot/DijetRootTreeAnalyzer).




                              

