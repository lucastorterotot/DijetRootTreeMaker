# Instructions for producing Gamma + jet step 1 root file 

**Target** CMS Gamma + Jet  - JEC calculation

**Author** Lucas Torterotot

**Last update** 20 Mar 2020

## Installation recipe
More informations on [this twiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2)
### CMSSW release
You may want to install this in your work space, for example
```
cd /afs/cern.ch/work/${USER:0:1}/$USER/
mkdir -p JEC && cd JEC
```

Now, you have to follow different instructions based on the era you want to work on:
- 2016: switch to [`JEC_JER_2016_master` branch](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/JEC_JER_2016_master/instructions/GammaJetTree_Instruction.md);
- 2017: switch to [`JEC_JER_2017_CMSSW_9_4_10_master` branch](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/JEC_JER_2017_CMSSW_9_4_10_master/instructions/GammaJetTree_Instruction.md);
- 2017UL: switch to [`JEC_JER_2017UL_CMSSW_10_6_3_master`
branch](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/JEC_JER_2017UL_CMSSW_10_6_3_master/instructions/GammaJetTree_Instruction.md);
- 2018: switch to [`JEC_JER_2018_CMSSW_10_2_5_master` branch](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/JEC_JER_2018_CMSSW_10_2_5_master/instructions/GammaJetTree_Instruction.md).

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
