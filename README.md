# Gamma + jet JERC analysis framework

This framework has 3 modules that have to be run in this order :
1. [DijetRootTreeMaker](https://github.com/lucastorterotot/DijetRootTreeMaker) which produce the Ntuples.
2. [DijetRootTreeAnalyzer](https://github.com/lucastorterotot/DijetRootTreeAnalyzer) which reduce the Ntuple by applying the selection and allow to run JEC on the fly.
3. [NewGammaJet](https://github.com/lucastorterotot/NewGammaJet) which is the plotting step.

Instruction can be found in each step repository.

This module is the first step of the gamma + jet analysis.
- To run on 2016 data, stay on the `master` branch.
- To run on 2017 data switch to `9_2_X` branch.
- To run on 2017 data switch to `JEC_JER_2018_CMSSW_10_2_5_master` branch.

Pick up your branch and then check instruction [here](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md)
and for more information on the original framework
see instruction [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/ExoDijet13TeV#Instructions_to_Create_RootTuple).
