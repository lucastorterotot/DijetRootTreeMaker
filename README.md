# Gamma + jet JERC analysis framework

This framework has 3 modules that have to be run in this order :
1. [DijetRootTreeMaker](https://github.com/lucastorterotot/DijetRootTreeMaker) which produce the Ntuples.
2. [DijetRootTreeAnalyzer](https://github.com/lucastorterotot/DijetRootTreeAnalyzer) which reduce the Ntuple by applying the selection and allow to run JEC on the fly.
3. [NewGammaJet](https://github.com/lucastorterotot/NewGammaJet) which is the plotting step.

Instruction can be found in each step repository.

This module is the first step of the gamma + jet analysis. Instructions can be found [here](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md)
and for more information on the original framework
see instruction [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/ExoDijet13TeV#Instructions_to_Create_RootTuple).
