Gamma + jet JERC analysis framework
==================

This framework has 3 modules that have to be run in this order :

1- DijetRootTreeMaker which produce the Ntuples                                                            : https://github.com/lattaud/DijetRootTreeMaker.git 

2- DijetRootTreeAnalyzer which reduce the Ntuple by applying the selection and allow to run JEC on the fly : https://github.com/lattaud/DijetRootTreeAnalyzer.git

3- NewGammaJet which is the plotting step                                                                  : https://github.com/lattaud/NewGammaJet.git

Instruction can be found in each step repository.

This module is the first step of the gamma + jet analysis : 

to run on 2016 data stay on the master branch
to run on 2017 data switch to 9_2_X branch

pick up your branch and then check instruction here:
https://github.com/lattaud/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.txt


and for more information on the original framework
See instruction here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/ExoDijet13TeV#Instructions_to_Create_RootTuple
