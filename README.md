Z Counting
==========
The project has to be set up inside the `src/` folder of a CMSSW environment. The latest version was tested for `CMSSW_12_4_2`. 

The project is divided into 4 packages that are explained in the following.

ZCountAnalyze
----------------
The `ZCountAnalyze` package is used to produce ntuples from cmssw files in MiniAOD or AOD format
and is run over MC or data.
These ntuples contain truth level information and can be used to analyze the method based on simulation, or re-analyze the data in a more flexible way.    
    
ZHarvester
-------------
This is the "core" package of the project. 
It can be used with the histograms from the [cmssw/DQMOffline/Lumi/](https://github.com/cms-sw/cmssw/tree/master/DQMOffline/Lumi) module or the ones produced from the ntuples obtained from the `ZCountAnalyze` package.
* Here, the measurement of the tag-and-probe muon efficiencies and calculation of e.g. the Z Rate can be done.
    * based on RooFit
    * HTCondor job submission is supported
* Several scipts are available to provide additional inputs and produce plots

ZUtils
-------------
Contains code that is used by multiple packages in this project in common.

BLUE
-------------
Combinations of Z counts and reference luminosity with the [BLUE](https://arxiv.org/abs/2001.10310) tool.

Useful Links
-------------
* Theory paper from Jakob Salfeld-Nebgen & Daniel Marlow: https://arxiv.org/abs/1806.02184
* ATLAS public note: https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PUBNOTES/ATL-DAPR-PUB-2021-001/
* Mattermost chat: https://mattermost.web.cern.ch/cms-lumi-pog/channels/zbosoncounting
* Histogram files: https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2022/SingleMuon/
* Run-2 results: https://dwalter.web.cern.ch/dwalter/ZCounting/
* Run-2 twiki: https://twiki.cern.ch/twiki/bin/view/CMS/LUM-21-001
