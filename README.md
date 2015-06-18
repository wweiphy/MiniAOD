MiniAOD
=======

Tools for miniAOD exploration

## Installation
Follow These Steps:

    cmsrel CMSSW_7_4_4_patch4
    cd CMSSW_7_4_4_patch4/src
    cmsenv

    voms-proxy-init --voms cms

    git clone -b run2mc git@github.com:cms-ttH/MiniAOD.git

    scram b -j 32

    cd MiniAOD/MiniAODHelper/

## Testing
To test using the full CMSSW framework

    cmsRun test/miniAOD_analyzer_cfg.py

To test using FWLite

    root -b -q test/head.C test/miniAOD_investigate.C+'(1000)'


