# RphiAnalysis

This is the ntuplizer for R(phi) measurement (Bs->phi ll). The ntuplizer includes new electron ID (low pT electron) and generated particle matching.

##### To work with CMSSW_10_2_15 and head version, you do :
```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git clone https://github.com/ottolau/RphiAnalysis.git
scram b clean; scram b -j20
```

To run the analyzer, in RphiAnalysis/RphiAnalyzer/test/, you do:
```
cmsRun run_data_102X_aod.py
```
