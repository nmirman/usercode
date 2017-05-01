# Usercode
Package for data skims and ntuples, interfaced with CMS Software (CMSSW).  Runs in batch mode over datasets collected by the CMS experiment.

## JetMETAnalysis
For MET significance variable development.
 - `METSigNtuple` selects events containing top quarks, W/Z bosons, and multiple jets for optimization and performance studies.
 - `METSignificance` is a producer for the MET significance variable, to be inserted into the CMS Software.
 
 ## TopMassAnalysis
 Skims datasets to select events containing top quarks, creates ntuples for analysis.
