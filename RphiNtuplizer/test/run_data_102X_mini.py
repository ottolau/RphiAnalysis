import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018ABC_v2', '')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/50000/00105BAD-63E0-E711-8640-02163E0146C5.root'
        '/store/data/Run2018A/ParkingBPH3/MINIAOD/14May2018-v1/70000/002DFC31-605A-E811-AA95-0CC47ABAC0CE.root'
        )
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
#runOnData( process, outputModules = [] )
removeMCMatching(process, names=['All'], outputModules=[])

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

process.load("RphiAnalysis.RphiNtuplizer.RphiNtuplizer_miniAOD_cfi")
process.RphiNtuplizer.development=cms.bool(True)
process.RphiNtuplizer.isAOD=cms.bool(False)
process.RphiNtuplizer.doGenParticles=cms.bool(False)
process.RphiNtuplizer.dumpMuons=cms.bool(True)
process.RphiNtuplizer.dumpElectrons=cms.bool(True)
process.RphiNtuplizer.dumpLowPtElectrons=cms.bool(True)
process.RphiNtuplizer.dumpMixElectrons=cms.bool(True)
process.RphiNtuplizer.removeDupEle=cms.bool(True)

process.p = cms.Path(
    process.ggNtuplizer
)

