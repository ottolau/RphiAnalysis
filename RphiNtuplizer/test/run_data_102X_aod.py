import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask

process = cms.Process('ggKit')

patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018ABC_v2', '')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/data/Run2018A/ParkingBPH1/AOD/05May2019-v1/00000/04879D85-8B10-E747-9EAA-A2324D784F44.root',
        )
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
patAlgosToolsTask.add(process.patCandidatesTask)
#process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
#patAlgosToolsTask.add(process.triggerProducerTask)
process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
patAlgosToolsTask.add(process.selectedPatCandidatesTask)
process.load('PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff')
patAlgosToolsTask.add(process.cleanPatCandidatesTask)

### add trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, None, None, None, None, '' )

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
#runOnData( process, outputModules = [] )
removeMCMatching(process, names=['All'], outputModules=[])

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))


process.load("RphiAnalysis.RphiNtuplizer.RphiNtuplizer_cfi")
process.RphiNtuplizer.development=cms.bool(True)

process.RphiNtuplizer.isAOD=cms.bool(True)
process.RphiNtuplizer.doGenParticles=cms.bool(False)
process.RphiNtuplizer.dumpMuons=cms.bool(True)
process.RphiNtuplizer.dumpElectrons=cms.bool(True)
process.RphiNtuplizer.dumpLowPtElectrons=cms.bool(True)
process.RphiNtuplizer.dumpMixElectrons=cms.bool(True)
process.RphiNtuplizer.removeDupEle=cms.bool(True)

process.p = cms.Path(
    process.RphiNtuplizer,
    patAlgosToolsTask
  )


