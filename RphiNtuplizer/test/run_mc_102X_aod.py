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
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/mc/RunIIAutumn18DR/BsToPhiJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/10000/90738B70-0856-1E4E-8902-86B3C339FF4F.root'
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

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

process.load("RphiAnalysis.RphiNtuplizer.RphiNtuplizer_cfi")
process.RphiNtuplizer.development=cms.bool(True)

process.RphiNtuplizer.isAOD=cms.bool(True)
process.RphiNtuplizer.doGenParticles=cms.bool(True)
process.RphiNtuplizer.dumpMuons=cms.bool(True)
process.RphiNtuplizer.dumpElectrons=cms.bool(True)
process.RphiNtuplizer.dumpLowPtElectrons=cms.bool(False)
process.RphiNtuplizer.dumpMixElectrons=cms.bool(False)
process.RphiNtuplizer.removeDupEle=cms.bool(True)

process.p = cms.Path(
    process.RphiNtuplizer,
    patAlgosToolsTask
  )


