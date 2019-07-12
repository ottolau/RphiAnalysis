import FWCore.ParameterSet.Config as cms

RphiNtuplizer = cms.EDAnalyzer("RphiNtuplizer",
                             doGenParticles       = cms.bool(False),
                             isAOD                = cms.bool(True), #### actually configured through run_data_74x.py
                             development          = cms.bool(False),
                             dumpElectrons        = cms.bool(True),
                             dumpMuons            = cms.bool(True),
                             dumpLowPtElectrons   = cms.bool(True),
                             dumpMixElectrons     = cms.bool(True),
                             removeDupEle         = cms.bool(True),

                             trgFilterDeltaPtCut  = cms.double(0.4),
                             trgFilterDeltaRCut   = cms.double(0.03),

                             triggerEvent         = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                             #triggerEvent         = cms.InputTag("selectedPatTrigger", "", ""),
                             #triggerEvent         = cms.InputTag("selectedPatTrigger", "", "PAT"),
                             triggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
                             patTriggerResults    = cms.InputTag("TriggerResults", "", "PAT"),
                             #patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),
                             genParticleSrc       = cms.InputTag("genParticles"),
                             generatorLabel       = cms.InputTag("generator"),
                             LHEEventLabel        = cms.InputTag("externalLHEProducer"),
                             pileupCollection     = cms.InputTag("addPileupInfo"),
                             VtxLabel             = cms.InputTag("offlinePrimaryVertices"),
                             VtxBSLabel           = cms.InputTag("offlinePrimaryVerticesWithBS"),
                             rhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
                             rhoCentralLabel      = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                             pfMETLabel           = cms.InputTag("patMETs"),
                             electronSrc          = cms.InputTag("selectedPatElectrons"),
                             #calibelectronSrc     = cms.InputTag("calibratedPatElectrons"),
                             calibelectronSrc     = cms.InputTag("selectedPatElectrons"),
                             photonSrc            = cms.InputTag("selectedPatPhotons"),
                             #calibphotonSrc       = cms.InputTag("calibratedPatPhotons"),
                             calibphotonSrc       = cms.InputTag("selectedPatPhotons"),
                             muonSrc              = cms.InputTag("selectedPatMuons"),
                             gsfTrackSrc          = cms.InputTag("electronGsfTracks"),
                             ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                             eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                             esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
                             recoPhotonSrc             = cms.InputTag("reducedEgamma", "reducedGedPhotonCores"),
                             TrackLabel                = cms.InputTag("generalTracks"),
                             gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                             PFAllCandidates           = cms.InputTag("particleFlow"),
                             #ak4JetSrc                 = cms.InputTag("updatedJets"),
                             #ak4JetSrc                 = cms.InputTag("slimmedJets"),
                             ak4JetSrc                 = cms.InputTag("patJets"),

                             #ak4JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
                             ak8JetSrc                 = cms.InputTag("slimmedJetsAK8"),
                             #ak8JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJECAK8"),
                             #boostedDoubleSVLabel      = cms.InputTag("pfBoostedDoubleSecondaryVertexAK8BJetTags"),
                             tauSrc                    = cms.InputTag("selectedPatTaus"),
                             #pfLooseId                 = pfJetIDSelector.clone(),

                             packedPFCands   = cms.InputTag("particleFlow::RECO"),
                             lostTracks      = cms.InputTag("lostTracks"),
                             lowpTelectrons  = cms.InputTag("lowPtGsfElectrons"),
                             eleBiasedWP = cms.InputTag("lowPtGsfElectronSeedValueMaps","ptbiased","RECO"),
                             eleUnbiasedWP = cms.InputTag("lowPtGsfElectronSeedValueMaps","unbiased","RECO"),
                             conversions = cms.InputTag("allConversions"),
                             PackedGenParticleSrc       = cms.InputTag("packedGenParticles"),
                             muPtCut = cms.double(2.0),
                             muEtaCut = cms.double(2.1),
                             muDzCut = cms.double(0.5),
                             elePtCut = cms.double(0.5),
                             eleEtaCut = cms.double(2.5),
                             eleDzCut = cms.double(0.5),
                             mixPfPtCut = cms.double(0.5),
                             mixPfEtaCut = cms.double(2.5),
                             mixPfDzCut = cms.double(0.5),           
                             mixLowPtPtLowCut = cms.double(0.5),     
                             mixLowPtPtUpCut = cms.double(5.0),
                             mixLowPtEtaCut = cms.double(2.5),
                             mixLowPtDzCut = cms.double(0.5),
                             mixLowPtUnBWPCut = cms.double(0.19),
                             lowPtPtLowCut = cms.double(0.5),
                             lowPtPtUpCut = cms.double(5.0),
                             lowPtEtaCut = cms.double(2.5),
                             lowPtDzCut = cms.double(0.5),
                             lowPtUnBWPLeadCut = cms.double(6.0),
                             lowPtUnBWPSubleadCut = cms.double(0.19),
                             kaonPtCut = cms.double(0.4),
                             kaonEtaCut = cms.double(2.5),
                             kaonDzCut = cms.double(0.5),
                             dilepMLowCut = cms.double(2.6),
                             dilepMUpCut = cms.double(3.6),
                             phiMLowCut = cms.double(0.98),
                             phiMUpCut = cms.double(1.06),
                             bsMLowCut = cms.double(4.5),
                             bsMUpCut = cms.double(6.0),
                             svProbCut = cms.double(0.001),
                             cosAngleCut = cms.double(0.0),
                             
)
