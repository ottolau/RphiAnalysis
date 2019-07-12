#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RphiAnalysis/RphiNtuplizer/interface/RphiNtuplizer.h"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
Int_t          nEle_;
vector<int>    eleCharge_lead_;
vector<int>    eleChargeConsistent_lead_;
vector<float>  eleEn_lead_;
vector<float>  eleSCEn_lead_;
vector<float>  eleEcalEn_lead_;
vector<float>  eleD0_lead_;
vector<float>  eleDz_lead_;
vector<float>  eleD0Error_lead_;
vector<float>  eleDzError_lead_;
vector<float>  eleSIP_lead_;
vector<float>  elePt_lead_;
vector<float>  eleEta_lead_;
vector<float>  elePhi_lead_;
vector<float>  eleGsfPt_lead_;
vector<float>  eleGsfEta_lead_;
vector<float>  eleGsfPhi_lead_;
vector<float>  eleCalibPt_lead_;
vector<float>  eleCalibEn_lead_;
vector<int>    eleConvVeto_lead_;

vector<int>    eleCharge_sublead_;
vector<int>    eleChargeConsistent_sublead_;
vector<float>  eleEn_sublead_;
vector<float>  eleSCEn_sublead_;
vector<float>  eleEcalEn_sublead_;
vector<float>  eleD0_sublead_;
vector<float>  eleDz_sublead_;
vector<float>  eleD0Error_sublead_;
vector<float>  eleDzError_sublead_;
vector<float>  eleSIP_sublead_;
vector<float>  elePt_sublead_;
vector<float>  eleEta_sublead_;
vector<float>  elePhi_sublead_;
vector<float>  eleGsfPt_sublead_;
vector<float>  eleGsfEta_sublead_;
vector<float>  eleGsfPhi_sublead_;
vector<float>  eleCalibPt_sublead_;
vector<float>  eleCalibEn_sublead_;
vector<int>    eleConvVeto_sublead_;

vector<float> eleSvChi2_;
vector<float> eleSvNDOF_;
vector<float> eleSvProb_;
vector<float> eleSvX_;
vector<float> eleSvY_;
vector<float> eleSvZ_;
vector<float> eleSvXError_;
vector<float> eleSvYError_;
vector<float> eleSvZError_;
vector<float> eleSvMass_;
vector<float> eleSvCtxy_;
vector<float> eleSvCosAngle_;
vector<float> eleSvLxy_;
vector<float> eleSvLxyError_;

vector<int>    kaonEECharge_lead_;
vector<float>  kaonEED0_lead_;
vector<float>  kaonEEDz_lead_;
vector<float>  kaonEED0Error_lead_;
vector<float>  kaonEEDzError_lead_;
vector<float>  kaonEEPt_lead_;
vector<float>  kaonEEEta_lead_;
vector<float>  kaonEEPhi_lead_;
vector<float>  kaonEEVx_lead_;
vector<float>  kaonEEVy_lead_;
vector<float>  kaonEEVz_lead_;
vector<float>  kaonEETrkChi2_lead_;
vector<float>  kaonEETrkNDOF_lead_;
vector<float>  kaonEETrkNormChi2_lead_;
vector<float>  kaonEEJPsiMass_lead_;
vector<float>  kaonEEPhiMass_lead_;
//vector<float>  kaonEETrkdEdx_lead_;

vector<int>    kaonEECharge_sublead_;
vector<float>  kaonEED0_sublead_;
vector<float>  kaonEEDz_sublead_;
vector<float>  kaonEED0Error_sublead_;
vector<float>  kaonEEDzError_sublead_;
vector<float>  kaonEEPt_sublead_;
vector<float>  kaonEEEta_sublead_;
vector<float>  kaonEEPhi_sublead_;
vector<float>  kaonEEVx_sublead_;
vector<float>  kaonEEVy_sublead_;
vector<float>  kaonEEVz_sublead_;
vector<float>  kaonEETrkChi2_sublead_;
vector<float>  kaonEETrkNDOF_sublead_;
vector<float>  kaonEETrkNormChi2_sublead_;
//vector<float>  kaonEETrkdEdx_sublead_;

vector<float>  bsEEdRele_;
vector<float>  bsEEdRkaon_;
vector<float>  bsEEdRJpsiPhi_;
vector<float>  bsEEJpsiMass_;
vector<float>  bsEEPhiMass_;
vector<float>  bsEEBsMass_;

void RphiNtuplizer::branchesElectrons(TTree* tree) {

  tree->Branch("nEle",                    &nEle_);
  tree->Branch("eleCharge_lead",               &eleCharge_lead_);
  tree->Branch("eleChargeConsistent_lead",     &eleChargeConsistent_lead_);
  tree->Branch("eleEn_lead",                   &eleEn_lead_);
  tree->Branch("eleSCEn_lead",                 &eleSCEn_lead_);
  tree->Branch("eleEcalEn_lead",               &eleEcalEn_lead_);
  tree->Branch("eleD0_lead",                   &eleD0_lead_);
  tree->Branch("eleDz_lead",                   &eleDz_lead_);
  tree->Branch("eleD0Error_lead",              &eleD0Error_lead_);
  tree->Branch("eleDzError_lead",              &eleDzError_lead_);
  tree->Branch("eleSIP_lead",                  &eleSIP_lead_);
  tree->Branch("elePt_lead",                   &elePt_lead_);
  tree->Branch("eleEta_lead",                  &eleEta_lead_);
  tree->Branch("elePhi_lead",                  &elePhi_lead_);
  tree->Branch("eleGsfPt_lead",                &eleGsfPt_lead_);
  tree->Branch("eleGsfEta_lead",               &eleGsfEta_lead_);
  tree->Branch("eleGsfPhi_lead",               &eleGsfPhi_lead_);
  tree->Branch("eleCalibPt_lead",              &eleCalibPt_lead_);
  tree->Branch("eleCalibEn_lead",              &eleCalibEn_lead_);
  tree->Branch("eleConvVeto_lead",             &eleConvVeto_lead_);

  tree->Branch("eleCharge_sublead",               &eleCharge_sublead_);
  tree->Branch("eleChargeConsistent_sublead",     &eleChargeConsistent_sublead_);
  tree->Branch("eleEn_sublead",                   &eleEn_sublead_);
  tree->Branch("eleSCEn_sublead",                 &eleSCEn_sublead_);
  tree->Branch("eleEcalEn_sublead",               &eleEcalEn_sublead_);
  tree->Branch("eleD0_sublead",                   &eleD0_sublead_);
  tree->Branch("eleDz_sublead",                   &eleDz_sublead_);
  tree->Branch("eleD0Error_sublead",              &eleD0Error_sublead_);
  tree->Branch("eleDzError_sublead",              &eleDzError_sublead_);
  tree->Branch("eleSIP_sublead",                  &eleSIP_sublead_);
  tree->Branch("elePt_sublead",                   &elePt_sublead_);
  tree->Branch("eleEta_sublead",                  &eleEta_sublead_);
  tree->Branch("elePhi_sublead",                  &elePhi_sublead_);
  tree->Branch("eleGsfPt_sublead",                &eleGsfPt_sublead_);
  tree->Branch("eleGsfEta_sublead",               &eleGsfEta_sublead_);
  tree->Branch("eleGsfPhi_sublead",               &eleGsfPhi_sublead_);
  tree->Branch("eleCalibPt_sublead",              &eleCalibPt_sublead_);
  tree->Branch("eleCalibEn_sublead",              &eleCalibEn_sublead_);
  tree->Branch("eleConvVeto_sublead",             &eleConvVeto_sublead_);

  tree->Branch("eleSvChi2",                    &eleSvChi2_);
  tree->Branch("eleSvNDOF",                    &eleSvNDOF_);
  tree->Branch("eleSvProb",                    &eleSvProb_);
  tree->Branch("eleSvX",                       &eleSvX_);
  tree->Branch("eleSvY",                       &eleSvY_);
  tree->Branch("eleSvZ",                       &eleSvZ_);
  tree->Branch("eleSvXError",                  &eleSvXError_);
  tree->Branch("eleSvYError",                  &eleSvYError_);
  tree->Branch("eleSvZError",                  &eleSvZError_);
  tree->Branch("eleSvMass",                    &eleSvMass_);
  tree->Branch("eleSvCtxy",                    &eleSvCtxy_);
  tree->Branch("eleSvCosAngle",                    &eleSvCosAngle_);
  tree->Branch("eleSvLxy",                    	   &eleSvLxy_);
  tree->Branch("eleSvLxyError",                    &eleSvLxyError_);

  tree->Branch("kaonEECharge_lead",               &kaonEECharge_lead_);
  tree->Branch("kaonEED0_lead",                   &kaonEED0_lead_);
  tree->Branch("kaonEEDz_lead",                   &kaonEEDz_lead_);
  tree->Branch("kaonEED0Error_lead",              &kaonEED0Error_lead_);
  tree->Branch("kaonEEDzError_lead",              &kaonEEDzError_lead_);
  tree->Branch("kaonEEPt_lead",                   &kaonEEPt_lead_);
  tree->Branch("kaonEEEta_lead",                  &kaonEEEta_lead_);
  tree->Branch("kaonEEPhi_lead",                  &kaonEEPhi_lead_);
  tree->Branch("kaonEEVx_lead",                   &kaonEEVx_lead_);
  tree->Branch("kaonEEVy_lead",                   &kaonEEVy_lead_);
  tree->Branch("kaonEEVz_lead",                   &kaonEEVz_lead_);
  tree->Branch("kaonEETrkChi2_lead",              &kaonEETrkChi2_lead_);
  tree->Branch("kaonEETrkNDOF_lead",              &kaonEETrkNDOF_lead_);
  tree->Branch("kaonEETrkNormChi2_lead",          &kaonEETrkNormChi2_lead_);
//  tree->Branch("kaonEETrkdEdx_lead",              &kaonEETrkdEdx_lead_);

  tree->Branch("kaonEECharge_sublead",               &kaonEECharge_sublead_);
  tree->Branch("kaonEED0_sublead",                   &kaonEED0_sublead_);
  tree->Branch("kaonEEDz_sublead",                   &kaonEEDz_sublead_);
  tree->Branch("kaonEED0Error_sublead",              &kaonEED0Error_sublead_);
  tree->Branch("kaonEEDzError_sublead",              &kaonEEDzError_sublead_);
  tree->Branch("kaonEEPt_sublead",                   &kaonEEPt_sublead_);
  tree->Branch("kaonEEEta_sublead",                  &kaonEEEta_sublead_);
  tree->Branch("kaonEEPhi_sublead",                  &kaonEEPhi_sublead_);
  tree->Branch("kaonEEVx_sublead",                   &kaonEEVx_sublead_);
  tree->Branch("kaonEEVy_sublead",                   &kaonEEVy_sublead_);
  tree->Branch("kaonEEVz_sublead",                   &kaonEEVz_sublead_);
  tree->Branch("kaonEETrkChi2_sublead",              &kaonEETrkChi2_sublead_);
  tree->Branch("kaonEETrkNDOF_sublead",              &kaonEETrkNDOF_sublead_);
  tree->Branch("kaonEETrkNormChi2_sublead",          &kaonEETrkNormChi2_sublead_);
//  tree->Branch("kaonEETrkdEdx_sublead",              &kaonEETrkdEdx_sublead_);

  tree->Branch("bsEEdRele",                     &bsEEdRele_);
  tree->Branch("bsEEdRkaon",                    &bsEEdRkaon_);
  tree->Branch("bsEEdRJpsiPhi",                 &bsEEdRJpsiPhi_);
  tree->Branch("bsEEJpsiMass",                  &bsEEJpsiMass_);
  tree->Branch("bsEEPhiMass",                   &bsEEPhiMass_);
  tree->Branch("bsEEBsMass",                    &bsEEBsMass_);
  
}

void RphiNtuplizer::fillElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  // cleanup from previous execution
  eleCharge_lead_                  .clear();
  eleChargeConsistent_lead_        .clear();
  eleEn_lead_                      .clear();
  eleSCEn_lead_                    .clear();
  eleEcalEn_lead_                  .clear();
  eleD0_lead_                      .clear();
  eleDz_lead_                      .clear();
  eleD0Error_lead_                 .clear();
  eleDzError_lead_                 .clear();
  eleSIP_lead_                     .clear();
  elePt_lead_                      .clear();
  eleEta_lead_                     .clear();
  elePhi_lead_                     .clear();
  eleGsfPt_lead_                   .clear();
  eleGsfEta_lead_                  .clear();
  eleGsfPhi_lead_                  .clear();
  eleCalibPt_lead_                 .clear();
  eleCalibEn_lead_                 .clear();
  eleConvVeto_lead_                .clear();

  eleCharge_sublead_                  .clear();
  eleChargeConsistent_sublead_        .clear();
  eleEn_sublead_                      .clear();
  eleSCEn_sublead_                    .clear();
  eleEcalEn_sublead_                  .clear();
  eleD0_sublead_                      .clear();
  eleDz_sublead_                      .clear();
  eleD0Error_sublead_                 .clear();
  eleDzError_sublead_                 .clear();
  eleSIP_sublead_                     .clear();
  elePt_sublead_                      .clear();
  eleEta_sublead_                     .clear();
  elePhi_sublead_                     .clear();
  eleGsfPt_sublead_                   .clear();
  eleGsfEta_sublead_                  .clear();
  eleGsfPhi_sublead_                  .clear();
  eleCalibPt_sublead_                 .clear();
  eleCalibEn_sublead_                 .clear();
  eleConvVeto_sublead_                .clear();

  eleSvChi2_.clear();
  eleSvNDOF_.clear();
  eleSvProb_.clear();
  eleSvX_.clear();
  eleSvY_.clear();
  eleSvZ_.clear();
  eleSvXError_.clear();
  eleSvYError_.clear();
  eleSvZError_.clear();
  eleSvMass_.clear();
  eleSvCtxy_.clear();
  eleSvCosAngle_.clear();
  eleSvLxy_.clear();
  eleSvLxyError_.clear();

  kaonEECharge_lead_                  .clear();
  kaonEED0_lead_                      .clear();
  kaonEEDz_lead_                      .clear();
  kaonEED0Error_lead_                 .clear();
  kaonEEDzError_lead_                 .clear();
  kaonEEPt_lead_                      .clear();
  kaonEEEta_lead_                     .clear();
  kaonEEPhi_lead_                     .clear();
  kaonEEVx_lead_                      .clear();
  kaonEEVy_lead_                      .clear();
  kaonEEVz_lead_                      .clear();
  kaonEETrkChi2_lead_                 .clear();
  kaonEETrkNDOF_lead_                 .clear();
  kaonEETrkNormChi2_lead_             .clear();
//  kaonEETrkdEdx_lead_	              .clear();

  kaonEECharge_sublead_                  .clear();
  kaonEED0_sublead_                      .clear();
  kaonEEDz_sublead_                      .clear();
  kaonEED0Error_sublead_                 .clear();
  kaonEEDzError_sublead_                 .clear();
  kaonEEPt_sublead_                      .clear();
  kaonEEEta_sublead_                     .clear();
  kaonEEPhi_sublead_                     .clear();
  kaonEEVx_sublead_                      .clear();
  kaonEEVy_sublead_                      .clear();
  kaonEEVz_sublead_                      .clear();
  kaonEETrkChi2_sublead_                 .clear();
  kaonEETrkNDOF_sublead_                 .clear();
  kaonEETrkNormChi2_sublead_             .clear();
//  kaonEETrkdEdx_sublead_	         .clear();

  bsEEdRele_		      .clear();
  bsEEdRkaon_                 .clear();
  bsEEdRJpsiPhi_              .clear();
  bsEEJpsiMass_		      .clear();
  bsEEPhiMass_		      .clear();
  bsEEBsMass_		      .clear();

  nEle_ = 0;

  if (isAOD_) {

    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<edm::View<pat::Electron> > calibelectronHandle;
    e.getByToken(calibelectronCollection_, calibelectronHandle);

    edm::Handle<reco::TrackCollection> tracksHandle;
    e.getByToken(tracklabel_, tracksHandle);

  //  edm::Handle<reco::DeDxDataValueMap> dEdxObjectHandle;
  //  e.getByToken(deDxProducer_, dEdxObjectHandle );
  //  const edm::ValueMap<reco::DeDxData> dEdxColl = *dEdxObjectHandle.product();


    if (!electronHandle.isValid()) {
      edm::LogWarning("RphiNtuplizer") << "no pat::Electrons in event";
      return;
    }

    if (!calibelectronHandle.isValid()) {
      edm::LogWarning("RphiNtuplizer") << "no calibrated pat::Electrons in event";
      return;
    }


    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      if (iEle->pt() < elePtCut_) continue;
      if (fabs(iEle->eta()) > eleEtaCut_) continue;
      if (fabs(iEle->vz() - pv.z()) > eleDzCut_) continue;

      for (edm::View<pat::Electron>::const_iterator jEle = iEle+1; jEle != electronHandle->end(); ++jEle) {
	if (jEle->pt() < elePtCut_) continue;
        if (fabs(iEle->eta()) > eleEtaCut_) continue;
	if (fabs(jEle->vz() - pv.z()) > eleDzCut_) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
	if ((iele_lv + jele_lv).M() < dilepMLowCut_ || (iele_lv + jele_lv).M() > dilepMUpCut_) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	//std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;
	//reco::Track ieletrk = eletrks[(iEle-electronHandle->begin())].second;
	//reco::Track jeletrk = eletrks[(jEle-electronHandle->begin())].second;
	//const reco::TransientTrack ielettk = getTransientTrack( ieletrk );
	//const reco::TransientTrack jelettk = getTransientTrack( jeletrk );

	//XParticles.push_back(pFactory.particle(getTransientTrack( ieletrk ), pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(getTransientTrack( jeletrk ), pmass, 0.0, 0, pmasse));

	//XParticles.push_back(pFactory.particle(ielettk, pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(jelettk, pmass, 0.0, 0, pmasse));

	//KinematicConstrainedVertexFitter kvFitter;
	//RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 || KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;
	//KinVtx->movePointerToTheTop();
	//RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();

	//for (pat::PackedCandidateCollection::const_iterator iHad = pfcands->begin(); iHad != pfcands->end(); ++iHad) {
	for (reco::TrackCollection::const_iterator iHad = tracksHandle->begin(); iHad != tracksHandle->end(); ++iHad) {
	  if (iHad->pt() < kaonPtCut_) continue;
	  if (fabs(iHad->eta()) > kaonEtaCut_) continue;
          if (fabs(iHad->vz() - pv.z()) > kaonDzCut_) continue;
	  if (iHad->normalizedChi2() < 0.0) continue;
	  if (iHad->normalizedChi2() > 20.0) continue;

	  //for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != pfcands->end(); ++jHad) {
	  for (reco::TrackCollection::const_iterator jHad = iHad+1; jHad != tracksHandle->end(); ++jHad) {
	    //if (iHad->charge()*jHad->charge() > 0.0) continue;
	    if (jHad->pt() <  kaonPtCut_) continue;
	    if (fabs(jHad->eta()) > kaonEtaCut_) continue;
	    if (fabs(jHad->vz() - pv.z()) > kaonDzCut_) continue;
	    if (jHad->normalizedChi2() < 0.0) continue;
	    if (jHad->normalizedChi2() > 20.0) continue;

	    // Phi mass window
	    float kpmass = 0.493677;
	    TLorentzVector iHad_lv, jHad_lv, bs_lv;
	    iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), kpmass);
	    jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), kpmass);      
	    bs_lv = iele_lv + jele_lv + iHad_lv + jHad_lv;
	    if ((iHad_lv+jHad_lv).M() < phiMLowCut_ || (iHad_lv+jHad_lv).M() > phiMUpCut_) continue; 
	    if (bs_lv.M() < bsMLowCut_ || bs_lv.M() > bsMUpCut_) continue;

	    std::vector<RefCountedKinematicParticle> BsParticles;
	    float kpmasse = 1.e-6 * pmass;
	    //float bsM = 5.3663;

	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	    //BsParticles.push_back(jpsi_part);

	    KinematicConstrainedVertexFitter BsKvFitter;
	    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
	    if (!(BsKinVtx->isValid())) continue;

	    RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

	    if (DecayVtx->chiSquared() < 0.0) continue;
	    if (TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()) < svProbCut_) continue;

	    // Accept these 4 tracks as a Bs candidate, fill ntuple

	    auto leadEle = iEle->pt() > jEle->pt() ? iEle : jEle;
	    auto subleadEle = iEle->pt() > jEle->pt() ? jEle : iEle;
	    auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
	    auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

	    double ctxy = ((DecayVtx->position().x() - pv.x())*bs_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_lv.Py())/(pow(bs_lv.Pt(),2))*bs_lv.M();
	    
	    math::XYZVector perp(bs_lv.Px(), bs_lv.Py(), 0.);
	    math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	    math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	    double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

	    if (cosAngle < cosAngleCut_) continue;

	    eleSvChi2_.push_back(DecayVtx->chiSquared());
	    eleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    eleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    eleSvX_.push_back(DecayVtx->position().x());
	    eleSvY_.push_back(DecayVtx->position().y());
	    eleSvZ_.push_back(DecayVtx->position().z());
	    eleSvXError_.push_back(DecayVtx->error().cxx());
	    eleSvYError_.push_back(DecayVtx->error().cyy());
	    eleSvZError_.push_back(DecayVtx->error().czz());
	    eleSvMass_.push_back(bs_lv.M());
	    eleSvCtxy_.push_back(ctxy);
	    eleSvCosAngle_.push_back(cosAngle);
	    eleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    eleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonEECharge_lead_            .push_back(leadHad->charge());
	    kaonEED0_lead_                .push_back(leadHad->dxy(pv));
	    kaonEEDz_lead_                .push_back(leadHad->dz(pv));
	    kaonEED0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonEEDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonEEPt_lead_                .push_back(leadHad->pt());
	    kaonEEEta_lead_               .push_back(leadHad->eta());
	    kaonEEPhi_lead_               .push_back(leadHad->phi());
	    kaonEEVx_lead_ 		.push_back(leadHad->vx());
	    kaonEEVy_lead_ 		.push_back(leadHad->vy());
	    kaonEEVz_lead_ 		.push_back(leadHad->vz());
	    kaonEETrkChi2_lead_ 		.push_back(leadHad->chi2());
	    kaonEETrkNDOF_lead_ 		.push_back(leadHad->ndof());
	    kaonEETrkNormChi2_lead_ 	.push_back(leadHad->normalizedChi2());
	    //kaonEETrkdEdx_lead_		.push_back(dEdxColl[reco::TrackRef(tracksHandle, leadHad-tracksHandle->begin())].dEdx());

	    kaonEECharge_sublead_         .push_back(subleadHad->charge());
	    kaonEED0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonEEDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonEED0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonEEDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonEEPt_sublead_             .push_back(subleadHad->pt());
	    kaonEEEta_sublead_            .push_back(subleadHad->eta());
	    kaonEEPhi_sublead_            .push_back(subleadHad->phi());
	    kaonEEVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonEEVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonEEVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonEETrkChi2_sublead_ 	.push_back(subleadHad->chi2());
	    kaonEETrkNDOF_sublead_ 	.push_back(subleadHad->ndof());
	    kaonEETrkNormChi2_sublead_ 	.push_back(subleadHad->normalizedChi2());
	    //kaonEETrkdEdx_sublead_	.push_back(dEdxColl[reco::TrackRef(tracksHandle, subleadHad-tracksHandle->begin())].dEdx());

	    bsEEdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsEEdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsEEdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsEEJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsEEPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsEEBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    Float_t corrPt_lead = -1;
	    Float_t corrEn_lead = -1;
	    Float_t corrPt_sublead = -1;
	    Float_t corrEn_sublead = -1;

	    for (edm::View<pat::Electron>::const_iterator iCEle = calibelectronHandle->begin(); iCEle != calibelectronHandle->end(); ++iCEle) {

	      if (fabs(leadEle->eta() - iCEle->eta()) < 0.001 && fabs(leadEle->phi() - iCEle->phi()) < 0.001) {
		corrPt_lead = iCEle->pt();
		corrEn_lead = iCEle->energy();
	      }
	      if (fabs(subleadEle->eta() - iCEle->eta()) < 0.001 && fabs(subleadEle->phi() - iCEle->phi()) < 0.001) {
		corrPt_sublead = iCEle->pt();
		corrEn_sublead = iCEle->energy();
	      }
	    }

	    eleCalibPt_lead_        .push_back(corrPt_lead);
	    eleCalibEn_lead_        .push_back(corrEn_lead);
	    eleCalibPt_sublead_     .push_back(corrPt_sublead);
	    eleCalibEn_sublead_     .push_back(corrEn_sublead);

	    eleCharge_lead_          .push_back(leadEle->charge());
	    eleChargeConsistent_lead_.push_back((Int_t)leadEle->isGsfCtfScPixChargeConsistent());
	    eleEn_lead_              .push_back(leadEle->energy());
	    eleD0_lead_              .push_back(leadEle->gsfTrack()->dxy(pv));
	    eleDz_lead_              .push_back(leadEle->gsfTrack()->dz(pv));
	    eleD0Error_lead_         .push_back(leadEle->gsfTrack()->dxyError());
	    eleDzError_lead_         .push_back(leadEle->gsfTrack()->dzError());
	    eleSIP_lead_             .push_back(fabs(leadEle->dB(pat::Electron::PV3D))/leadEle->edB(pat::Electron::PV3D));
	    elePt_lead_              .push_back(leadEle->pt());
	    eleEta_lead_             .push_back(leadEle->eta());
	    elePhi_lead_             .push_back(leadEle->phi());
	    eleGsfPt_lead_           .push_back(leadEle->gsfTrack()->ptMode());
	    eleGsfEta_lead_          .push_back(leadEle->gsfTrack()->etaMode());
	    eleGsfPhi_lead_          .push_back(leadEle->gsfTrack()->phiMode());
	    eleSCEn_lead_            .push_back(leadEle->superCluster()->energy());
	    eleEcalEn_lead_          .push_back(leadEle->ecalEnergy());
	    eleConvVeto_lead_        .push_back((Int_t)leadEle->passConversionVeto()); // ConvVtxFit || missHit == 0

	    eleCharge_sublead_          .push_back(subleadEle->charge());
	    eleChargeConsistent_sublead_.push_back((Int_t)subleadEle->isGsfCtfScPixChargeConsistent());
	    eleEn_sublead_              .push_back(subleadEle->energy());
	    eleD0_sublead_              .push_back(subleadEle->gsfTrack()->dxy(pv));
	    eleDz_sublead_              .push_back(subleadEle->gsfTrack()->dz(pv));
	    eleD0Error_sublead_         .push_back(subleadEle->gsfTrack()->dxyError());
	    eleDzError_sublead_         .push_back(subleadEle->gsfTrack()->dzError());
	    eleSIP_sublead_             .push_back(fabs(subleadEle->dB(pat::Electron::PV3D))/subleadEle->edB(pat::Electron::PV3D));
	    elePt_sublead_              .push_back(subleadEle->pt());
	    eleEta_sublead_             .push_back(subleadEle->eta());
	    elePhi_sublead_             .push_back(subleadEle->phi());
	    eleGsfPt_sublead_           .push_back(subleadEle->gsfTrack()->ptMode());
	    eleGsfEta_sublead_          .push_back(subleadEle->gsfTrack()->etaMode());
	    eleGsfPhi_sublead_          .push_back(subleadEle->gsfTrack()->phiMode());
	    eleSCEn_sublead_            .push_back(subleadEle->superCluster()->energy());
	    eleEcalEn_sublead_          .push_back(subleadEle->ecalEnergy());
	    eleConvVeto_sublead_        .push_back((Int_t)subleadEle->passConversionVeto()); // ConvVtxFit || missHit == 0

	    nEle_++;
	  }
	}
      }
    }
  } else {

    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<edm::View<pat::Electron> > calibelectronHandle;
    e.getByToken(calibelectronCollection_, calibelectronHandle);

    edm::Handle<pat::PackedCandidateCollection> pfcands;
    e.getByToken(pckPFCandidateCollection_, pfcands);

    edm::Handle<pat::PackedCandidateCollection> losttracks;
    e.getByToken(lostTracksLabel_, losttracks);

    std::vector<pat::PackedCandidate> alltracks;
    alltracks.reserve(pfcands->size() + losttracks->size());
    alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
    alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

  //  edm::Handle<reco::DeDxDataValueMap> dEdxObjectHandle;
  //  e.getByToken(deDxProducer_, dEdxObjectHandle );
  //  const edm::ValueMap<reco::DeDxData> dEdxColl = *dEdxObjectHandle.product();


    if (!electronHandle.isValid()) {
      edm::LogWarning("RphiNtuplizer") << "no pat::Electrons in event";
      return;
    }

    if (!calibelectronHandle.isValid()) {
      edm::LogWarning("RphiNtuplizer") << "no calibrated pat::Electrons in event";
      return;
    }


    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    VertexDistanceXY vertTool;

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      if (iEle->pt() < elePtCut_) continue;
      if (fabs(iEle->eta()) > eleEtaCut_) continue;
      if (fabs(iEle->vz() - pv.z()) > eleDzCut_) continue;

      for (edm::View<pat::Electron>::const_iterator jEle = iEle+1; jEle != electronHandle->end(); ++jEle) {
	if (jEle->pt() < elePtCut_) continue;
	if (fabs(jEle->eta()) > eleEtaCut_) continue;
	if (fabs(jEle->vz() - pv.z()) > eleDzCut_) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
	if ((iele_lv + jele_lv).M() < dilepMLowCut_ || (iele_lv + jele_lv).M() > dilepMUpCut_) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	//std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;
	//reco::Track ieletrk = eletrks[(iEle-electronHandle->begin())].second;
	//reco::Track jeletrk = eletrks[(jEle-electronHandle->begin())].second;
	//const reco::TransientTrack ielettk = getTransientTrack( ieletrk );
	//const reco::TransientTrack jelettk = getTransientTrack( jeletrk );

	//XParticles.push_back(pFactory.particle(getTransientTrack( ieletrk ), pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(getTransientTrack( jeletrk ), pmass, 0.0, 0, pmasse));

	//XParticles.push_back(pFactory.particle(ielettk, pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(jelettk, pmass, 0.0, 0, pmasse));

	//KinematicConstrainedVertexFitter kvFitter;
	//RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0 ||  KinVtx->currentDecayVertex()->chiSquared() > 30.0) continue;
	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;
	//KinVtx->movePointerToTheTop();
	//RefCountedKinematicParticle jpsi_part = KinVtx->currentParticle();

	for (pat::PackedCandidateCollection::const_iterator iHad = alltracks.begin(); iHad != alltracks.end(); ++iHad) {
	  if (iHad->pt() <= kaonPtCut_) continue;
          if (iHad->charge() == 0) continue;
          if (abs(iHad->pdgId()) != 211) continue;
          if (iHad->bestTrack() == nullptr) continue;
	  if (fabs(iHad->eta()) > kaonEtaCut_) continue;
	  if (fabs(iHad->vz() - pv.z()) > kaonDzCut_) continue;

	  for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != alltracks.end(); ++jHad) {
	    if (jHad->pt() <= kaonPtCut_) continue;
	    if (fabs(jHad->eta()) > kaonEtaCut_) continue;
            if (jHad->charge() == 0) continue;
            if (abs(jHad->pdgId()) != 211) continue;
            if (jHad->bestTrack() == nullptr) continue;
	    if (fabs(jHad->vz() - pv.z()) > kaonDzCut_) continue;
	    //if (iHad->charge()*jHad->charge() > 0.0) continue;

	    // Phi mass window
	    float kpmass = 0.493677;
	    TLorentzVector iHad_lv, jHad_lv, bs_lv;
	    iHad_lv.SetPtEtaPhiM(iHad->pt(), iHad->eta(), iHad->phi(), kpmass);
	    jHad_lv.SetPtEtaPhiM(jHad->pt(), jHad->eta(), jHad->phi(), kpmass);      
	    bs_lv = iele_lv + jele_lv + iHad_lv + jHad_lv;
	    if ((iHad_lv+jHad_lv).M() < phiMLowCut_ || (iHad_lv+jHad_lv).M() > phiMUpCut_) continue; 
	    if (bs_lv.M() < bsMLowCut_ || bs_lv.M() > bsMUpCut_) continue;

	    std::vector<RefCountedKinematicParticle> BsParticles;
	    float kpmasse = 1.e-6 * pmass;
	    //float bsM = 5.3663;

	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad->bestTrack()) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	    //BsParticles.push_back(jpsi_part);

	    KinematicConstrainedVertexFitter BsKvFitter;
	    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
	    if (!(BsKinVtx->isValid())) continue;

	    RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

	    if (DecayVtx->chiSquared() < 0.0) continue;
	    if (TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()) < svProbCut_) continue;

	    // Accept these 4 tracks as a Bs candidate, fill ntuple

	    auto leadEle = iEle->pt() > jEle->pt() ? iEle : jEle;
	    auto subleadEle = iEle->pt() > jEle->pt() ? jEle : iEle;
	    auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
	    auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

	    double ctxy = ((DecayVtx->position().x() - pv.x())*bs_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_lv.Py())/(pow(bs_lv.Pt(),2))*bs_lv.M();
	    
	    math::XYZVector perp(bs_lv.Px(), bs_lv.Py(), 0.);
	    math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	    math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	    double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

	    if (cosAngle < cosAngleCut_) continue;

	    eleSvChi2_.push_back(DecayVtx->chiSquared());
	    eleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    eleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    eleSvX_.push_back(DecayVtx->position().x());
	    eleSvY_.push_back(DecayVtx->position().y());
	    eleSvZ_.push_back(DecayVtx->position().z());
	    eleSvXError_.push_back(DecayVtx->error().cxx());
	    eleSvYError_.push_back(DecayVtx->error().cyy());
	    eleSvZError_.push_back(DecayVtx->error().czz());
	    eleSvMass_.push_back(bs_lv.M());
	    eleSvCtxy_.push_back(ctxy);
	    eleSvCosAngle_.push_back(cosAngle);
	    eleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    eleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonEECharge_lead_            .push_back(leadHad->charge());
	    kaonEED0_lead_                .push_back(leadHad->dxy(pv));
	    kaonEEDz_lead_                .push_back(leadHad->dz(pv));
	    kaonEED0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonEEDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonEEPt_lead_                .push_back(leadHad->pt());
	    kaonEEEta_lead_               .push_back(leadHad->eta());
	    kaonEEPhi_lead_               .push_back(leadHad->phi());
	    kaonEEVx_lead_ 		.push_back(leadHad->vx());
	    kaonEEVy_lead_ 		.push_back(leadHad->vy());
	    kaonEEVz_lead_ 		.push_back(leadHad->vz());
	    kaonEETrkChi2_lead_ 		.push_back(leadHad->bestTrack()->chi2());
	    kaonEETrkNDOF_lead_ 		.push_back(leadHad->bestTrack()->ndof());
	    kaonEETrkNormChi2_lead_ 	.push_back(leadHad->bestTrack()->normalizedChi2());
	    //kaonEETrkdEdx_lead_		.push_back(dEdxColl[reco::TrackRef(tracksHandle, leadHad-tracksHandle->begin())].dEdx());

	    kaonEECharge_sublead_         .push_back(subleadHad->charge());
	    kaonEED0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonEEDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonEED0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonEEDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonEEPt_sublead_             .push_back(subleadHad->pt());
	    kaonEEEta_sublead_            .push_back(subleadHad->eta());
	    kaonEEPhi_sublead_            .push_back(subleadHad->phi());
	    kaonEEVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonEEVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonEEVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonEETrkChi2_sublead_ 	.push_back(subleadHad->bestTrack()->chi2());
	    kaonEETrkNDOF_sublead_ 	.push_back(subleadHad->bestTrack()->ndof());
	    kaonEETrkNormChi2_sublead_ 	.push_back(subleadHad->bestTrack()->normalizedChi2());
	    //kaonEETrkdEdx_sublead_	.push_back(dEdxColl[reco::TrackRef(tracksHandle, subleadHad-tracksHandle->begin())].dEdx());

	    bsEEdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsEEdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsEEdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsEEJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsEEPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsEEBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    Float_t corrPt_lead = -1;
	    Float_t corrEn_lead = -1;
	    Float_t corrPt_sublead = -1;
	    Float_t corrEn_sublead = -1;

	    for (edm::View<pat::Electron>::const_iterator iCEle = calibelectronHandle->begin(); iCEle != calibelectronHandle->end(); ++iCEle) {

	      if (fabs(leadEle->eta() - iCEle->eta()) < 0.001 && fabs(leadEle->phi() - iCEle->phi()) < 0.001) {
		corrPt_lead = iCEle->pt();
		corrEn_lead = iCEle->energy();
	      }
	      if (fabs(subleadEle->eta() - iCEle->eta()) < 0.001 && fabs(subleadEle->phi() - iCEle->phi()) < 0.001) {
		corrPt_sublead = iCEle->pt();
		corrEn_sublead = iCEle->energy();
	      }
	    }

	    eleCalibPt_lead_        .push_back(corrPt_lead);
	    eleCalibEn_lead_        .push_back(corrEn_lead);
	    eleCalibPt_sublead_     .push_back(corrPt_sublead);
	    eleCalibEn_sublead_     .push_back(corrEn_sublead);

	    eleCharge_lead_          .push_back(leadEle->charge());
	    eleChargeConsistent_lead_.push_back((Int_t)leadEle->isGsfCtfScPixChargeConsistent());
	    eleEn_lead_              .push_back(leadEle->energy());
	    eleD0_lead_              .push_back(leadEle->gsfTrack()->dxy(pv));
	    eleDz_lead_              .push_back(leadEle->gsfTrack()->dz(pv));
	    eleD0Error_lead_         .push_back(leadEle->gsfTrack()->dxyError());
	    eleDzError_lead_         .push_back(leadEle->gsfTrack()->dzError());
	    eleSIP_lead_             .push_back(fabs(leadEle->dB(pat::Electron::PV3D))/leadEle->edB(pat::Electron::PV3D));
	    elePt_lead_              .push_back(leadEle->pt());
	    eleEta_lead_             .push_back(leadEle->eta());
	    elePhi_lead_             .push_back(leadEle->phi());
	    eleGsfPt_lead_           .push_back(leadEle->gsfTrack()->ptMode());
	    eleGsfEta_lead_          .push_back(leadEle->gsfTrack()->etaMode());
	    eleGsfPhi_lead_          .push_back(leadEle->gsfTrack()->phiMode());
	    eleSCEn_lead_            .push_back(leadEle->superCluster()->energy());
	    eleEcalEn_lead_          .push_back(leadEle->ecalEnergy());
	    eleConvVeto_lead_        .push_back((Int_t)leadEle->passConversionVeto()); // ConvVtxFit || missHit == 0

	    eleCharge_sublead_          .push_back(subleadEle->charge());
	    eleChargeConsistent_sublead_.push_back((Int_t)subleadEle->isGsfCtfScPixChargeConsistent());
	    eleEn_sublead_              .push_back(subleadEle->energy());
	    eleD0_sublead_              .push_back(subleadEle->gsfTrack()->dxy(pv));
	    eleDz_sublead_              .push_back(subleadEle->gsfTrack()->dz(pv));
	    eleD0Error_sublead_         .push_back(subleadEle->gsfTrack()->dxyError());
	    eleDzError_sublead_         .push_back(subleadEle->gsfTrack()->dzError());
	    eleSIP_sublead_             .push_back(fabs(subleadEle->dB(pat::Electron::PV3D))/subleadEle->edB(pat::Electron::PV3D));
	    elePt_sublead_              .push_back(subleadEle->pt());
	    eleEta_sublead_             .push_back(subleadEle->eta());
	    elePhi_sublead_             .push_back(subleadEle->phi());
	    eleGsfPt_sublead_           .push_back(subleadEle->gsfTrack()->ptMode());
	    eleGsfEta_sublead_          .push_back(subleadEle->gsfTrack()->etaMode());
	    eleGsfPhi_sublead_          .push_back(subleadEle->gsfTrack()->phiMode());
	    eleSCEn_sublead_            .push_back(subleadEle->superCluster()->energy());
	    eleEcalEn_sublead_          .push_back(subleadEle->ecalEnergy());
	    eleConvVeto_sublead_        .push_back((Int_t)subleadEle->passConversionVeto()); // ConvVtxFit || missHit == 0

	    nEle_++;
	  }
	}
      }
    }

  }

}


