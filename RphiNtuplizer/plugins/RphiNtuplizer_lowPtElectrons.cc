#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
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
Int_t          nLowPt_;
vector<int>    lowPtCharge_lead_;
vector<float>  lowPtD0_lead_;
vector<float>  lowPtDz_lead_;
vector<float>  lowPtD0Error_lead_;
vector<float>  lowPtDzError_lead_;
vector<float>  lowPtPt_lead_;
vector<float>  lowPtEta_lead_;
vector<float>  lowPtPhi_lead_;
vector<float>  lowPtPtMean_lead_;
vector<float>  lowPtEtaMean_lead_;
vector<float>  lowPtPhiMean_lead_;
vector<float>  lowPtMVABWP_lead_;
vector<float>  lowPtMVAUnBWP_lead_;

vector<int>    lowPtCharge_sublead_;
vector<float>  lowPtD0_sublead_;
vector<float>  lowPtDz_sublead_;
vector<float>  lowPtD0Error_sublead_;
vector<float>  lowPtDzError_sublead_;
vector<float>  lowPtPt_sublead_;
vector<float>  lowPtEta_sublead_;
vector<float>  lowPtPhi_sublead_;
vector<float>  lowPtPtMean_sublead_;
vector<float>  lowPtEtaMean_sublead_;
vector<float>  lowPtPhiMean_sublead_;
vector<float>  lowPtMVABWP_sublead_;
vector<float>  lowPtMVAUnBWP_sublead_;

vector<float> lowPtSvChi2_;
vector<float> lowPtSvNDOF_;
vector<float> lowPtSvProb_;
vector<float> lowPtSvX_;
vector<float> lowPtSvY_;
vector<float> lowPtSvZ_;
vector<float> lowPtSvXError_;
vector<float> lowPtSvYError_;
vector<float> lowPtSvZError_;
vector<float> lowPtSvMass_;
vector<float> lowPtSvCtxy_;
vector<float> lowPtSvCosAngle_;
vector<float> lowPtSvLxy_;
vector<float> lowPtSvLxyError_;

vector<int>    kaonLowPtCharge_lead_;
vector<float>  kaonLowPtD0_lead_;
vector<float>  kaonLowPtDz_lead_;
vector<float>  kaonLowPtD0Error_lead_;
vector<float>  kaonLowPtDzError_lead_;
vector<float>  kaonLowPtPt_lead_;
vector<float>  kaonLowPtEta_lead_;
vector<float>  kaonLowPtPhi_lead_;
vector<float>  kaonLowPtVx_lead_;
vector<float>  kaonLowPtVy_lead_;
vector<float>  kaonLowPtVz_lead_;
vector<float>  kaonLowPtTrkChi2_lead_;
vector<float>  kaonLowPtTrkNDOF_lead_;
vector<float>  kaonLowPtTrkNormChi2_lead_;
vector<float>  kaonLowPtJPsiMass_lead_;
vector<float>  kaonLowPtPhiMass_lead_;

vector<int>    kaonLowPtCharge_sublead_;
vector<float>  kaonLowPtD0_sublead_;
vector<float>  kaonLowPtDz_sublead_;
vector<float>  kaonLowPtD0Error_sublead_;
vector<float>  kaonLowPtDzError_sublead_;
vector<float>  kaonLowPtPt_sublead_;
vector<float>  kaonLowPtEta_sublead_;
vector<float>  kaonLowPtPhi_sublead_;
vector<float>  kaonLowPtVx_sublead_;
vector<float>  kaonLowPtVy_sublead_;
vector<float>  kaonLowPtVz_sublead_;
vector<float>  kaonLowPtTrkChi2_sublead_;
vector<float>  kaonLowPtTrkNDOF_sublead_;
vector<float>  kaonLowPtTrkNormChi2_sublead_;

vector<float>  bsLowPtdRele_;
vector<float>  bsLowPtdRkaon_;
vector<float>  bsLowPtdRJpsiPhi_;
vector<float>  bsLowPtJpsiMass_;
vector<float>  bsLowPtPhiMass_;
vector<float>  bsLowPtBsMass_;

void RphiNtuplizer::branchesLowPtElectrons(TTree* tree) {

  tree->Branch("nLowPt",                    &nLowPt_);
  tree->Branch("lowPtCharge_lead",               &lowPtCharge_lead_);
  tree->Branch("lowPtD0_lead",                   &lowPtD0_lead_);
  tree->Branch("lowPtDz_lead",                   &lowPtDz_lead_);
  tree->Branch("lowPtD0Error_lead",              &lowPtD0Error_lead_);
  tree->Branch("lowPtDzError_lead",              &lowPtDzError_lead_);
  tree->Branch("lowPtPt_lead",                   &lowPtPt_lead_);
  tree->Branch("lowPtEta_lead",                  &lowPtEta_lead_);
  tree->Branch("lowPtPhi_lead",                  &lowPtPhi_lead_);
  tree->Branch("lowPtPtMean_lead",               &lowPtPtMean_lead_);
  tree->Branch("lowPtEtaMean_lead",              &lowPtEtaMean_lead_);
  tree->Branch("lowPtPhiMean_lead",              &lowPtPhiMean_lead_);
  tree->Branch("lowPtMVABWP_lead",               &lowPtMVABWP_lead_);
  tree->Branch("lowPtMVAUnBWP_lead",             &lowPtMVAUnBWP_lead_);

  tree->Branch("lowPtCharge_sublead",               &lowPtCharge_sublead_);
  tree->Branch("lowPtD0_sublead",                   &lowPtD0_sublead_);
  tree->Branch("lowPtDz_sublead",                   &lowPtDz_sublead_);
  tree->Branch("lowPtD0Error_sublead",              &lowPtD0Error_sublead_);
  tree->Branch("lowPtDzError_sublead",              &lowPtDzError_sublead_);
  tree->Branch("lowPtPt_sublead",                   &lowPtPt_sublead_);
  tree->Branch("lowPtEta_sublead",                  &lowPtEta_sublead_);
  tree->Branch("lowPtPhi_sublead",                  &lowPtPhi_sublead_);
  tree->Branch("lowPtPtMean_sublead",               &lowPtPtMean_sublead_);
  tree->Branch("lowPtEtaMean_sublead",              &lowPtEtaMean_sublead_);
  tree->Branch("lowPtPhiMean_sublead",              &lowPtPhiMean_sublead_);
  tree->Branch("lowPtMVABWP_sublead",               &lowPtMVABWP_sublead_);
  tree->Branch("lowPtMVAUnBWP_sublead",             &lowPtMVAUnBWP_sublead_);

  tree->Branch("lowPtSvChi2",                    &lowPtSvChi2_);
  tree->Branch("lowPtSvNDOF",                    &lowPtSvNDOF_);
  tree->Branch("lowPtSvProb",                    &lowPtSvProb_);
  tree->Branch("lowPtSvX",                       &lowPtSvX_);
  tree->Branch("lowPtSvY",                       &lowPtSvY_);
  tree->Branch("lowPtSvZ",                       &lowPtSvZ_);
  tree->Branch("lowPtSvXError",                  &lowPtSvXError_);
  tree->Branch("lowPtSvYError",                  &lowPtSvYError_);
  tree->Branch("lowPtSvZError",                  &lowPtSvZError_);
  tree->Branch("lowPtSvMass",                    &lowPtSvMass_);
  tree->Branch("lowPtSvCtxy",                    &lowPtSvCtxy_);
  tree->Branch("lowPtSvCosAngle",                    &lowPtSvCosAngle_);
  tree->Branch("lowPtSvLxy",                    	   &lowPtSvLxy_);
  tree->Branch("lowPtSvLxyError",                    &lowPtSvLxyError_);

  tree->Branch("kaonLowPtCharge_lead",               &kaonLowPtCharge_lead_);
  tree->Branch("kaonLowPtD0_lead",                   &kaonLowPtD0_lead_);
  tree->Branch("kaonLowPtDz_lead",                   &kaonLowPtDz_lead_);
  tree->Branch("kaonLowPtD0Error_lead",              &kaonLowPtD0Error_lead_);
  tree->Branch("kaonLowPtDzError_lead",              &kaonLowPtDzError_lead_);
  tree->Branch("kaonLowPtPt_lead",                   &kaonLowPtPt_lead_);
  tree->Branch("kaonLowPtEta_lead",                  &kaonLowPtEta_lead_);
  tree->Branch("kaonLowPtPhi_lead",                  &kaonLowPtPhi_lead_);
  tree->Branch("kaonLowPtVx_lead",                   &kaonLowPtVx_lead_);
  tree->Branch("kaonLowPtVy_lead",                   &kaonLowPtVy_lead_);
  tree->Branch("kaonLowPtVz_lead",                   &kaonLowPtVz_lead_);
  tree->Branch("kaonLowPtTrkChi2_lead",              &kaonLowPtTrkChi2_lead_);
  tree->Branch("kaonLowPtTrkNDOF_lead",              &kaonLowPtTrkNDOF_lead_);
  tree->Branch("kaonLowPtTrkNormChi2_lead",          &kaonLowPtTrkNormChi2_lead_);

  tree->Branch("kaonLowPtCharge_sublead",               &kaonLowPtCharge_sublead_);
  tree->Branch("kaonLowPtD0_sublead",                   &kaonLowPtD0_sublead_);
  tree->Branch("kaonLowPtDz_sublead",                   &kaonLowPtDz_sublead_);
  tree->Branch("kaonLowPtD0Error_sublead",              &kaonLowPtD0Error_sublead_);
  tree->Branch("kaonLowPtDzError_sublead",              &kaonLowPtDzError_sublead_);
  tree->Branch("kaonLowPtPt_sublead",                   &kaonLowPtPt_sublead_);
  tree->Branch("kaonLowPtEta_sublead",                  &kaonLowPtEta_sublead_);
  tree->Branch("kaonLowPtPhi_sublead",                  &kaonLowPtPhi_sublead_);
  tree->Branch("kaonLowPtVx_sublead",                   &kaonLowPtVx_sublead_);
  tree->Branch("kaonLowPtVy_sublead",                   &kaonLowPtVy_sublead_);
  tree->Branch("kaonLowPtVz_sublead",                   &kaonLowPtVz_sublead_);
  tree->Branch("kaonLowPtTrkChi2_sublead",              &kaonLowPtTrkChi2_sublead_);
  tree->Branch("kaonLowPtTrkNDOF_sublead",              &kaonLowPtTrkNDOF_sublead_);
  tree->Branch("kaonLowPtTrkNormChi2_sublead",          &kaonLowPtTrkNormChi2_sublead_);

  tree->Branch("bsLowPtdRele",                     &bsLowPtdRele_);
  tree->Branch("bsLowPtdRkaon",                    &bsLowPtdRkaon_);
  tree->Branch("bsLowPtdRJpsiPhi",                 &bsLowPtdRJpsiPhi_);
  tree->Branch("bsLowPtJpsiMass",                  &bsLowPtJpsiMass_);
  tree->Branch("bsLowPtPhiMass",                   &bsLowPtPhiMass_);
  tree->Branch("bsLowPtBsMass",                    &bsLowPtBsMass_);
  
}

void RphiNtuplizer::fillLowPtElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  // cleanup from previous execution
  lowPtCharge_lead_                  .clear();
  lowPtD0_lead_                      .clear();
  lowPtDz_lead_                      .clear();
  lowPtD0Error_lead_                 .clear();
  lowPtDzError_lead_                 .clear();
  lowPtPt_lead_                      .clear();
  lowPtEta_lead_                     .clear();
  lowPtPhi_lead_                     .clear();
  lowPtPtMean_lead_                  .clear();
  lowPtEtaMean_lead_                 .clear();
  lowPtPhiMean_lead_                 .clear();
  lowPtMVABWP_lead_                  .clear();
  lowPtMVAUnBWP_lead_                .clear();

  lowPtCharge_sublead_                  .clear();
  lowPtD0_sublead_                      .clear();
  lowPtDz_sublead_                      .clear();
  lowPtD0Error_sublead_                 .clear();
  lowPtDzError_sublead_                 .clear();
  lowPtPt_sublead_                      .clear();
  lowPtEta_sublead_                     .clear();
  lowPtPhi_sublead_                     .clear();
  lowPtPtMean_sublead_                  .clear();
  lowPtEtaMean_sublead_                 .clear();
  lowPtPhiMean_sublead_                 .clear();
  lowPtMVABWP_sublead_                  .clear();
  lowPtMVAUnBWP_sublead_                .clear();

  lowPtSvChi2_.clear();
  lowPtSvNDOF_.clear();
  lowPtSvProb_.clear();
  lowPtSvX_.clear();
  lowPtSvY_.clear();
  lowPtSvZ_.clear();
  lowPtSvXError_.clear();
  lowPtSvYError_.clear();
  lowPtSvZError_.clear();
  lowPtSvMass_.clear();
  lowPtSvCtxy_.clear();
  lowPtSvCosAngle_.clear();
  lowPtSvLxy_.clear();
  lowPtSvLxyError_.clear();

  kaonLowPtCharge_lead_                  .clear();
  kaonLowPtD0_lead_                      .clear();
  kaonLowPtDz_lead_                      .clear();
  kaonLowPtD0Error_lead_                 .clear();
  kaonLowPtDzError_lead_                 .clear();
  kaonLowPtPt_lead_                      .clear();
  kaonLowPtEta_lead_                     .clear();
  kaonLowPtPhi_lead_                     .clear();
  kaonLowPtVx_lead_                      .clear();
  kaonLowPtVy_lead_                      .clear();
  kaonLowPtVz_lead_                      .clear();
  kaonLowPtTrkChi2_lead_                 .clear();
  kaonLowPtTrkNDOF_lead_                 .clear();
  kaonLowPtTrkNormChi2_lead_             .clear();

  kaonLowPtCharge_sublead_                  .clear();
  kaonLowPtD0_sublead_                      .clear();
  kaonLowPtDz_sublead_                      .clear();
  kaonLowPtD0Error_sublead_                 .clear();
  kaonLowPtDzError_sublead_                 .clear();
  kaonLowPtPt_sublead_                      .clear();
  kaonLowPtEta_sublead_                     .clear();
  kaonLowPtPhi_sublead_                     .clear();
  kaonLowPtVx_sublead_                      .clear();
  kaonLowPtVy_sublead_                      .clear();
  kaonLowPtVz_sublead_                      .clear();
  kaonLowPtTrkChi2_sublead_                 .clear();
  kaonLowPtTrkNDOF_sublead_                 .clear();
  kaonLowPtTrkNormChi2_sublead_             .clear();

  bsLowPtdRele_		      .clear();
  bsLowPtdRkaon_                 .clear();
  bsLowPtdRJpsiPhi_              .clear();
  bsLowPtJpsiMass_		      .clear();
  bsLowPtPhiMass_		      .clear();
  bsLowPtBsMass_		      .clear();

  nLowPt_ = 0;

  if (isAOD_) {

    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<std::vector<reco::GsfElectron> > lowpTelectronHandle;
    e.getByToken(lowpTelectronlabel_, lowpTelectronHandle);

    edm::Handle<reco::TrackCollection> tracksHandle;
    e.getByToken(tracklabel_, tracksHandle);

    edm::Handle<edm::ValueMap<float> > ele_mva_wp_biased;
    e.getByToken( eleBWPToken_ ,ele_mva_wp_biased); 
    edm::Handle<edm::ValueMap<float> > ele_mva_wp_unbiased;
    e.getByToken( eleUnBWPToken_ ,ele_mva_wp_unbiased);

    edm::Handle<reco::ConversionCollection> conversions;
    e.getByToken(conversionsToken_, conversions);  

    if (!lowpTelectronHandle.isValid()) {
      edm::LogWarning("RphiNtuplizer") << "no low pT electrons in event";
      return;
    }

    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    VertexDistanceXY vertTool;

    for (reco::GsfElectronCollection::const_iterator iEle = lowpTelectronHandle->begin(); iEle != lowpTelectronHandle->end(); ++iEle) {
      if (iEle->gsfTrack()->ptMode() < lowPtPtLowCut_ || iEle->gsfTrack()->ptMode() > lowPtPtUpCut_) continue;
      if (fabs(iEle->eta()) > lowPtEtaCut_) continue;
      if (fabs(iEle->gsfTrack()->vz() - pv.z()) > lowPtDzCut_) continue;
      if (ConversionTools::hasMatchedConversion(*iEle, conversions, pv)) continue;

      
      // remove duplicated electrons
      bool duplicatedEle = false;
      if (removeDupEle_) {
	for (edm::View<pat::Electron>::const_iterator pfEle = electronHandle->begin(); pfEle != electronHandle->end(); ++pfEle) {
	  if (deltaR(iEle->eta(), iEle->phi(), pfEle->eta(), pfEle->phi()) < 0.001) {
	    duplicatedEle = true;
	    break;
	  }
	}
      }

      if (duplicatedEle) continue;
      

      for (reco::GsfElectronCollection::const_iterator jEle = iEle+1; jEle != lowpTelectronHandle->end(); ++jEle) {
	if (jEle->gsfTrack()->ptMode() < lowPtPtLowCut_ || jEle->gsfTrack()->ptMode() > lowPtPtUpCut_) continue;
	if (fabs(jEle->eta()) > lowPtEtaCut_) continue;
	if (fabs(jEle->gsfTrack()->vz() - pv.z()) > lowPtDzCut_) continue;
	if (ConversionTools::hasMatchedConversion(*jEle, conversions, pv)) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	//iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	//jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
      	iele_lv.SetPtEtaPhiM(iEle->gsfTrack()->ptMode(), iEle->gsfTrack()->etaMode(), iEle->gsfTrack()->phiMode(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->gsfTrack()->ptMode(), jEle->gsfTrack()->etaMode(), jEle->gsfTrack()->phiMode(), pmass);

	if ((iele_lv + jele_lv).M() < dilepMLowCut_ || (iele_lv + jele_lv).M() > dilepMUpCut_) continue;

	auto leadEle = iele_lv.Pt() > jele_lv.Pt() ? iEle : jEle;
	auto subleadEle = iele_lv.Pt() > jele_lv.Pt() ? jEle : iEle;

	float leadUnBDT = (*ele_mva_wp_unbiased)[leadEle->gsfTrack()] > (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()] ? (*ele_mva_wp_unbiased)[leadEle->gsfTrack()] : (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()];
	float subleadUnBDT = (*ele_mva_wp_unbiased)[leadEle->gsfTrack()] > (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()] ? (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()] : (*ele_mva_wp_unbiased)[leadEle->gsfTrack()];

	if (leadUnBDT < lowPtUnBWPLeadCut_) continue;
	if (subleadUnBDT < lowPtUnBWPSubleadCut_) continue;

	
	// remove duplicated electrons
	if (removeDupEle_) {
	  for (edm::View<pat::Electron>::const_iterator pfEle = electronHandle->begin(); pfEle != electronHandle->end(); ++pfEle) {
	    if (deltaR(jEle->eta(), jEle->phi(), pfEle->eta(), pfEle->phi()) < 0.001) {
	      duplicatedEle = true;
	      break;
	    }
	  }
	}

	if (duplicatedEle) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	//std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;

	//XParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));


	//KinematicConstrainedVertexFitter kvFitter;
	//RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;

	for (reco::TrackCollection::const_iterator iHad = tracksHandle->begin(); iHad != tracksHandle->end(); ++iHad) {
	  if (iHad->pt() < kaonPtCut_) continue;
	  if (fabs(iHad->eta()) > kaonEtaCut_) continue;
          if (fabs(iHad->vz() - pv.z()) > kaonDzCut_) continue;
	  if (iHad->normalizedChi2() < 0.0) continue;
	  if (iHad->normalizedChi2() > 20.0) continue;

	  for (reco::TrackCollection::const_iterator jHad = iHad+1; jHad != tracksHandle->end(); ++jHad) {
	    if (jHad->pt() < kaonPtCut_) continue;
	    if (fabs(jHad->eta()) > kaonEtaCut_) continue;
	    if (fabs(jHad->vz() - pv.z()) > kaonDzCut_) continue;
	    if (jHad->normalizedChi2() < 0.0) continue;
	    if (jHad->normalizedChi2() > 20.0) continue;
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

	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jHad) ), kpmass, 0.0, 0, kpmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));


	    KinematicConstrainedVertexFitter BsKvFitter;
	    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
	    if (!(BsKinVtx->isValid())) continue;

	    RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

	    if (DecayVtx->chiSquared() < 0.0) continue;
	    if (TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()) < svProbCut_) continue;

	    // Accept these 4 tracks as a Bs candidate, fill ntuple

	    auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
	    auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

	    double ctxy = ((DecayVtx->position().x() - pv.x())*bs_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_lv.Py())/(pow(bs_lv.Pt(),2))*bs_lv.M();
	    
	    math::XYZVector perp(bs_lv.Px(), bs_lv.Py(), 0.);
	    math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	    math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	    double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

	    if (cosAngle < cosAngleCut_) continue;

	    lowPtSvChi2_.push_back(DecayVtx->chiSquared());
	    lowPtSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    lowPtSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    lowPtSvX_.push_back(DecayVtx->position().x());
	    lowPtSvY_.push_back(DecayVtx->position().y());
	    lowPtSvZ_.push_back(DecayVtx->position().z());
	    lowPtSvXError_.push_back(DecayVtx->error().cxx());
	    lowPtSvYError_.push_back(DecayVtx->error().cyy());
	    lowPtSvZError_.push_back(DecayVtx->error().czz());
	    lowPtSvMass_.push_back(bs_lv.M());
	    lowPtSvCtxy_.push_back(ctxy);
	    lowPtSvCosAngle_.push_back(cosAngle);
	    lowPtSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    lowPtSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonLowPtCharge_lead_            .push_back(leadHad->charge());
	    kaonLowPtD0_lead_                .push_back(leadHad->dxy(pv));
	    kaonLowPtDz_lead_                .push_back(leadHad->dz(pv));
	    kaonLowPtD0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonLowPtDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonLowPtPt_lead_                .push_back(leadHad->pt());
	    kaonLowPtEta_lead_               .push_back(leadHad->eta());
	    kaonLowPtPhi_lead_               .push_back(leadHad->phi());
	    kaonLowPtVx_lead_ 		.push_back(leadHad->vx());
	    kaonLowPtVy_lead_ 		.push_back(leadHad->vy());
	    kaonLowPtVz_lead_ 		.push_back(leadHad->vz());
	    kaonLowPtTrkChi2_lead_ 		.push_back(leadHad->chi2());
	    kaonLowPtTrkNDOF_lead_ 		.push_back(leadHad->ndof());
	    kaonLowPtTrkNormChi2_lead_ 	.push_back(leadHad->normalizedChi2());

	    kaonLowPtCharge_sublead_         .push_back(subleadHad->charge());
	    kaonLowPtD0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonLowPtDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonLowPtD0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonLowPtDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonLowPtPt_sublead_             .push_back(subleadHad->pt());
	    kaonLowPtEta_sublead_            .push_back(subleadHad->eta());
	    kaonLowPtPhi_sublead_            .push_back(subleadHad->phi());
	    kaonLowPtVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonLowPtVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonLowPtVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonLowPtTrkChi2_sublead_ 	.push_back(subleadHad->chi2());
	    kaonLowPtTrkNDOF_sublead_ 	.push_back(subleadHad->ndof());
	    kaonLowPtTrkNormChi2_sublead_ 	.push_back(subleadHad->normalizedChi2());

	    bsLowPtdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsLowPtdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsLowPtdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsLowPtJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsLowPtPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsLowPtBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    lowPtCharge_lead_          .push_back(leadEle->gsfTrack()->charge());
	    lowPtD0_lead_              .push_back(leadEle->gsfTrack()->dxy(pv));
	    lowPtDz_lead_              .push_back(leadEle->gsfTrack()->dz(pv));
	    lowPtD0Error_lead_         .push_back(leadEle->gsfTrack()->dxyError());
	    lowPtDzError_lead_         .push_back(leadEle->gsfTrack()->dzError());
	    lowPtPt_lead_              .push_back(leadEle->gsfTrack()->ptMode());
	    lowPtEta_lead_             .push_back(leadEle->gsfTrack()->etaMode());
	    lowPtPhi_lead_             .push_back(leadEle->gsfTrack()->phiMode());
	    lowPtPtMean_lead_              .push_back(leadEle->gsfTrack()->pt());
	    lowPtEtaMean_lead_             .push_back(leadEle->gsfTrack()->eta());
	    lowPtPhiMean_lead_             .push_back(leadEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_lead = leadEle->gsfTrack();
	    lowPtMVABWP_lead_          .push_back((*ele_mva_wp_biased)[mvaSeed_lead]);
	    lowPtMVAUnBWP_lead_     .push_back((*ele_mva_wp_unbiased)[mvaSeed_lead]);

	    lowPtCharge_sublead_          .push_back(subleadEle->gsfTrack()->charge());
	    lowPtD0_sublead_              .push_back(subleadEle->gsfTrack()->dxy(pv));
	    lowPtDz_sublead_              .push_back(subleadEle->gsfTrack()->dz(pv));
	    lowPtD0Error_sublead_         .push_back(subleadEle->gsfTrack()->dxyError());
	    lowPtDzError_sublead_         .push_back(subleadEle->gsfTrack()->dzError());
	    lowPtPt_sublead_              .push_back(subleadEle->gsfTrack()->ptMode());
	    lowPtEta_sublead_             .push_back(subleadEle->gsfTrack()->etaMode());
	    lowPtPhi_sublead_             .push_back(subleadEle->gsfTrack()->phiMode());
	    lowPtPtMean_sublead_              .push_back(subleadEle->gsfTrack()->pt());
	    lowPtEtaMean_sublead_             .push_back(subleadEle->gsfTrack()->eta());
	    lowPtPhiMean_sublead_             .push_back(subleadEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_sublead = subleadEle->gsfTrack();
	    lowPtMVABWP_sublead_          .push_back((*ele_mva_wp_biased)[mvaSeed_sublead]);
	    lowPtMVAUnBWP_sublead_        .push_back((*ele_mva_wp_unbiased)[mvaSeed_sublead]);

	    nLowPt_++;
	  }
	}
      }
    }
  } else {

    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    edm::Handle<std::vector<reco::GsfElectron> > lowpTelectronHandle;
    e.getByToken(lowpTelectronlabel_, lowpTelectronHandle);

    edm::Handle<edm::ValueMap<float> > ele_mva_wp_biased;
    e.getByToken( eleBWPToken_ ,ele_mva_wp_biased); 
    edm::Handle<edm::ValueMap<float> > ele_mva_wp_unbiased;
    e.getByToken( eleUnBWPToken_ ,ele_mva_wp_unbiased);

    edm::Handle<reco::ConversionCollection> conversions;
    e.getByToken(conversionsToken_, conversions);  

    if (!lowpTelectronHandle.isValid()) {
      edm::LogWarning("RphiNtuplizer") << "no low pT electrons in event";
      return;
    }

    edm::Handle<pat::PackedCandidateCollection> pfcands;
    e.getByToken(pckPFCandidateCollection_, pfcands);

    edm::Handle<pat::PackedCandidateCollection> losttracks;
    e.getByToken(lostTracksLabel_, losttracks);

    std::vector<pat::PackedCandidate> alltracks;
    alltracks.reserve(pfcands->size() + losttracks->size());
    alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
    alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

    edm::Handle<reco::VertexCollection> recVtxs;
    e.getByToken(vtxLabel_, recVtxs);

    VertexDistanceXY vertTool;

    for (reco::GsfElectronCollection::const_iterator iEle = lowpTelectronHandle->begin(); iEle != lowpTelectronHandle->end(); ++iEle) {
      if (iEle->gsfTrack()->ptMode() < lowPtPtLowCut_ || iEle->gsfTrack()->ptMode() > lowPtPtUpCut_) continue;
      if (fabs(iEle->eta()) > lowPtEtaCut_) continue;
      if (fabs(iEle->gsfTrack()->vz() - pv.z()) > lowPtDzCut_) continue;
      if (ConversionTools::hasMatchedConversion(*iEle, conversions, pv)) continue;

      // remove duplicated electrons
      bool duplicatedEle = false;
      if (removeDupEle_) {
	for (edm::View<pat::Electron>::const_iterator pfEle = electronHandle->begin(); pfEle != electronHandle->end(); ++pfEle) {
	  if (deltaR(iEle->eta(), iEle->phi(), pfEle->eta(), pfEle->phi()) < 0.001) {
	    duplicatedEle = true;
	    break;
	  }
	}
      }
      if (duplicatedEle) continue;

      for (reco::GsfElectronCollection::const_iterator jEle = iEle+1; jEle != lowpTelectronHandle->end(); ++jEle) {
	if (jEle->gsfTrack()->ptMode() < lowPtPtLowCut_ || jEle->gsfTrack()->ptMode() > lowPtPtUpCut_) continue;
	if (fabs(jEle->eta()) > lowPtEtaCut_) continue;
	if (fabs(jEle->gsfTrack()->vz() - pv.z()) > lowPtDzCut_) continue;
	if (ConversionTools::hasMatchedConversion(*jEle, conversions, pv)) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
	if ((iele_lv + jele_lv).M() < dilepMLowCut_ || (iele_lv + jele_lv).M() > dilepMUpCut_) continue;

	auto leadEle = iEle->pt() > jEle->pt() ? iEle : jEle;
	auto subleadEle = iEle->pt() > jEle->pt() ? jEle : iEle;

	float leadUnBDT = (*ele_mva_wp_unbiased)[leadEle->gsfTrack()] > (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()] ? (*ele_mva_wp_unbiased)[leadEle->gsfTrack()] : (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()];
	float subleadUnBDT = (*ele_mva_wp_unbiased)[leadEle->gsfTrack()] > (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()] ? (*ele_mva_wp_unbiased)[subleadEle->gsfTrack()] : (*ele_mva_wp_unbiased)[leadEle->gsfTrack()];

	if (leadUnBDT < lowPtUnBWPLeadCut_) continue;
	if (subleadUnBDT < lowPtUnBWPSubleadCut_) continue;

	// remove duplicated electrons
	if (removeDupEle_) {
	  for (edm::View<pat::Electron>::const_iterator pfEle = electronHandle->begin(); pfEle != electronHandle->end(); ++pfEle) {
	    if (deltaR(jEle->eta(), jEle->phi(), pfEle->eta(), pfEle->phi()) < 0.001) {
	      duplicatedEle = true;
	      break;
	    }
	  }
	}
	if (duplicatedEle) continue;

	KinematicParticleFactoryFromTransientTrack pFactory;  
	//std::vector<RefCountedKinematicParticle> XParticles;
	float pmasse = 1.e-6 * pmass;

	//XParticles.push_back(pFactory.particle(getTransientTrack( *(iEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));
	//XParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	//KinematicConstrainedVertexFitter kvFitter;
	//RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles);

	//if (!(KinVtx->isValid()) || KinVtx->currentDecayVertex()->chiSquared() < 0.0) continue;

	for (pat::PackedCandidateCollection::const_iterator iHad = alltracks.begin(); iHad != alltracks.end(); ++iHad) {
	  if (iHad->pt() <= kaonPtCut_) continue;
          if (iHad->charge() == 0) continue;
          if (abs(iHad->pdgId()) != 211) continue;
          if (iHad->bestTrack() == nullptr) continue;
	  if (fabs(iHad->eta()) > kaonEtaCut_) continue;
	  if (fabs(iHad->vz() - pv.z()) > kaonDzCut_) continue;

	  for (pat::PackedCandidateCollection::const_iterator jHad = iHad+1; jHad != alltracks.end(); ++jHad) {
	    if (jHad->pt() <= kaonPtCut_) continue;
            if (jHad->charge() == 0) continue;
            if (abs(jHad->pdgId()) != 211) continue;
            if (jHad->bestTrack() == nullptr) continue;
	    if (fabs(jHad->eta()) > kaonEtaCut_) continue;
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
	    BsParticles.push_back(pFactory.particle(getTransientTrack( *(jEle->gsfTrack()) ), pmass, 0.0, 0, pmasse));

	    KinematicConstrainedVertexFitter BsKvFitter;
	    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
	    if (!(BsKinVtx->isValid())) continue;

	    RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();

	    if (DecayVtx->chiSquared() < 0.0) continue;
	    if (TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()) < svProbCut_) continue;

	    // Accept these 4 tracks as a Bs candidate, fill ntuple

	    auto leadHad = iHad->pt() > jHad->pt() ? iHad : jHad;
	    auto subleadHad = iHad->pt() > jHad->pt() ? jHad : iHad;

	    double ctxy = ((DecayVtx->position().x() - pv.x())*bs_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_lv.Py())/(pow(bs_lv.Pt(),2))*bs_lv.M();
	    
	    math::XYZVector perp(bs_lv.Px(), bs_lv.Py(), 0.);
	    math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
	    math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
	    double cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

	    if (cosAngle < cosAngleCut_) continue;

	    lowPtSvChi2_.push_back(DecayVtx->chiSquared());
	    lowPtSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    lowPtSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    lowPtSvX_.push_back(DecayVtx->position().x());
	    lowPtSvY_.push_back(DecayVtx->position().y());
	    lowPtSvZ_.push_back(DecayVtx->position().z());
	    lowPtSvXError_.push_back(DecayVtx->error().cxx());
	    lowPtSvYError_.push_back(DecayVtx->error().cyy());
	    lowPtSvZError_.push_back(DecayVtx->error().czz());
	    lowPtSvMass_.push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());
	    lowPtSvCtxy_.push_back(ctxy);
	    lowPtSvCosAngle_.push_back(cosAngle);
	    lowPtSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    lowPtSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonLowPtCharge_lead_            .push_back(leadHad->charge());
	    kaonLowPtD0_lead_                .push_back(leadHad->dxy(pv));
	    kaonLowPtDz_lead_                .push_back(leadHad->dz(pv));
	    kaonLowPtD0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonLowPtDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonLowPtPt_lead_                .push_back(leadHad->pt());
	    kaonLowPtEta_lead_               .push_back(leadHad->eta());
	    kaonLowPtPhi_lead_               .push_back(leadHad->phi());
	    kaonLowPtVx_lead_ 		.push_back(leadHad->vx());
	    kaonLowPtVy_lead_ 		.push_back(leadHad->vy());
	    kaonLowPtVz_lead_ 		.push_back(leadHad->vz());
	    kaonLowPtTrkChi2_lead_ 		.push_back(leadHad->bestTrack()->chi2());
	    kaonLowPtTrkNDOF_lead_ 		.push_back(leadHad->bestTrack()->ndof());
	    kaonLowPtTrkNormChi2_lead_ 	.push_back(leadHad->bestTrack()->normalizedChi2());

	    kaonLowPtCharge_sublead_         .push_back(subleadHad->charge());
	    kaonLowPtD0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonLowPtDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonLowPtD0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonLowPtDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonLowPtPt_sublead_             .push_back(subleadHad->pt());
	    kaonLowPtEta_sublead_            .push_back(subleadHad->eta());
	    kaonLowPtPhi_sublead_            .push_back(subleadHad->phi());
	    kaonLowPtVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonLowPtVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonLowPtVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonLowPtTrkChi2_sublead_ 	.push_back(subleadHad->bestTrack()->chi2());
	    kaonLowPtTrkNDOF_sublead_ 	.push_back(subleadHad->bestTrack()->ndof());
	    kaonLowPtTrkNormChi2_sublead_ 	.push_back(subleadHad->bestTrack()->normalizedChi2());

	    bsLowPtdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsLowPtdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsLowPtdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsLowPtJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsLowPtPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsLowPtBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    lowPtCharge_lead_          .push_back(leadEle->gsfTrack()->charge());
	    lowPtD0_lead_              .push_back(leadEle->gsfTrack()->dxy(pv));
	    lowPtDz_lead_              .push_back(leadEle->gsfTrack()->dz(pv));
	    lowPtD0Error_lead_         .push_back(leadEle->gsfTrack()->dxyError());
	    lowPtDzError_lead_         .push_back(leadEle->gsfTrack()->dzError());
	    lowPtPt_lead_              .push_back(leadEle->gsfTrack()->ptMode());
	    lowPtEta_lead_             .push_back(leadEle->gsfTrack()->etaMode());
	    lowPtPhi_lead_             .push_back(leadEle->gsfTrack()->phiMode());
	    lowPtPtMean_lead_              .push_back(leadEle->gsfTrack()->pt());
	    lowPtEtaMean_lead_             .push_back(leadEle->gsfTrack()->eta());
	    lowPtPhiMean_lead_             .push_back(leadEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_lead = leadEle->gsfTrack();
	    lowPtMVABWP_lead_          .push_back((*ele_mva_wp_biased)[mvaSeed_lead]);
	    lowPtMVAUnBWP_lead_     .push_back((*ele_mva_wp_unbiased)[mvaSeed_lead]);

	    lowPtCharge_sublead_          .push_back(subleadEle->gsfTrack()->charge());
	    lowPtD0_sublead_              .push_back(subleadEle->gsfTrack()->dxy(pv));
	    lowPtDz_sublead_              .push_back(subleadEle->gsfTrack()->dz(pv));
	    lowPtD0Error_sublead_         .push_back(subleadEle->gsfTrack()->dxyError());
	    lowPtDzError_sublead_         .push_back(subleadEle->gsfTrack()->dzError());
	    lowPtPt_sublead_              .push_back(subleadEle->gsfTrack()->ptMode());
	    lowPtEta_sublead_             .push_back(subleadEle->gsfTrack()->etaMode());
	    lowPtPhi_sublead_             .push_back(subleadEle->gsfTrack()->phiMode());
	    lowPtPtMean_sublead_              .push_back(subleadEle->gsfTrack()->pt());
	    lowPtEtaMean_sublead_             .push_back(subleadEle->gsfTrack()->eta());
	    lowPtPhiMean_sublead_             .push_back(subleadEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_sublead = subleadEle->gsfTrack();
	    lowPtMVABWP_sublead_          .push_back((*ele_mva_wp_biased)[mvaSeed_sublead]);
	    lowPtMVAUnBWP_sublead_        .push_back((*ele_mva_wp_unbiased)[mvaSeed_sublead]);


	    nLowPt_++;
	  }
	}
      }
    }

  }

}


