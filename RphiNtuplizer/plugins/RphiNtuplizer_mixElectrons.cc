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
Int_t          nMix_;
vector<int>    mixEleCharge_pf_;
vector<float>  mixEleD0_pf_;
vector<float>  mixEleDz_pf_;
vector<float>  mixEleD0Error_pf_;
vector<float>  mixEleDzError_pf_;
vector<float>  mixElePt_pf_;
vector<float>  mixEleEta_pf_;
vector<float>  mixElePhi_pf_;
vector<float>  mixEleGsfPt_pf_;
vector<float>  mixEleGsfEta_pf_;
vector<float>  mixEleGsfPhi_pf_;
vector<float>  mixEleMVABWP_pf_;
vector<float>  mixEleMVAUnBWP_pf_;

vector<int>    mixEleCharge_lowPt_;
vector<float>  mixEleD0_lowPt_;
vector<float>  mixEleDz_lowPt_;
vector<float>  mixEleD0Error_lowPt_;
vector<float>  mixEleDzError_lowPt_;
vector<float>  mixElePt_lowPt_;
vector<float>  mixEleEta_lowPt_;
vector<float>  mixElePhi_lowPt_;
vector<float>  mixElePtMean_lowPt_;
vector<float>  mixEleEtaMean_lowPt_;
vector<float>  mixElePhiMean_lowPt_;
vector<float>  mixEleMVABWP_lowPt_;
vector<float>  mixEleMVAUnBWP_lowPt_;

vector<float> mixEleSvChi2_;
vector<float> mixEleSvNDOF_;
vector<float> mixEleSvProb_;
vector<float> mixEleSvX_;
vector<float> mixEleSvY_;
vector<float> mixEleSvZ_;
vector<float> mixEleSvXError_;
vector<float> mixEleSvYError_;
vector<float> mixEleSvZError_;
vector<float> mixEleSvMass_;
vector<float> mixEleSvCtxy_;
vector<float> mixEleSvCosAngle_;
vector<float> mixEleSvLxy_;
vector<float> mixEleSvLxyError_;

vector<int>    kaonMixCharge_lead_;
vector<float>  kaonMixD0_lead_;
vector<float>  kaonMixDz_lead_;
vector<float>  kaonMixD0Error_lead_;
vector<float>  kaonMixDzError_lead_;
vector<float>  kaonMixPt_lead_;
vector<float>  kaonMixEta_lead_;
vector<float>  kaonMixPhi_lead_;
vector<float>  kaonMixVx_lead_;
vector<float>  kaonMixVy_lead_;
vector<float>  kaonMixVz_lead_;
vector<float>  kaonMixTrkChi2_lead_;
vector<float>  kaonMixTrkNDOF_lead_;
vector<float>  kaonMixTrkNormChi2_lead_;
vector<float>  kaonMixJPsiMass_lead_;
vector<float>  kaonMixPhiMass_lead_;

vector<int>    kaonMixCharge_sublead_;
vector<float>  kaonMixD0_sublead_;
vector<float>  kaonMixDz_sublead_;
vector<float>  kaonMixD0Error_sublead_;
vector<float>  kaonMixDzError_sublead_;
vector<float>  kaonMixPt_sublead_;
vector<float>  kaonMixEta_sublead_;
vector<float>  kaonMixPhi_sublead_;
vector<float>  kaonMixVx_sublead_;
vector<float>  kaonMixVy_sublead_;
vector<float>  kaonMixVz_sublead_;
vector<float>  kaonMixTrkChi2_sublead_;
vector<float>  kaonMixTrkNDOF_sublead_;
vector<float>  kaonMixTrkNormChi2_sublead_;

vector<float>  bsMixdRele_;
vector<float>  bsMixdRkaon_;
vector<float>  bsMixdRJpsiPhi_;
vector<float>  bsMixJpsiMass_;
vector<float>  bsMixPhiMass_;
vector<float>  bsMixBsMass_;

void RphiNtuplizer::branchesMixElectrons(TTree* tree) {

  tree->Branch("nMix",                    &nMix_);
  tree->Branch("mixEleCharge_pf",               &mixEleCharge_pf_);
  tree->Branch("mixEleD0_pf",                   &mixEleD0_pf_);
  tree->Branch("mixEleDz_pf",                   &mixEleDz_pf_);
  tree->Branch("mixEleD0Error_pf",              &mixEleD0Error_pf_);
  tree->Branch("mixEleDzError_pf",              &mixEleDzError_pf_);
  tree->Branch("mixElePt_pf",                   &mixElePt_pf_);
  tree->Branch("mixEleEta_pf",                  &mixEleEta_pf_);
  tree->Branch("mixElePhi_pf",                  &mixElePhi_pf_);
  tree->Branch("mixEleGsfPt_pf",                &mixEleGsfPt_pf_);
  tree->Branch("mixEleGsfEta_pf",               &mixEleGsfEta_pf_);
  tree->Branch("mixEleGsfPhi_pf",               &mixEleGsfPhi_pf_);

  tree->Branch("mixEleCharge_lowPt",               &mixEleCharge_lowPt_);
  tree->Branch("mixEleD0_lowPt",                   &mixEleD0_lowPt_);
  tree->Branch("mixEleDz_lowPt",                   &mixEleDz_lowPt_);
  tree->Branch("mixEleD0Error_lowPt",              &mixEleD0Error_lowPt_);
  tree->Branch("mixEleDzError_lowPt",              &mixEleDzError_lowPt_);
  tree->Branch("mixElePt_lowPt",                   &mixElePt_lowPt_);
  tree->Branch("mixEleEta_lowPt",                  &mixEleEta_lowPt_);
  tree->Branch("mixElePhi_lowPt",                  &mixElePhi_lowPt_);
  tree->Branch("mixElePtMean_lowPt",               &mixElePtMean_lowPt_);
  tree->Branch("mixEleEtaMean_lowPt",              &mixEleEtaMean_lowPt_);
  tree->Branch("mixElePhiMean_lowPt",              &mixElePhiMean_lowPt_);
  tree->Branch("mixEleMVABWP_lowPt",               &mixEleMVABWP_lowPt_);
  tree->Branch("mixEleMVAUnBWP_lowPt",             &mixEleMVAUnBWP_lowPt_);

  tree->Branch("mixEleSvChi2",                    &mixEleSvChi2_);
  tree->Branch("mixEleSvNDOF",                    &mixEleSvNDOF_);
  tree->Branch("mixEleSvProb",                    &mixEleSvProb_);
  tree->Branch("mixEleSvX",                       &mixEleSvX_);
  tree->Branch("mixEleSvY",                       &mixEleSvY_);
  tree->Branch("mixEleSvZ",                       &mixEleSvZ_);
  tree->Branch("mixEleSvXError",                  &mixEleSvXError_);
  tree->Branch("mixEleSvYError",                  &mixEleSvYError_);
  tree->Branch("mixEleSvZError",                  &mixEleSvZError_);
  tree->Branch("mixEleSvMass",                    &mixEleSvMass_);
  tree->Branch("mixEleSvCtxy",                    &mixEleSvCtxy_);
  tree->Branch("mixEleSvCosAngle",                    &mixEleSvCosAngle_);
  tree->Branch("mixEleSvLxy",                    	   &mixEleSvLxy_);
  tree->Branch("mixEleSvLxyError",                    &mixEleSvLxyError_);

  tree->Branch("kaonMixCharge_lead",               &kaonMixCharge_lead_);
  tree->Branch("kaonMixD0_lead",                   &kaonMixD0_lead_);
  tree->Branch("kaonMixDz_lead",                   &kaonMixDz_lead_);
  tree->Branch("kaonMixD0Error_lead",              &kaonMixD0Error_lead_);
  tree->Branch("kaonMixDzError_lead",              &kaonMixDzError_lead_);
  tree->Branch("kaonMixPt_lead",                   &kaonMixPt_lead_);
  tree->Branch("kaonMixEta_lead",                  &kaonMixEta_lead_);
  tree->Branch("kaonMixPhi_lead",                  &kaonMixPhi_lead_);
  tree->Branch("kaonMixVx_lead",                   &kaonMixVx_lead_);
  tree->Branch("kaonMixVy_lead",                   &kaonMixVy_lead_);
  tree->Branch("kaonMixVz_lead",                   &kaonMixVz_lead_);
  tree->Branch("kaonMixTrkChi2_lead",              &kaonMixTrkChi2_lead_);
  tree->Branch("kaonMixTrkNDOF_lead",              &kaonMixTrkNDOF_lead_);
  tree->Branch("kaonMixTrkNormChi2_lead",          &kaonMixTrkNormChi2_lead_);

  tree->Branch("kaonMixCharge_sublead",               &kaonMixCharge_sublead_);
  tree->Branch("kaonMixD0_sublead",                   &kaonMixD0_sublead_);
  tree->Branch("kaonMixDz_sublead",                   &kaonMixDz_sublead_);
  tree->Branch("kaonMixD0Error_sublead",              &kaonMixD0Error_sublead_);
  tree->Branch("kaonMixDzError_sublead",              &kaonMixDzError_sublead_);
  tree->Branch("kaonMixPt_sublead",                   &kaonMixPt_sublead_);
  tree->Branch("kaonMixEta_sublead",                  &kaonMixEta_sublead_);
  tree->Branch("kaonMixPhi_sublead",                  &kaonMixPhi_sublead_);
  tree->Branch("kaonMixVx_sublead",                   &kaonMixVx_sublead_);
  tree->Branch("kaonMixVy_sublead",                   &kaonMixVy_sublead_);
  tree->Branch("kaonMixVz_sublead",                   &kaonMixVz_sublead_);
  tree->Branch("kaonMixTrkChi2_sublead",              &kaonMixTrkChi2_sublead_);
  tree->Branch("kaonMixTrkNDOF_sublead",              &kaonMixTrkNDOF_sublead_);
  tree->Branch("kaonMixTrkNormChi2_sublead",          &kaonMixTrkNormChi2_sublead_);

  tree->Branch("bsMixdRele",                     &bsMixdRele_);
  tree->Branch("bsMixdRkaon",                    &bsMixdRkaon_);
  tree->Branch("bsMixdRJpsiPhi",                 &bsMixdRJpsiPhi_);
  tree->Branch("bsMixJpsiMass",                  &bsMixJpsiMass_);
  tree->Branch("bsMixPhiMass",                   &bsMixPhiMass_);
  tree->Branch("bsMixBsMass",                    &bsMixBsMass_);
  
}

void RphiNtuplizer::fillMixElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  // cleanup from previous execution
  mixEleCharge_pf_                  .clear();
  mixEleD0_pf_                      .clear();
  mixEleDz_pf_                      .clear();
  mixEleD0Error_pf_                 .clear();
  mixEleDzError_pf_                 .clear();
  mixElePt_pf_                      .clear();
  mixEleEta_pf_                     .clear();
  mixElePhi_pf_                     .clear();
  mixEleGsfPt_pf_                   .clear();
  mixEleGsfEta_pf_                  .clear();
  mixEleGsfPhi_pf_                  .clear();

  mixEleCharge_lowPt_                  .clear();
  mixEleD0_lowPt_                      .clear();
  mixEleDz_lowPt_                      .clear();
  mixEleD0Error_lowPt_                 .clear();
  mixEleDzError_lowPt_                 .clear();
  mixElePt_lowPt_                      .clear();
  mixEleEta_lowPt_                     .clear();
  mixElePhi_lowPt_                     .clear();
  mixElePtMean_lowPt_                  .clear();
  mixEleEtaMean_lowPt_                 .clear();
  mixElePhiMean_lowPt_                 .clear();
  mixEleMVABWP_lowPt_                  .clear();
  mixEleMVAUnBWP_lowPt_                .clear();

  mixEleSvChi2_.clear();
  mixEleSvNDOF_.clear();
  mixEleSvProb_.clear();
  mixEleSvX_.clear();
  mixEleSvY_.clear();
  mixEleSvZ_.clear();
  mixEleSvXError_.clear();
  mixEleSvYError_.clear();
  mixEleSvZError_.clear();
  mixEleSvMass_.clear();
  mixEleSvCtxy_.clear();
  mixEleSvCosAngle_.clear();
  mixEleSvLxy_.clear();
  mixEleSvLxyError_.clear();

  kaonMixCharge_lead_                  .clear();
  kaonMixD0_lead_                      .clear();
  kaonMixDz_lead_                      .clear();
  kaonMixD0Error_lead_                 .clear();
  kaonMixDzError_lead_                 .clear();
  kaonMixPt_lead_                      .clear();
  kaonMixEta_lead_                     .clear();
  kaonMixPhi_lead_                     .clear();
  kaonMixVx_lead_                      .clear();
  kaonMixVy_lead_                      .clear();
  kaonMixVz_lead_                      .clear();
  kaonMixTrkChi2_lead_                 .clear();
  kaonMixTrkNDOF_lead_                 .clear();
  kaonMixTrkNormChi2_lead_             .clear();

  kaonMixCharge_sublead_                  .clear();
  kaonMixD0_sublead_                      .clear();
  kaonMixDz_sublead_                      .clear();
  kaonMixD0Error_sublead_                 .clear();
  kaonMixDzError_sublead_                 .clear();
  kaonMixPt_sublead_                      .clear();
  kaonMixEta_sublead_                     .clear();
  kaonMixPhi_sublead_                     .clear();
  kaonMixVx_sublead_                      .clear();
  kaonMixVy_sublead_                      .clear();
  kaonMixVz_sublead_                      .clear();
  kaonMixTrkChi2_sublead_                 .clear();
  kaonMixTrkNDOF_sublead_                 .clear();
  kaonMixTrkNormChi2_sublead_             .clear();

  bsMixdRele_		      .clear();
  bsMixdRkaon_                 .clear();
  bsMixdRJpsiPhi_              .clear();
  bsMixJpsiMass_		      .clear();
  bsMixPhiMass_		      .clear();
  bsMixBsMass_		      .clear();

  nMix_ = 0;

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

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      if (iEle->pt() < mixPfPtCut_) continue;
      if (fabs(iEle->eta()) > mixPfEtaCut_) continue;
      if (fabs(iEle->vz() - pv.z()) > mixPfDzCut_) continue;
      if (ConversionTools::hasMatchedConversion(*iEle, conversions, pv)) continue;

      for (reco::GsfElectronCollection::const_iterator jEle = lowpTelectronHandle->begin(); jEle != lowpTelectronHandle->end(); ++jEle) {
	if (jEle->gsfTrack()->ptMode() < mixLowPtPtLowCut_ || jEle->gsfTrack()->ptMode() > mixLowPtPtUpCut_) continue;
	if (fabs(jEle->eta()) > mixLowPtEtaCut_) continue;
	if (fabs(jEle->gsfTrack()->vz() - pv.z()) > mixLowPtDzCut_) continue;
	if (ConversionTools::hasMatchedConversion(*jEle, conversions, pv)) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->gsfTrack()->ptMode(), jEle->gsfTrack()->etaMode(), jEle->gsfTrack()->phiMode(), pmass);

	if ((iele_lv + jele_lv).M() < dilepMLowCut_ || (iele_lv + jele_lv).M() > dilepMUpCut_) continue;

	//auto leadEle = iele_lv.Pt() > jele_lv.Pt() ? iEle : jEle;
	//auto subleadEle = iele_lv.Pt() > jele_lv.Pt() ? jEle : iEle;

	if ((*ele_mva_wp_unbiased)[jEle->gsfTrack()] < mixLowPtUnBWPCut_) continue;

	// remove duplicated electrons
	bool duplicatedEle = false;
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
	    if (jHad->pt() <  kaonPtCut_) continue;
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

	    mixEleSvChi2_.push_back(DecayVtx->chiSquared());
	    mixEleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    mixEleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    mixEleSvX_.push_back(DecayVtx->position().x());
	    mixEleSvY_.push_back(DecayVtx->position().y());
	    mixEleSvZ_.push_back(DecayVtx->position().z());
	    mixEleSvXError_.push_back(DecayVtx->error().cxx());
	    mixEleSvYError_.push_back(DecayVtx->error().cyy());
	    mixEleSvZError_.push_back(DecayVtx->error().czz());
	    mixEleSvMass_.push_back(bs_lv.M());
	    mixEleSvCtxy_.push_back(ctxy);
	    mixEleSvCosAngle_.push_back(cosAngle);
	    mixEleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    mixEleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonMixCharge_lead_            .push_back(leadHad->charge());
	    kaonMixD0_lead_                .push_back(leadHad->dxy(pv));
	    kaonMixDz_lead_                .push_back(leadHad->dz(pv));
	    kaonMixD0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonMixDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonMixPt_lead_                .push_back(leadHad->pt());
	    kaonMixEta_lead_               .push_back(leadHad->eta());
	    kaonMixPhi_lead_               .push_back(leadHad->phi());
	    kaonMixVx_lead_ 		.push_back(leadHad->vx());
	    kaonMixVy_lead_ 		.push_back(leadHad->vy());
	    kaonMixVz_lead_ 		.push_back(leadHad->vz());
	    kaonMixTrkChi2_lead_ 		.push_back(leadHad->chi2());
	    kaonMixTrkNDOF_lead_ 		.push_back(leadHad->ndof());
	    kaonMixTrkNormChi2_lead_ 	.push_back(leadHad->normalizedChi2());

	    kaonMixCharge_sublead_         .push_back(subleadHad->charge());
	    kaonMixD0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonMixDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonMixD0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonMixDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonMixPt_sublead_             .push_back(subleadHad->pt());
	    kaonMixEta_sublead_            .push_back(subleadHad->eta());
	    kaonMixPhi_sublead_            .push_back(subleadHad->phi());
	    kaonMixVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonMixVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonMixVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonMixTrkChi2_sublead_ 	.push_back(subleadHad->chi2());
	    kaonMixTrkNDOF_sublead_ 	.push_back(subleadHad->ndof());
	    kaonMixTrkNormChi2_sublead_ 	.push_back(subleadHad->normalizedChi2());

	    bsMixdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsMixdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsMixdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsMixJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsMixPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsMixBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    mixEleCharge_pf_          .push_back(iEle->charge());
	    mixEleD0_pf_              .push_back(iEle->gsfTrack()->dxy(pv));
	    mixEleDz_pf_              .push_back(iEle->gsfTrack()->dz(pv));
	    mixEleD0Error_pf_         .push_back(iEle->gsfTrack()->dxyError());
	    mixEleDzError_pf_         .push_back(iEle->gsfTrack()->dzError());
	    mixElePt_pf_              .push_back(iEle->pt());
	    mixEleEta_pf_             .push_back(iEle->eta());
	    mixElePhi_pf_             .push_back(iEle->phi());
	    mixEleGsfPt_pf_           .push_back(iEle->gsfTrack()->ptMode());
	    mixEleGsfEta_pf_          .push_back(iEle->gsfTrack()->etaMode());
	    mixEleGsfPhi_pf_          .push_back(iEle->gsfTrack()->phiMode());

	    mixEleCharge_lowPt_          .push_back(jEle->gsfTrack()->charge());
	    mixEleD0_lowPt_              .push_back(jEle->gsfTrack()->dxy(pv));
	    mixEleDz_lowPt_              .push_back(jEle->gsfTrack()->dz(pv));
	    mixEleD0Error_lowPt_         .push_back(jEle->gsfTrack()->dxyError());
	    mixEleDzError_lowPt_         .push_back(jEle->gsfTrack()->dzError());
	    mixElePt_lowPt_              .push_back(jEle->gsfTrack()->ptMode());
	    mixEleEta_lowPt_             .push_back(jEle->gsfTrack()->etaMode());
	    mixElePhi_lowPt_             .push_back(jEle->gsfTrack()->phiMode());
	    mixElePtMean_lowPt_              .push_back(jEle->gsfTrack()->pt());
	    mixEleEtaMean_lowPt_             .push_back(jEle->gsfTrack()->eta());
	    mixElePhiMean_lowPt_             .push_back(jEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_lowPt = jEle->gsfTrack();
	    mixEleMVABWP_lowPt_          .push_back((*ele_mva_wp_biased)[mvaSeed_lowPt]);
	    mixEleMVAUnBWP_lowPt_        .push_back((*ele_mva_wp_unbiased)[mvaSeed_lowPt]);

	    nMix_++;
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

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
      if (iEle->pt() < mixPfPtCut_) continue;
      if (fabs(iEle->eta()) > mixPfEtaCut_) continue;
      if (fabs(iEle->vz() - pv.z()) > mixPfDzCut_) continue;
      if (ConversionTools::hasMatchedConversion(*iEle, conversions, pv)) continue;

      for (reco::GsfElectronCollection::const_iterator jEle = lowpTelectronHandle->begin(); jEle != lowpTelectronHandle->end(); ++jEle) {
	if (jEle->gsfTrack()->ptMode() < mixLowPtPtLowCut_ || jEle->gsfTrack()->ptMode() > mixLowPtPtUpCut_) continue;
	if (fabs(jEle->eta()) > mixLowPtEtaCut_) continue;
	if (fabs(jEle->gsfTrack()->vz() - pv.z()) > mixLowPtDzCut_) continue;
	if (ConversionTools::hasMatchedConversion(*jEle, conversions, pv)) continue;
	//if (iEle->charge()*jEle->charge() > 0.0) continue;

	float pmass  = 0.0005109989461;
	TLorentzVector iele_lv, jele_lv;
	iele_lv.SetPtEtaPhiM(iEle->pt(), iEle->eta(), iEle->phi(), pmass);
	jele_lv.SetPtEtaPhiM(jEle->pt(), jEle->eta(), jEle->phi(), pmass);
	if ((iele_lv + jele_lv).M() < dilepMLowCut_ || (iele_lv + jele_lv).M() > dilepMUpCut_) continue;

	//auto leadEle = iEle->pt() > jEle->pt() ? iEle : jEle;
	//auto subleadEle = iEle->pt() > jEle->pt() ? jEle : iEle;

	if ((*ele_mva_wp_unbiased)[jEle->gsfTrack()] < mixLowPtUnBWPCut_) continue;

	// remove duplicated electrons
	bool duplicatedEle = false;
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

	    mixEleSvChi2_.push_back(DecayVtx->chiSquared());
	    mixEleSvNDOF_.push_back(DecayVtx->degreesOfFreedom());
	    mixEleSvProb_.push_back(TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom()));
	    mixEleSvX_.push_back(DecayVtx->position().x());
	    mixEleSvY_.push_back(DecayVtx->position().y());
	    mixEleSvZ_.push_back(DecayVtx->position().z());
	    mixEleSvXError_.push_back(DecayVtx->error().cxx());
	    mixEleSvYError_.push_back(DecayVtx->error().cyy());
	    mixEleSvZError_.push_back(DecayVtx->error().czz());
	    mixEleSvMass_.push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());
	    mixEleSvCtxy_.push_back(ctxy);
	    mixEleSvCosAngle_.push_back(cosAngle);
	    mixEleSvLxy_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).value());
	    mixEleSvLxyError_.push_back(vertTool.distance(vtx, DecayVtx.get()->vertexState()).error());

	    kaonMixCharge_lead_            .push_back(leadHad->charge());
	    kaonMixD0_lead_                .push_back(leadHad->dxy(pv));
	    kaonMixDz_lead_                .push_back(leadHad->dz(pv));
	    kaonMixD0Error_lead_ 		.push_back(leadHad->dxyError());
	    kaonMixDzError_lead_ 		.push_back(leadHad->dzError());
	    kaonMixPt_lead_                .push_back(leadHad->pt());
	    kaonMixEta_lead_               .push_back(leadHad->eta());
	    kaonMixPhi_lead_               .push_back(leadHad->phi());
	    kaonMixVx_lead_ 		.push_back(leadHad->vx());
	    kaonMixVy_lead_ 		.push_back(leadHad->vy());
	    kaonMixVz_lead_ 		.push_back(leadHad->vz());
	    kaonMixTrkChi2_lead_ 		.push_back(leadHad->bestTrack()->chi2());
	    kaonMixTrkNDOF_lead_ 		.push_back(leadHad->bestTrack()->ndof());
	    kaonMixTrkNormChi2_lead_ 	.push_back(leadHad->bestTrack()->normalizedChi2());

	    kaonMixCharge_sublead_         .push_back(subleadHad->charge());
	    kaonMixD0_sublead_             .push_back(subleadHad->dxy(pv));
	    kaonMixDz_sublead_             .push_back(subleadHad->dz(pv));
	    kaonMixD0Error_sublead_ 	.push_back(subleadHad->dxyError());
	    kaonMixDzError_sublead_ 	.push_back(subleadHad->dzError());
	    kaonMixPt_sublead_             .push_back(subleadHad->pt());
	    kaonMixEta_sublead_            .push_back(subleadHad->eta());
	    kaonMixPhi_sublead_            .push_back(subleadHad->phi());
	    kaonMixVx_sublead_ 		.push_back(subleadHad->vx());
	    kaonMixVy_sublead_ 		.push_back(subleadHad->vy());
	    kaonMixVz_sublead_ 		.push_back(subleadHad->vz());
	    kaonMixTrkChi2_sublead_ 	.push_back(subleadHad->bestTrack()->chi2());
	    kaonMixTrkNDOF_sublead_ 	.push_back(subleadHad->bestTrack()->ndof());
	    kaonMixTrkNormChi2_sublead_ 	.push_back(subleadHad->bestTrack()->normalizedChi2());

	    bsMixdRele_             .push_back(iele_lv.DeltaR(jele_lv));
	    bsMixdRkaon_            .push_back(iHad_lv.DeltaR(jHad_lv));
	    bsMixdRJpsiPhi_         .push_back((iele_lv+jele_lv).DeltaR(iHad_lv+jHad_lv));
	    bsMixJpsiMass_          .push_back((iele_lv+jele_lv).M());
	    bsMixPhiMass_           .push_back((iHad_lv+jHad_lv).M());
	    bsMixBsMass_            .push_back((iele_lv+jele_lv+iHad_lv+jHad_lv).M());

	    mixEleCharge_pf_          .push_back(iEle->charge());
	    mixEleD0_pf_              .push_back(iEle->gsfTrack()->dxy(pv));
	    mixEleDz_pf_              .push_back(iEle->gsfTrack()->dz(pv));
	    mixEleD0Error_pf_         .push_back(iEle->gsfTrack()->dxyError());
	    mixEleDzError_pf_         .push_back(iEle->gsfTrack()->dzError());
	    mixElePt_pf_              .push_back(iEle->pt());
	    mixEleEta_pf_             .push_back(iEle->eta());
	    mixElePhi_pf_             .push_back(iEle->phi());
	    mixEleGsfPt_pf_           .push_back(iEle->gsfTrack()->ptMode());
	    mixEleGsfEta_pf_          .push_back(iEle->gsfTrack()->etaMode());
	    mixEleGsfPhi_pf_          .push_back(iEle->gsfTrack()->phiMode());

	    mixEleCharge_lowPt_          .push_back(jEle->gsfTrack()->charge());
	    mixEleD0_lowPt_              .push_back(jEle->gsfTrack()->dxy(pv));
	    mixEleDz_lowPt_              .push_back(jEle->gsfTrack()->dz(pv));
	    mixEleD0Error_lowPt_         .push_back(jEle->gsfTrack()->dxyError());
	    mixEleDzError_lowPt_         .push_back(jEle->gsfTrack()->dzError());
	    mixElePt_lowPt_              .push_back(jEle->gsfTrack()->ptMode());
	    mixEleEta_lowPt_             .push_back(jEle->gsfTrack()->etaMode());
	    mixElePhi_lowPt_             .push_back(jEle->gsfTrack()->phiMode());
	    mixElePtMean_lowPt_              .push_back(jEle->gsfTrack()->pt());
	    mixEleEtaMean_lowPt_             .push_back(jEle->gsfTrack()->eta());
	    mixElePhiMean_lowPt_             .push_back(jEle->gsfTrack()->phi());

	    reco::GsfTrackRef mvaSeed_lowPt = jEle->gsfTrack();
	    mixEleMVABWP_lowPt_          .push_back((*ele_mva_wp_biased)[mvaSeed_lowPt]);
	    mixEleMVAUnBWP_lowPt_        .push_back((*ele_mva_wp_unbiased)[mvaSeed_lowPt]);


	    nMix_++;
	  }
	}
      }
    }

  }

}


