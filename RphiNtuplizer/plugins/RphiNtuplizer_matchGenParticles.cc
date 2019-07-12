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

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"

#include "RphiAnalysis/RphiNtuplizer/interface/RphiNtuplizer.h"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
float          gen_elep_pt_;
float	       gen_elep_eta_;
float          gen_elep_phi_;
float          gen_elem_pt_;
float	       gen_elem_eta_;
float          gen_elem_phi_;
float          gen_kp_pt_;
float	       gen_kp_eta_;
float          gen_kp_phi_;
float          gen_km_pt_;
float	       gen_km_eta_;
float          gen_km_phi_;
float          gen_ee_pt_;
float          gen_ee_eta_;
float          gen_ee_phi_;
float          gen_phi_pt_;
float          gen_phi_eta_;
float          gen_phi_phi_;
float          gen_bs_pt_;
float          gen_bs_eta_;
float          gen_bs_phi_;
float          gen_ee_m_;
float          gen_phi_m_;
float          gen_bs_m_;

bool           found_elep_reco_;
float          reco_elep_pt_;
float          reco_elep_eta_;
float          reco_elep_phi_;
float          reco_elep_d0_;
float          reco_elep_dz_;
float          reco_elep_d0error_;
float          reco_elep_dzerror_;
float          reco_elep_gsfpt_;
float          reco_elep_gsfeta_;
float          reco_elep_gsfphi_;

bool           found_elem_reco_;
float          reco_elem_pt_;
float          reco_elem_eta_;
float          reco_elem_phi_;
float          reco_elem_d0_;
float          reco_elem_dz_;
float          reco_elem_d0error_;
float          reco_elem_dzerror_;
float          reco_elem_gsfpt_;
float          reco_elem_gsfeta_;
float          reco_elem_gsfphi_;

bool           found_kp_reco_;
float          reco_kp_pt_;
float          reco_kp_eta_;
float          reco_kp_phi_;
float          reco_kp_d0_;
float          reco_kp_dz_;
float          reco_kp_d0error_;
float          reco_kp_dzerror_;
float          reco_kp_chi2_;
float          reco_kp_ndof_;
float          reco_kp_normchi2_;

bool           found_km_reco_;
float          reco_km_pt_;
float          reco_km_eta_;
float          reco_km_phi_;
float          reco_km_d0_;
float          reco_km_dz_;
float          reco_km_d0error_;
float          reco_km_dzerror_;
float          reco_km_chi2_;
float          reco_km_ndof_;
float          reco_km_normchi2_;

bool           found_sv_reco_;
float          reco_sv_chi2_;
float          reco_sv_ndof_;
float          reco_sv_prob_;
float          reco_sv_ctxy_;
float          reco_sv_cosangle_;
float          reco_sv_lxy_;
float          reco_sv_lxyerror_;

bool           found_kalmansv_reco_;
float          reco_kalmansv_chi2_;
float          reco_kalmansv_ndof_;
float          reco_kalmansv_prob_;
float          reco_kalmansv_ctxy_;
float          reco_kalmansv_cosangle_;
float          reco_kalmansv_lxy_;
float          reco_kalmansv_lxyerror_;

float          reco_ee_m_;
float          reco_phi_m_;
float          reco_bs_m_;


void RphiNtuplizer::branchesMatchGenParticles(TTree* tree) {

  tree->Branch("gen_elep_pt",               &gen_elep_pt_);
  tree->Branch("gen_elep_eta",              &gen_elep_eta_);
  tree->Branch("gen_elep_phi",              &gen_elep_phi_);
  tree->Branch("gen_elem_pt",               &gen_elem_pt_);
  tree->Branch("gen_elem_eta",              &gen_elem_eta_);
  tree->Branch("gen_elem_phi",              &gen_elem_phi_);
  tree->Branch("gen_kp_pt",                 &gen_kp_pt_);
  tree->Branch("gen_kp_eta",                &gen_kp_eta_);
  tree->Branch("gen_kp_phi",                &gen_kp_phi_);
  tree->Branch("gen_km_pt",                 &gen_km_pt_);
  tree->Branch("gen_km_eta",                &gen_km_eta_);
  tree->Branch("gen_km_phi",                &gen_km_phi_);
  tree->Branch("gen_ee_pt",                 &gen_ee_pt_);
  tree->Branch("gen_ee_eta",                &gen_ee_eta_);
  tree->Branch("gen_ee_phi",                &gen_ee_phi_);
  tree->Branch("gen_phi_pt",                &gen_phi_pt_);
  tree->Branch("gen_phi_eta",               &gen_phi_eta_);
  tree->Branch("gen_phi_phi",               &gen_phi_phi_);
  tree->Branch("gen_bs_pt",                 &gen_bs_pt_);
  tree->Branch("gen_bs_eta",                &gen_bs_eta_);
  tree->Branch("gen_bs_phi",                &gen_bs_phi_);
  tree->Branch("gen_ee_m",                  &gen_ee_m_);
  tree->Branch("gen_phi_m",                 &gen_phi_m_);
  tree->Branch("gen_bs_m",                  &gen_bs_m_);

  tree->Branch("found_elep_reco",           &found_elep_reco_);
  tree->Branch("reco_elep_pt",              &reco_elep_pt_);
  tree->Branch("reco_elep_eta",             &reco_elep_eta_);
  tree->Branch("reco_elep_phi",             &reco_elep_phi_);
  tree->Branch("reco_elep_d0",              &reco_elep_d0_);
  tree->Branch("reco_elep_dz",              &reco_elep_dz_);
  tree->Branch("reco_elep_d0error",         &reco_elep_d0error_);
  tree->Branch("reco_elep_dzerror",         &reco_elep_dzerror_);
  tree->Branch("reco_elep_gsfpt",           &reco_elep_gsfpt_);
  tree->Branch("reco_elep_gsfeta",          &reco_elep_gsfeta_);
  tree->Branch("reco_elep_gsfphi",          &reco_elep_gsfphi_);

  tree->Branch("found_elem_reco",           &found_elem_reco_);
  tree->Branch("reco_elem_pt",              &reco_elem_pt_);
  tree->Branch("reco_elem_eta",             &reco_elem_eta_);
  tree->Branch("reco_elem_phi",             &reco_elem_phi_);
  tree->Branch("reco_elem_d0",              &reco_elem_d0_);
  tree->Branch("reco_elem_dz",              &reco_elem_dz_);
  tree->Branch("reco_elem_d0error",         &reco_elem_d0error_);
  tree->Branch("reco_elem_dzerror",         &reco_elem_dzerror_);
  tree->Branch("reco_elem_gsfpt",           &reco_elem_gsfpt_);
  tree->Branch("reco_elem_gsfeta",          &reco_elem_gsfeta_);
  tree->Branch("reco_elem_gsfphi",          &reco_elem_gsfphi_);

  tree->Branch("found_kp_reco",             &found_kp_reco_);
  tree->Branch("reco_kp_pt",                &reco_kp_pt_);
  tree->Branch("reco_kp_eta",               &reco_kp_eta_);
  tree->Branch("reco_kp_phi",               &reco_kp_phi_);
  tree->Branch("reco_kp_d0",                &reco_kp_d0_);
  tree->Branch("reco_kp_dz",                &reco_kp_dz_);
  tree->Branch("reco_kp_d0error",           &reco_kp_d0error_);
  tree->Branch("reco_kp_dzerror",           &reco_kp_dzerror_);
  tree->Branch("reco_kp_chi2",              &reco_kp_chi2_);
  tree->Branch("reco_kp_ndof",              &reco_kp_ndof_);
  tree->Branch("reco_kp_normchi2",          &reco_kp_normchi2_);

  tree->Branch("found_km_reco",             &found_km_reco_);
  tree->Branch("reco_km_pt",                &reco_km_pt_);
  tree->Branch("reco_km_eta",               &reco_km_eta_);
  tree->Branch("reco_km_phi",               &reco_km_phi_);
  tree->Branch("reco_km_d0",                &reco_km_d0_);
  tree->Branch("reco_km_dz",                &reco_km_dz_);
  tree->Branch("reco_km_d0error",           &reco_km_d0error_);
  tree->Branch("reco_km_dzerror",           &reco_km_dzerror_);
  tree->Branch("reco_km_chi2",              &reco_km_chi2_);
  tree->Branch("reco_km_ndof",              &reco_km_ndof_);
  tree->Branch("reco_km_normchi2",          &reco_km_normchi2_);

  tree->Branch("found_sv_reco",             &found_sv_reco_);
  tree->Branch("reco_sv_chi2",              &reco_sv_chi2_);
  tree->Branch("reco_sv_ndof",              &reco_sv_ndof_);
  tree->Branch("reco_sv_prob",              &reco_sv_prob_);
  tree->Branch("reco_sv_ctxy",              &reco_sv_ctxy_);
  tree->Branch("reco_sv_cosangle",          &reco_sv_cosangle_);
  tree->Branch("reco_sv_lxy",               &reco_sv_lxy_);
  tree->Branch("reco_sv_lxyerror",          &reco_sv_lxyerror_);

  tree->Branch("found_kalmansv_reco",       &found_kalmansv_reco_);
  tree->Branch("reco_kalmansv_chi2",        &reco_kalmansv_chi2_);
  tree->Branch("reco_kalmansv_ndof",        &reco_kalmansv_ndof_);
  tree->Branch("reco_kalmansv_prob",        &reco_kalmansv_prob_);
  tree->Branch("reco_kalmansv_ctxy",        &reco_kalmansv_ctxy_);
  tree->Branch("reco_kalmansv_cosangle",    &reco_kalmansv_cosangle_);
  tree->Branch("reco_kalmansv_lxy",         &reco_kalmansv_lxy_);
  tree->Branch("reco_kalmansv_lxyerror",    &reco_kalmansv_lxyerror_);

  tree->Branch("reco_ee_m",                 &reco_ee_m_);
  tree->Branch("reco_phi_m",                &reco_phi_m_);
  tree->Branch("reco_bs_m",                 &reco_bs_m_);

  
}

void RphiNtuplizer::fillMatchGenParticles(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  gen_elep_pt_          = -99;
  gen_elep_eta_         = -99;
  gen_elep_phi_         = -99;
  gen_elem_pt_          = -99;
  gen_elem_eta_         = -99;
  gen_elem_phi_         = -99;
  gen_kp_pt_            = -99;
  gen_kp_eta_           = -99;
  gen_kp_phi_           = -99;
  gen_km_pt_            = -99;
  gen_km_eta_           = -99;
  gen_km_phi_           = -99;
  gen_ee_pt_            = -99;
  gen_ee_eta_           = -99;
  gen_ee_phi_           = -99;
  gen_phi_pt_           = -99;
  gen_phi_eta_          = -99;
  gen_phi_phi_          = -99;
  gen_bs_pt_            = -99;
  gen_bs_eta_           = -99;
  gen_bs_phi_           = -99;
  gen_ee_m_             = -99;
  gen_phi_m_            = -99;
  gen_bs_m_             = -99;

  found_elep_reco_      = false;
  reco_elep_pt_         = -99;
  reco_elep_eta_        = -99;
  reco_elep_phi_        = -99;
  reco_elep_d0_         = -99;
  reco_elep_dz_         = -99;
  reco_elep_d0error_    = -99;
  reco_elep_dzerror_    = -99;
  reco_elep_gsfpt_      = -99;
  reco_elep_gsfeta_     = -99;
  reco_elep_gsfphi_     = -99;

  found_elem_reco_      = false;
  reco_elem_pt_         = -99;
  reco_elem_eta_        = -99;
  reco_elem_phi_        = -99;
  reco_elem_d0_         = -99;
  reco_elem_dz_         = -99;
  reco_elem_d0error_    = -99;
  reco_elem_dzerror_    = -99;
  reco_elem_gsfpt_      = -99;
  reco_elem_gsfeta_     = -99;
  reco_elem_gsfphi_     = -99;

  found_kp_reco_        = false;
  reco_kp_pt_           = -99;
  reco_kp_eta_          = -99;
  reco_kp_phi_          = -99;
  reco_kp_d0_           = -99;
  reco_kp_dz_           = -99;
  reco_kp_d0error_      = -99;
  reco_kp_dzerror_      = -99;
  reco_kp_chi2_         = -99;
  reco_kp_ndof_         = -99;
  reco_kp_normchi2_     = -99;

  found_km_reco_        = false;
  reco_km_pt_           = -99;
  reco_km_eta_          = -99;
  reco_km_phi_          = -99;
  reco_km_d0_           = -99;
  reco_km_dz_           = -99;
  reco_km_d0error_      = -99;
  reco_km_dzerror_      = -99;
  reco_km_chi2_         = -99;
  reco_km_ndof_         = -99;
  reco_km_normchi2_     = -99;

  found_sv_reco_        = false;
  reco_sv_chi2_         = -99;
  reco_sv_ndof_         = -99;
  reco_sv_prob_         = -99;
  reco_sv_ctxy_         = -99;
  reco_sv_cosangle_     = -99;
  reco_sv_lxy_          = -99;
  reco_sv_lxyerror_     = -99;

  found_kalmansv_reco_        = false;
  reco_kalmansv_chi2_         = -99;
  reco_kalmansv_ndof_         = -99;
  reco_kalmansv_prob_         = -99;
  reco_kalmansv_ctxy_         = -99;
  reco_kalmansv_cosangle_     = -99;
  reco_kalmansv_lxy_          = -99;
  reco_kalmansv_lxyerror_     = -99;

  reco_ee_m_            = -99;
  reco_phi_m_           = -99;
  reco_bs_m_            = -99;


  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  //edm::Handle<edm::View<pat::PackedGenParticle> > packedGenParticlesHandle;
  //e.getByToken(packedGenParticlesCollection_, packedGenParticlesHandle);

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  edm::Handle<reco::TrackCollection> tracksHandle;
  e.getByToken(tracklabel_, tracksHandle);

  if (!electronHandle.isValid()) {
    edm::LogWarning("RphiNtuplizer") << "no pat::Electrons in event";
    return;
  }

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  VertexDistanceXY vertTool;

  bool found_elep_gen = false, found_elem_gen = false, found_kp_gen = false, found_km_gen = false;
  int elep_gen_id = -1, elem_gen_id = -1, kp_gen_id = -1, km_gen_id = -1;

  for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    // get Bs meson (pdgid=531)
    if (abs(ip->pdgId()) != 531) continue;
    found_elep_gen = false, found_elem_gen = false; found_kp_gen = false, found_km_gen = false;
    const Candidate* bMeson = &(*ip);

    // loop over all generated final state particles, get 
    //for (edm::View<pat::PackedGenParticle>::const_iterator ipp = packedGenParticlesHandle->begin(); ipp != packedGenParticlesHandle->end(); ++ipp) {
    for (vector<reco::GenParticle>::const_iterator ipp = genParticlesHandle->begin(); ipp != genParticlesHandle->end(); ++ipp) {
      // check if the generated final state particle comes from the Bs
      if (ipp->status() != 1) continue;
      if (! (ipp->mother(0) && isAncestor(bMeson, ipp->mother(0)))) continue;
      if (ipp->pdgId() == -abs(11)) {
	found_elep_gen = true;
	elep_gen_id = ipp - genParticlesHandle->begin();
      }
      if (ipp->pdgId() == abs(11)) {
	found_elem_gen = true;
	elem_gen_id = ipp - genParticlesHandle->begin();
      }
      if (ipp->pdgId() == 321) {
	found_kp_gen = true;
	kp_gen_id = ipp - genParticlesHandle->begin();
      }
      if (ipp->pdgId() == -321) {
	found_km_gen = true;
	km_gen_id = ipp - genParticlesHandle->begin();
      }
      if (found_elep_gen && found_elem_gen && found_kp_gen && found_km_gen) break;

    }
    if (found_elep_gen && found_elem_gen && found_kp_gen && found_km_gen) break;

  }

  if (! (found_elep_gen && found_elem_gen && found_kp_gen && found_km_gen)) return;

  double temp_dR = 9999999999.0, bestdR_elep = 9999999999.0, bestdR_elem = 9999999999.0, bestdR_kp = 9999999999.0, bestdR_km = 9999999999.0;
  bool found_elep_reco = false, found_elem_reco = false, found_kp_reco = false, found_km_reco =  false;
  int elep_reco_id = -1, elem_reco_id = -1, kp_reco_id = -1, km_reco_id = -1;

  // match the reco electrons to the gen electrons
  for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
    if (iEle->charge() > 0) {
      temp_dR = deltaR(iEle->eta(), iEle->phi(), genParticlesHandle->at(elep_gen_id).eta(), genParticlesHandle->at(elep_gen_id).phi());
      if (temp_dR < bestdR_elep) {
	bestdR_elep = temp_dR;
	elep_reco_id = iEle - electronHandle->begin();
      }
    } else {
      temp_dR = deltaR(iEle->eta(), iEle->phi(), genParticlesHandle->at(elem_gen_id).eta(), genParticlesHandle->at(elem_gen_id).phi());
      if (temp_dR < bestdR_elem) {
        bestdR_elem = temp_dR;
        elem_reco_id = iEle - electronHandle->begin();
      }
    }
  }
  
  for (reco::TrackCollection::const_iterator iHad = tracksHandle->begin(); iHad != tracksHandle->end(); ++iHad) {
     if (iHad->charge() > 0) {
      temp_dR = deltaR(iHad->eta(), iHad->phi(), genParticlesHandle->at(kp_gen_id).eta(), genParticlesHandle->at(kp_gen_id).phi());
      if (temp_dR < bestdR_kp) {
	bestdR_kp = temp_dR;
	kp_reco_id = iHad - tracksHandle->begin();
      }
    } else {
      temp_dR = deltaR(iHad->eta(), iHad->phi(), genParticlesHandle->at(km_gen_id).eta(), genParticlesHandle->at(km_gen_id).phi());
      if (temp_dR < bestdR_km) {
        bestdR_km = temp_dR;
        km_reco_id = iHad - tracksHandle->begin();
      }
    }
  }

  if (bestdR_elep < 0.1) found_elep_reco = true;
  if (bestdR_elem < 0.1) found_elem_reco = true;
  if (bestdR_kp < 0.1) found_kp_reco = true;
  if (bestdR_km < 0.1) found_km_reco = true;

  float eleM = 0.0005109989461;
  float eleMe = 1.e-6 * eleM;
  float kaonM = 0.493677;
  float kaonMe = 1.e-6 * kaonM;

  TLorentzVector elep_gen_lv, elem_gen_lv, kp_gen_lv, km_gen_lv, elep_reco_lv, elem_reco_lv, kp_reco_lv, km_reco_lv, ee_reco_lv, phi_reco_lv, bs_reco_lv;
  elep_gen_lv.SetPtEtaPhiM(genParticlesHandle->at(elep_gen_id).pt(), genParticlesHandle->at(elep_gen_id).eta(), genParticlesHandle->at(elep_gen_id).phi(), eleM);
  elem_gen_lv.SetPtEtaPhiM(genParticlesHandle->at(elem_gen_id).pt(), genParticlesHandle->at(elem_gen_id).eta(), genParticlesHandle->at(elem_gen_id).phi(), eleM);
  kp_gen_lv.SetPtEtaPhiM(genParticlesHandle->at(kp_gen_id).pt(), genParticlesHandle->at(kp_gen_id).eta(), genParticlesHandle->at(kp_gen_id).phi(), kaonM);
  km_gen_lv.SetPtEtaPhiM(genParticlesHandle->at(km_gen_id).pt(), genParticlesHandle->at(km_gen_id).eta(), genParticlesHandle->at(km_gen_id).phi(), kaonM);

  if (found_elep_reco && found_elem_reco && found_kp_reco && found_km_reco) {

    // kinematic vertex fit
    KinematicParticleFactoryFromTransientTrack pFactory;  
    std::vector<RefCountedKinematicParticle> BsParticles;

    BsParticles.push_back(pFactory.particle(getTransientTrack( tracksHandle->at(kp_reco_id) ), kaonM, 0.0, 0, kaonMe));
    BsParticles.push_back(pFactory.particle(getTransientTrack( tracksHandle->at(km_reco_id) ), kaonM, 0.0, 0, kaonMe));
    BsParticles.push_back(pFactory.particle(getTransientTrack( *(electronHandle->at(elep_reco_id).gsfTrack()) ), eleM, 0.0, 0, eleMe));
    BsParticles.push_back(pFactory.particle(getTransientTrack( *(electronHandle->at(elem_reco_id).gsfTrack()) ), eleM, 0.0, 0, eleMe));

    KinematicConstrainedVertexFitter BsKvFitter;
    RefCountedKinematicTree BsKinVtx = BsKvFitter.fit(BsParticles);
    if (BsKinVtx->isValid() && BsKinVtx->currentDecayVertex()->chiSquared() > 0.0) {
      RefCountedKinematicVertex DecayVtx = BsKinVtx->currentDecayVertex();
      elep_reco_lv.SetPtEtaPhiM(electronHandle->at(elep_reco_id).pt(), electronHandle->at(elep_reco_id).eta(), electronHandle->at(elep_reco_id).phi(), eleM);
      elem_reco_lv.SetPtEtaPhiM(electronHandle->at(elem_reco_id).pt(), electronHandle->at(elem_reco_id).eta(), electronHandle->at(elem_reco_id).phi(), eleM);
      kp_reco_lv.SetPtEtaPhiM(tracksHandle->at(kp_reco_id).pt(), tracksHandle->at(kp_reco_id).eta(), tracksHandle->at(kp_reco_id).phi(), kaonM);
      km_reco_lv.SetPtEtaPhiM(tracksHandle->at(km_reco_id).pt(), tracksHandle->at(km_reco_id).eta(), tracksHandle->at(km_reco_id).phi(), kaonM);
      ee_reco_lv = elep_reco_lv + elem_reco_lv;
      phi_reco_lv = kp_reco_lv + km_reco_lv;
      bs_reco_lv = ee_reco_lv + phi_reco_lv; 

      float ctxy = ((DecayVtx->position().x() - pv.x())*bs_reco_lv.Px() + (DecayVtx->position().y() - pv.y())*bs_reco_lv.Py())/(pow(bs_reco_lv.Pt(),2))*bs_reco_lv.M();
      
      math::XYZVector perp(bs_reco_lv.Px(), bs_reco_lv.Py(), 0.);
      math::XYZPoint dxybs(-1*(pv.x() - DecayVtx->position().x()), -1*(pv.y() - DecayVtx->position().y()), 0.);
      math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
      float cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

      found_sv_reco_ = true;
      reco_sv_chi2_ = DecayVtx->chiSquared();
      reco_sv_ndof_ = DecayVtx->degreesOfFreedom();
      reco_sv_prob_ = TMath::Prob(DecayVtx->chiSquared(), DecayVtx->degreesOfFreedom());
      reco_sv_ctxy_ = ctxy;
      reco_sv_cosangle_ = cosAngle;
      reco_sv_lxy_ = vertTool.distance(vtx, DecayVtx.get()->vertexState()).value();
      reco_sv_lxyerror_ = vertTool.distance(vtx, DecayVtx.get()->vertexState()).error();
      reco_ee_m_ = ee_reco_lv.M();
      reco_phi_m_ = phi_reco_lv.M();
      reco_bs_m_ = bs_reco_lv.M();
    }

    // kalman vertex fit
    KalmanVertexFitter theKalmanFitter(false);
    TransientVertex BsKalmanVertex;
    std::vector<reco::TransientTrack> tempTracks;
    tempTracks.push_back(getTransientTrack( tracksHandle->at(kp_reco_id) ));
    tempTracks.push_back(getTransientTrack( tracksHandle->at(km_reco_id) ));
    tempTracks.push_back(getTransientTrack( *(electronHandle->at(elep_reco_id).gsfTrack()) ));
    tempTracks.push_back(getTransientTrack( *(electronHandle->at(elem_reco_id).gsfTrack()) ));

    BsKalmanVertex = theKalmanFitter.vertex(tempTracks);
    if (BsKalmanVertex.isValid() && BsKalmanVertex.totalChiSquared() > 0.0) {
      elep_reco_lv.SetPtEtaPhiM(electronHandle->at(elep_reco_id).pt(), electronHandle->at(elep_reco_id).eta(), electronHandle->at(elep_reco_id).phi(), eleM);
      elem_reco_lv.SetPtEtaPhiM(electronHandle->at(elem_reco_id).pt(), electronHandle->at(elem_reco_id).eta(), electronHandle->at(elem_reco_id).phi(), eleM);
      kp_reco_lv.SetPtEtaPhiM(tracksHandle->at(kp_reco_id).pt(), tracksHandle->at(kp_reco_id).eta(), tracksHandle->at(kp_reco_id).phi(), kaonM);
      km_reco_lv.SetPtEtaPhiM(tracksHandle->at(km_reco_id).pt(), tracksHandle->at(km_reco_id).eta(), tracksHandle->at(km_reco_id).phi(), kaonM);
      ee_reco_lv = elep_reco_lv + elem_reco_lv;
      phi_reco_lv = kp_reco_lv + km_reco_lv;
      bs_reco_lv = ee_reco_lv + phi_reco_lv; 
   
      float ctxy = ((BsKalmanVertex.position().x() - pv.x())*bs_reco_lv.Px() + (BsKalmanVertex.position().y() - pv.y())*bs_reco_lv.Py())/(pow(bs_reco_lv.Pt(),2))*bs_reco_lv.M();
      
      math::XYZVector perp(bs_reco_lv.Px(), bs_reco_lv.Py(), 0.);
      math::XYZPoint dxybs(-1*(pv.x() - BsKalmanVertex.position().x()), -1*(pv.y() - BsKalmanVertex.position().y()), 0.);
      math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
      float cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

      found_kalmansv_reco_ = true;
      reco_kalmansv_chi2_ = BsKalmanVertex.totalChiSquared();
      reco_kalmansv_ndof_ = BsKalmanVertex.degreesOfFreedom();
      reco_kalmansv_prob_ = TMath::Prob(BsKalmanVertex.totalChiSquared(), BsKalmanVertex.degreesOfFreedom());
      reco_kalmansv_ctxy_ = ctxy;
      reco_kalmansv_cosangle_ = cosAngle;
      reco_kalmansv_lxy_ = vertTool.distance(vtx, BsKalmanVertex).value();
      reco_kalmansv_lxyerror_ = vertTool.distance(vtx, BsKalmanVertex).error();
      reco_ee_m_ = ee_reco_lv.M();
      reco_phi_m_ = phi_reco_lv.M();
      reco_bs_m_ = bs_reco_lv.M();

    }

  }

  // fill tree for gen particles
  auto elep_gen_it = genParticlesHandle->at(elep_gen_id);
  auto elem_gen_it = genParticlesHandle->at(elem_gen_id);
  auto kp_gen_it = genParticlesHandle->at(kp_gen_id);
  auto km_gen_it = genParticlesHandle->at(km_gen_id);

  gen_elep_pt_ = elep_gen_it.pt();
  gen_elep_eta_ = elep_gen_it.eta();
  gen_elep_phi_ = elep_gen_it.phi();
  gen_elem_pt_ = elem_gen_it.pt();
  gen_elem_eta_ = elem_gen_it.eta();
  gen_elem_phi_ = elem_gen_it.phi();
  gen_kp_pt_ = kp_gen_it.pt();
  gen_kp_eta_ = kp_gen_it.eta();
  gen_kp_phi_ = kp_gen_it.phi();
  gen_km_pt_ = km_gen_it.pt();
  gen_km_eta_ = km_gen_it.eta();
  gen_km_phi_ = km_gen_it.phi();
  gen_ee_pt_ = (elep_gen_lv + elem_gen_lv).Pt();
  gen_ee_eta_ = (elep_gen_lv + elem_gen_lv).Eta();
  gen_ee_phi_ = (elep_gen_lv + elem_gen_lv).Phi();
  gen_phi_pt_ = (kp_gen_lv + km_gen_lv).Pt();
  gen_phi_eta_ = (kp_gen_lv + km_gen_lv).Eta();
  gen_phi_phi_ = (kp_gen_lv + km_gen_lv).Phi();
  gen_bs_pt_ = (elep_gen_lv + elem_gen_lv + kp_gen_lv + km_gen_lv).Pt();
  gen_bs_eta_ = (elep_gen_lv + elem_gen_lv + kp_gen_lv + km_gen_lv).Eta();
  gen_bs_phi_ = (elep_gen_lv + elem_gen_lv + kp_gen_lv + km_gen_lv).Phi();
  gen_ee_m_ = (elep_gen_lv + elem_gen_lv).M();
  gen_phi_m_ = (kp_gen_lv + km_gen_lv).M();
  gen_bs_m_ = (elep_gen_lv + elem_gen_lv + kp_gen_lv + km_gen_lv).M();

  if (found_elep_reco) {
    auto elep_reco_it = electronHandle->at(elep_reco_id);
    found_elep_reco_ = true;
    reco_elep_pt_ = elep_reco_it.pt();
    reco_elep_eta_ = elep_reco_it.eta();
    reco_elep_phi_ = elep_reco_it.phi();
    reco_elep_d0_ = elep_reco_it.gsfTrack()->dxy(pv);
    reco_elep_dz_ = elep_reco_it.gsfTrack()->dz(pv);
    reco_elep_d0error_ = elep_reco_it.gsfTrack()->dxyError();
    reco_elep_dzerror_ = elep_reco_it.gsfTrack()->dzError();
    reco_elep_gsfpt_ = elep_reco_it.gsfTrack()->ptMode();
    reco_elep_gsfeta_ = elep_reco_it.gsfTrack()->etaMode();
    reco_elep_gsfphi_ = elep_reco_it.gsfTrack()->phiMode();

  }

  if (found_elem_reco) {
    auto elem_reco_it = electronHandle->at(elem_reco_id);
    found_elem_reco_ = true;
    reco_elem_pt_ = elem_reco_it.pt();
    reco_elem_eta_ = elem_reco_it.eta();
    reco_elem_phi_ = elem_reco_it.phi();
    reco_elem_d0_ = elem_reco_it.gsfTrack()->dxy(pv);
    reco_elem_dz_ = elem_reco_it.gsfTrack()->dz(pv);
    reco_elem_d0error_ = elem_reco_it.gsfTrack()->dxyError();
    reco_elem_dzerror_ = elem_reco_it.gsfTrack()->dzError();
    reco_elem_gsfpt_ = elem_reco_it.gsfTrack()->ptMode();
    reco_elem_gsfeta_ = elem_reco_it.gsfTrack()->etaMode();
    reco_elem_gsfphi_ = elem_reco_it.gsfTrack()->phiMode();

  }

  if (found_kp_reco) {
    auto kp_reco_it = tracksHandle->at(kp_reco_id);
    found_kp_reco_ = true;
    reco_kp_pt_ = kp_reco_it.pt();
    reco_kp_eta_ = kp_reco_it.eta();
    reco_kp_phi_ = kp_reco_it.phi();
    reco_kp_d0_ = kp_reco_it.dxy(pv);
    reco_kp_dz_ = kp_reco_it.dz(pv);
    reco_kp_d0error_ = kp_reco_it.dxyError();
    reco_kp_dzerror_ = kp_reco_it.dzError();
    reco_kp_chi2_ = kp_reco_it.chi2();
    reco_kp_ndof_ = kp_reco_it.ndof();
    reco_kp_normchi2_ = kp_reco_it.normalizedChi2();

  }

  if (found_km_reco) {
    auto km_reco_it = tracksHandle->at(km_reco_id);
    found_km_reco_ = true;
    reco_km_pt_ = km_reco_it.pt();
    reco_km_eta_ = km_reco_it.eta();
    reco_km_phi_ = km_reco_it.phi();
    reco_km_d0_ = km_reco_it.dxy(pv);
    reco_km_dz_ = km_reco_it.dz(pv);
    reco_km_d0error_ = km_reco_it.dxyError();
    reco_km_dzerror_ = km_reco_it.dzError();
    reco_km_chi2_ = km_reco_it.chi2();
    reco_km_ndof_ = km_reco_it.ndof();
    reco_km_normchi2_ = km_reco_it.normalizedChi2();

  }

}

bool RphiNtuplizer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++) {
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
    return false;
}



