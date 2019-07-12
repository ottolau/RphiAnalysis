#include "RphiAnalysis/RphiNtuplizer/interface/RphiNtuplizer.h"

using namespace std;
using namespace edm;

void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}

RphiNtuplizer::RphiNtuplizer(const edm::ParameterSet& ps) :
  hltPrescaleProvider_(ps, consumesCollector(), *this)
{

  development_               = ps.getParameter<bool>("development");
  doGenParticles_            = ps.getParameter<bool>("doGenParticles");
  isAOD_                     = ps.getParameter<bool>("isAOD");
  dumpElectrons_             = ps.getParameter<bool>("dumpElectrons");
  dumpMuons_                 = ps.getParameter<bool>("dumpMuons");
  dumpLowPtElectrons_        = ps.getParameter<bool>("dumpLowPtElectrons");
  dumpMixElectrons_          = ps.getParameter<bool>("dumpMixElectrons");
  removeDupEle_              = ps.getParameter<bool>("removeDupEle");

  trgFilterDeltaPtCut_       = ps.getParameter<double>("trgFilterDeltaPtCut");
  trgFilterDeltaRCut_        = ps.getParameter<double>("trgFilterDeltaRCut");

  vtxLabel_                  = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxLabel"));
  vtxBSLabel_                = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxBSLabel"));
  rhoLabel_                  = consumes<double>                        (ps.getParameter<InputTag>("rhoLabel"));
  rhoCentralLabel_           = consumes<double>                        (ps.getParameter<InputTag>("rhoCentralLabel"));
  trgEventLabel_             = consumes<trigger::TriggerEvent>         (ps.getParameter<InputTag>("triggerEvent"));
  triggerObjectsLabel_       = consumes<pat::TriggerObjectStandAloneCollection>(ps.getParameter<edm::InputTag>("triggerEvent"));
  trgResultsLabel_           = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("triggerResults"));
  patTrgResultsLabel_        = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("patTriggerResults"));
  trgResultsProcess_         =                                          ps.getParameter<InputTag>("triggerResults").process();
  generatorLabel_            = consumes<GenEventInfoProduct>           (ps.getParameter<InputTag>("generatorLabel"));
  lheEventLabel_             = consumes<LHEEventProduct>               (ps.getParameter<InputTag>("LHEEventLabel"));
  puCollection_              = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("pileupCollection"));
  genParticlesCollection_    = consumes<vector<reco::GenParticle> >    (ps.getParameter<InputTag>("genParticleSrc"));
  pfMETlabel_                = consumes<View<pat::MET> >               (ps.getParameter<InputTag>("pfMETLabel"));
  electronCollection_        = consumes<View<pat::Electron> >          (ps.getParameter<InputTag>("electronSrc"));
  calibelectronCollection_   = consumes<View<pat::Electron> >          (ps.getParameter<InputTag>("calibelectronSrc"));
  gsfTracks_                 = consumes<View<reco::GsfTrack>>          (ps.getParameter<InputTag>("gsfTrackSrc"));

  photonCollection_          = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("photonSrc"));
  calibphotonCollection_     = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("calibphotonSrc"));
  muonCollection_            = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("muonSrc"));
  ebReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("ebReducedRecHitCollection"));
  eeReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("eeReducedRecHitCollection"));
  esReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("esReducedRecHitCollection")); 
  recophotonCollection_      = consumes<reco::PhotonCollection>        (ps.getParameter<InputTag>("recoPhotonSrc"));
  tracklabel_                = consumes<reco::TrackCollection>         (ps.getParameter<InputTag>("TrackLabel"));
  gsfElectronlabel_          = consumes<reco::GsfElectronCollection>   (ps.getParameter<InputTag>("gsfElectronLabel"));
  tauCollection_             = consumes<vector<pat::Tau> >             (ps.getParameter<InputTag>("tauSrc"));
  pfAllParticles_            = consumes<reco::PFCandidateCollection>   (ps.getParameter<InputTag>("PFAllCandidates"));
  pckPFCandidateCollection_  = consumes<pat::PackedCandidateCollection>(ps.getParameter<InputTag>("packedPFCands"));
  pckPFCdsLabel_             = consumes<vector<pat::PackedCandidate>>  (ps.getParameter<InputTag>("packedPFCands"));
  recoCdsLabel_              = consumes<View<reco::Candidate>>         (ps.getParameter<InputTag>("packedPFCands"));
  lostTracksLabel_           = consumes<pat::PackedCandidateCollection>(ps.getParameter<InputTag>("lostTracks"));
  packedGenParticlesCollection_    = consumes<edm::View<pat::PackedGenParticle> >    (ps.getParameter<InputTag>("PackedGenParticleSrc"));

  jetsAK4Label_              = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("ak4JetSrc"));
  jetsAK8Label_              = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("ak8JetSrc"));
  //boostedDoubleSVLabel_      = consumes<reco::JetTagCollection>        (ps.getParameter<InputTag>("boostedDoubleSVLabel"));
  //jecAK8PayloadNames_        =                                          ps.getParameter<std::vector<std::string> >("jecAK8PayloadNames"); 

  //pfLooseId_                 = ps.getParameter<ParameterSet>("pfLooseId");

  deDxProducer_ = consumes<reco::DeDxDataValueMap>(edm::InputTag("dedxHarmonic2"));

  lowpTelectronlabel_        = consumes<std::vector<reco::GsfElectron> >(ps.getParameter<edm::InputTag>("lowpTelectrons"));
  eleBWPToken_               = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("eleBiasedWP"));
  eleUnBWPToken_             = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("eleUnbiasedWP"));
  conversionsToken_          = consumes<reco::ConversionCollection >(ps.getParameter<edm::InputTag>("conversions"));

  muPtCut_                   = ps.getParameter<double>("muPtCut");
  muEtaCut_                  = ps.getParameter<double>("muEtaCut");
  muDzCut_                   = ps.getParameter<double>("muDzCut");
  elePtCut_                  = ps.getParameter<double>("elePtCut");
  eleEtaCut_                 = ps.getParameter<double>("eleEtaCut");
  eleDzCut_                  = ps.getParameter<double>("eleDzCut");
  mixPfPtCut_                = ps.getParameter<double>("mixPfPtCut");
  mixPfEtaCut_               = ps.getParameter<double>("mixPfEtaCut");
  mixPfDzCut_                = ps.getParameter<double>("mixPfDzCut");
  mixLowPtPtLowCut_          = ps.getParameter<double>("mixLowPtPtLowCut");
  mixLowPtPtUpCut_           = ps.getParameter<double>("mixLowPtPtUpCut");
  mixLowPtEtaCut_            = ps.getParameter<double>("mixLowPtEtaCut");
  mixLowPtDzCut_             = ps.getParameter<double>("mixLowPtDzCut");
  mixLowPtUnBWPCut_          = ps.getParameter<double>("mixLowPtUnBWPCut");
  lowPtPtLowCut_             = ps.getParameter<double>("lowPtPtLowCut");
  lowPtPtUpCut_              = ps.getParameter<double>("lowPtPtUpCut");
  lowPtEtaCut_               = ps.getParameter<double>("lowPtEtaCut");
  lowPtDzCut_                = ps.getParameter<double>("lowPtDzCut");
  lowPtUnBWPLeadCut_         = ps.getParameter<double>("lowPtUnBWPLeadCut");
  lowPtUnBWPSubleadCut_      = ps.getParameter<double>("lowPtUnBWPSubleadCut");

  kaonPtCut_                 = ps.getParameter<double>("kaonPtCut");
  kaonEtaCut_                = ps.getParameter<double>("kaonEtaCut");
  kaonDzCut_                 = ps.getParameter<double>("kaonDzCut");
  dilepMLowCut_              = ps.getParameter<double>("dilepMLowCut");
  dilepMUpCut_               = ps.getParameter<double>("dilepMUpCut");
  phiMLowCut_                = ps.getParameter<double>("phiMLowCut");
  phiMUpCut_                 = ps.getParameter<double>("phiMUpCut");
  bsMLowCut_                 = ps.getParameter<double>("bsMLowCut");
  bsMUpCut_                  = ps.getParameter<double>("bsMUpCut");
  svProbCut_                 = ps.getParameter<double>("svProbCut");
  cosAngleCut_               = ps.getParameter<double>("cosAngleCut");


  Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "Event data (tag V08_00_26_03)");
  hEvents_ = fs->make<TH1F>("hEvents",    "total processed and skimmed events",   2,  0,   2);

  htrgMudpT_ = fs->make<TH1F>("htrgMudpT",    "htrgMudpT",   100,  0,   1);
  htrgMudR_  = fs->make<TH1F>("htrgMudR",     "htrgMudR",    100,  0, 0.2);

  branchesTriggers(tree_);

  if (doGenParticles_) {
    branchesMatchGenParticles(tree_);
  }

  if (dumpElectrons_) branchesElectrons(tree_);
  if (dumpMuons_) branchesMuons(tree_);
  if (dumpLowPtElectrons_) branchesLowPtElectrons(tree_);
  if (dumpMixElectrons_) branchesMixElectrons(tree_);

}

RphiNtuplizer::~RphiNtuplizer() {

}

void RphiNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es) {

  hEvents_->Fill(0.5);

  // Update when necessary the trigger table
  
  bool changedConfig = true;
  if (!hltConfig_.init(e.getRun(), es, trgResultsProcess_, changedConfig)) {
    edm::LogError("RphiNtuplizer") << "Initialization of HLTConfigProvider failed!";
    return ;
  }
  
  bool cfg_changed = true;
  hltPrescaleProvider_.init(e.getRun(), es, trgResultsProcess_, cfg_changed);

  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);

  initTriggerFilters(e);

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);
  double thetagmu_vz = 0.0;

  if (! doGenParticles_) {
    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      if (matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()) == 1) {
        thetagmu_vz = iMu->vz();
        break;
      }
    }
  } else {
    edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
    e.getByToken(genParticlesCollection_, genParticlesHandle);
    float thetagmuPt = -99.0;
    int thetagmu_gen_id = -1;
    for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
      if (abs(ip->pdgId()) != 13 || ip->status() != 1) continue;
      if (ip->pt() > thetagmuPt) {
        thetagmuPt = ip->pt();
        thetagmu_gen_id = ip - genParticlesHandle->begin();
      }
    }
    if (thetagmu_gen_id != -1) {
      double temp_dR = 9999999999.0, bestdR_thetagmu = 999999999.0;
      for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
        temp_dR = deltaR(iMu->eta(), iMu->phi(), genParticlesHandle->at(thetagmu_gen_id).eta(), genParticlesHandle->at(thetagmu_gen_id).phi());
        if (temp_dR < bestdR_thetagmu) {
              bestdR_thetagmu = temp_dR;
              thetagmu_vz =  iMu->vz();
        }
      }
    }
  }

  reco::Vertex vtx;
  double thetagmuPvDz = 1.e+10;
  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0);
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    // replace isFake() for miniAOD since it requires tracks while miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = isAOD_ ? v->isFake() : (v->chi2() == 0 && v->ndof() == 0);

    if (!isFake) {
      if (fabs(v->z() - thetagmu_vz) < thetagmuPvDz) {
        pv.SetXYZ(v->x(), v->y(), v->z());
        vtx = *v;
        thetagmuPvDz = fabs(v->z() - thetagmu_vz);
      }
    }
  }

  
  if (!e.isRealData()) {
    if (doGenParticles_)
      fillMatchGenParticles(e, es, pv, vtx);
  }
  

  if (dumpElectrons_) fillElectrons(e, es, pv, vtx);
  if (dumpMuons_) fillMuons(e, pv, vtx);
  if (dumpLowPtElectrons_) fillLowPtElectrons(e, es, pv, vtx);
  if (dumpMixElectrons_) fillMixElectrons(e, es, pv, vtx);

  hEvents_->Fill(1.5);
  tree_->Fill();

}

const reco::TransientTrack RphiNtuplizer::getTransientTrack(const reco::Track& track) {   
    //OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T"); 
    const reco::TransientTrack transientTrack(track, paramField);
    return transientTrack;
}

const reco::TransientTrack RphiNtuplizer::getTransientTrack(const reco::GsfTrack& track) {    
    reco::GsfTransientTrack gsfTransientTrack(track, paramField);
    reco::TransientTrack transientTrack(gsfTransientTrack.track(), paramField);
    return transientTrack;
}

Double_t RphiNtuplizer::deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}

Double_t RphiNtuplizer::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {

  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}

// void RphiNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
// {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }

DEFINE_FWK_MODULE(RphiNtuplizer);
