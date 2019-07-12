#ifndef RphiNtuplizer_h
#define RphiNtuplizer_h

#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

using namespace std;

void setbit(UShort_t& x, UShort_t bit);

class RphiNtuplizer : public edm::EDAnalyzer {
 public:

  explicit RphiNtuplizer(const edm::ParameterSet&);
  ~RphiNtuplizer();
  
  //   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  
  //   virtual void beginJob() {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //   virtual void endJob() {};
  
  void initTriggerFilters(const edm::Event&);
  ULong64_t matchMuonTriggerFilters(double pt, double eta, double phi);
  vector<bool> muonTriggerMap(const edm::Event&);
  Double_t deltaPhi(Double_t phi1, Double_t phi2);
  Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
  Double_t getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands, const reco::Candidate* ptcl,  
			    double r_iso_min, double r_iso_max, double kt_scale, bool charged_only);
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle);

  void branchesTriggers   (TTree*);
  void branchesMatchGenParticles    (TTree*);
  void branchesElectrons  (TTree*);
  void branchesMuons      (TTree*);
  void branchesLowPtElectrons(TTree*);
  void branchesMixElectrons(TTree*);

  void fillMatchGenParticles  (const edm::Event&, const edm::EventSetup&, math::XYZPoint&, const reco::Vertex);
  void fillElectrons  (const edm::Event&, const edm::EventSetup&, math::XYZPoint&, const reco::Vertex);
  void fillMuons      (const edm::Event&, math::XYZPoint&, const reco::Vertex);
  void fillLowPtElectrons(const edm::Event&, const edm::EventSetup&, math::XYZPoint&, const reco::Vertex);
  void fillMixElectrons(const edm::Event&, const edm::EventSetup&, math::XYZPoint&, const reco::Vertex);
  const reco::TransientTrack getTransientTrack(const reco::Track& track);
  const reco::TransientTrack getTransientTrack(const reco::GsfTrack& track);
  OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T"); 

  bool development_;
  bool isAOD_;
  bool doGenParticles_;
  bool dumpElectrons_;
  bool dumpMuons_;
  bool dumpLowPtElectrons_;
  bool dumpMixElectrons_;
  bool removeDupEle_;

  double trgFilterDeltaPtCut_;
  double trgFilterDeltaRCut_;

  edm::EDGetTokenT<reco::VertexCollection>         vtxLabel_;
  edm::EDGetTokenT<reco::VertexCollection>         vtxBSLabel_;
  edm::EDGetTokenT<double>                         rhoLabel_;
  edm::EDGetTokenT<double>                         rhoCentralLabel_;
  edm::EDGetTokenT<trigger::TriggerEvent>          trgEventLabel_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsLabel_;
  edm::EDGetTokenT<edm::TriggerResults>            trgResultsLabel_;
  string                                           trgResultsProcess_;
  edm::EDGetTokenT<edm::TriggerResults>            patTrgResultsLabel_;
  edm::EDGetTokenT<GenEventInfoProduct>            generatorLabel_;
  edm::EDGetTokenT<LHEEventProduct>                lheEventLabel_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<vector<reco::GenParticle> >     genParticlesCollection_;
  edm::EDGetTokenT<edm::View<pat::MET> >           pfMETlabel_;
  edm::EDGetTokenT<edm::View<pat::Electron> >      electronCollection_;
  edm::EDGetTokenT<edm::View<pat::Electron> >      calibelectronCollection_;
  edm::EDGetTokenT<edm::View<pat::Photon> >        photonCollection_;
  edm::EDGetTokenT<edm::View<pat::Photon> >        calibphotonCollection_;
  edm::EDGetTokenT<edm::View<pat::Muon> >          muonCollection_;
  edm::EDGetTokenT<vector<pat::Tau> >              tauCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           esReducedRecHitCollection_; 
  edm::EDGetTokenT<reco::PhotonCollection>         recophotonCollection_;
  edm::EDGetTokenT<reco::TrackCollection>          tracklabel_;
  edm::EDGetTokenT<reco::GsfElectronCollection>    gsfElectronlabel_;
  edm::EDGetTokenT<edm::View<reco::GsfTrack> >     gsfTracks_;
  edm::EDGetTokenT<reco::PFCandidateCollection>    pfAllParticles_;
  edm::EDGetTokenT<vector<pat::PackedCandidate> >  pckPFCdsLabel_;
  edm::EDGetTokenT<edm::View<reco::Candidate> >    recoCdsLabel_;
  edm::EDGetTokenT<edm::View<pat::Jet> >           jetsAK4Label_;
  edm::EDGetTokenT<edm::View<pat::Jet> >           jetsAK8Label_;
  edm::EDGetTokenT<reco::JetTagCollection>         boostedDoubleSVLabel_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pckPFCandidateCollection_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksLabel_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> >packedGenParticlesCollection_;
  edm::EDGetTokenT<reco::PFCandidateCollection>    pfCandidateCollection_;

  //check
  edm::EDGetTokenT< reco::DeDxDataValueMap > deDxProducer_;  
  edm::EDGetToken lowpTelectronlabel_;
  edm::EDGetTokenT<edm::ValueMap<float> > eleBWPToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > eleUnBWPToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

  double muPtCut_;
  double muEtaCut_;
  double muDzCut_;
  double elePtCut_;             
  double eleEtaCut_;            
  double eleDzCut_;             
  double mixPfPtCut_;           
  double mixPfEtaCut_;          
  double mixPfDzCut_;           
  double mixLowPtPtLowCut_;     
  double mixLowPtPtUpCut_;      
  double mixLowPtEtaCut_;       
  double mixLowPtDzCut_;        
  double mixLowPtUnBWPCut_;        
  double lowPtPtLowCut_;        
  double lowPtPtUpCut_;         
  double lowPtEtaCut_;          
  double lowPtDzCut_;           
  double lowPtUnBWPLeadCut_;    
  double lowPtUnBWPSubleadCut_; 
  double kaonPtCut_;
  double kaonEtaCut_;
  double kaonDzCut_;
  double dilepMLowCut_;
  double dilepMUpCut_;
  double phiMLowCut_;
  double phiMUpCut_;
  double bsMLowCut_;
  double bsMUpCut_;
  double svProbCut_;
  double cosAngleCut_;

  TTree   *tree_;
  TH1F    *hEvents_;
  TH1F    *htrgMudpT_;
  TH1F    *htrgMudR_;

  string processName_;
  HLTConfigProvider hltConfig_;
  HLTPrescaleProvider hltPrescaleProvider_;
};

#endif
