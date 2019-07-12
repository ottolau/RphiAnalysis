#include "map"
#include "FWCore/Common/interface/TriggerNames.h"
#include "RphiAnalysis/RphiNtuplizer/interface/RphiNtuplizer.h"

using namespace std;

// local variables: per-filter per-electron/muon/photon/jet arrays of matched trigger objects
// NOTE: number of elements in the arrays equals sizeof(ULong64_t)
vector<float> trgMuPt[64],  trgMuEta[64],  trgMuPhi[64];

// (local) variables associated with tree branches
Int_t          nTrg_;
vector<float>  trgPt_;
vector<float>  trgEta_;
vector<float>  trgPhi_;
Int_t          hltMu9IP6_;
vector<std::string> trgPath_;

void RphiNtuplizer::branchesTriggers(TTree* tree) {

  tree->Branch("nTrg",                    &nTrg_);
  tree->Branch("trgPt",                   &trgPt_);
  tree->Branch("trgEta",                  &trgEta_);
  tree->Branch("trgPhi",                  &trgPhi_);
  tree->Branch("hltMu9IP6",               &hltMu9IP6_);
  tree->Branch("trgPath",                 &trgPath_);

}

void RphiNtuplizer::initTriggerFilters(const edm::Event &e) {
  // Fills the arrays above.

  // cleanup from previous execution
  for (size_t i = 0; i < 64; ++i) {
    trgMuPt  [i].clear();
    trgMuEta [i].clear();
    trgMuPhi [i].clear();
  }

  // filter => index (in trg*[] arrays) mappings
  static std::map<string,size_t> muFilters;

 
  // cleanup from previous execution
  trgPt_                      .clear();
  trgEta_                     .clear();
  trgPhi_                     .clear();
  trgPath_                    .clear();
  nTrg_ = 0;
  hltMu9IP6_ = 0;
  std::string parkingPathName = "youcantmatchthisstring234235";

  if (isAOD_) {
    edm::Handle<edm::TriggerResults> trgResultsHandle;
    e.getByToken(trgResultsLabel_, trgResultsHandle);

    // Get the HLT trigger path name, desired path name: HLT_Mu9_IP6_part2_v1

    const edm::TriggerNames& trigNames = e.triggerNames(*trgResultsHandle);   
    for (unsigned int iTrig=0; iTrig<trigNames.size(); ++iTrig) {
      if (trgResultsHandle->accept(iTrig)) {
        std::string pathName = trigNames.triggerName(iTrig);
        //std::cout << "[PATH]: " << pathName << std::endl;
        trgPath_.push_back(pathName);
        //if (pathName.find("HLT_Mu9_IP") !=std::string::npos || pathName.find("HLT_Mu8p5_IP3p5") !=std::string::npos || pathName.find("HLT_Mu8p5_IP3p5") !=std::string::npos) {
        if (pathName.find("HLT_") != std::string::npos && pathName.find("Mu") != std::string::npos && pathName.find("IP") != std::string::npos) {
          hltMu9IP6_ = 1;
          parkingPathName = pathName;
        //std::cout << "[PATH]: " << pathName << std::endl; 
        }
      }
    }

    if (hltMu9IP6_ == 0) return;

    edm::Handle<trigger::TriggerEvent> triggerHandle;
    e.getByToken(trgEventLabel_, triggerHandle);

    const trigger::TriggerObjectCollection& trgObjects = triggerHandle->getObjects();
    const std::vector<std::string> tagName = hltConfig_.moduleLabels(parkingPathName);

    // loop over particular filters (and not over full HLTs)
    for (trigger::size_type iF = 0; iF != triggerHandle->sizeFilters(); ++iF) {
      // full filter name and its keys each corresponding to a matched (pt, eta, phi, ...) object
      string const&        label = triggerHandle->filterTag(iF).label();

      // Matching between path and filter. Only select filters which belong to our desired path
      bool isParkingTrig = false;
      for (std::vector<std::string>::const_iterator tag = tagName.begin(); tag != tagName.end(); ++tag) {
        if (tag->find(label) != std::string::npos || label.find(*tag) != std::string::npos) {
          isParkingTrig = true;
        }
      }
      
      if (!isParkingTrig) continue;

      const trigger::Keys& keys  = triggerHandle->filterKeys(iF);

      size_t idx = 0;

      for (size_t iK = 0; iK < keys.size(); ++iK) {
        const trigger::TriggerObject& trgV = trgObjects.at(keys[iK]);
        if (abs(trgV.id()) != 13) continue;
        //cout<<trgV.id()<<endl;
        trgMuPt [idx].push_back(trgV.pt());
        trgMuEta[idx].push_back(trgV.eta());
        trgMuPhi[idx].push_back(trgV.phi());

        trgPt_.push_back(trgV.pt());
        trgEta_.push_back(trgV.eta());
        trgPhi_.push_back(trgV.phi());
        nTrg_++;
      }
    } // HLT filter loop

  } else {

    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerHandleMiniAOD;
    e.getByToken(triggerObjectsLabel_, triggerHandleMiniAOD);

    edm::Handle<edm::TriggerResults> trgResultsHandle;
    e.getByToken(trgResultsLabel_, trgResultsHandle);

    edm::Handle<std::string> filterLabels_;
    e.getByLabel("slimmedPatTrigger:filterLabels", filterLabels_);

    const edm::TriggerNames &names = e.triggerNames(*trgResultsHandle);

    for (pat::TriggerObjectStandAlone obj : *triggerHandleMiniAOD) {
      obj.unpackPathNames(names);
      //obj.unpackPathNames(e);
      obj.unpackFilterLabels(e, *trgResultsHandle);

      // loop over filters    
      for (size_t iF = 0; iF < obj.pathNames().size(); ++iF) {
        string pathName = obj.pathNames()[iF];
        trgPath_.push_back(pathName);
        if (pathName.find("HLT_") != std::string::npos && pathName.find("Mu") != std::string::npos && pathName.find("IP") != std::string::npos) {
          hltMu9IP6_ = 1;
          parkingPathName = pathName;
          //for (size_t iD = 0; iD < obj.filterIds().size(); ++iD) {
          //  std::cout<<"type id: "<<obj.filterIds()[iD]<<std::endl;
          //}
          //std::cout << "[PATH]: " << pathName << std::endl; 
          size_t idx = 0;
          trgMuPt [idx].push_back(obj.pt());
          trgMuEta[idx].push_back(obj.eta());
          trgMuPhi[idx].push_back(obj.phi());

          trgPt_.push_back(obj.pt());
          trgEta_.push_back(obj.eta());
          trgPhi_.push_back(obj.phi());
          nTrg_++;

        }
      }
    }


  }


  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);

  // Muon trigger matching

  for (size_t v=0; v<trgMuPt[0].size(); ++v){
    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      float pt = iMu->pt();
      float eta = iMu->eta();
      float phi = iMu->phi();
      if (fabs(pt - trgMuPt[0][v])/trgMuPt[0][v] < trgFilterDeltaPtCut_ &&
        deltaR(eta, phi, trgMuEta[0][v], trgMuPhi[0][v]) < trgFilterDeltaRCut_) {
        htrgMudpT_->Fill(fabs(pt - trgMuPt[0][v])/trgMuPt[0][v]);
        htrgMudR_->Fill(deltaR(eta, phi, trgMuEta[0][v], trgMuPhi[0][v]));
      }
    }
  }
  
  return;

 
}

ULong64_t RphiNtuplizer::matchMuonTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgMuPt[f].size(); ++v)
      if (fabs(pt - trgMuPt[f][v])/trgMuPt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgMuEta[f][v], trgMuPhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

vector<bool> RphiNtuplizer::muonTriggerMap(const edm::Event &e) {

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);
  vector<bool> muTrkMap(muonHandle->size(), false);

  size_t f = 0;
  for (size_t v = 0; v < trgMuPt[f].size(); ++v) {
    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      if (fabs(iMu->pt() - trgMuPt[f][v])/trgMuPt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(iMu->eta(), iMu->phi(), trgMuEta[f][v], trgMuPhi[f][v]) < trgFilterDeltaRCut_ && muTrkMap[iMu - muonHandle->begin()] != true) {
          muTrkMap[iMu - muonHandle->begin()] = true;
          break;
      }
    }
  }
  return muTrkMap;
}


