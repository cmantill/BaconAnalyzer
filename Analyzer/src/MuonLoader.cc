#include "../include/MuonLoader.hh"

#include "TFile.h"
#include <cmath>
#include <iostream> 
#include <sstream> 

using namespace baconhep;

MuonLoader::MuonLoader(TTree *iTree) { 
  fMuons  = new TClonesArray("baconhep::TMuon");
  iTree->SetBranchAddress("Muon",       &fMuons);
  fMuonBr  = iTree->GetBranch("Muon");
  fN = 3;

  for(int i0 = 0; i0 < fN*3.; i0++) {double pVar = 0; fVars.push_back(pVar);}
  for(int i0 = 0; i0 <     4; i0++) {double pVar = 0; fVars.push_back(pVar);}
}
MuonLoader::~MuonLoader() { 
  delete fMuons;
  delete fMuonBr;
}
void MuonLoader::reset() { 
  fNMuonsLoose = 0; 
  fNMuonsTight = 0;
  fismu0Tight  = 0; 
  fismu1Tight  = 0;
  fLooseMuons.clear();
  fMediumMuons.clear();
  fTightMuons.clear();
  fHighPtMuons.clear();
  looseMuons.clear();
  mediumMuons.clear();
  tightMuons.clear();
  highptMuons.clear();
  for(int i0 = 0; i0 < int(fmuId.size()); i0++){ fmuId[i0] = -1; fmuSel[i0] = -1;}
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
}
void MuonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("nmuLoose"  ,&fNMuonsLoose ,"fNMuonsLoose/I"); 
  fTree->Branch("nmuMedium" ,&fNMuonsMedium,"fNMuonsMedium/I");
  fTree->Branch("nmuTight"  ,&fNMuonsTight ,"fNMuonsTight/I");
  fTree->Branch("nmuHighPt" ,&fNMuonsHighPt,"fNMuonsHighPt/I");
  for(int i0 = 0; i0 < fN; i0++) {
    fmuId.push_back(-1);
    fmuSel.push_back(-1);
  }
  for(int i0 = 0; i0 < int(fmuId.size()); i0++) {
    std::stringstream muId; muId << "vmuoid_" << i0;
    std::stringstream muSel; muSel << "vmuosel_" << i0;
    fTree->Branch(muId.str().c_str(),  &fmuId[i0], (muId.str()+"/I").c_str());
    fTree->Branch(muSel.str().c_str(),  &fmuSel[i0], (muSel.str()+"/I").c_str());
  }
  setupNtuple("vmuoLoose"   ,iTree,fN,fVars);       // add leading 2 muons: pt,eta,phi,mass (2*4=8)
}
void MuonLoader::load(int iEvent) { 
  fMuons  ->Clear();
  fMuonBr ->GetEntry(iEvent);
}
void MuonLoader::selectMuons(std::vector<TLorentzVector> &iMuons, float met, float metPhi) {
  reset(); 
  int lCount = 0,lMCount = 0, lHPCount =0,lTCount =0; 
  fvMetNoMu.SetMagPhi(met,metPhi);

  // Muon selection
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) { 
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(pMuon->pt        <=  10)                      continue;
    if(fabs(pMuon->eta) >=  2.4)                     continue;
    if(!passMuonLooseId(pMuon))                     continue;
    // if(!passMuonLooseSel(pMuon))                     continue;
    lCount++;

    if(!passMuonMediumSel(pMuon))                   lMCount++;
    if(!passMuonHighPtSel(pMuon))                   lHPCount++;

    TVector2 vMu; vMu.SetMagPhi(pMuon->pt, pMuon->phi);
    fvMetNoMu = fvMetNoMu + vMu;

    addMuon(pMuon,fLooseMuons);
    if(pMuon->pt>20 && fabs(pMuon->eta)< 2.4 && passMuonTightSel(pMuon)){
      if(lCount==1) fismu0Tight = 1;
      if(lCount==2) fismu1Tight = 1;
      lTCount++;
    }
    addMuon(pMuon,fTightMuons);

  }
  addVMuon(fLooseMuons,looseMuons,MUON_MASS);
  addVMuon(fMediumMuons,mediumMuons,MUON_MASS);
  addVMuon(fLooseMuons,tightMuons,MUON_MASS);
  addVMuon(fHighPtMuons,highptMuons,MUON_MASS);

  fNMuonsLoose = lCount;
  fNMuonsMedium = lMCount; 
  fNMuonsTight = lTCount;
  fNMuonsHighPt = lHPCount;

  // Cleaning iMuons
  if(fTightMuons.size()>1){
    iMuons.push_back(tightMuons[0]); // save first tight muon
    iMuons.push_back(looseMuons[1]); // if leading lepton is tight, save the subleading loose one
  }

  // Muon id and sel:
  // 1: loose
  // 2: medium
  // 3: tight
  // 4: soft
  // 5: high-pt
  for(int i0 = 0; i0 < int(fLooseMuons.size()); i0++){
    TMuon *pMuon = (TMuon*)(fLooseMuons[i0]);
    if(passMuonLooseId(pMuon)) fmuId[i0] = 1; 
    if(passMuonMediumId(pMuon)) fmuId[i0] = 2;
    if(passMuonTightId(pMuon)) fmuId[i0] = 3;
    if(passMuonSoftId(pMuon)) fmuId[i0] = 4;
    if(passMuonHighPtId(pMuon)) fmuId[i0] = 5;
    if(passMuonLooseSel(pMuon)) fmuSel[i0] = 1;
    if(passMuonMediumSel(pMuon)) fmuSel[i0] = 2;
    if(passMuonTightSel(pMuon)) fmuSel[i0] = 3;
    if(passMuonSoftSel(pMuon)) fmuSel[i0] = 4;
    if(passMuonHighPtSel(pMuon)) fmuSel[i0] = 5;
  }

  if(fVars.size() > 0) fillMuon(fN,fLooseMuons,fVars);
}
