#include <iostream>
#include <assert.h> 
#include <string> 
#include <sstream>

#include "../include/GenLoader.hh"

using namespace baconhep;

GenLoader::GenLoader(TTree *iTree, bool isPs) { 
  fGenInfo  = new TGenEventInfo();
  iTree->SetBranchAddress("GenEvtInfo",       &fGenInfo);
  fGenInfoBr  = iTree->GetBranch("GenEvtInfo");

  fGens  = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("GenParticle",       &fGens);
  fGenBr  = iTree->GetBranch("GenParticle");

  if(isPs) {
    fPSWeights  = new TClonesArray("baconhep::TPSWeight");
    iTree->SetBranchAddress("PSWeight", &fPSWeights);
    fPSWeightBr  = iTree->GetBranch("PSWeight");
  }
}
GenLoader::~GenLoader() { 
  delete fGenInfo;
  delete fGenInfoBr;

  delete fGens;
  delete fGenBr;

  delete fPSWeights;
  delete fPSWeightBr;
}
void GenLoader::reset() { 
  fBosonPt  = -1;
  fBosonPhi = -999;
  fBosonEta = -999;
  fBosonMass = -1;
  fBosonPdgId = -1;
  genEleFromW = -1;
  genMuFromW = -1;
  genTauFromW = -1;
  fgenbPt = -999;
  fgenbEta = -999;
  fgenbPhi = -999;
  fgenbMass = -999;
  fgenlPt = -999;
  fgenlEta = -999;
  fgenlPhi = -999;
  fgenlId = -999;
  fTopPt = -1;
  fAntitopPt = -1;
  fTopPtWeight = -1;
}
void GenLoader::setupTree(TTree *iTree,float iXSIn) { 
  reset();
  fTree = iTree;
  fTree->Branch("genVPt"     ,&fBosonPt   ,"fBosonPt/F");
  fTree->Branch("genVPhi"    ,&fBosonPhi  ,"fBosonPhi/F");
  fTree->Branch("genVMass"    ,&fBosonMass  ,"fBosonMass/F");
  fTree->Branch("genVEta"    ,&fBosonEta  ,"fBosonEta/F");
  fTree->Branch("genVPdfId"    ,&fBosonPdgId  ,"fBosonPdgId/I");
  fTree->Branch("genEleFromW"    ,&genEleFromW  ,"genEleFromW/I");
  fTree->Branch("genMuFromW"    ,&genMuFromW  ,"genMuFromW/I");
  fTree->Branch("genTauFromW"    ,&genTauFromW  ,"genTauFromW/I");
  fTree->Branch("genbPt"     ,&fgenbPt   ,"fgenbPt/F");
  fTree->Branch("genbPhi"    ,&fgenbPhi  ,"fgenbPhi/F");
  fTree->Branch("genbMass"    ,&fgenbMass  ,"fgenbMass/F");
  fTree->Branch("genbEta"    ,&fgenbEta  ,"fgenbEta/F");
  fTree->Branch("genlPt"     ,&fgenlPt   ,"fgenlPt/F");
  fTree->Branch("genlPhi"    ,&fgenlPhi  ,"fgenlPhi/F");
  fTree->Branch("genlId"    ,&fgenlId  ,"fgenlId/I");
  fTree->Branch("genlEta"    ,&fgenlEta  ,"fgenlEta/F");
  fTree->Branch("topPt"     ,&fTopPt   ,"fTopPt/F");
  fTree->Branch("antitopPt"     ,&fAntitopPt   ,"fAntitopPt/F");
  fTree->Branch("topPtWeight"     ,&fTopPtWeight   ,"fTopPtWeight/F");
}
void GenLoader::resetHiggs() {
  for(int i0 = 0; i0 < int(fgenHPt.size()); i0++) {
    fgenHPt[i0] = -99;
    fgenHEta[i0] = -99;
    fgenHPhi[i0] = -99;
    fgenHMass[i0] = -99;
    for(int i1 = 1; i1 < int(fgenHDauPt[i0].size()); i1++) {
      fgenHDauPt[i0][i1] = -99;
      fgenHDauEta[i0][i1] = -99;
      fgenHDauPhi[i0][i1] = -99;
      fgenHDauM[i0][i1] = -99;
      fgenHDauId[i0][i1] = -99;
      fgenHDauDecay[i0][i1] = -99;
    }
  }
}
void GenLoader::setupTreeHiggs(TTree *iTree) {
  resetHiggs();
  fTree = iTree;
  fgenHPt.clear(); fgenHEta.clear(); fgenHPhi.clear(); fgenHMass.clear(); 
  fgenHDauPt.clear(); fgenHDauEta.clear(); fgenHDauPhi.clear(); fgenHDauM.clear(); fgenHDauId.clear(); fgenHDauDecay.clear();
  for(int i0 = 0; i0 < 2; i0++) {
    fgenHPt.push_back(-999);
    fgenHEta.push_back(-999);
    fgenHPhi.push_back(-999);
    fgenHMass.push_back(-999);
    std::vector<float> lPt,lEta,lPhi,lM,lId,lDecay;
    for(int i1 = 0; i1 < 4; i1++) {
      lPt.push_back(-999);
      lEta.push_back(-999);
      lPhi.push_back(-999);
      lM.push_back(-999);
      lId.push_back(-999);
      lDecay.push_back(-999);
    }
    fgenHDauPt.push_back(lPt);
    fgenHDauEta.push_back(lEta);
    fgenHDauPhi.push_back(lPhi);
    fgenHDauM.push_back(lM);
    fgenHDauId.push_back(lId);
    fgenHDauDecay.push_back(lDecay);
  }
  for(int i0 = 0; i0 <  int(fgenHPt.size()); i0++) {
    std::stringstream pSPt;    pSPt << "genHPt" << i0;
    std::stringstream pSEta;   pSEta << "genHEta" << i0;
    std::stringstream pSPhi;   pSPhi << "genHPhi" << i0;
    std::stringstream pSMass;  pSMass << "genHMass" << i0;
    fTree->Branch(pSPt.str().c_str()   ,&fgenHPt[i0]   ,(pSPt.str()+"/F").c_str());
    fTree->Branch(pSEta.str().c_str()  ,&fgenHEta[i0]  ,(pSEta.str()+"/F").c_str());
    fTree->Branch(pSMass.str().c_str() ,&fgenHMass[i0] ,(pSMass.str()+"/F").c_str());
    fTree->Branch(pSPhi.str().c_str()  ,&fgenHPhi[i0]  ,(pSPhi.str()+"/F").c_str());
    for(int i1 = 0; i1 < int(fgenHDauPt[i0].size()); i1++) {
      std::stringstream pSPt;   pSPt << "genHDau" << i0 << "Pt" << i1;
      std::stringstream pSEta;  pSEta << "genHDau" << i0 << "Eta" << i1;
      std::stringstream pSPhi;  pSPhi << "genHDau" << i0 << "Phi" << i1;
      std::stringstream pSM;    pSM << "genHDau" << i0 << "Mass" << i1;
      std::stringstream pSId;   pSId << "genHDau" << i0 << "Id" << i1;
      std::stringstream pSDecay;   pSDecay << "genHDau" << i0 << "Decay" << i1;
      fTree->Branch(pSPt.str().c_str()   ,&fgenHDauPt[i0][i1]   ,(pSPt.str()+"/F").c_str());
      fTree->Branch(pSEta.str().c_str()  ,&fgenHDauEta[i0][i1]  ,(pSEta.str()+"/F").c_str());
      fTree->Branch(pSPhi.str().c_str()  ,&fgenHDauPhi[i0][i1]  ,(pSPhi.str()+"/F").c_str());
      fTree->Branch(pSM.str().c_str()    ,&fgenHDauM[i0][i1]    ,(pSM.str()+"/F").c_str());
      fTree->Branch(pSId.str().c_str()   ,&fgenHDauId[i0][i1]   ,(pSId.str()+"/F").c_str());
      fTree->Branch(pSDecay.str().c_str()  ,&fgenHDauDecay[i0][i1]     ,(pSDecay.str()+"/F").c_str());
    }
  }
}
void GenLoader::load(int iEvent) { 
  reset();
  fGens     ->Clear();
  fGenBr    ->GetEntry(iEvent);
  fGenInfoBr->GetEntry(iEvent);
  fWeight = fGenInfo->weight;
}

TGenParticle* GenLoader::getParticle(int i) {
  if ( i >= 0 && i < fGens->GetEntriesFast() ) {
    return (TGenParticle*) (*fGens)[i];
  }
  return nullptr;
}

bool GenLoader::isType(std::string boson,std::string mode)
{
  int iPDGID,iId;
  if (boson.find("Z")==0) iPDGID = 23;
  if (boson.find("W")==0) iPDGID = 24;
  if (boson.find("H")==0) iPDGID = 25;
  if (boson.find("Zprime")==0) iPDGID = 10031;
  if (boson.find("DMSpin0")==0) iPDGID = 9900032;  
  if (mode.find("bb")==0) iId = 5;
  if (mode.find("cc")==0) iId = 4;

  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *pGen = getParticle(i0);
    if (mode.find("bb")==0 || mode.find("cc")==0) {
      if(abs(pGen->pdgId)==iId) {
	if(pGen->parent<0) continue;
	TGenParticle *parent = getParticle(pGen->parent);
	if(abs(parent->pdgId)==iPDGID) return true;
      }
    }
    if (mode.find("cs")==0) {
      if(abs(pGen->pdgId)==4 || abs(pGen->pdgId)==3) {
	if(pGen->parent<0) continue;
	TGenParticle *parent = getParticle(pGen->parent);
	if(abs(parent->pdgId)==iPDGID) return true;
      }
    }
  }
  
  return false;
}

// utils
bool GenLoader::hard(int &iP)
{
  TGenParticle *p = getParticle(iP);
  if(p->status>20 && p->status<30) return true;
  return false;
}
bool GenLoader::hasChild(int &iparent, bool beHard)
{
  TGenParticle *parent = getParticle(iparent);
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *child = getParticle(i0);
    if (child->pdgId != parent->pdgId)
      continue;
    if (beHard && !hard(i0))
      continue;
    if (child->parent !=-2 &&
        child->parent == iparent) {
      return true;
    }
  }
  return false;
}
TGenParticle* GenLoader::findDaughter(int iparent, int dauId)
{
  for(int k = iparent+1; k < fGens->GetEntriesFast(); k++) {
    TGenParticle *genp = getParticle(k);
    if(genp->parent == iparent) {
      if(abs(genp->pdgId) == dauId) {
        return genp;
      }
    }
  }
  return 0;
}
int GenLoader::findDaughterId(int iparent, int dauId)
{
  for(int k = iparent+1; k < fGens->GetEntriesFast(); k++) {
    const baconhep::TGenParticle *genp = getParticle(k);;
    if(genp->parent == iparent) {
      if(abs(genp->pdgId) == dauId) {
        return k;
      }
    }
  }
  return -1;
}
int GenLoader::findLastParent(int iparent,int iId){
  Bool_t foundLast = kFALSE;
  int iLast = iparent;
  while (!foundLast) {
    int tmpId = findDaughterId(iLast,iId);
    if (tmpId>=0) iLast = tmpId;
    else foundLast = kTRUE;
  }
  return iLast;
}
void GenLoader::findBoson(int iId, int lOption){
  reset();
  float pbosonPt(-1),pbosonPhi(-999),pbosonMass(-1),pbosonEta(-999);
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *genp0 = getParticle(i0);

    // find highest Pt boson G(22)                                                                                                                                                                                                                            
    if(lOption == 0){
      if(fabs(genp0->pdgId)==iId && genp0->pt > pbosonPt){
        pbosonPt = genp0->pt;
        pbosonPhi = genp0->phi;
        pbosonEta = genp0->eta;
        pbosonMass = genp0->mass;
      }
    }

    // find last boson Z(23),W(24),Z'(10031),H(25),DMSpin(9900032)                                                                                                                                                                                            
    if(lOption == 1){
      if(fabs(genp0->pdgId)==iId){
        int iL0 = findLastParent(i0,iId);
        for(int k0 = 0; k0 < fGens->GetEntriesFast(); k0++) {
          TGenParticle *genp1 = getParticle(k0);
          if(k0==iL0){
            pbosonPt = genp1->pt;
            pbosonPhi = genp1->phi;
            pbosonEta = genp1->eta;
            pbosonMass = genp1->mass;
            break;
          }
        }
      }
    }

    // find W(24) for ttbar semilep                                                                                                                                                                                                                           
    if(lOption == 2){
      if(fabs(genp0->pdgId)==iId) {
        TGenParticle *dau1 = findDaughter(i0, 11);
        TGenParticle *dau2 = findDaughter(i0, 13);
        if(dau1 || dau2){
          pbosonPt = genp0->pt;
          pbosonPhi = genp0->phi;
          pbosonEta = genp0->eta;
          pbosonMass = genp0->mass;
        }
      }
    }

    // find W for ttbar dileptonic (6)                                                                                                                                                                                                                        
    if(lOption == 3){
      if(fabs(genp0->pdgId)==iId) {
        int iW0 = findLastParent(i0,24);
        for(int k0 = 0; k0 < fGens->GetEntriesFast(); k0++) {
          TGenParticle *dau0 = getParticle(k0);
          TGenParticle *ele0 = findDaughter(iW0, 11); TGenParticle *muo0 = findDaughter(iW0, 13);
          if(k0==iW0 && (ele0 || muo0)){
            for(int i1=iW0+1; i0 < fGens->GetEntriesFast(); i1++) {
              TGenParticle *genp1 = getParticle(i1);
              if(genp1->pdgId == iId) {
                int iW1 = findLastParent(i1,24);
                for(int k1 = 0; k1 < fGens->GetEntriesFast(); k1++) {
                  TGenParticle *dau1 = getParticle(k1);
                  TGenParticle *ele1 = findDaughter(iW1, 11); TGenParticle *muo1 = findDaughter(iW1, 13);
                  if(k1==iW1 && (ele1 || muo1)){
                    TLorentzVector vDau0; vDau0.SetPtEtaPhiM(dau0->pt, dau0->eta, dau0->phi, dau0->mass);
                    TLorentzVector vDau1; vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
                    pbosonPt = (vDau0 + vDau1).Pt();
                    pbosonPhi = (vDau0 + vDau1).Phi();
                  }
                }
              }
            }
          }
        }
      }
    }

    // find hadronic W(24) for ttbar semilep                                                                                                                                                                                                                  
    if(lOption == 4){
      if(fabs(genp0->pdgId)==iId) {
        TGenParticle *dau1 = findDaughter(i0, 2);
        TGenParticle *dau2 = findDaughter(i0, 3);
        if(dau1 || dau2){
          pbosonPt = genp0->pt;
          pbosonPhi = genp0->phi;
          pbosonEta = genp0->eta;
          pbosonMass = genp0->mass;
        }
      }
    }


  }
  fBosonPt = pbosonPt;
  fBosonPhi = pbosonPhi;
  fBosonMass = pbosonMass;
  fBosonEta = pbosonEta;
  fBosonPdgId = iId;
}

// is W(qq) in Top->Wb?                                                                                                                                                                                                                                       
// Return 1 if W(light)
// Return 2 if W(cs)
// If b quark inside dR, add 8
//  i.e. (x&8)==8 lets you know if true
int GenLoader::getHadronicWInTopFlavor(TGenParticle *genp,int iW,TLorentzVector jet,double dR,double &wMatching, double &wSize)
{
  TLorentzVector vW,vDau1,vDau2,b;
  TGenParticle *dau1{nullptr};
  TGenParticle *dau2{nullptr};
  TGenParticle *glep{nullptr};

  vW.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
  int iWlast = findLastParent(iW, 24);

  int iQ=0, jQ=0;
  for (; iQ<fGens->GetEntriesFast(); ++iQ) {
    dau1 = getParticle(iQ);
    if( dau1->parent==iWlast && (abs(dau1->pdgId)<6 || abs(dau1->pdgId)<15)) {
      vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
      wMatching = jet.DeltaR(vDau1);
      wSize     = vW.DeltaR(vDau1);
      break; // found the first quark                                                                                                                                            
    }
  }
  for (jQ=iQ+1; jQ<fGens->GetEntriesFast(); ++jQ) {
    dau2 = getParticle(jQ);
    if(dau2->parent==iWlast && (abs(dau2->pdgId)<6 || abs(dau1->pdgId)<15)) {
      vDau2.SetPtEtaPhiM(dau2->pt, dau2->eta, dau2->phi, dau2->mass);
      wMatching = TMath::Max(wMatching,jet.DeltaR(vDau2));
      wSize     = TMath::Max(wSize,vW.DeltaR(vDau2));
      break;
    }
  }
  
  // Missing gen particles from the list?
  if ( dau1 == nullptr || dau2 == nullptr ) return -2;

  int wType = 0;
  if ( std::abs(dau1->pdgId) <= 3 && std::abs(dau2->pdgId) <= 3 ) wType = 1;
  else if ( std::abs(dau1->pdgId) == 4 || std::abs(dau2->pdgId) == 4 ) wType = 2;
  else if ( std::abs(dau1->pdgId) == 11 || std::abs(dau2->pdgId) == 11 ) wType = 3;
  else if ( std::abs(dau1->pdgId) == 13 || std::abs(dau2->pdgId) == 13 ) wType = 4;

  // save gen info
  if ( std::abs(dau1->pdgId) == 11 || std::abs(dau1->pdgId) == 13 ) glep = dau1;
  if ( std::abs(dau2->pdgId) == 11 || std::abs(dau2->pdgId) == 13 ) glep = dau2;
  if ( glep != nullptr ) {
    fgenlPt = glep->pt; fgenlPhi = glep->phi; fgenlId = glep->pdgId; fgenlEta =glep->eta;
  }

  // Check if b is in jet cone
  int iTop = genp->parent;
  while ( std::abs(getParticle(iTop)->pdgId) == 24 ) {
    iTop = getParticle(iTop)->parent;
  }
  if ( std::abs(getParticle(iTop)->pdgId) == 6 || std::abs(getParticle(iTop)->pdgId) == 624 ) {
    TGenParticle* genB = nullptr;
    for(int iB=0; iB<fGens->GetEntriesFast(); ++iB) {
      TGenParticle* p = getParticle(iB);
      if ( p->parent == iTop && std::abs(p->pdgId) == 5 ) {
        genB = p;
        break;
      }
    }
    if ( genB == nullptr ) return -1;
    if ( glep != nullptr ) { //save b only for lep W
      fgenbPt = genB->pt;fgenbPhi = genB->phi; fgenbMass =genB->mass; fgenbEta =genB->eta;
    }
    b.SetPtEtaPhiM(genB->pt, genB->eta, genB->phi, genB->mass);
    if ( b.DeltaR(jet) < dR ) wType += 8;
  }
  // One could in principle check the other W in the Z' case as well
  return wType;
}


// jet matched to V(qq) with quark flavor
// 0: error
// 1: light quark, e.g. W(ud), Z(ss)
// 2: c quark, e.g. W(cs), Z(cc), H(cc)
// 3: b quark, e.g. Z(bb), H(bb)
int GenLoader::getHadronicVflav(TGenParticle *genp,int j,int iId, TLorentzVector jet,double dR,double &vMatching, double &vSize)
{
  TLorentzVector vV,vDau1,vDau2;
  TGenParticle *dau1{nullptr}, *dau2{nullptr};

  vV.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
  int iV = findLastParent(j, iId);

  int iQ;
  for (iQ=0; iQ<fGens->GetEntriesFast(); ++iQ) {
    dau1 = getParticle(iQ);
    if(dau1->parent==iV) {
      vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
      vMatching = jet.DeltaR(vDau1);
      vSize     = vV.DeltaR(vDau1);
      break;
    }
  }
  for (int jQ=iQ+1; jQ<fGens->GetEntriesFast(); ++jQ) {
    dau2 = getParticle(jQ);
    if(dau2->parent==iV) {
      vDau2.SetPtEtaPhiM(dau2->pt, dau2->eta, dau2->phi, dau2->mass);
      vMatching = TMath::Max(vMatching,jet.DeltaR(vDau2));
      vSize     = TMath::Max(vSize,vV.DeltaR(vDau2));
      break;
    }
  }

  // Missing gen particles from the list?
  if ( dau1 == nullptr || dau2 == nullptr ) return 0;

  if ( std::abs(dau1->pdgId) <= 3 && std::abs(dau2->pdgId) <= 3 ) return 1;
  else if ( std::abs(dau1->pdgId) == 4 || std::abs(dau2->pdgId) == 4 ) return 2;
  else if ( std::abs(dau1->pdgId) == 5 || std::abs(dau2->pdgId) == 5 ) return 3;
  // W'(tb) ?
  else return 0;
}

// jet matched to something,
// For tops (6, or 624 in 2016 sample?) result = flavor of matched W
//   if (result&8)==8, then the b quark is within dR of the jet
// Z' to WW (9000001) result = flavor of matched W
// For bosons(24,23,10031,25,55) result = flavor:
// 1: light quark, e.g. W(ud), Z(ss)
// 2: c quark, e.g. W(cs), Z(cc), H(cc)
// 3: b quark, e.g. Z(bb), H(bb)
// In all cases, 0 = unmatched if matching or size are -999, otherwise an error occurred
int GenLoader::ismatchedJet(TLorentzVector jet0, double dR,double &matching, double &size, int iId){
  int result=0;
  matching = -999;
  size = -999;
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *genp0 = getParticle(i0);
    TLorentzVector mcMom; mcMom.SetPtEtaPhiM(genp0->pt,genp0->eta,genp0->phi,genp0->mass);
    if (mcMom.DeltaR(jet0) < dR) {
      if( (iId == 624 || iId == 6 || iId == 9000001) && std::abs(genp0->pdgId)==24 ) {
        result = getHadronicWInTopFlavor(genp0,i0,jet0,dR,matching,size);
	break;
      }
      if( iId == std::abs(genp0->pdgId) && (iId == 24 || iId == 23 || iId == 10031 || iId == 25 || iId == 55) ){
        result = getHadronicVflav(genp0,i0,iId,jet0,dR,matching,size);
        break;
      }
    }
  }
  return result;
}

int GenLoader::ismatchedSubJet(TLorentzVector subjet0){
  int lOption =0;
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *genp0 = getParticle(i0);
    TLorentzVector vq;
    if(abs(genp0->pdgId)==5) {
      vq.SetPtEtaPhiM(genp0->pt, genp0->eta, genp0->phi, 5);
      if(vq.DeltaR(subjet0) < 0.4) lOption = 1; //isB                                                                                                                                                                                                         
    }
    else if(abs(genp0->pdgId)==4) {
      vq.SetPtEtaPhiM(genp0->pt, genp0->eta, genp0->phi, 1.29);
      if(vq.DeltaR(subjet0) < 0.4) lOption = 2; //isC                                                                                                                                                                                                         
    }
    else lOption = 3; //isLF                                                                                                                                                                                                                                  
  }
  return lOption;
}

int GenLoader::isHWWsemilepBoson(int iH,int iId,int iIdDau,float &genSize){
  TGenParticle *genH = getParticle(iH);
  TLorentzVector vH,vDau1,vDau2,vDau3,vDau4;
  std::vector<TGenParticle*> lDaus;
  std::vector<int> lDausIndex;
  if(abs(genH->pdgId)==iId) {
    vH.SetPtEtaPhiM(genH->pt, genH->eta, genH->phi, genH->mass);
    for (int iD=0; iD<fGens->GetEntriesFast(); ++iD) {
      TGenParticle *dau = getParticle(iD);
      if(dau->parent==iH && abs(dau->pdgId)!=23) {
	lDaus.push_back(dau); lDausIndex.push_back(iD);
      }
    }
    //std::cout << "dau size "<< lDaus.size() << std::endl;
    unsigned int lMax =4;
    if(lDaus.size()<lMax) lMax = lDaus.size();
    std::vector<TGenParticle*> lLeps,lQuarks;
    // 2 daughters
    if(lMax == 2){
      if( abs(lDaus[0]->pdgId)==iIdDau && abs(lDaus[1]->pdgId)==iIdDau) {
	for (int i0=0; i0<fGens->GetEntriesFast(); ++i0) {
	  TGenParticle *pGen = getParticle(i0);
	  if(pGen->parent==lDausIndex[0] || pGen->parent==lDausIndex[1]) {
	    if ( abs(pGen->pdgId)>=11 && abs(pGen->pdgId)<=14 ) lLeps.push_back(pGen);
            if ( abs(pGen->pdgId)>=14 && abs(pGen->pdgId)<=5 ) lQuarks.push_back(pGen);
	  }
	}
      }
      else return 0;
    }
    // 3 daughters
    if(lMax == 3){
      for(unsigned int i0 = 0; i0 < lMax; i0++){
	if( abs(lDaus[i0]->pdgId)==iIdDau ) {
	  for (int i1=0; i1<fGens->GetEntriesFast(); ++i1) {
	    TGenParticle *pGen = getParticle(i1);
	    if ( pGen->parent==lDausIndex[i0] ){
	      if ( abs(pGen->pdgId)>=11 && abs(pGen->pdgId)<=14 ) lLeps.push_back(pGen);
	      if ( abs(pGen->pdgId)>=14 && abs(pGen->pdgId)<=5 ) lQuarks.push_back(pGen);
	    }
          }
        }
	else if ( abs(lDaus[i0]->pdgId)>=11 && abs(lDaus[i0]->pdgId)<=14 ) lLeps.push_back(lDaus[i0]);
	else if ( abs(lDaus[i0]->pdgId)>=14 && abs(lDaus[i0]->pdgId)<=5 ) lQuarks.push_back(lDaus[i0]);
	else {
	  std::cout << "3 daus return " << abs(lDaus[i0]->pdgId) << std::endl; 
	  return 0;
	}
      }
    }
    // 4 daugthers
    if(lMax == 4){
      for(unsigned int i0 = 0; i0 < lMax; i0++){
	if ( abs(lDaus[i0]->pdgId)>=11 && abs(lDaus[i0]->pdgId)<=14 ) lLeps.push_back(lDaus[i0]);
	else if ( abs(lDaus[i0]->pdgId)>=14 && abs(lDaus[i0]->pdgId)<=5 ) lQuarks.push_back(lDaus[i0]);
	else {
	  std::cout << "4 daus return "<< abs(lDaus[i0]->pdgId) << std::endl;
          return 0;
        }
      }
    }
    // WW->lvqq                                                                                                                                           
    std::cout << "leps size " << lLeps.size() << " quarks size " << lQuarks.size() << std::endl;
    if (lLeps.size() == 2 && lQuarks.size() == 2){
      vDau1.SetPtEtaPhiM(lLeps[0]->pt, lLeps[0]->eta, lLeps[0]->phi, lLeps[0]->mass);
      vDau2.SetPtEtaPhiM(lLeps[1]->pt, lLeps[1]->eta, lLeps[1]->phi, lLeps[1]->mass);
      vDau3.SetPtEtaPhiM(lQuarks[0]->pt, lQuarks[0]->eta, lQuarks[0]->phi, lQuarks[0]->mass);
      vDau4.SetPtEtaPhiM(lQuarks[1]->pt, lQuarks[1]->eta, lQuarks[1]->phi, lQuarks[1]->mass);
      genSize = TMath::Max(TMath::Max(TMath::Max(vH.DeltaR(vDau1),vH.DeltaR(vDau2)),vH.DeltaR(vDau3)),vH.DeltaR(vDau4));
      return 1;
    }
  }
  return 0;
}
int GenLoader::isHDau(int iId, int iDauId, TLorentzVector jet, int iHiggs){
  int iH(-1),  iR(-1);
  std::vector<TLorentzVector> lParts;
  std::vector<int> lHiggs;
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *pGen = getParticle(i0);
    if(fabs(pGen->pdgId == iId)) {
      lHiggs.push_back(findLastParent(i0,iId));
    }
  }
  unsigned int lMax =2;
  float lMaxR = 0.8;
  if(lHiggs.size()<lMax) lMax = lHiggs.size();
  for(unsigned int i0 = 0; i0 < lMax; i0++)
    {
      int iD(-1);
      iH = lHiggs[i0];
      TGenParticle *pH = getParticle(iH);
      TLorentzVector ivH; ivH.SetPtEtaPhiM(pH->pt, pH->eta, pH->phi, pH->mass);
      fgenHPt[i0] = pH->pt;
      fgenHEta[i0] = pH->eta;
      fgenHPhi[i0] = pH->phi;
      fgenHMass[i0] = pH->mass;
      for(int i1 = 0; i1 < fGens->GetEntriesFast(); i1++) {
        TGenParticle *pGen = getParticle(i1);
        if(pGen->parent == iH) {
          //if(fabs(pGen->pdgId) == iDauId) {
	  iD += 1;
	  if(iD > 3) break;
	  int iP = findLastParent(i1,iDauId);
	  TGenParticle *pP = getParticle(iP);
	  fgenHDauPt[i0][iD] = pP->pt;
	  fgenHDauEta[i0][iD] = pP->eta;
	  fgenHDauPhi[i0][iD] = pP->phi;
	  fgenHDauM[i0][iD] = pP->mass;
	  fgenHDauId[i0][iD] = pP->pdgId;
	  for(int i2 = 0; i2 < fGens->GetEntriesFast(); i2++) {
	    TGenParticle *pPart = getParticle(i2);
	    TLorentzVector vPart; vPart.SetPtEtaPhiM(pPart->pt, pPart->eta, pPart->phi, pPart->mass);
	    if(fabs(pPart->parent) ==iP) {
	      fgenHDauDecay[i0][iD] = pPart->pdgId;
	      break;
	    }
	  }
	  if(fabs(pGen->pdgId) == iDauId) {
	    TLorentzVector ivP; ivP.SetPtEtaPhiM(pP->pt, pP->eta, pP->phi, pP->mass);
	    lParts.push_back(ivP);
	  }
	}
      }

      if(jet.DeltaR(ivH)<lMaxR && iD > -1) {
	lMaxR = jet.DeltaR(ivH);
        for(unsigned int i1 = 0; i1<lParts.size(); i1++) {
          if(jet.DeltaR(lParts[i1])<lMaxR){ iR = i0;}
        }
      }
    }
  return iR;
}
int GenLoader::isttbarType(int lepPdgId)
{
  assert(fGens);
  int nlep=0;
  for(int i0=0; i0<fGens->GetEntriesFast(); i0++) {
    TGenParticle *pGen = getParticle(i0);
    if(abs(pGen->pdgId)==abs(lepPdgId)) {
      if(pGen->parent<0) continue;
      TGenParticle *lparent = getParticle(pGen->parent);
      if(abs(lparent->pdgId)==24) nlep++;
    }
  }
  return nlep;
}

void GenLoader::saveTTbarType() {
  genEleFromW = isttbarType(11);
  genMuFromW = isttbarType(13);
  genTauFromW = isttbarType(15);
}

// tt mc correction: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#MC_SFs_Reweighting     
float GenLoader::computeTTbarCorr() {
  const int TOP_PDGID = 6;
  double pt1=0, pt2=0;
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *p = getParticle(i0);
    if(p->pdgId ==  TOP_PDGID) {
      pt1 = p->pt;
      fTopPt = pt1;
    }
    if(p->pdgId == -TOP_PDGID) {
      pt2 = p->pt;
      fAntitopPt = pt2;
    }
  }

  double w1 = exp(0.0615 - 0.0005*pt1);
  double w2 = exp(0.0615 - 0.0005*pt2);
  fTopPtWeight = sqrt(w1*w2);
  return sqrt(w1*w2);
}
void GenLoader::setPSWeights(TTree *iTree) {
  fTree = iTree;
  fgenPSWeight.clear();
  for(int i1 = 0; i1 < fNWeights; i1++ ) {
    fgenPSWeight.push_back(-99);
  }
  for(int i0 = 0; i0 < fNWeights; i0++) {
    std::stringstream pSPSWeight; pSPSWeight << "psWeight" << i0;
    fTree->Branch(pSPSWeight.str().c_str()   ,&fgenPSWeight[i0]   ,(pSPSWeight.str()+"/F").c_str());
  }
}
void GenLoader::loadPSWeights(int iEvent) {
  fPSWeights    ->Clear();
  fPSWeightBr   ->GetEntry(iEvent);
}
void GenLoader::fillPSWeights(){
  for(int i1 = 0; i1 < fNWeights; i1++) {
    TPSWeight *psweight = (TPSWeight*)((*fPSWeights)[i1]);
    fgenPSWeight[i1] = psweight->weight;
  }
}
