#include <iostream>
#include <assert.h> 
#include <string> 
#include "../include/GenScaledPartLoader.hh"

using namespace baconhep;

#define MINPT 0.01
#define NPART 5


GenScaledPartLoader::GenScaledPartLoader(TTree *iTree,bool iHadrons/*=true*/) { 
  fGenInfo  = new TGenEventInfo();
  iTree->SetBranchStatus("GenEvtInfo*", 1);
  iTree->SetBranchAddress("GenEvtInfo",       &fGenInfo);
  fGenInfoBr  = iTree->GetBranch("GenEvtInfo");

  fGens  = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchStatus("GenParticle*", 1);
  iTree->SetBranchAddress("GenParticle",       &fGens);
  fGenBr  = iTree->GetBranch("GenParticle");

  fDRHeavy = 0.5;
  //  fGenVtx   = new TClonesArray("baconhep::TGenVtx");
  //  iTree->SetBranchAddress("GenVtx",       &fGenVtx);
  //  fGenVtxBr = iTree->GetBranch("GenVtx");

  int activeAreaRepeats = 1;
  double ghostArea = MINPT;
  double ghostEtaMax = 7.0;
  activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);

  jetDefCA = new fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.8);
  jetDefAK = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8);

  sd = new fastjet::contrib::SoftDrop(0,0.1,0.8);
  fastjet::contrib::OnePass_KT_Axes onepass;
  tau = new fastjet::contrib::Njettiness(onepass, fastjet::contrib::NormalizedMeasure(1., 0.8));
}

GenScaledPartLoader::~GenScaledPartLoader() { 
  delete fGenInfo;
  delete fGenInfoBr;

  delete fGens;
  delete fGenBr;

  delete activeArea; 
  delete areaDef; 

  delete jetDefCA;
  delete jetDefAK;

  delete sd; 
  delete tau;
  //delete fGenVtx;
  //delete fGenVtxBr;
}

void GenScaledPartLoader::reset() { 
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0].clear();
}

void GenScaledPartLoader::setupTree(TTree *iTree,float iXSIn,float iRadius) { 
  reset();
  fRadius = iRadius;
  fTree   = iTree;
  fLabels.push_back("weight");
  fLabels.push_back("jet_pt");
  fLabels.push_back("jet_eta");
  fLabels.push_back("jet_phi");
  fLabels.push_back("jet_mass");
  fLabels.push_back("jet_t32");
  fLabels.push_back("jet_t21");
  fLabels.push_back("jet_msd");
  fLabels.push_back("jet_t32sd");
  fLabels.push_back("jet_t21sd");
  fLabels.push_back("jet_parton_t32");
  fLabels.push_back("jet_parton_t21");
  fLabels.push_back("jet_parton_mass");
  fLabels.push_back("jet_parton_N");
  fLabels.push_back("jet_nProngs");
  fLabels.push_back("jet_prongM");
  fPartonBase = fLabels.size();
  fLabels.push_back("parton_pt");
  fLabels.push_back("parton_eta");
  fLabels.push_back("parton_phi");
  fLabels.push_back("parton_mass");
  fLabels.push_back("parton_pdgId");
  fParticleBase = fLabels.size();
  fLabels.push_back("part_pt");
  fLabels.push_back("part_eta");
  fLabels.push_back("part_phi");
  fLabels.push_back("part_pdgId");
  fLabels.push_back("part_d0");
  fLabels.push_back("part_mthrIndx");
  for(unsigned int i0 = 0; i0 < fLabels.size(); i0++) { 
    std::vector<float> lVars; 
    fVars.push_back(lVars); 
  }
 setupNtupleArr(iTree,fVars,fLabels);
 fXS = iXSIn;
}

void GenScaledPartLoader::load(int iEvent) { 
  fGens     ->Clear();
  //fGenVtx   ->Clear();
  fGenBr    ->GetEntry(iEvent);
  fGenInfoBr->GetEntry(iEvent);
  //fGenVtxBr ->GetEntry(iEvent);
  fWeight = fXS*(fGenInfo->weight);
}

void GenScaledPartLoader::chain(float iZCut) {


  // first define the potential "partons" in the event by fixing an energy scale
  unsigned nG = fGens->GetEntriesFast();
  
  // now we need to recluster the event
  std::vector<fastjet::PseudoJet> vpj;
  TLorentzVector v;
  for (unsigned iG = 0; iG != nG; ++iG) {
    TGenParticle *part = (TGenParticle*)((*fGens)[iG]);
    fprintf(stderr,"%3u %5i %-3.3f %3u\n", iG, part->pdgId, part->pt, part->parent);
    if (part->status != 1)
      continue;
    if (part->pt < MINPT)
      continue;
    if (isNeutrino(part->pdgId))
      continue;
    v.SetPtEtaPhiM(part->pt, part->eta, part->phi, part->mass);
    bool duplicate = false;
    for (auto &_pj : vpj) {
      if (fabs(_pj.perp() - part->pt) < 0.01 * part->pt && 
          deltaR2(_pj.eta(), _pj.phi(), part->eta, part->phi) < 0.01) 
      {
        fprintf(stderr,"Found duplicates with pdgid=%i, pT=%.3f, %.3f\n", part->pdgId, part->pt, _pj.perp());
        duplicate = true; 
        break;
      }
    }
    if (!duplicate) {
      vpj.emplace_back(v.Px(), v.Py(), v.Pz(), v.E());
      vpj.back().set_user_index(iG);
    }
  }

  fastjet::ClusterSequenceArea seq(vpj, *jetDefAK, *areaDef);
  std::vector<fastjet::PseudoJet> alljets(seq.inclusive_jets(200.));

  auto start = alljets.begin() + 0;
  auto end = alljets.begin() + std::min(alljets.size(), (long unsigned)2);
  std::vector<fastjet::PseudoJet> myjets(start,end); // reduce number of jets to run on 

  for (auto &jet : myjets) {
    if (jet.perp() < 200)
      continue;
    reset();
    fVars[0].push_back(fWeight);
    fVars[1].push_back(jet.perp());
    fVars[2].push_back(jet.eta());
    fVars[3].push_back(jet.phi());
    fVars[4].push_back(jet.m());

    float t3 = tau->getTau(3,jet.constituents()); 
    float t2 = tau->getTau(2,jet.constituents()); 
    float t1 = tau->getTau(1,jet.constituents()); 
    fVars[5].push_back(t2 ? t3/t2 : 1); 
    fVars[6].push_back(t1 ? t2/t1 : 1); 

    fastjet::PseudoJet sdjet = (*sd)(jet);
    t3 = tau->getTau(3,sdjet.constituents()); 
    t2 = tau->getTau(2,sdjet.constituents()); 
    t1 = tau->getTau(1,sdjet.constituents()); 
    fVars[7].push_back(sdjet.m());
    fVars[8].push_back(t2 ? t3/t2 : 1); 
    fVars[9].push_back(t1 ? t2/t1 : 1); 


    std::map<TGenParticle*, TGenParticle*> c2p;
    std::vector<TGenParticle*> jetConstituents;
    std::unordered_set<TGenParticle*> partons; 
    std::set<unsigned> partonIDs = {1, 2, 3, 4, 5, 21, 511, 513};
    for (auto &c : jet.constituents()) {
      if (c.perp() < MINPT)
        continue; // must have been a fastjet ghost
      // convert back to bacon format
      TGenParticle *bc = (TGenParticle*)((*fGens)[c.user_index()]);
      jetConstituents.push_back(bc);

      TGenParticle *parent = bc;
      TGenParticle *parton = NULL;
      fprintf(stderr,"particle pdgid=%i; parents=",bc->pdgId);
      while (parent) {
        unsigned apdgid = abs(parent->pdgId);
        fprintf(stderr,"%i, ",parent->pdgId);
        if (((apdgid < 6 || apdgid == 21) && parent->pt > 5) || // non-color singlets
            ((apdgid > 99 && apdgid < 9999) && parent->pt > 30)) {  // hadrons
          parton = parent; 
          if (partons.find(parent) == partons.end()) {
            partons.insert(parton);
          }
          break;
        }
        if (parent->parent <= 0) {
          break;
        }
        parent = (TGenParticle*)((*fGens)[parent->parent]);
      }

      c2p[bc] = parton;
      fprintf(stderr,"\n");
      if (parton)
        fprintf(stderr, "  particle pT=%.3f, parton pT=%.3f, pdgid=%i\n", bc->pt, parton->pt, parton->pdgId);
      else
        fprintf(stderr, "  particle pT=%.3f, parton IS NULL\n", bc->pt);
    }

    fprintf(stderr,"\n\n");

    std::vector<TGenParticle*> sortedPartons;
    std::map<TGenParticle*, TLorentzVector> partonKinematics;
    for (auto &iter : c2p) {
      if (partonKinematics.find(iter.second) == partonKinematics.end()) {
        if (iter.second != NULL)
          sortedPartons.push_back(iter.second);
        partonKinematics[iter.second] = TLorentzVector(0,0,0,0);
      }
      v.SetPtEtaPhiM(iter.first->pt, iter.first->eta, iter.first->phi, iter.first->mass);
      partonKinematics[iter.second] += v;
    }
    std::sort(sortedPartons.begin(),sortedPartons.end(),
              [](TGenParticle* x,TGenParticle* y) -> bool { return x->pt > y->pt; });
    std::map<TGenParticle*, int> p2i;
    for (unsigned i = 0; i!= sortedPartons.size(); ++i) {
      p2i[sortedPartons[i]] = i + 1;
    }

    unsigned lBase = fPartonBase;
    for(auto *p : sortedPartons) { 
      fVars[lBase+0].push_back(p->pt);
      fVars[lBase+1].push_back(p->eta-jet.eta());
      fVars[lBase+2].push_back(phi(p->phi,jet.phi()));
      fVars[lBase+3].push_back(p->mass);
      fVars[lBase+4].push_back(float(p->pdgId));
    }

    lBase=fParticleBase;
    for (auto *c : jetConstituents) {
      fVars[lBase+0].push_back(c->pt); 
      fVars[lBase+1].push_back(c->eta-jet.eta()); 
      fVars[lBase+2].push_back(phi(c->phi,jet.phi())); 
      fVars[lBase+3].push_back(float(simplifiedPdg(c->pdgId))); 
      fVars[lBase+4].push_back(float(c->d0)); 
      
      auto *parton = c2p[c];
      if (parton == NULL) {
        fVars[lBase+5].push_back(0.);
      } else {
        fVars[lBase+5].push_back(float(p2i[parton]));
      }
    }

    // convert partons to pseudojets and run tauN
    std::vector<fastjet::PseudoJet> vpj_partons;
    TLorentzVector vSum(0,0,0,0);
    for (auto &iter : partonKinematics) {
      TLorentzVector &vParton = iter.second;
      vSum += vParton;
      vpj_partons.emplace_back(vParton.Px(), vParton.Py(), vParton.Pz(), vParton.E());
    }
    t3 = tau->getTau(3,vpj_partons);
    t2 = tau->getTau(2,vpj_partons);
    t1 = tau->getTau(1,vpj_partons);
    fVars[10].push_back(t2 ? t3/t2 : 1); 
    fVars[11].push_back(t1 ? t2/t1 : 1); 
    fVars[12].push_back(vSum.M()); 
    fVars[13].push_back(vpj_partons.size()); 


    // count prongs 
    std::unordered_set<TGenParticle*> prongs; // avoid double-counting
    double threshold = 0.2 * jet.perp();
    for (unsigned iG = 0; iG != nG; ++iG) {
      TGenParticle *part = (TGenParticle*)((*fGens)[iG]);
      unsigned apdgid = abs(part->pdgId);
      if (apdgid > 5 &&
          apdgid != 21 &&
          apdgid != 15 &&
          apdgid != 11 &&
          apdgid != 13) 
        continue;

      if (part->pt < threshold)
        continue;

      if (partons.find(part) == partons.end()) // a prong must be a parton
        continue;

      TGenParticle *parent = part;
      TGenParticle *foundParent = NULL;
      while (parent->parent > 0) {
        parent = (TGenParticle*)((*fGens)[parent->parent]);
        if (prongs.find(parent) != prongs.end()) {
          foundParent = parent;
          break;
        }
      }

      // check if the particle has a 1->2 splitting where the daughters satisfy
      // the z-cut condition and match the jet cone
      TGenParticle *dau1 = NULL, *dau2 = NULL;
      for (unsigned jG = 0; jG != nG; ++jG) {
        TGenParticle *child = (TGenParticle*)(*fGens)[jG];
        if (child->parent != (int)iG)
          continue; // only consider splittings from the current particle 

        if (dau1)
          dau2 = child; 
        else 
          dau1 = child; 
        if (dau1 && dau2)
          break; // ! assume it's not possible to have 1->N for N>2
      }
      if (dau1 && dau1->pt > threshold && (partons.find(dau1) != partons.end()) && 
          dau2 && dau2->pt > threshold && (partons.find(dau2) != partons.end())) {
        if (foundParent) {
          prongs.erase(prongs.find(foundParent)); // remove the ancestor
        }
        prongs.insert(dau1);
        prongs.insert(dau2);
      } else if (foundParent) {
        // this means we found an ancestor parton, but this isn't the vertex that gives
        // a large 1->2 splitting, so we can skip it as an intermediary
        continue; 
      } else {
        // it passes all the checks and doesn't have an ancestor - keep it!
        prongs.insert(part);
      }

    }  
    unsigned nP = std::min((int)prongs.size(), NPART);
    fVars[14].push_back(nP);
    TLorentzVector vProngSum(0,0,0,0);
    for (auto &p : prongs) {
      v.SetPtEtaPhiM(p->pt, p->eta, p->phi, p->mass);
      vProngSum += v;
    }
    fVars[15].push_back(vProngSum.M());

    fTree->Fill();
  }
 
}

float GenScaledPartLoader::deltaR(float iPhi0,float iEta0,float iPhi1,float iEta1) { 
  float pDEta = fabs(iEta0-iEta1);
  float pDPhi = fabs(iPhi0-iPhi1); if(pDPhi > 2.*TMath::Pi()-pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
  return sqrt(pDPhi*pDPhi+pDEta*pDEta);
}


bool GenScaledPartLoader::isNeutrino(int iPdgId) { 
  return (abs(iPdgId) == 12 || abs(iPdgId) == 14 || abs(iPdgId) == 16);
}

int GenScaledPartLoader::simplifiedPdg(int iPdgId) { 
  if(abs(iPdgId) == 130)  return 2112; //kLong  (neutral)
  if(abs(iPdgId) == 310)  return 2112; //kShort (neutral)
  if(abs(iPdgId) == 2112) return 2112; //Neutron
  if(abs(iPdgId) == 3322 || abs(iPdgId) == 3212) return 2112; //Neutral Strange Baryons
  if(abs(iPdgId) == 321)  return (iPdgId/abs(iPdgId))*211; //Charged Kaon
  if(abs(iPdgId) == 2212) return (iPdgId/abs(iPdgId))*211; //Proton
  if(abs(iPdgId) > 3000 && abs(iPdgId) < 4000 && abs(iPdgId) != 3322 && abs(iPdgId) != 3212) return (iPdgId/abs(iPdgId))*211; //Strange Baryons
  return iPdgId;
}

float  GenScaledPartLoader::phi(float iPhi0,float iPhi1) { 
  double pDPhi = iPhi0-iPhi1;
  if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi) && pDPhi < 0) pDPhi =   2.*TMath::Pi()+pDPhi;
  if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi) && pDPhi > 0) pDPhi =  -2.*TMath::Pi()+pDPhi;
  return pDPhi;
}

bool GenScaledPartLoader::leptonVeto() { 
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) {   
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(fabs(pGen->pdgId) != 11 && fabs(pGen->pdgId) != 13) continue;
    if(pGen->pt < 15) continue;
    TGenParticle *pParent = (TGenParticle*)((*fGens)[pGen->parent]);
    if(fabs(pParent->pdgId) == 24 || fabs(pParent->pdgId) == 23) return true;
  }
  return false;
}
