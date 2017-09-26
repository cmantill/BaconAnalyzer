#ifndef GenScaledPartLoader_H
#define GenScaledPartLoader_H
#include "Utils.hh"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/MeasureDefinition.hh"

using namespace baconhep;

class GenScaledPartLoader { 
public:
  GenScaledPartLoader(TTree *iTree,bool iHadrons=true);
  ~GenScaledPartLoader();
  void reset();
  void setupTree(TTree *iTree,float iXSIn,float iRadius);
  void load (int iEvent);
  void chain(float iZCut);
  float deltaR(float iPhi0,float iEta0,float iPhi1,float iEta1);
  bool isNeutrino(int iPdgId);
  float  phi(float iPhi0,float iPhi1);
  int  simplifiedPdg(int iPdgId);
  bool leptonVeto();

  TClonesArray  *fGens;
  TBranch       *fGenBr;
  TGenEventInfo *fGenInfo;
  TBranch       *fGenInfoBr;

  float fWeight;
  std::vector<std::string> fLabels;
  std::vector<std::vector<float> > fVars;

  fastjet::AreaDefinition *areaDef=0;
  fastjet::GhostedAreaSpec *activeArea=0;
  fastjet::JetDefinition *jetDefCA=0;
  fastjet::JetDefinition *jetDefAK=0;
  fastjet::contrib::SoftDrop *sd=0;
  fastjet::contrib::Njettiness *tau=0;

protected: 
  TTree         *fTree;
  int fPartonBase;
  int fParticleBase;
  float fXS;
  float fRadius;
  float fDRHeavy;
};
#endif
