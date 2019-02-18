// run:
// ./runPFJetsToJUNIPR input_directory input_file output_directory recluster_def

#include "../include/PerJetLoader.hh"
#include "../include/fastjet_JUNIPR_utilities.hh"

#include "fastjet/ClusterSequence.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace fastjet;
using namespace std;

#define PI 3.141592654

PerJetLoader      *fVJet8      = 0;

TTree* load(std::string iName) {
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}

double SignedDeltaPhi(double phi1, double phi2) {
  double dPhi = phi1-phi2;
  if (dPhi<-PI)
    dPhi = 2*PI+dPhi;
  else if (dPhi>PI)
    dPhi = -2*PI+dPhi;
  return dPhi;
}

int main(int argc, char *argv[]) {

  if (argc != 5) {    
    cout << "Wrong number of arguments" << endl;
    return 0;
  }
  
  // User input
  string input_file = argv[1];
  string output_file = argv[2];
  double recluster_def = atof(argv[3]);
  string isdata = argv[4];

  bool isData;
  if(isdata.compare("data")!=0) isData = false;
  else isData = true;

  // Input, output
  stringstream s_infile;
  s_infile << input_file;

  stringstream s_outfile;
  s_outfile << "JUNIPR_format_" << output_file ;
  ofstream outfile;
  outfile.open( s_outfile.str().c_str());
  
  // Read Tree
  const std::string lName = s_infile.str().c_str();
  TTree *lTree = load( lName );
  fVJet8     = new PerJetLoader  (lTree,"AK8Puppi","AddAK8Puppi","AK8CHS","AddAK8CHS",3, isData);

  // Select jets from Tree
  int nsubjets;
  int jet_counter = 0;
  for(int i0 = 0; i0 < int(lTree->GetEntriesFast()); i0++) {
    fVJet8->load(i0);
    int ifill = 0;
    for  (int i1 = 0; i1 < fVJet8->fVJets->GetEntriesFast(); i1++) {
      TJet *pVJet = (TJet*)((*fVJet8->fVJets)[i1]);
      if(pVJet->pt   <=  200) continue;
      if(fabs(pVJet->eta) >  2.4) continue;
      // select pt bin
      if(pVJet->pt   <  500 || pVJet->pt > 600) continue;
      if(ifill > 1) continue;
      ifill += 1;
      std::vector<TPFPart*> jetPFs;
      for (auto idx : pVJet->pfCands) {
        jetPFs.push_back( (TPFPart*)(fVJet8->fPFs->At(idx)) );
      }
      vector<PseudoJet> particles;
      nsubjets = jetPFs.size();
      // Particle loop
      for (auto *pf : jetPFs) {
	double pPt, pEta, pPhi, pM;
	pPt = pf->pup * pf->pt / pVJet->pt;
	pEta = pf->eta - pVJet->eta;
	pPhi = SignedDeltaPhi(pf->phi, pVJet->phi);
	pM = pf->m;
	TLorentzVector jetPF; jetPF.SetPtEtaPhiM(pPt,pEta,pPhi,pM);
	PseudoJet particle(jetPF.Px(),jetPF.Py(),jetPF.Pz(),jetPF.E());
	particles.push_back(particle);
      } // end loop (particles)
      // Recluster jet constituents
      JetDefinition reclust_def(ee_genkt_algorithm, 7, recluster_def); // radius unused
      ClusterSequence reclust_seq(particles, reclust_def);
      vector<PseudoJet> reclust_jets = reclust_seq.exclusive_jets(1);
      if (jet_counter == 0)
	cout << "Reclustered with " << reclust_def.description() << endl;
      if(int(reclust_jets[0].constituents().size()) != nsubjets){
	cout << "Mismatch between nsubjets and reclustered jet" << endl;
      }
      outfile << "J " << jet_counter << endl;
      convert_cluster_sequence_to_JUNIPR(reclust_seq, outfile);

      jet_counter++;
    } // end loop (jets)
  } // end loop (events)

  outfile.close();
  return 0;

} // end function (main)

