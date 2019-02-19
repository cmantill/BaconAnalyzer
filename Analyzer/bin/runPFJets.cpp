// run:
// ./runPFJetsToJUNIPR input_directory input_file output_directory recluster_def

#include "../include/PerJetLoader.hh"

#include "fastjet/ClusterSequence.hh"
#include "../include/LundGenerator.hh"
#include "../include/LundJSON.hh"

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
  //double recluster_def = atof(argv[3]);
  string isdata = argv[4];

  bool isData;
  if(isdata.compare("data")!=0) isData = false;
  else isData = true;

  // Input, output
  stringstream s_infile;
  s_infile << input_file;

  stringstream s_outfile;
  s_outfile << "EFN_format_" << output_file ;
  ofstream outfile;
  outfile.open( s_outfile.str().c_str());
  
  // Read Tree
  const std::string lName = s_infile.str().c_str();
  TTree *lTree = load( lName );
  fVJet8     = new PerJetLoader  (lTree,"AK8Puppi","AddAK8Puppi","AK8CHS","AddAK8CHS",3, isData);

  // Select jets from Tree
  //int nsubjets;
  //int jet_counter = 0;
  for(int i0 = 0; i0 < int(lTree->GetEntriesFast()); i0++) {
    fVJet8->load(i0);
    int ifill = 0; // only leading jet
    for  (int i1 = 0; i1 < fVJet8->fVJets->GetEntriesFast(); i1++) {
      TJet *pVJet = (TJet*)((*fVJet8->fVJets)[i1]);
      if(pVJet->pt   <=  200) continue;
      if(fabs(pVJet->eta) >  2.4) continue;
      // select pt bin
      //if(pVJet->pt <  500 || pVJet->pt > 600) continue;
      //if(pVJet->pt <  600) continue;
      if(ifill > 1) continue;
      ifill += 1;
      std::vector<TPFPart*> jetPFs;
      for (auto idx : pVJet->pfCands) {
        jetPFs.push_back( (TPFPart*)(fVJet8->fPFs->At(idx)) );
      }
      //nsubjets = jetPFs.size();
      //outfile << "J " << jet_counter << " " << nsubjets << endl;
      //outfile << "J " << jet_counter << endl;
      //outfile << "N " << nsubjets << endl;
      // Particle loop
      vector<PseudoJet> particles;
      for (auto *pf : jetPFs) {
	double pPt, pEta, pPhi, pM;
	pPt = pf->pup * pf->pt / pVJet->pt;
	pEta = pf->eta - pVJet->eta;
	pPhi = SignedDeltaPhi(pf->phi, pVJet->phi);
	pM = pf->m;
	TLorentzVector pVec; pVec.SetPtEtaPhiM(pPt,pEta,pPhi,pM);
	//std::cout << "C " << pf->pup << " " << pVec.Px() << " " << pVec.Py() << " " << pVec.Pz() << " " << pf->pfType << std::endl;   
	//if( pVec.E() == 0.139526){
	if( pf->pup == 0) {
	  continue;
	}
	//outfile << "C " << pPt << " " << pEta << " " << pPhi << " " << pf->pfType << endl;
	//outfile << pVec.Px() << " " << pVec.Py() << " " << pVec.Pz() << " " << pVec.E() << endl;
	// get some pseudojets
	PseudoJet particle(pVec.Px(),pVec.Py(),pVec.Pz(),pVec.E());
	particles.push_back(particle);
      } // end loop (particles)                                                                                                                                                      
      JetDefinition jet_def(antikt_algorithm, 1.0);
      ClusterSequence cs(particles, jet_def);
      //vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(100.));
      if(particles.size()>0){
	vector<PseudoJet> jets = cs.exclusive_jets(1);
	
	LundGenerator lund;
	//std::cout << lund.description() << endl;
	for (unsigned int ijet = 0; ijet < jets.size(); ijet++) {
	  //std::cout << std::endl << "Lund coordinates ( ln 1/Delta, ln kt ) of declusterings of jet " 
	  //<< ijet << " are:" << std::endl;
	  vector<LundDeclustering> declusts = lund(jets[ijet]);
	  
	  //for (unsigned int idecl = 0; idecl < declusts.size(); idecl++) {
	  //  pair<double,double> coords = declusts[idecl].lund_coordinates();
	    //std::cout << "(" << coords.first << ", " << coords.second << ")";
	    //if (idecl < declusts.size() - 1) std::cout << "; ";
	  //}
	  //std::cout << endl;
    
	  // outputs the primary Lund plane
	  lund_to_json(outfile, declusts);
	  outfile << endl;
	  // outputs the full Lund tree
	  //to_json(cout, lund_gen, jets[ijet]); cout << endl;
	} // end loop (jets of particles)
      }
    } // end loop (vjets)
  } // end loop (events)

  outfile.close();
  return 0;

} // end function (main)

