// run:
// ./runPFJetsToJUNIPR input_directory input_file output_directory recluster_def

#include "../include/PerJetLoader.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"

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

  //TFile *lFile = TFile::Open("Output.root","RECREATE");
  TFile *lFile = TFile::Open(output_file.c_str(),"RECREATE");
  TTree *lOut  = new TTree("Events","Events");

  stringstream s_outfile;
  s_outfile << "EFN_format_" << output_file ;
  ofstream outfile;
  outfile.open( s_outfile.str().c_str());

  int jetnum(-1);
  float mass(-1),softdropmass(-1),pT(-1),N2(-1),tau32(-1),tau21(-1),D2(-1);
  float p_pt(-1),p_m(-1),h_pt(-1),s_pt(-1),z(-1),delta(-1),kt(-1),psi(-1);

  lOut->Branch("jetnum",&jetnum,"jetnum/I");
  lOut->Branch("mass",&mass,"mass/F");
  lOut->Branch("softdropmass",&softdropmass,"softdropmass/F");
  lOut->Branch("pT",&pT,"pT/F");
  lOut->Branch("N2",&N2,"N2/F");
  lOut->Branch("D2",&D2,"D2/F");
  lOut->Branch("tau32",&tau32,"tau32/F");
  lOut->Branch("tau21",&tau21,"tau21/F");
  lOut->Branch("p_pt",&p_pt,"p_pt/F");
  lOut->Branch("p_m",&p_m,"p_m/F");
  lOut->Branch("h_pt",&h_pt,"h_pt/F");
  lOut->Branch("s_pt",&s_pt,"s_pt/F");
  lOut->Branch("z",&z,"z/F");
  lOut->Branch("delta",&delta,"delta/F");
  lOut->Branch("kt",&kt,"kt/F");
  lOut->Branch("psi",&psi,"psi/F");

  // Read Tree
  const std::string lName = s_infile.str().c_str();
  TTree *lTree = load( lName );
  fVJet8     = new PerJetLoader  (lTree,"AK8Puppi","AddAK8Puppi","AK8CHS","AddAK8CHS",3, isData);

  // Select jets from Tree
  //int nsubjets;
  int jet_counter = 0;
  for(int i0 = 0; i0 < int(lTree->GetEntriesFast()); i0++) {
    fVJet8->load(i0);
    int ifill = 0; // only leading jet
    for  (int i1 = 0; i1 < fVJet8->fVJets->GetEntriesFast(); i1++) {
      TJet *pVJet = (TJet*)((*fVJet8->fVJets)[i1]);
      if(pVJet->pt   <=  200) continue;
      if(fabs(pVJet->eta) >  2.4) continue;
      if(ifill > 1) continue;
      ifill += 1;
      std::vector<TPFPart*> jetPFs;
      for (auto idx : pVJet->pfCands) {
        jetPFs.push_back( (TPFPart*)(fVJet8->fPFs->At(idx)) );
      }
      jetnum = jet_counter;

      // addjet info
      TAddJet *pAddJet = fVJet8->getAddJet(pVJet);
      mass = pVJet->mass;
      softdropmass = pAddJet->mass_sd0;
      pT = pVJet->pt;
      N2 = pAddJet->e3_v2_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1);
      D2 = pAddJet->e3_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1*pAddJet->e2_sdb1);
      tau32 = (pAddJet->tau3/pAddJet->tau2);
      tau21 = (pAddJet->tau2/pAddJet->tau1);

      // jet
      TLorentzVector pVec; pVec.SetPtEtaPhiM(pVJet->pt,pVJet->eta,pVJet->phi,pVJet->mass);
      PseudoJet jet(pVec.Px(),pVec.Py(),pVec.Pz(),pVec.E());
      
      // trying to do lund on that jet
      // LundGenerator lund;                                                                                                                                                      
      // vector<LundDeclustering> declusts = lund(jet);
      // for (unsigned int idecl = 0; idecl < declusts.size(); idecl++) {                                                                                                           
      // 	p_pt = declusts[idecl].pair().pt();                                                                                                                                
      // 	p_m = declusts[idecl].pair().m();                                                                                                                                 
      // 	h_pt =  declusts[idecl].harder().pt();                                                                                                                           
      // 	s_pt = declusts[idecl].softer().pt();                                                                                                                           
      // 	z = declusts[idecl].z();                                                                                                                                         
      // 	delta = declusts[idecl].Delta();                                                                                                                                 
      // 	kt =  declusts[idecl].kt();                                                                                                                                      
      // 	psi = declusts[idecl].psi();                                                                                                                                     
      // }                                                                                                                                                                          
      // lOut->Fill(); 

      // Vector of psuedojets for jet particles
      vector<PseudoJet> particles;
      for (auto *pf : jetPFs) {
	double pPt, pEta, pPhi, pM;
	pPt = pf->pup * pf->pt / pVJet->pt;
	pEta = pf->eta - pVJet->eta;
	pPhi = SignedDeltaPhi(pf->phi, pVJet->phi);
	pM = pf->m;
	TLorentzVector pVec; pVec.SetPtEtaPhiM(pPt,pEta,pPhi,pM);
	if( pf->pup == 0) {
	  continue;
	}
	// get some pseudojets
	PseudoJet particle(pVec.Px(),pVec.Py(),pVec.Pz(),pVec.E());
	particles.push_back(particle);
      } // end loop (particles)                                                                                                                                                      
      // reclustering
      JetDefinition jet_def(antikt_algorithm, 1.0);
      ClusterSequence cs(particles, jet_def);
      //vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
      //vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(100.));

      if(particles.size()>0){
	vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));

	LundGenerator lund;

	//std::cout << "size of recluster jets " << jets.size() << std::endl;
	for (unsigned int ijet = 0; ijet < jets.size(); ijet++) {
	  vector<LundDeclustering> declusts = lund(jets[ijet]);
	  
	  //std::cout << "size of lund plane " <<declusts.size() << std::endl;

	  //for (unsigned int ijet = 0; ijet < particles.size(); ijet++) {
	  //vector<LundDeclustering> declusts = lund(particles[ijet]);
	  
	  for (unsigned int idecl = 0; idecl < declusts.size(); idecl++) {
	    p_pt = declusts[idecl].pair().pt();
	    p_m = declusts[idecl].pair().m();
	    h_pt =  declusts[idecl].harder().pt();
	    s_pt = declusts[idecl].softer().pt();
	    z = declusts[idecl].z();
	    delta = declusts[idecl].Delta();
	    kt =  declusts[idecl].kt();
	    psi = declusts[idecl].psi();
	    lOut->Fill();
	  }

	  // outputs the primary Lund plane
	  //lund_to_json(outfile, declusts);
	  //outfile << endl;
	  // outputs the full Lund tree
	  //to_json(cout, lund_gen, particles[ijet]); cout << endl;

	} // end loop (jets of particles)
      }

      jet_counter+=1;
      //lOut->Fill();
    } // end loop (vjets)
  } // end loop (events)

  outfile.close();
  lFile->cd();
  lOut->Write();  
  lFile->Close();
  return 0;

} // end function (main)

