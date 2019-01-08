#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
using namespace fastjet;
using namespace std;

vector<double> PJ_3mom(PseudoJet PJ);
double dot3(vector<double> v1, vector<double> v2);
vector<double> hat(vector<double> v);
double get_theta(vector<double> v1, vector<double> v2);
double get_thetaPJ(PseudoJet PJ1, PseudoJet PJ2);
vector<double> cross(vector<double> v1, vector<double> v2);
double get_phiPJ(PseudoJet p, PseudoJet ref);
vector<double> get_branching(PseudoJet mother, PseudoJet daughter1, PseudoJet daughter2);
vector<double> get_mother_momenta(PseudoJet mother, PseudoJet ref);
vector<double> get_daughter_momenta(PseudoJet daughter1, PseudoJet daughter2, PseudoJet ref);
void convert_cluster_sequence_to_JUNIPR(ClusterSequence cs, ostream& outfile);
