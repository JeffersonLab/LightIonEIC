// Made by Oz Amram
// June 2014

#include <stdio.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TMath.h>
#include <utility>
#include <iostream>

void Tran_to_fixed_target (TLorentzVector, TLorentzVector, 
                        TLorentzVector &, TLorentzVector &);

void Tran_to_cm(TLorentzVector e_vec, TLorentzVector d_vec,
                TLorentzVector &e_cm, TLorentzVector &d_cm);

void print_vec(TLorentzVector);

TLorentzVector rotate(TLorentzVector v, Float_t theta, Float_t phi);

void print_se(TLorentzVector, TLorentzVector);
const double_t deut_mass = 1.875628;
const double_t pro_mass = 0.938272;
const double_t neu_mass = 0.939565;

using namespace std;
