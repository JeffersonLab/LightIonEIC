/* for the fixed target configuration
 *  SpectatorMC.h
 *  
 *
 *  Created by Charles Hyde on 3/4/14.
 *  Copyright 2014. All rights reserved.
 *
 */
#include <stdio.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iomanip>


using std::cout;
using std::endl;
using std::scientific;
using std::fixed;
using std::ios;

 const int NEvts = 110000;


// Asymmetry and g1 neutron spin structure 
// function study for various x, and Q2




const double pSMax=  0.3;


	//  Physical Constants
const double MProton   = 0.93827203;
const double MNeutron  = 0.93956536;
const double mElectron = 0.5110e-3;
const double MDeut     = MProton+MNeutron-0.0022;
const double MBindA    = -0.008; // Average Binding Energy
const double MPSq      = MProton*MProton;
const double alphaQED  = 1./137.03;
const double pi        = acos(-1.0);

	// Initialize Beam
const double PBeam = 100.0;  // ion Beam momentum/Z, GeV/c
const double kBeam =  5.;  // Electron Beam Momentum
//  electron and ion beam polarization
const double eBeamPol = 1.;
const double DBeamPol = 1.;



const double ZBeam =   1.;
const double ABeam =   2.;
//const double CrossingTheta = 0.050;
const double CrossingTheta = -0.050;
//0.050; // 0.010 for eRHIC/ePHENIX
const double CrossingPhi   = 0.0;   // pi/2.0 for eRHIC/ePHENIX

const double eBetaStarX = 0.10; // electron, ion beta at IP (m)
const double eBetaStarY = 0.02;
const double iBetaStarX = 0.10;
const double iBetaStarY = 0.02;

const double eEpsNX     = 54.e-6; // Normalized emittance values (m*radian) nominal setup : Transverse
const double eEpsNY     = 11.e-6;

const double iEpsNX     = 0.35e-6;
const double iEpsNY     = 0.07e-6;



const double eDkOverk   = 7.1e-4; // Fractional energy spread  Normalized emittance values : Longitudinal
const double iDPoverP   = 3.0e-4;



	// Global event-by-event Invariant Variable Structure

