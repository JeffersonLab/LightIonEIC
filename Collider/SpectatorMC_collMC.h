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

	// Initialize Monte Carlo
// const double xMin = 0.01;
// const double xMax = 0.05;
// const double Q2Min= 15.0;
// const double Q2Max= 20.0;

	// Initialize Monte Carlo : fixed target 
//const double xMin = 0.02;
//const double xMax = 0.04;
//const double Q2Min= 15.0;
//const double Q2Max= 20.0;

// slight expand to avoid bin migration
 // const double xMin = 0.015;
 // const double xMax = 0.045;
//const double xMin = 0.06;
//const double xMax = 0.08;
 const double xMin = 0.1;
 const double xMax = 0.9;
// const double Q2Min= 15.;
// const double Q2Max= 20.;

const double Q2Min= 5.;
const double Q2Max= 10.;



// const double pSMax=  0.3;
const double pSMax=  0.3;


	//  Physical Constants
// const double MProton   = 0.93827;
// const double MNeutron  = 0.93957;
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
//const double kBeam =   20.;  // Electron Beam Momentum
const double kBeam =  5.;  // Electron Beam Momentum
	// Initialize Beam : fixed target case
// const double PBeam = 0.0;  // ion Beam momentum/Z, GeV/c
// const double kBeam =   10.;  // Electron Beam Momentum
const double ZBeam =   1.;
const double ABeam =   2.;
const double CrossingTheta = 0.050;
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


// const double eEpsNX     = 54.e-6/5; // Normalized emittance values (m*radian) 1/5x times : Transverse
// const double eEpsNY     = 11.e-6/5;
// const double iEpsNX     = 0.35e-6;
// const double iEpsNY     = 0.07e-6/5;


// const double eEpsNX     = 54.e-6*5; // Normalized emittance values (m*radian) 5x times : Transverse
// const double eEpsNY     = 11.e-6*5;
// const double iEpsNX     = 0.35e-6*0.3;
//  const double iEpsNY     = 0.07e-6*0.3;

// const double eEpsNX     = 0.; // Normalized emittance values (m*radian) 0x times : Transverse
// const double eEpsNY     = 0.;
// const double iEpsNX     = 0.;
//const double iEpsNY     = 0.;



const double eDkOverk   = 7.1e-4; // Fractional energy spread  Normalized emittance values : Longitudinal
const double iDPoverP   = 3.0e-4;


// const double eDkOverk   = 7.1e-4/5; // Fractional energy spread : Longitudinal
// const double iDPoverP   = 3.0e-4/5;
// const double eDkOverk   = 7.1e-4*5; // Fractional energy spread : Longitudinal
// const double iDPoverP   = 3.0e-4/2;

// const double eDkOverk   = 0.; // Fractional energy spread : Longitudinal
// const double iDPoverP   = 0.;



	// Global event-by-event Invariant Variable Structure
