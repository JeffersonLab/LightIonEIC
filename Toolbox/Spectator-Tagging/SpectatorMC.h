/*
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

	// Initialize Monte Carlo
const int NEvts = 100000;
const double xMin = 0.1;
const double xMax = 0.5;
const double Q2Min= 0.5;
const double Q2Max= 5.0;
const double pSMax=  0.25;

	//  Physical Constants
const double MProton   = 0.93827;
const double MNeutron  = 0.93957;
const double mElectron = 0.5110e-3;
const double MDeut     = MProton+MNeutron-0.0022;
const double MBindA    = -0.008; // Average Binding Energy
const double MPSq      = MProton*MProton;
const double alphaQED  = 1./137.03;
const double pi        = acos(-1.0);

	// Initialize Beam
const double PBeam = 100.0;  // ion Beam momentum/Z, GeV/c
const double kBeam =   5.;  // Electron Beam Momentum
const double ZBeam =   1.;
const double ABeam =   2.;
const double CrossingTheta = 0.050; // 0.010 for eRHIC/ePHENIX
const double CrossingPhi   = 0.0;   // pi/2.0 for eRHIC/ePHENIX
const double eEpsNX     = 54.e-6; // Normalized emittance values (m*radian)
const double eEpsNY     = 11.e-6; //  11.e-6;
const double iEpsNX     = 0.35e-6; //  0.35e-6;
const double iEpsNY     = 0.07e-6; //  0.07e-6;
const double eBetaStarX = 0.10; // electron, ion beta at IP (m)
const double eBetaStarY = 0.02;
const double iBetaStarX = 0.10;
const double iBetaStarY = 0.02;
const double eDkOverk   = 7.1e-4; // Fractional energy spread
const double iDPoverP   = 3.0e-4; //   3.0e-4;



	// Global event-by-event Invariant Variable Structure

