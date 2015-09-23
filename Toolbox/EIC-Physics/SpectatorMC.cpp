/*
 *  SpectatorMC.cpp
 *  
 *  Copyright (C) 2014 by Charles Hyde
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 */

#include "SpectatorMC.h"

	// Define function calls
double sigma_th(double pInc, double mInc, double NormEmit, double betaSt);
	//
	// ********************************************* main()
int main(){
	TRandom3 ran3;
	ran3.SetSeed(13579);
	Bool_t iran=kTRUE;      // TRUE==include incident beam emittance
	Bool_t iProton = kTRUE; // TRUE == Spectator Proton
	double MSpectator;
	if (iProton) {
		MSpectator=MProton;
	}else {
		MSpectator=MNeutron;
	}

	char tTitle[80], tName[18], rName[32];
	sprintf(tTitle,"D(e,e'N)X Event Generation %3.0f GeV/c x%4.0f GeV/c",
			kBeam, PBeam);
	sprintf(rName,"DeepS-MC%02.0fx%03.0f.root", kBeam, PBeam);
		//	TFile *fRoot = TFile::Open(rName,"Recreate");
	TFile fRoot(rName,"Recreate", tTitle);
	sprintf(tName,"DeepS_MC%02.0fx%03.0f", kBeam, PBeam);
	TTree *tree = new TTree("Evnts",tTitle);
	TLorentzVector kIncident_Vertex, PIncident_Vertex, qVirtual_Vertex;
	TLorentzVector pSpectator_Vertex, kScattered_Vertex, PX_Vertex;
	typedef struct{
		Double_t s_e, Q2, xBj, x_D, y_D, tSpectator, tPrime, TwoPdotk, tPrime0,
		         MX2, alphaS, pPerpS;
	} Invariants;
	static Invariants invts;
	Double_t Jacob;

	const Int_t bufsize=32000;
		//	TBranch *branch = tree->Branch("e_inc.", &kIncident_Vertex,bufsize,1);
	tree->Branch("e_inc.", &kIncident_Vertex,bufsize,1);
	tree->Branch("P_inc.",  &PIncident_Vertex, bufsize, 1);
	tree->Branch("e_scat.", &kScattered_Vertex, bufsize, 1);
	tree->Branch("p_Sp.", &pSpectator_Vertex, bufsize, 1);
	tree->Branch("q_V.", &qVirtual_Vertex, bufsize, 1);
	tree->Branch("P_X.", &PX_Vertex, bufsize, 1);
	tree->Branch("invts",&invts,
	    "s_e/D:Q2/D:xBj/D:x_D/D:y_D/D:tSpectator/D:tPrime/D:TwoPdotk/D:tPrime0/D:MX2/D:alphaS/D:pPerpS/D");
		//	tree->Branch("psf",&psf,"psf/D");
	tree->Branch("Jacob",&Jacob,"Jacob/D");

	double sig_eThx=0.0, sig_eThy=0.0, sig_iThx=0.0, sig_iThy=0.0;
		//	double sig_k=0.0, sig_P=0.0;
	double MIon;
		// Global Phase Space Factor;
		//s_e - M_ion^2 - m_electron^2
	double uu, norm;
	if (ZBeam==1.0&&ABeam==2.0) {
		MIon = MDeut;
	} else {
		MIon = ZBeam*MProton+(ABeam-ZBeam)*MNeutron+ABeam*MBindA; // Average Binding Energy
	}
	double PIon = PBeam*ZBeam;
	double s_e0 = 2.*kBeam*(sqrt(PIon*PIon+MIon*MIon)+PIon*cos(CrossingTheta))
	+ MIon*MIon + mElectron*mElectron;
	printf("Incident Ion Mass %6.2f GeV \n", MIon);
	printf("Incident Electron, Ion Momenta:  %8.4f, %8.2f GeV/c;  s_0 = %10.4f GeV^2 \n", 
		   kBeam, PIon, s_e0);
	TH1D *hse = new TH1D("hse","Invariant s_e =(k+P)^{2} (GeV^{2})",
						 100, s_e0*0.99,s_e0*1.01);
		//	tree->Branch("hse","TH1D",&hse,bufsize,0);
	int nQ2 = 20;
	if (Q2Min>0.0) {
		nQ2 = Q2Max/Q2Min+1;
	}
		// Set up log(x) bins
	const int nxBins=50+1;
	double xLowBinEdge[nxBins], xBinFact;
	xBinFact = pow(xMax/xMin,1./(nxBins-1.0));
	for (int ix=0; ix<nxBins; ix++) {
		xLowBinEdge[ix] = xMin*pow(xBinFact,ix);
	}
	TH2D *hQ2vsX = new TH2D("hQ2vsX",
							"Q^{2} vs. x_{Bj};  x_{Bj}  ; Q^{2} (GeV^{2})  ",
							nxBins-1,xLowBinEdge, nQ2, 0.0,Q2Max);
		//	tree->Branch("hQ2vsX","TH2D",&hQ2vsX,bufsize,0);
	TH2D *hRejects = new TH2D("hRejects",
							"Q^{2} vs. x_{Bj};  x_{Bj}  ; Q^{2} (GeV^{2})  ",
							nxBins-1,xLowBinEdge, nQ2, 0.0,Q2Max);
		//	tree->Branch("hRejects","TH2D",&hRejects,bufsize,0);
	TH1D *htS = new TH1D("htS", 
		"t'=M_{S}^{2} - (P_{D}-p_{S})^{2} ; t' (GeV^{2})  ; Jacobian weighted events ",
						 100,0.0,1.0);
		//	tree->Branch("htS","TH1D",&htS,bufsize,0);
	TH1D *hnu = new TH1D("hnu", "#nu_{Rest}", 100.,-100.,100.);
	
	if (iran) {
		sig_eThx = sigma_th(kBeam, mElectron, eEpsNX, eBetaStarX);
		sig_eThy = sigma_th(kBeam, mElectron, eEpsNY, eBetaStarY);
		sig_iThx = sigma_th(PBeam, MIon, iEpsNX, iBetaStarX);
		sig_iThy = sigma_th(PBeam, MIon, iEpsNY, iBetaStarY);
	}
	double kBeamMC, kBeamMCx, kBeamMCy, kBeamMCz;
	double PBeamMC, PBeamMCx, PBeamMCy, PBeamMCz;
	double EScatRest, kScatRest, csTheRest, PhiScatRest;
	TVector3       UnitXLab, UnitYLab, UnitZLab;
	UnitXLab.SetXYZ(1.0,0.0,0.0);
	UnitYLab.SetXYZ(0.0,1.0,0.0);
	UnitZLab.SetXYZ(0.0,0.0,1.0);
	TVector3       UnitXRest, UnitYRest, UnitZRest;
	TVector3       UnitXqCM,  UnitYqCM,  UnitZqCM;
	TVector3       BoostCM, BoostRest;
	TVector3       kScat3vec, pS_3vec;
	TLorentzVector kIncident_Rest, PIncident_Rest, qVirtual_Rest;
	TLorentzVector kScattered_Rest, pSpectator_Rest;
	TLorentzVector kIncident_0, PIncident_0; // unsmeared
	kIncident_0.SetXYZM(0.0,0.0,-kBeam,mElectron);
	PIncident_0.SetXYZM(PIon*sin(CrossingTheta)*cos(CrossingPhi),
						PIon*sin(CrossingTheta)*sin(CrossingPhi),
						PIon*cos(CrossingTheta),MIon);
	TTree *itree = new TTree("Init",tTitle);
	itree->Branch("e0.", &kIncident_0,bufsize,1);
	itree->Branch("P0.",  &PIncident_0, bufsize, 1);
	typedef struct {
		Int_t nEvts;
		Float_t PhSpFct;
	} MonteCarlo;
	static	MonteCarlo mc;
	mc.nEvts = NEvts;
	mc.PhSpFct   = (Q2Max-Q2Min)*(log(xMax)-log(xMin))*4.*pi*pow(pSMax,3)/3.;
	itree->Branch("MC", &mc, "nEvts/I:PhSpFct/F");

	double pS_rest, csThRecoil, phiRecoil;
	Int_t MEvts=0;
	for (int iEvt=0; iEvt<NEvts; iEvt++) {
		Jacob = 1.0;
		if (iran) {
			kBeamMC = kBeam*ran3.Gaus(1.0,eDkOverk);
			kBeamMCx= kBeamMC*ran3.Gaus(0.0,sig_eThx);
			kBeamMCy= kBeamMC*ran3.Gaus(0.0,sig_eThy);
			kBeamMCz=-kBeamMC;  //  Angles are really Tangents
			PBeamMC = PIon*ran3.Gaus(1.0,iDPoverP);
			PBeamMCx= PBeamMC*ran3.Gaus(0.0,sig_iThx);
			PBeamMCy= PBeamMC*ran3.Gaus(0.0,sig_iThy);
			PBeamMCz= PBeamMC;
			
		} else {
			kBeamMC = kBeam;
			kBeamMCx= 0.0;
			kBeamMCy= 0.0;
			kBeamMCz= kBeam;
			PBeamMC = PIon;
			PBeamMCx= 0.0;
			PBeamMCy= 0.0;
			PBeamMCz= PBeamMC;
		}
		kIncident_Vertex.SetXYZM(kBeamMCx, kBeamMCy, kBeamMCz, mElectron);
		PIncident_Vertex.SetXYZM(PBeamMCx, PBeamMCy, PBeamMCz, MIon);
			//  Crossing angle MEIC or eRHIC
		PIncident_Vertex.RotateY(CrossingTheta);
		PIncident_Vertex.RotateZ(CrossingPhi);
			//
		invts.TwoPdotk = 2.*PIncident_Vertex.Dot(kIncident_Vertex);
		invts.s_e = MIon*MIon + mElectron*mElectron + invts.TwoPdotk;
			// Boost vectors from electron+Ion CM frame
			//		BoostCM = (kIncident_Vertex+PIncident_Vertex).BoostVector();
			// Boost vector from Ion rest frame
		BoostRest = PIncident_Vertex.BoostVector();
		PIncident_Rest = PIncident_Vertex;
		PIncident_Rest.Boost(-BoostRest); // should result in (0.,0.,0.,MIon)
		kIncident_Rest = kIncident_Vertex;
		kIncident_Rest.Boost(-BoostRest);
			// Generate Q2, ln(xBj) uniformly
			//  psf *= (Q2Max-Q2Min)*(log(xMax)-log(xMin));
		uu   = ran3.Uniform();
		invts.Q2   = Q2Max*uu + Q2Min*(1.-uu);
		uu   = ran3.Uniform();
		invts.xBj  = pow(xMin,1.-uu)*pow(xMax,uu);
		invts.x_D  = invts.xBj*MProton/MIon;
		invts.y_D  = invts.Q2/(invts.x_D*invts.TwoPdotk);
		if (invts.y_D>=(1.0-2.*mElectron*MIon/invts.TwoPdotk) ) {
				// Unphysical kinematics
			hRejects->Fill(invts.xBj,invts.Q2);
			continue;
		}
		if (invts.y_D>=1./(1.+invts.x_D*MIon*MIon/invts.TwoPdotk) ) {
				// Unphysical kinematics
			hRejects->Fill(invts.xBj,invts.Q2);
			continue;
		}

		Jacob *=invts.xBj;   // Jacobian  dx/dlnx
							 //		psf *=2.*pi;
		PhiScatRest = pi*(2.*ran3.Uniform()-1.0);
			//  mElectron-->0 approximation
		EScatRest = invts.TwoPdotk*(1.-invts.y_D)/(2.*MIon);
		if (EScatRest<mElectron) {
				// should never happen
			printf("illegal Rest frame scattered electron energy =%10.6f \n",
				   EScatRest);
			continue;
		}
		kScatRest = sqrt(EScatRest*EScatRest-mElectron*mElectron);
		csTheRest = (2.*EScatRest*kIncident_Rest.E() - invts.Q2 - 2.*mElectron*mElectron) / 
					(2.*kScatRest*kIncident_Rest.P());
		if (csTheRest*csTheRest>1.0) {
				// should never happen
			printf("illegal Rest frame cos(the) = %6.2f \n", csTheRest);
			printf(" (k_Rest, k'_Rest, Q2, xBj) = (%8.3f,  %8.4f, %6.2f, %5.3f) \n",
				   kIncident_Rest.E(),kScatRest,invts.Q2,invts.xBj);
			printf(" (2k.P, invts.s_e, y_D, x_D) = (%6.2f, %6.2f,%10.6f, %8.4f) \n",
				   invts.TwoPdotk, invts.s_e, invts.y_D, invts.x_D);
			continue;
		}
		UnitZRest = -kIncident_Rest.Vect();
		norm      = 1./UnitZRest.Mag();
		UnitZRest*= norm;
		UnitYRest = UnitZRest.Cross(UnitXLab);
		norm      = 1./UnitYRest.Mag();
		UnitYRest*= norm;
		UnitXRest = UnitYRest.Cross(UnitZRest);

		kScat3vec = kScatRest*(
					sin(acos(csTheRest))*(cos(PhiScatRest)*UnitXRest
										  +sin(PhiScatRest)*UnitYRest)
		           -csTheRest*UnitZRest);
		kScattered_Rest.SetVectM(kScat3vec, mElectron);
		qVirtual_Rest = kIncident_Rest-kScattered_Rest;
		kScattered_Vertex = kScattered_Rest;
		// Back to Lab frame
		kScattered_Vertex.Boost(BoostRest);
		qVirtual_Vertex = kIncident_Vertex-kScattered_Vertex;
			// Hadronic Unit vectors relative to q
		UnitZqCM = -qVirtual_Rest.Vect();
		norm     = 1./UnitZqCM.Mag();
		UnitZqCM*= norm;
		UnitYqCM = -kScat3vec.Cross(kIncident_Rest.Vect());
		norm     = 1./UnitYqCM.Mag();
		UnitYqCM*=norm;
		UnitXqCM = UnitYqCM.Cross(UnitZqCM);
			// Generate pSpectator Spherically Symmetric in rest frame 
			//		psf     *= 4.*pi*pSMax^3/3;
		uu       = ran3.Uniform();
		pS_rest   = pSMax*pow(uu,1./3.); // uniform in 3p^2 dp = d(p^3)
		uu       = ran3.Uniform();
		csThRecoil = (2.*uu-1.);
		uu       = ran3.Uniform();
		phiRecoil  = pi*(2.*uu-1.);
		pS_3vec  = pS_rest*sin(acos(csThRecoil))*
						(cos(phiRecoil)*UnitXqCM + sin(phiRecoil)*UnitYqCM);
		pS_3vec += pS_rest*csThRecoil*UnitZqCM;
		pSpectator_Rest.SetVectM(pS_3vec,MSpectator);
		Jacob     *= 1./(2.*pSpectator_Rest.E());
		pSpectator_Vertex = pSpectator_Rest;
		pSpectator_Vertex.Boost(BoostRest);
		PX_Vertex = kIncident_Vertex+PIncident_Vertex
		            -kScattered_Vertex-pSpectator_Vertex;
		invts.MX2 = PX_Vertex.M2();
		invts.alphaS = ABeam*(pS_rest*csThRecoil+pSpectator_Rest.E())/MIon;
		invts.pPerpS = pS_rest*sqrt(1.-csThRecoil*csThRecoil);
		hse->Fill(invts.s_e,Jacob);
		hQ2vsX->Fill(invts.xBj,invts.Q2, Jacob);
		
		invts.tSpectator = MIon*MIon+MSpectator*MSpectator
		           - 2.*pSpectator_Vertex.Dot(PIncident_Vertex);
		invts.tPrime     = 2.*pSpectator_Vertex.Dot(PIncident_Vertex) - MIon*MIon;
		invts.tPrime0    = 2.*pSpectator_Vertex.Dot(PIncident_0) - MIon*MIon;
		htS->Fill(invts.tPrime,Jacob);
		hnu->Fill(qVirtual_Rest.E());
		tree->Fill();
		MEvts++;
	}
	itree->Fill();
	fRoot.Write();
	fRoot.Close();
	printf("Total of %d events out of %d Trials \n",MEvts,NEvts);
	return MEvts;
}
double sigma_th(double pInc, double mInc, double NormEmit, double betaSt){
	double gamma = sqrt(pInc*pInc+mInc*mInc);
	double sig   = sqrt(NormEmit/(gamma*betaSt));
	return sig;
}
