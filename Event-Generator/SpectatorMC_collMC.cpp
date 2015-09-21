/*
 *  SpectatorMC.cpp
 *  
 *
 *  Created by Charles Hyde on 3/4/14.
 *  Copyright 2014. All rights reserved.
 *  
 *  Major modification is done 2/20/15 @ K.Park
 */

#include "SpectatorMC_collMC.h"
#include <time.h>

	// Define function calls
double sigma_th(double pInc, double mInc, double NormEmit, double betaSt);
	//
	// ********************************************* main()
int mainx(double xMin,double xMax, double Q2Min,double Q2Max){

    srand (time(NULL));
    double rnumber = rand() ;
    TRandom3 ran3;
    // ran3.SetSeed(13579);
    ran3.SetSeed(rnumber);

    cout << "test random number = "<< rnumber << endl;

    Bool_t iran=kTRUE;      // TRUE==include incident beam emittance
    //Bool_t iran=kFALSE;      // Falut NO incident beam emittance

    Bool_t iProton = kTRUE; // TRUE == Spectator Proton
    //  Bool_t iProton = kFALSE; // FALSE == Spectator Neutron

    double MSpectator;
    double W2;
    int NumPtls;
    int charge,e_charge;

    // spectator particle information
    double_t psx_Lab,psy_Lab,psz_Lab,EsE_Lab ;
    double_t psx_Res,psy_Res,psz_Res,EsE_Res ;
    double_t vsx_Lab,vsy_Lab,vsz_Lab ;
    double_t vsx_Res,vsy_Res,vsz_Res ;
    int sp_particle_id;
    double spmass;
    //    double alpha_r;


    // scattered electron information
    double pex_Res,pey_Res,pez_Res,EeE_Res;
    double pex_Lab,pey_Lab,pez_Lab,EeE_Lab ;
    double_t vex_Lab,vey_Lab,vez_Lab ;
    double_t vex_Res,vey_Res,vez_Res ;
    int e_particle_id;
    double emass;
    emass =  mElectron;
    e_particle_id = 11;
    e_charge = -1;

    // missing particle (X) information
    double pXx_Res,pXy_Res,pXz_Res,EXE_Res ;
    double_t vXx_Res,vXy_Res,vXz_Res ;
    double pXx_Lab,pXy_Lab,pXz_Lab,EXE_Lab ;
    double_t vXx_Lab,vXy_Lab,vXz_Lab ;
    double dummass,zeros;
    int dum_particle_id;
    dum_particle_id = 71;

	if (iProton) {
		MSpectator=MProton;
		charge =1;
	        spmass =MSpectator;
		sp_particle_id = 2212;
	}else {
		MSpectator=MNeutron;
		charge =0;
	        spmass =MSpectator;
		sp_particle_id = 2112;
	}

	char tTitle[80], tName[18], rName[32];
	sprintf(tTitle,"D(e,e'N)X Event Generation %3.0f GeV/c x%4.0f GeV/c",
			kBeam, PBeam);
	sprintf(rName,"DeepS-MC%02.0fx%03.0f.root", kBeam, PBeam);
		//	TFile *fRoot = TFile::Open(rName,"Recreate");
	TFile fRoot(rName,"Recreate", tTitle);
	sprintf(tName,"DeepS_MC%02.0fx%03.0f", kBeam, PBeam);
	TTree *tree = new TTree("Evnts",tTitle);
	TLorentzVector kIncident_Vertex, PIncident_Vertex, qVirtual_Vertex, kScattered_Vertex;
	TLorentzVector kIncident_Vertex0, PIncident_Vertex0, qVirtual_Vertex0, kScattered_Vertex0;
	TLorentzVector pSpectator_Vertex, pSpectator_Vertex0, PX_Vertex, PX_Vertex_Rest;
	typedef struct Invariants {
	    Double_t s_e, s_e0, s_q, s_q0, Q2, xBj, nu, W, x_D, y_D, x_D0, y_D0, 
                 tSpectator, tPrime, TwoPdotk, TwoPdotq, TwoPdotk0, TwoPdotq0, tPrime0, p_RT, pDrest, pDrest0, tempVar,
		MX2, alphaS, pPerpS, pPerpZ;
	};
	Invariants invts;

        Double_t Jacob;


	Double_t qMag, pDotq;
	Double_t qMag0, pDotq0;

	const Int_t bufsize=32000;
		//	TBranch *branch = tree->Branch("e_inc.", &kIncident_Vertex,bufsize,1);
	tree->Branch("e_inc.", &kIncident_Vertex,bufsize,1);
	tree->Branch("P_inc.",  &PIncident_Vertex, bufsize, 1);
	tree->Branch("P0_inc.",  &PIncident_Vertex0, bufsize, 1);
	tree->Branch("e_scat.", &kScattered_Vertex, bufsize, 1);
	tree->Branch("p_Sp.", &pSpectator_Vertex, bufsize, 1);
	tree->Branch("q_V.", &qVirtual_Vertex, bufsize, 1);
	tree->Branch("P_X.", &PX_Vertex, bufsize, 1);
/*	tree->Branch("invts",&invts,
	    "s_e/D:s_e0:s_q/s_q0/D:Q2/D:xBj/D:nu/D:W/D:x_D/D:y_D/D:x_D0/D:y_D0/D:tSpectator/D:tPrime/TwoPdotk/TwoPdotq/TwoPdotk0/TwoPdotq0/tPrime0/D:p_RT/D:pDrest/D:pDrest0/tempVar/D:MX2/D:alphaS/D:pPerpS/D:pPerpZ");
		//	tree->Branch("psf",&psf,"psf/D");
*/
	tree->Branch("invts",&invts,
	    "s_e/D:s_e0:s_q:s_q0:Q2/D:xBj/D:nu/D:W/D:x_D/D:y_D/D:x_D0/D:y_D0/D:tSpectator/D:tPrime:TwoPdotk:TwoPdotq:TwoPdotk0:TwoPdotq0:tPrime0/D:p_RT/D:pDrest/D:pDrest0:tempVar/D:MX2/D:alphaS/D:pPerpS/D:pPerpZ/D");
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
	printf("Your kinematics: [xBj_min:xBj_max] = [%9.6f:%9.6f] \n", xMin,xMax);
	printf("Your kinematics: [Q2_min:Q2_max] = [%9.6f:%9.6f] \n", Q2Min,Q2Max);

	printf("Incident Ion Mass %6.2f GeV \n", MIon);
	printf("Incident Electron, Ion Momenta:  %8.4f, %8.2f GeV/c;  s_0 = %10.4f GeV^2 \n", 
		   kBeam, PIon, s_e0);
	 TH1D *hse = new TH1D("hse","Invariant s_e =(k+P)^{2} (GeV^{2})",
	 					 100, s_e0*0.99,s_e0*1.01);
	// 	//	tree->Branch("hse","TH1D",&hse,bufsize,0);
	 int nQ2 = 50;
	 if (Q2Min>0.0) {
	     nQ2 = Q2Max/Q2Min+1;
	 }
	// Set up log(y) bins
	// const int nyBins=10+1;
	// double yLowBinEdge[nyBins], yBinFact;
	// yBinFact = pow(yMax/yMin,1./(nyBins-1.0));
        // for(int iy=0; iy<nyBins; iy++){
	//     yLowBinEdge[iy] = yMin*pow(yBinFact,iy);
	// }
	//	int s_eN = 2002.4442;

 	// Set up log(x) bins
	const int nxBins=50+1;
	double xLowBinEdge[nxBins], xBinFact;
	xBinFact = pow(xMax/xMin,1./(nxBins-1.0));
	for (int ix=0; ix<nxBins; ix++) {
		xLowBinEdge[ix] = xMin*pow(xBinFact,ix);
	 }
	//const int nQ2 = 100;
	//	double Q2Max = xMax*yMax*s_eN;
	//double Q2Min = xMin*yMin*s_eN;

	TH2D *hQ2vsX = new TH2D("hQ2vsX",
				"Q^{2} vs. x_{Bj};  x_{Bj}  ; Q^{2} (GeV^{2})  ",
				// 200,0.049,0.051,200,34.5,35.5);
				// 200,0.005,0.11,200,0.,55.);
				200,0.0001,1.0,500,0.1,200.);
        
				// nxBins-1,xLowBinEdge, nQ2, 0.01,Q2Max);
        

	// TH1D *htS = new TH1D("htS", 
	// 		     "#delta p_{Dx} Collider ; #delta p_{Dx} (GeV) ; No.EVNT",
	// 		     200,-0.1,0.1);
	// TH1D *htSY = new TH1D("htSY", 
	// 		      "#delta p_{Dy} Collider; #delta p_{Dy} (GeV) ; No.EVNT",
	// 		      200,-0.1,0.1);
	// TH1D *htSZ = new TH1D("htSZ", 
	// 		      "#delta p_{Dz} Collider; #delta p_{Dz} (GeV) ; No.EVNT",
	// 		      200,99.5,100.5);
	TH1D *htS = new TH1D("htS", 
			     "#delta p_{Dx} Collider ; #delta p_{Dx} (GeV) ; No.EVNT",
			     200,4.9,5.1);
	TH1D *htSY = new TH1D("htSY", 
			      "#delta p_{Dy} Collider; #delta p_{Dy} (GeV) ; No.EVNT",
			      200,-0.1,0.1);
	TH1D *htSZ = new TH1D("htSZ", 
			      "#delta p_{Dz} Collider; #delta p_{Dz} (GeV) ; No.EVNT",
			      200,99.5,100.5);





	 TH1D *htSXr = new TH1D("htSXr",
	 		     "#delta p_{Dx} Collinear; #delta p_{Dx'} (GeV) ; No.EVNT",
	 		       200,-0.1,0.1);
	 TH1D *htSYr = new TH1D("htSYr", 
	 		      "#delta p_{Dy} Collinear; #delta p_{Dy'} (GeV) ; No.EVNT",
	 		      200,-0.1,0.1);
	 TH1D *htSZr = new TH1D("htSZr", 
	 		      "#delta p_{Dz} Collinear; #delta p_{Dz'} (GeV) ; No.EVNT",
	 		       200,-0.1,0.1);



	 TH1D *htSXs = new TH1D("htSXs",
	 		     " Collinear; pSpectator(#delta p_{x}) (GeV) ; No.EVNT",
	 		       200,-0.1,0.1);
	 TH1D *htSYs = new TH1D("htSYs", 
	 		      " Collinear; pSpectator(#delta p_{y}) (GeV) ; No.EVNT",
	 		      200,-0.1,0.1);
	 TH1D *htSZs = new TH1D("htSZs", 
	 		      " Collinear; pSpectator(#delta p_{z}) (GeV) ; No.EVNT",
	 		       200,-0.1,0.1);


	//TH1D *htSXr = new TH1D("htSXr",
	//		     "#delta p_{Dx} Collinear; #delta p_{Dx'} (GeV) ; No.EVNT",
	//		       200,-1.,0.);
	//TH1D *htSYr = new TH1D("htSYr", 
	//		      "#delta p_{Dy} Collinear; #delta p_{Dy'} (GeV) ; No.EVNT",
	//		      200,-0.02,0.02);
	//TH1D *htSZr = new TH1D("htSZr", 
	//			      "#delta p_{Dz} Collinear; #delta p_{Dz'} (GeV) ; No.EVNT",
	//		       200,-0.02,0.02);



	 TH1D *hnu = new TH1D("hnu", "#nu_{Rest}", 100.,-100.,100.);
	 TH1D *halphaS = new TH1D("halphaS", "#alpha_{S}", 100.,0.,2.);
	
	if (iran) {
		sig_eThx = sigma_th(kBeam, mElectron, eEpsNX, eBetaStarX);
		sig_eThy = sigma_th(kBeam, mElectron, eEpsNY, eBetaStarY);
		sig_iThx = sigma_th(PBeam, MIon, iEpsNX, iBetaStarX);
		sig_iThy = sigma_th(PBeam, MIon, iEpsNY, iBetaStarY);
	}
	double kBeamMC, kBeamMCx, kBeamMCy, kBeamMCz;
	double PBeamMC, PBeamMCx, PBeamMCy, PBeamMCz;
//	double p_RT, tempVar, pDrest, pDrest0;  DO NOT DOUBLE DECLARE VARIABLES
	double EScatRest, kScatRest, csTheRest, PhiScatRest;
	double EScatRest0, kScatRest0, csTheRest0, PhiScatRest0;

//	double eBeamPol, DBeamPol;  Declared as constants in .h

	double_t tchnX2, tchnXX2, tchnXXd;
	TVector3       UnitXLab, UnitYLab, UnitZLab;
	UnitXLab.SetXYZ(1.0,0.0,0.0);
	UnitYLab.SetXYZ(0.0,1.0,0.0);
	UnitZLab.SetXYZ(0.0,0.0,1.0);
	TVector3       UnitXRest, UnitYRest, UnitZRest;
	TVector3       UnitXqCM,  UnitYqCM,  UnitZqCM;
	TVector3       BoostCM, BoostRest, BoostRest0;
	TVector3       kScat3vec, pS_3vec;
	TVector3       kScat3vec0, pS_3vec0;
	TLorentzVector kIncident_Rest, PIncident_Rest, qVirtual_Rest, p_DT;
	TLorentzVector kIncident_Rest0, PIncident_Rest0, qVirtual_Rest0, p_DT0, p_ST0, p_ST;
	TLorentzVector kScattered_Rest, pSpectator_Rest, PX_Vector_Rest;
	TLorentzVector kScattered_Rest0, pSpectator_Rest0, PX_Vector_Rest0;
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
	//  mc.nEvts = 999;
	 mc.PhSpFct   = (Q2Max-Q2Min)*(log(xMax)-log(xMin))*4.*pi*pow((pSMax),3)/3.;
	 // Tel Aviv's high x configuration
	 //	 mc.PhSpFct   = (Q2Max-Q2Min)*(log(xMax)-log(xMin))*4.*pi*pow((pSMax+0.3),3)/3.;
	itree->Branch("MC", &mc, "nEvts/I:PhSpFct/F");

	double pS_rest, csThRecoil, phiRecoil;
	//	Int_t MEvts=0;
	Int_t MEvts=1;

	//name of output file : = "testout.txt";

	ofstream OUT ("testout.txt", ios::app);

	for (int iEvt=0; iEvt<NEvts; iEvt++) {
		Jacob = 1.0;
		if (iran) {
			kBeamMC = kBeam*ran3.Gaus(1.0,eDkOverk);
			kBeamMCx= kBeamMC*ran3.Gaus(0.0,sig_eThx);
			kBeamMCy= kBeamMC*ran3.Gaus(0.0,sig_eThy);
			kBeamMCz=-kBeamMC;  //  Angles are really Tangents
			PBeamMC = PIon*ran3.Gaus(1.0,iDPoverP);
			//PBeamMCx= PBeamMC*ran3.Gaus(0.0,sig_iThx); // IP2 configuration
			PBeamMCx= PBeamMC*ran3.Gaus(0.0,sig_iThx); // IP1 configuration
			PBeamMCy= PBeamMC*ran3.Gaus(0.0,sig_iThy);
			PBeamMCz= PBeamMC;
			
		} else {
			kBeamMC = kBeam;
			kBeamMCx= 0.0;
			kBeamMCy= 0.0;
			kBeamMCz= -kBeam;
			PBeamMC = PIon;
			PBeamMCx= 0.0;
			PBeamMCy= 0.0;
			PBeamMCz= PBeamMC;
		}

		kIncident_Vertex.SetXYZM(kBeamMCx, kBeamMCy, kBeamMCz, mElectron);
		PIncident_Vertex.SetXYZM(PBeamMCx, PBeamMCy, PBeamMCz, MIon);

		kIncident_Vertex0.SetXYZM(0.0, 0.0, kBeamMCz, mElectron);
		PIncident_Vertex0.SetXYZM(0.0, 0.0, PBeamMCz, MIon);

			//  Crossing angle MEIC or eRHIC
		PIncident_Vertex.RotateY(CrossingTheta);
		PIncident_Vertex.RotateZ(CrossingPhi);

		PIncident_Vertex0.RotateY(CrossingTheta);
		PIncident_Vertex0.RotateZ(CrossingPhi);
			//
		invts.TwoPdotk = 2.*PIncident_Vertex.Dot(kIncident_Vertex);
		invts.s_e = MIon*MIon + mElectron*mElectron + invts.TwoPdotk;
		invts.TwoPdotq = 2.*PIncident_Vertex.Dot(qVirtual_Vertex); 
                invts.s_q = MIon*MIon + invts.TwoPdotq;   


		invts.TwoPdotk0 = 2.*PIncident_Vertex0.Dot(kIncident_Vertex0);
		invts.s_e0 = MIon*MIon + mElectron*mElectron + invts.TwoPdotk0;
		invts.TwoPdotq0 = 2.*PIncident_Vertex0.Dot(qVirtual_Vertex0); 
                invts.s_q0 = MIon*MIon + invts.TwoPdotq0;   



                                                                 
			// Boost vectors from electron+Ion CM frame
		//		BoostCM = (kIncident_Vertex+PIncident_Vertex).BoostVector();
			// Boost vector from Ion rest frame
		BoostRest = PIncident_Vertex.BoostVector();
		BoostRest0 = PIncident_Vertex0.BoostVector();

		PIncident_Rest = PIncident_Vertex;
		PIncident_Rest0 = PIncident_Vertex0;

		PIncident_Rest.Boost(-BoostRest); // should result in (0.,0.,0.,MIon)
		PIncident_Rest0.Boost(-BoostRest0); // should result in (0.,0.,0.,MIon)

		kIncident_Rest = kIncident_Vertex;
		kIncident_Rest.Boost(-BoostRest);

		kIncident_Rest0 = kIncident_Vertex0;
		kIncident_Rest0.Boost(-BoostRest0);




			// Generate Q2, ln(xBj) uniformly
			//  psf *= (Q2Max-Q2Min)*(log(xMax)-log(xMin));
		uu   = ran3.Uniform();
		invts.Q2   = Q2Max*uu + Q2Min*(1.-uu);
		uu   = ran3.Uniform();
		invts.xBj  = pow(xMin,1.-uu)*pow(xMax,uu);
		invts.x_D  = invts.xBj*MProton/MIon;
		invts.y_D  = invts.Q2/(invts.x_D*invts.TwoPdotk);
		invts.x_D0  = invts.xBj*MProton/MIon;
		invts.y_D0  = invts.Q2/(invts.x_D*invts.TwoPdotk0);


		if (invts.y_D>=(1.0-2.*mElectron*MIon/invts.TwoPdotk) ) {
				// Unphysical kinematics
		    //			hRejects->Fill(invts.xBj,invts.Q2);
			continue;
		}
		if (invts.y_D>=1./(1.+invts.x_D*MIon*MIon/invts.TwoPdotk) ) {
				// Unphysical kinematics
		    //			hRejects->Fill(invts.xBj,invts.Q2);
			continue;
		}

		Jacob *=invts.xBj;   // Jacobian  dx/dlnx
							 //		psf *=2.*pi;
		PhiScatRest = pi*(2.*ran3.Uniform()-1.0);
		PhiScatRest0 = pi*(2.*ran3.Uniform()-1.0);
			//  mElectron-->0 approximation
		EScatRest = invts.TwoPdotk*(1.-invts.y_D)/(2.*MIon);
		EScatRest0 = invts.TwoPdotk0*(1.-invts.y_D0)/(2.*MIon);
		if (EScatRest<mElectron) {
				// should never happen
			printf("illegal Rest frame scattered electron energy =%10.6f \n",
				   EScatRest);
			continue;
		}
		kScatRest = sqrt(EScatRest*EScatRest-mElectron*mElectron);
		kScatRest0 = sqrt(EScatRest0*EScatRest0-mElectron*mElectron);
		csTheRest = (2.*EScatRest*kIncident_Rest.E() - invts.Q2 - 2.*mElectron*mElectron) / 
					(2.*kScatRest*kIncident_Rest.P());
		csTheRest0 = (2.*EScatRest0*kIncident_Rest0.E() - invts.Q2 - 2.*mElectron*mElectron) / 
					(2.*kScatRest0*kIncident_Rest0.P());



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


		kScat3vec0 = kScatRest0*(
				       sin(acos(csTheRest0))*(cos(PhiScatRest0)*UnitXRest
							     +sin(PhiScatRest0)*UnitYRest)
				       -csTheRest0*UnitZRest);
		kScattered_Rest0.SetVectM(kScat3vec0, mElectron);
		qVirtual_Rest0 = kIncident_Rest0-kScattered_Rest0;
		kScattered_Vertex0 = kScattered_Rest0;


		invts.pDrest = sqrt(PIncident_Rest(0)*PIncident_Rest(0)+PIncident_Rest(1)*PIncident_Rest(1)+PIncident_Rest(2)*PIncident_Rest(2));
		invts.pDrest0 = sqrt(PIncident_Rest0(0)*PIncident_Rest0(0)+PIncident_Rest0(1)*PIncident_Rest0(1)+PIncident_Rest0(2)*PIncident_Rest0(2));

		//cout << "pDrest= "<< PIncident_Rest(1) << ", pDrest0= "<<  PBeamMCy << endl;
			//		cout << "pDrest= "<< pDrest  << ", pDrest0= "<<  pDrest0 << endl;
			//cout << "pDrest= "<< sqrt(kIncident_Rest(2)*kIncident_Rest(2))  << ", pDrest0= "<< kBeamMCz  << endl;


		// Back to Lab frame
		kScattered_Vertex.Boost(BoostRest);
		qVirtual_Vertex = kIncident_Vertex-kScattered_Vertex;
		kScattered_Vertex0.Boost(BoostRest0);
		qVirtual_Vertex0 = kIncident_Vertex0-kScattered_Vertex0;


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
		pS_rest   = (pSMax)*pow(uu,1./3.); // uniform in 3p^2 dp = d(p^3)
		// Tel Aviv's high pT configuration
		//		pS_rest   = (pSMax+0.3)*pow(uu,1./3.); // uniform in 3p^2 dp = d(p^3)
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

		pSpectator_Vertex0 = pSpectator_Rest;
		pSpectator_Vertex0.Boost(BoostRest0);



		PX_Vertex = kIncident_Vertex+PIncident_Vertex
		            -kScattered_Vertex-pSpectator_Vertex;
		invts.MX2 = PX_Vertex.M2();
		invts.alphaS = ABeam*(pS_rest*csThRecoil+pSpectator_Rest.E())/MIon;
		invts.pPerpS = pS_rest*sqrt(1.-csThRecoil*csThRecoil);
		invts.pPerpZ = pS_rest*sqrt(csThRecoil*csThRecoil);






		//qMag = qVirtual_Rest.P();
		//pDotq = PIncident_Rest(0)*qVirtual_Rest(0)+PIncident_Rest(1)*qVirtual_Rest(1)+PIncident_Rest(2)*qVirtual_Rest(2);
		//p_DT = PIncident_Rest - (pDotq/qMag)*qVirtual_Rest;
		//qMag0 = qVirtual_Rest0.P();
		//pDotq0 = PIncident_Rest0(0)*qVirtual_Rest0(0)+PIncident_Rest0(1)*qVirtual_Rest0(1)+PIncident_Rest0(2)*qVirtual_Rest0(2);
		//p_DT0 = PIncident_Rest0 - (pDotq0/qMag0/qMag0)*qVirtual_Rest0;




		//qMag = qVirtual_Rest.P();
		//pDotq = PIncident_Vertex(0)*qVirtual_Rest(0)+PIncident_Vertex(1)*qVirtual_Rest(1)+PIncident_Vertex(2)*qVirtual_Rest(2);
		//p_DT = PIncident_Vertex - (pDotq/qMag/qMag)*qVirtual_Rest;
		//qMag0 = qVirtual_Rest0.P();
		//pDotq0 = PIncident_Vertex0(0)*qVirtual_Rest0(0)+PIncident_Vertex0(1)*qVirtual_Rest0(1)+PIncident_Vertex0(2)*qVirtual_Rest0(2);
		//p_DT0 = PIncident_Vertex0 -(pDotq0/qMag0/qMag0)*qVirtual_Rest0; 

		qMag = qVirtual_Rest.P();
	        pDotq = BoostRest(0)*qVirtual_Rest(0)+BoostRest(1)*qVirtual_Rest(1)+BoostRest(2)*qVirtual_Rest(2);
		p_DT = PIncident_Vertex - (pDotq/qMag/qMag)*qVirtual_Rest;
		qMag0 = qVirtual_Rest0.P();
		pDotq0 = BoostRest0(0)*qVirtual_Rest0(0)+BoostRest0(1)*qVirtual_Rest0(1)+BoostRest0(2)*qVirtual_Rest0(2);
		p_DT0 = PIncident_Vertex0 -(pDotq0/qMag0/qMag0)*qVirtual_Rest0;
		


		p_ST = pSpectator_Vertex -(pDotq/qMag/qMag)*qVirtual_Rest;
		p_ST0 = pSpectator_Vertex0 -(pDotq0/qMag0/qMag0)*qVirtual_Rest0;


		//qMag = qVirtual_Vertex.P();
		//pDotq = PIncident_Vertex(0)*qVirtual_Vertex(0)+PIncident_Vertex(1)*qVirtual_Vertex(1)+PIncident_Vertex(2)*qVirtual_Vertex(2);
		//p_DT = PIncident_Vertex - (pDotq/qMag/qMag)*qVirtual_Vertex;
		//qMag0 = qVirtual_Vertex0.P();
		//pDotq0 = PIncident_Vertex0(0)*qVirtual_Vertex0(0)+PIncident_Vertex0(1)*qVirtual_Vertex0(1)+PIncident_Vertex0(2)*qVirtual_Vertex0(2);
		//p_DT0 = PIncident_Vertex0 -(pDotq0/qMag0/qMag0)*qVirtual_Vertex0;




	              if(pow(invts.alphaS-1,2)<0.0001){
			  //			  cout << "%delta# p_DTy= "<< p_DT(0) - p_DT0(0) << endl;

			  // cout << qMag << ""<< "p_DT_X= "<< qVirtual_Rest(0) << ",  p_DT_y= "<< qVirtual_Rest(1) << ",  p_DT_z= "<< qVirtual_Rest(2) <<  ",  p_DT_xxx = "<< qVirtual_Rest.P() << endl;





		    // cout << "Deuteron:(BOOST)"<< "x= "<< BoostRest(0) << " y= "<< BoostRest(1) << " z= "<< BoostRest(2) << endl;
		    // cout << "Deuteron:Rest   "<< "x= "<< PIncident_Rest(0) << " y= "<< PIncident_Rest(1) << " z= "<< PIncident_Rest(2) << endl;
		    // cout << "Deuteron:Vertex "<< "x= "<< PIncident_Vertex(0) << " y= "<< PIncident_Vertex(1) << " z= "<< PIncident_Vertex(2) << endl;
  
		    // cout << ""<< endl;


		    // cout << "Electron:(BOOST)"<< "x= "<< BoostCM(0) << " y= "<< BoostCM(1) << " z= "<< BoostCM(2) << endl;
		    // cout << "Electron:Rest   "<< "x= "<< kIncident_Rest(0) << " y= "<< kIncident_Rest(1) << " z= "<< kIncident_Rest(2) << endl;
		    // cout << "Electron:Vertex "<< "x= "<< kIncident_Vertex(0) << " y= "<< kIncident_Vertex(1) << " z= "<< kIncident_Vertex(2) << endl;
		    
		    // cout << ""<< endl;

		    // cout << "Virtual Vertex: "<< "x= " << qVirtual_Vertex(0) <<" y= " << qVirtual_Vertex(1) <<" z= " << qVirtual_Vertex(2) << endl;
		    // cout << "Virtual Rest:   "<< "x= " << qVirtual_Rest(0) <<" y= " << qVirtual_Rest(1) <<" z= " << qVirtual_Rest(2) << endl;
		    // cout << ""<< endl;


			  //			  cout << "pSpectRest"<< "x= "<<  p_ST(0)  << endl;


		    // htS->Fill( PBeamMCx);
		    // htSY->Fill( PBeamMCy);
		    // htSZ->Fill( PBeamMCz);
		  

		    // Deuteron in Detector frame....
		    htS->Fill(PIncident_Vertex(0));
		    htSY->Fill(PIncident_Vertex(1) );
		    htSZ->Fill(PIncident_Vertex(2) );
	    
		    // Deuteron in Boost Vector from the Deuteron Rest frame....
		    // htSXr->Fill( BoostRest(0));
		    // htSYr->Fill( BoostRest(1));
		    // htSZr->Fill( BoostRest(2));
		    htSXr->Fill( p_DT(0) - p_DT0(0) );
		    htSYr->Fill( p_DT(1) - p_DT0(1) );
		    htSZr->Fill( p_DT(2) - p_DT0(2) );

		    // These are perfectly ZERO !!!  <=== Deuteron Rest frame....
		    // htSXr->Fill( PIncident_Rest(0));
		    // htSYr->Fill( PIncident_Rest(1));
		    // htSZr->Fill( PIncident_Rest(2));
		    
		    // pSpectator Vector in Rest frame (part of Collinear frame)
		    htSXs->Fill(  p_ST(0) - p_ST0(0) );
		    htSYs->Fill(  p_ST(1) - p_ST0(1) );
		    htSZs->Fill(  p_ST(2) - p_ST0(2) );
		    // htSXs->Fill(  pSpectator_Vertex(0) -pSpectator_Rest(0) );
		    // htSYs->Fill(  pSpectator_Vertex(1) -pSpectator_Rest(1) );
		    // htSZs->Fill(  pSpectator_Vertex(2) -pSpectator_Rest(2) );


		}




		hse->Fill(invts.s_e,Jacob);
		hQ2vsX->Fill(invts.xBj,invts.Q2, 1.);

		//p_R = sqrt(pSpectator_Rest(0)*pSpectator_Rest(0)+pSpectator_Rest(1)*pSpectator_Rest(1)+pSpectator_Rest(2)*pSpectator_Rest(2));
		//  p_RT = sqrt(pSpectator_Rest(0)*pSpectator_Rest(0)+pSpectator_Rest(1)*pSpectator_Rest(1));
		// this is same as pPerpS

		invts.tSpectator = MIon*MIon+MSpectator*MSpectator
		           - 2.*pSpectator_Vertex.Dot(PIncident_Vertex);
		invts.tPrime     = 2.*pSpectator_Vertex.Dot(PIncident_Vertex) - MIon*MIon;
		invts.tPrime0    = 2.*pSpectator_Vertex.Dot(PIncident_0) - MIon*MIon;

		// This relation is from Ch. Weiss' note Eq.(35) , p_R**2
		invts.tempVar = (invts.tPrime/2)*(1+invts.tPrime/(2*MIon*MIon))+MIon*MIon/4-MSpectator*MSpectator;
		
		//p_RT = sqrt(invts.tempVar-MIon*MIon/4+MSpectator*MSpectator); //completly wrong
		
		//p_RT = sqrt((invts.tPrime)/2)*0.9; 
		invts.p_RT = sqrt((invts.tPrime)/2)*1.0; 


		// tempVar is exactly same as = p_R**2 = (invts.pPerpS**2+invts.pPerpZ**2)
		//		cout << (invts.tPrime/(2*MIon*MIon))*(invts.tPrime/(2*MIon*MIon))  << endl;



		//htS->Fill(invts.tPrime-invts.tPrime0,Jacob);
		    hnu->Fill(qVirtual_Rest.E());
		    halphaS->Fill(invts.alphaS,1.);

		PX_Vertex_Rest = kIncident_Rest+PIncident_Rest-kScattered_Rest-pSpectator_Rest;



		    pex_Res = kScat3vec(0) ;
		    pey_Res = kScat3vec(1) ;
		    pez_Res = kScat3vec(2) ;
		    EeE_Res = sqrt(pex_Res*pex_Res+pey_Res*pey_Res+pez_Res*pez_Res+emass*emass);
		    vex_Res = 0.0;
		    vey_Res = 0.0;
		    vez_Res = 0.0;

		    pex_Lab = kScattered_Vertex(0) ;
		    pey_Lab = kScattered_Vertex(1) ;
		    pez_Lab = kScattered_Vertex(2) ;
		    EeE_Lab = sqrt(pex_Lab*pex_Lab+pey_Lab*pey_Lab+pez_Lab*pez_Lab+emass*emass);
		    vex_Lab = 0.0;
		    vey_Lab = 0.0;
		    vez_Lab = 0.0;

		    psx_Res = pSpectator_Rest(0) ;
		    psy_Res = pSpectator_Rest(1) ;
		    psz_Res = pSpectator_Rest(2) ;
		    EsE_Res = pSpectator_Rest.E() ;
		    vsx_Res = 0.0;
		    vsy_Res = 0.0;
		    vsz_Res = 0.0;

		    psx_Lab = pSpectator_Vertex(0) ;
		    psy_Lab = pSpectator_Vertex(1) ;
		    psz_Lab = pSpectator_Vertex(2) ;
		    EsE_Lab = pSpectator_Vertex.E() ;
		    vsx_Lab = 0.0;
		    vsy_Lab = 0.0;
		    vsz_Lab = 0.0;


		    pXx_Res = PX_Vertex_Rest(0) ;
		    pXy_Res = PX_Vertex_Rest(1) ;
		    pXz_Res = PX_Vertex_Rest(2) ;
		    EXE_Res = PX_Vertex_Rest.E() ;
		    vXx_Res = 0.0;
		    vXy_Res = 0.0;
		    vXz_Res = 0.0;

		    pXx_Lab = PX_Vertex(0) ;
		    pXy_Lab = PX_Vertex(1) ;
		    pXz_Lab = PX_Vertex(2) ;
		    EXE_Lab = PX_Vertex.E() ;
		    vXx_Lab = 0.0;
		    vXy_Lab = 0.0;
		    vXz_Lab = 0.0;

		    dummass = sqrt(EXE_Res*EXE_Res -(pXx_Res*pXx_Res+pXy_Res*pXy_Res+pXz_Res*pXz_Res));

		    //alpha_r = 2*(EsE_Res - psz_Res)/(1.875 - 0) ;

		    invts.nu = invts.Q2 / (2 * MIon * invts.xBj) ;
		    invts.W = MIon*MIon + invts.Q2/invts.xBj -invts.Q2;
		    W2 = invts.W*invts.W;
		    NumPtls =3;
		    zeros = 0.0001;
//		    invts.Jacob = Jacob;
/*   Declared as constants (=1.0)  in .h file
		    eBeamPol = 0.;
		    DBeamPol = 0.;
*/



		    if( invts.alphaS<2.0 &&  invts.alphaS>0.0 ){

			//
			//  OUTPUT FOR GEMC 
			//
			//	OUT << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  << MEvts << " \t " << NumPtls << " \t " << setprecision(8)<< scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << invts.s_q << " \t" << invts.alphaS << " \t" << invts.pPerpS << " \t"  << invts.tPrime  << " \t"  <<  Jacob << " \t"  <<  eBeamPol <<  " \t"  <<  DBeamPol << endl;

			OUT << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  << MEvts << " \t " << NumPtls << " \t " << setprecision(8)<< scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << invts.s_q << " \t" << invts.alphaS << " \t" << invts.pPerpS << " \t"  << invts.tPrime  << " \t"  <<  Jacob << " \t"  <<  invts.p_RT <<  " \t"  <<  sqrt((invts.tPrime)/2.) << " \t" <<  kBeamMC << " \t" << PBeamMC << " \t" << eBeamPol << " \t" << DBeamPol << endl;

			//			OUT << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  << MEvts << " \t " << NumPtls << " \t " << setprecision(8)<< scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << invts.s_q << " \t" << invts.alphaS << " \t" << invts.pPerpS << " \t"  << invts.tPrime  << " \t"  <<  Jacob << " \t"  <<  invts.p_RT <<  " \t"  <<  sqrt((invts.tPrime)/2.) << endl;


       OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(8) << "\t" << "1" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< scientific <<pex_Lab << " \t " << pey_Lab << " \t " << pez_Lab << " \t " << EeE_Lab << " \t " << emass << " \t " << vex_Lab  << " \t " << vey_Lab << " \t " << vez_Lab << endl; 
       OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(8) << "\t" << "1" << " \t " << sp_particle_id << " \t " << "0" <<  " \t "<< scientific<< psx_Lab << " \t " << psy_Lab << " \t " << psz_Lab << " \t " << EsE_Lab << " \t " << spmass << " \t " << vsx_Lab  << " \t " << vsy_Lab << " \t " << vsz_Lab << endl; 



		    }

// p_n = (P_D - p_S)  as a four-vector, then compute  p_n^2  = t
// compare this event-by-event with




	 tchnX2  =  (PIncident_Vertex.E()-EsE_Lab)*(PIncident_Vertex.E()-EsE_Lab) -  (PIncident_Vertex(0)-psx_Lab)*(PIncident_Vertex(0)-psx_Lab) - (PIncident_Vertex(1)-psy_Lab)*(PIncident_Vertex(1)-psy_Lab) - (PIncident_Vertex(2)-psz_Lab)*(PIncident_Vertex(2)-psz_Lab);
	 // tprimXX2 = 2.*(PIncident_Rest(0)*psx_Res + PIncident_Rest(1)*psy_Res + PIncident_Rest(2)*psz_Res) - MDeut*MDeut + MProton*MProton ;
	 tchnXX2 =  MProton*MProton - invts.tPrime ;
	 tchnXXd = tchnX2 - tchnXX2;



	 //	invts.tPrime     = 2.*pSpectator_Vertex.Dot(PIncident_Vertex) - MIon*MIon;

	 // cout << setiosflags(ios::fixed) << setprecision(16) <<  " \t "<< tchnX2 << " \t "<< tchnXX2 << " \t "<< tchnXXd << endl;

	 MEvts++;

	 //	}  // constraint of alphaS

		tree->Fill();
		//	 MEvts++;
	}

       	OUT.close();

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
