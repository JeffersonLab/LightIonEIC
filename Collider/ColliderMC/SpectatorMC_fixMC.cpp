/*
 *  SpectatorMC.cpp
 *  
 *
 *  Created by Charles Hyde on 3/4/14.
 *  Copyright 2014. All rights reserved.
 *
 */

#include "SpectatorMC_fixMC.h"
#include <time.h>

	// Define function calls
double sigma_th(double pInc, double mInc, double NormEmit, double betaSt);
	//
	// ********************************************* main()
int mainx(){

    srand (time(NULL));
    double rnumber = rand()%6 + 12345.;
    TRandom3 ran3;
    // ran3.SetSeed(13579);
    ran3.SetSeed(rnumber);
    //	Bool_t iran=kTRUE;      // TRUE==include incident beam emittance
    Bool_t iran=kFALSE;      // Falut NO incident beam emittance
    Bool_t iProton = kTRUE; // TRUE == Spectator Proton
    double MSpectator;
    double W2;
    int NumPtls;
    int charge,e_charge;

    // spectator particle information
    double psx_Lab,psy_Lab,psz_Lab,EsE_Lab ;
    double psx_Res,psy_Res,psz_Res,EsE_Res ;
    double_t vsx_Lab,vsy_Lab,vsz_Lab ;
    double_t vsx_Res,vsy_Res,vsz_Res ;
    int sp_particle_id;
    double spmass;

    // scattered electron information
    double pex_Res,pey_Res,pez_Res,EeE_Res;
    double pex_Lab,pey_Lab,pez_Lab,EeE_Lab ;
    double_t vex_Lab,vey_Lab,vez_Lab ;
    double_t vex_Res,vey_Res,vez_Res ;
    int e_particle_id;
    double emass;
    emass =  0.000511;
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
	        spmass =0.93827;
		sp_particle_id = 2212;
	}else {
		MSpectator=MNeutron;
		charge =0;
	        spmass =0.93955;
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
	TLorentzVector kIncident_Vertex, PIncident_Vertex, qVirtual_Vertex;
	TLorentzVector pSpectator_Vertex, kScattered_Vertex, PX_Vertex, PX_Vertex_Rest;
	typedef struct{
	    Double_t s_e, s_q, Q2, xBj, W, x_D, y_D, tSpectator, tPrime, TwoPdotk, TwoPdotq, tPrime0, p_RT,
		MX2, alphaS, pPerpS, nu;
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
	    "s_e/s_q/D:Q2/D:xBj/D:nu/D:W/D:p_RT/D:x_D/D:y_D/D:tSpectator/D:tPrime/D:TwoPdotk/D:TwoPdotq/D:tPrime0/D:MX2/D:alphaS/D:pPerpS/D");
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
	// 	//	tree->Branch("hse","TH1D",&hse,bufsize,0);
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
	// 	//	tree->Branch("hQ2vsX","TH2D",&hQ2vsX,bufsize,0);
	// TH2D *hRejects = new TH2D("hRejects",
	// 						"Q^{2} vs. x_{Bj};  x_{Bj}  ; Q^{2} (GeV^{2})  ",
	// 						nxBins-1,xLowBinEdge, nQ2, 0.0,Q2Max);
		//	tree->Branch("hRejects","TH2D",&hRejects,bufsize,0);
	TH1D *htS = new TH1D("htS", 
		"t'=M_{S}^{2} - (P_{D}-p_{S})^{2} ; t' (GeV^{2})  ; Jacobian weighted events ",
			     //		     		 5000,-0.1,0.1);
							 1000,0.0,0.1);
	// 	//	tree->Branch("htS","TH1D",&htS,bufsize,0);
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
	double p_RT;
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
	TLorentzVector kScattered_Rest, pSpectator_Rest, PX_Vector_Rest;
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
	//mc.nEvts = NEvts;
	 mc.nEvts = 99999;
	mc.PhSpFct   = (Q2Max-Q2Min)*(log(xMax)-log(xMin))*4.*pi*pow(pSMax,3)/3.;
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
		invts.TwoPdotq = 2.*PIncident_Vertex.Dot(qVirtual_Vertex); 
                invts.s_q = MIon*MIon + invts.TwoPdotq;                                                                    
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
			
		p_RT = sqrt(pSpectator_Vertex(0)*pSpectator_Vertex(0)+pSpectator_Vertex(1)*pSpectator_Vertex(1));

		invts.tSpectator = MIon*MIon+MSpectator*MSpectator
		           - 2.*pSpectator_Vertex.Dot(PIncident_Vertex);
		invts.tPrime     = 2.*pSpectator_Vertex.Dot(PIncident_Vertex) - MIon*MIon;
		invts.tPrime0    = 2.*pSpectator_Vertex.Dot(PIncident_0) - MIon*MIon;
		//htS->Fill(invts.tPrime-invts.tPrime0,Jacob);
		    hnu->Fill(qVirtual_Rest.E());
		    halphaS->Fill(invts.alphaS,1.);

		PX_Vertex_Rest = kIncident_Rest+PIncident_Rest-kScattered_Rest-pSpectator_Rest;

		if(invts.alphaS>0.97&&invts.alphaS<1.03){
			htS->Fill(invts.tPrime,Jacob);

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


		    invts.nu = invts.Q2 / (2 * MIon * invts.xBj) ;
		    invts.W = MIon*MIon + invts.Q2/invts.xBj -invts.Q2;
		    W2 = invts.W*invts.W;
		    NumPtls =3;
		    zeros = 0.0001;

        // OUT << setiosflags(ios::left)  << setiosflags(ios::fixed) << setprecision(4) <<"                 "  << MEvts << " \t " << NumPtls << " \t " << scientific << W2 << " \t " << invts.Q2  << " \t " << invts.nu  << " \t " << zeros << endl;

		    OUT << setiosflags(ios::left)  << setiosflags(ios::fixed) << setprecision(3) <<"                 "  << MEvts << " \t " << NumPtls << " \t " << scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << invts.s_q << " \t" << invts.alphaS << " \t" << invts.pPerpS << " \t"  << invts.tPrime  << endl;
		    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(3) << "\t" << "1" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< scientific <<pex_Res << " \t " << pey_Res << " \t " << pez_Res << " \t " << EeE_Res << " \t " << emass << " \t " << vex_Res  << " \t " << vey_Res << " \t " << vez_Res << endl; 
		    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(3) << "\t" << "1" << " \t " << sp_particle_id << " \t " << "0" <<  " \t "<< scientific<< psx_Res << " \t " << psy_Res << " \t " << psz_Res << " \t " << EsE_Res << " \t " << spmass << " \t " << vsx_Res  << " \t " << vsy_Res << " \t " << vsz_Res << endl; 
		    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(3) << "\t" << "1" << " \t " << dum_particle_id << " \t " << "0" <<  " \t "<< scientific << pXx_Res << " \t " << pXy_Res << " \t " << pXz_Res << " \t " << EXE_Res << " \t " << dummass << " \t " << vXx_Res  << " \t " << vXy_Res << " \t " << vXz_Res << endl; 
		    


	 MEvts++;

	 	}  // constraint of alphaS

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
