//Written by Oz Amram July 2014
//amram@jlab.org
//
//This file reads cross section / pole data from christians model and does
//various methods of extracting F2N. It can do monte carlo methods with adding
//random error to look at spreads in F2Ns
//NOTE: You cannot run this program twice in the same session of root or you
//will get a seg-fault. If you need to rerun quit root and then run it again.
//NOTE:be careful when doing the model dependent extraction or blind test (mode=2 or mode=3)
//because it can take up a lot of memory if nruns is high (uses about 5GB ram for 500 runs)

#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <vector>
#include <cmath>

using namespace std;

//Define Globals
const int nfiles = 2; //global const needed to build arrays
bool do_mc = true; //Do monte carlo simulations by adding a random gaussian error to the data
int mode = 3; // 1: model independent extraction, 2: model dependent extraction, 3: blind model dependent extraction
//be careful when doing the model dependent extraction or blind test because it can take up a lot
//of memory if nruns is high (uses about 5GB ram for 500 runs)
float percent_error = 0.003;//amount of random error to add
bool single_run_graphs =false;
//make graphs corresponding to a single randomization (used to look for regions
//of good agreement). Uses the last randomization



void do_fits_mi(const int nfits, TGraph * grph, float tp_start, float tp_bin, const int nbin, 
                  Double_t  tp_max[nbin], Double_t  f2n[][nbin]){
   //MODEL INDEPENDENT VERSION
   //do the various fits to the data for different orders and maximum tp's
   //tp_max and d f2n are ouputs and the rest are inputs
   float tp;
   TF1 * f[nfits]; 
   TF1 * cur_f;
   for (int i = 0; i<nfits ; i++){ //set up the different order fit functions
       f[i] = new TF1(TString::Format("f%i", i), TString::Format("pol%i", i));
   }
   for (int i = 0; i<nbin; i++){
      tp = tp_start + (i * tp_bin);
      tp_max[i] = tp;
      for (int j = 0; j<nfits; j++){
          grph -> Fit(TString::Format("f%i", j), "QO", "", 0.0, tp); //fit the function to the data
          cur_f = f[j];
          Double_t temp = cur_f -> GetParameter(0); // get the y-intercept (F2N)
          f2n[j][i] = temp;
      }
   }
   return;
}

void do_fits_nom(const int nfits, TGraph * grph_nom, float tp_start, float tp_bin, 
        const int nbin, Double_t params[nbin][nfits][nfits]){
    //do the fits for the unrandomized data and output the value of the fitting parameters as a 3d array (params)
    //used to set shape of the plot to reduce fit to 1 parameter
    float tp;
    TF1 * f[nfits];
    TF1 * cur_f;
   for (int i = 0; i<nfits ; i++){ //set up the different order fit functions
       f[i] = new TF1(TString::Format("f_nom%i", i), TString::Format("pol%i", i));
   }
   for (int i = 0; i<nbin; i++){
      tp = tp_start + (i * tp_bin);
      for (int j = 0; j<nfits; j++){
          grph_nom -> Fit(TString::Format("f_nom%i", j), "QO", "", 0.0, tp);
          cur_f = f[j];
          for (int k = 0; k <= j; k++){//store all the coefficients from the fit
              Double_t temp = cur_f -> GetParameter(k);
              params[i][j][k] = temp;
         }
      }
   }
   return;
    
    

}
void do_fits_md(const int nfits, TGraph * grph_rand, float tp_start, float tp_bin, 
       const int nbin, Double_t params[][nfits][nfits], 
       Double_t  tp_max[], Double_t  f2n[nfits][nbin]){
   //MODEL DEPENDENT VERSION
   //do the various fits to the data for different orders and maximum tp's
   //uses params to constrain all the parameters of the fit but 1
   //tp_max and d f2n are ouputs and the rest are inputs
   float tp;
   vector<vector<TF1 *>>  f;
   f.resize(nbin);
   Double_t temp;
   char * f_str[7];
   //setup different order functions
   //must be manually specified to set the parameters but they are the same as
   //the pol2, pol3 ect functions used before
   f_str[0] = "[0]";
   f_str[1] = "[0] + [1] * x";
   f_str[2] = "[0] + [1] * x + [2] * pow(x,2)";
   f_str[3] = "[0] + [1] * x + [2] * pow(x,2) + [3] * pow(x,3)";
   f_str[4] = "[0] + [1] * x + [2] * pow(x,2) + [3] * pow(x,3) + [4] * pow(x,4)";
   f_str[5] = "[0] + [1] * x + [2] * pow(x,2) + [3] * pow(x, 3) + [4] * pow(x,4) + [5] * pow(x,5)";
   f_str[6] = "[0] + [1] * x + [2] * pow(x,2) + [3] * pow(x, 3) + [4] * pow(x,4) + [5] * pow(x,5) + [6] * pow(x,6)";
   for (int i = 0; i<nbin ; i++){ //set up the different order fit functions
       f[i].resize(nfits);
       for (int j = 0; j<nfits; j++){
            f[i][j] = new TF1(TString::Format("f_rand%i%i", i, j), f_str[j]);
            for (int k = 1; k<=j; k++){ //note that we start with k=1 because want 0 parameter free
                temp = params[i][j][k];
                f[i][j] -> FixParameter(k, temp);
            } 
        }
   }
   
   for (int i = 0; i<nbin; i++){
      tp = tp_start + (i * tp_bin);
      tp_max[i] = tp;
      for (int j = 0; j<nfits; j++){
          grph_rand -> Fit(TString::Format("f_rand%i%i", i,j), "QO", "", 0.0, tp);
          f2n[j][i] = f[i][j] -> GetParameter(0);
      }
   }

   return;
   
}

void do_fits_blind(const int nfits, TGraph * grph_rand, float tp_start, float tp_bin, 
       const int nbin, Double_t params[][nfits][nfits], 
       Double_t  tp_max[], Double_t  f2n[nfits][nbin]){
   //BLIND TEST VERSION
   //do the various fits to the data for different orders and maximum tp's
   //uses params to constrain the shape of the plot and only fit scaling ratio
   //(should be same as ratio of F2N's used)
   //tp_max and d f2n are ouputs and the rest are inputs
   float tp;
   vector<vector<TF1 *>>  f;
   f.resize(nbin);
   Double_t temp;
   char * f_str[7];
   Double_t scaling;
   Double_t f2n_val;
   //setup different order fit functions. param [0] is the scaling ratio
   //the other params will be fixed
   f_str[0] = "[0] * ([1])";
   f_str[1] = "[0] *([1] + [2] * x)";
   f_str[2] = "[0] * ([1] + [2] * x + [3] * pow(x,2))";
   f_str[3] = "[0] * ([1] + [2] * x + [3] * pow(x,2) + [4] * pow(x,3))";
   f_str[4] = "[0] * ([1] + [2] * x + [3] * pow(x,2) + [4] * pow(x,3) +  [5] * pow(x,4))";
   f_str[5] = "[0] * ([1] + [2] * x + [3] * pow(x,2) + [4] * pow(x, 3) + [5] * pow(x,4) + [6] * pow(x,5))";
   f_str[6] = "[0] * ([1] + [2] * x + [3] * pow(x,2) + [4] * pow(x, 3) + [5] * pow(x,4) + [6] * pow(x,5) + [7] * pow(x,6))";
   for (int i = 0; i<nbin ; i++){ //set up the different order fit functions
       f[i].resize(nfits);
       for (int j = 0; j<nfits; j++){
            f[i][j] = new TF1(TString::Format("f_rand%i%i", i, j), f_str[j]);
            for (int k = 0; k <= j; k++){ //fix all params but scaling ratio ([0])
                temp = params[i][j][k];
                f[i][j] -> FixParameter(k+1, temp);
                //the index of the correct param is now 1 different because we
                //have to use param [0] as the scaling ratio
            } 
       } 
   }
   
   for (int i = 0; i<nbin; i++){
      tp = tp_start + (i * tp_bin);
      tp_max[i] = tp;
      for (int j = 0; j<nfits; j++){
          grph_rand -> Fit(TString::Format("f_rand%i%i", i,j), "QO", "", 0.0, tp);
          scaling = f[i][j] -> GetParameter(0); //scaling variable is 0 (fitted value)
          f2n_val = scaling * (f[i][j] -> GetParameter(1)); // scaling * intial F2N (set) is extracted F2N
          f2n[j][i] = f2n_val;
      }
   }

   return;
   
}
bool do_slice(const int nfits, int bin_slice, const int nbin, Double_t f2n[][nbin],
              Double_t slice[nfits], Double_t n[nfits], 
              vector<Double_t> &extracted_f2n, vector<Double_t> &errors){
   //choose a particular tp max and compare the F2Ns of the different orders of extractions
   //slice, n are outputs. extracted_f2n and errors get modified to include this iteration 
   //and the rest are inputs
   Double_t last = 0;
   for (int i = 0; i<nfits; i++){
       slice[i] = f2n[i][bin_slice];
       n[i] = i;
       if ((slice[i] !=0) && i==4){
           Double_t error = 100.0 * abs((slice[i] - last)/slice[i]);
           extracted_f2n.push_back(slice[i-1]); //take the lower order in the convergence
           errors.push_back(error);
           if (error <= 1.0) continue; 
           else {
               //printf("Convergence Failure between %ith and %ith order. F2N is %f and %f. There was a difference of %1.3f%% \n", 
               //                        i-1, i, slice[i-1], slice[i], error);
               return false;
           }
       }
       last = slice[i];
   }
   return true;

}

void make_errors(const vector<Float_t> FSIG, 
                vector<Float_t> &XError, vector<Float_t> &YError){
    for(int i = 0; i < FSIG.size(); i++){
        //create error bars for data for graphs
        XError.push_back(0);
        YError.push_back(percent_error * FSIG[i]);
    }

}

void randomize(vector<Float_t> &FSIG, vector<Float_t> &FSIG_rand){
    //output a vector with random gaussian error added to the values of the
    //input
    TRandom3 rand;
    rand.SetSeed(0);
    for (int i = 0; i<FSIG.size(); i++){
         FSIG_rand.push_back(rand.Gaus(FSIG[i], percent_error*FSIG[i]));
     }
     return;
}

void HallA_style() {
  gROOT->SetStyle("Plain");
  gStyle->SetPaperSize(TStyle::kUSLetter);
  gStyle->SetPaperSize(18,22);
  gStyle->SetOptFit(1111);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505);
  
  gStyle->SetCanvasColor(10);
  gStyle->SetPadTopMargin(.15);
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadRightMargin(.1);
  gStyle->SetPadBottomMargin(.15);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");

  // prepare gStyle to be useful
  //   1 = solid
  //   2 = long dash (30 10)
  //   3 = dotted (4 8)
  //   4 = dash-dot (15 12 4 12)
  //   5 = short dash ( 15 15 )
  //   6 = dash-dot-dot 
  gStyle->SetLineStyleString(1,"[]");
  gStyle->SetLineStyleString(2,"[30 10]");
  gStyle->SetLineStyleString(3,"[4 8]");
  gStyle->SetLineStyleString(4,"[15 12 4 12]");
  gStyle->SetLineStyleString(5,"[15 15]");
  gStyle->SetLineStyleString(6,"[15 12 4 12 4 12]");
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");
  gStyle->SetNdivisions(505,"Y");
  gStyle->SetOptDate(0);
  gStyle->SetDateY(.98);
  gStyle->SetStripDecimals(kFALSE);
}

void analyze_data(vector<Float_t> TP_vec,  vector<Float_t> FSIG1_vec, vector<Float_t> FSIG2_vec,
                 vector<Float_t> XErrors, vector<Float_t> YErrors){
   HallA_style(); //make it pretty
   int size = TP_vec.size();




/*
 * Uncomment this if you want a graph showing the crossection/pole with the
 * last fit performed shown
 
   //all the data has the same TP values so we use the same
   TGraph * grph2 = new TGraph(size, &TP_vec[0], &FSIG1_vec[0]);
   grph2 -> SetMarkerStyle(21);
   grph2 -> SetMarkerColor(3);  
   grph2 -> SetTitle("Extrapolation");
   TCanvas *c1 = new TCanvas("c1", "Cross Sections", 700, 700);
   TVirtualPad  *pad1 = c1->cd();
   //pad1->SetLogy(); //make Y axis on log scale

   grph2 -> Draw("ALP");
   grph2->GetXaxis()->SetTitle("t' (GEV^2)");
   grph2->GetXaxis()->CenterTitle();
   grph2->GetYaxis()->SetTitle("Cross Section / Pole Factor");
   grph2->GetYaxis()->CenterTitle();
   grph2 -> SetMaximum(2.8);

*/


   //do the various fits
   const int nfits = 5; //highest order fit is 1 less because of 0th order fit
   const float tp_start = 0.050; //minimum tp range to start performing fits
   const float tp_bin = 0.005; //how often to perform a new fit
   const float tp_end = 0.15; //end point for doing fits
   const int nbin = ceil((tp_end - tp_start)/tp_bin); //how many fits will be done
   const int nbin_ = nbin; //CINT didnt think nbin was const for some reason so this is a work around

   float tp_slice = 0.075; //what tp range of fit to look at specifically with monte carlo simmulations
   int bin_slice = (tp_slice - tp_start)/tp_bin; //index for this range of git

   int nruns = 500; //how many iterations to perform


   Double_t slice[nfits], n[nfits];
   Double_t tp_max[nbin_], f2n[nfits][nbin_];


   if (do_mc){
       if (mode == 1){ //model independent extraction
           vector<Double_t> extracted_f2n, error;
           vector<Float_t> *  FSIG_rand = 0;
           int success = 0;
           int failure = 0;
           for (int i = 0; i<nruns; i++){
                FSIG_rand = new vector<Float_t>;
                randomize(FSIG1_vec, *FSIG_rand); //randomize the data
                TGraph * grph1 = new TGraphErrors(size, &TP_vec[0], &(*FSIG_rand)[0], &XErrors[0], &YErrors[0]);
                do_fits_mi(nfits, grph1, tp_start, tp_bin, nbin, tp_max, f2n);//do the fits
                bool worked = do_slice(nfits, bin_slice, nbin, f2n, slice, n, extracted_f2n, error);//look at the slice we want
                if (worked) success ++;
                else failure++;
                delete FSIG_rand;
           }
        }
       else if (mode == 2) {//model dependent extraction
           vector<Double_t> extracted_f2n, error;
           TGraph *grph_nom = new TGraph(size, &TP_vec[0], &FSIG1_vec[0]); 
           Double_t params[nbin_][nfits][nfits];
           do_fits_nom(nfits, grph_nom, tp_start, tp_bin, nbin, params); //do fits for "model" data
           vector<Float_t> *  FSIG_rand = 0;
           int success = 0;
           int failure = 0;
           for (int i = 0; i<nruns; i++){
                FSIG_rand = new vector<Float_t>;
                randomize(FSIG1_vec, *FSIG_rand); //randomize the data
                TGraph * grph_rand = new TGraphErrors(size, &TP_vec[0], &(*FSIG_rand)[0], &XErrors[0], &YErrors[0]);
                do_fits_md(nfits, grph_rand, tp_start, tp_bin, nbin_, params, tp_max, f2n); //do constrained fits
                bool worked = do_slice(nfits, bin_slice, nbin, f2n, slice, n, extracted_f2n, error); //look at the slice we want
                if (worked) success ++;
                else failure++;
                delete FSIG_rand;
           }

       }
       else if (mode ==3){//blind test
           vector<Double_t> extracted_f2n, error;
           TGraph *grph_nom = new TGraph(size, &TP_vec[0], &FSIG1_vec[0]);//use 1st data set for "model" data
           Double_t params[nbin_][nfits][nfits];
           do_fits_nom(nfits, grph_nom, tp_start, tp_bin, nbin, params);
           vector<Float_t> *  FSIG_rand = 0;
           int success = 0;
           int failure = 0;
           for (int i = 0; i<nruns; i++){
                FSIG_rand = new vector<Float_t>;
                randomize(FSIG2_vec, *FSIG_rand);//randomize the 2nd data set (the real data with unknown F2N)
                TGraph * grph_rand = new TGraphErrors(size, &TP_vec[0], &(*FSIG_rand)[0], &XErrors[0], &YErrors[0]);
                do_fits_blind(nfits, grph_rand, tp_start, tp_bin, nbin_, params, tp_max, f2n);//do the fits of scaling ratio
                bool worked = do_slice(nfits, bin_slice, nbin, f2n, slice, n, extracted_f2n, error); //look at the slice we want
                if (worked) success ++;
                else failure++;
                delete FSIG_rand;
           }
    
       }

       printf("We were under 1%% error between 3rd and 4th order  %2.1f%% of the time for t' max = %1.3f \n", success/10.0, tp_slice);
       //draw histograms of extraced F2N and error
       //I changed around the titles of these plots many times
       TCanvas *hc1 = new TCanvas("hc1", "Extracted F2N", 700, 700);
       TCanvas *hc2 = new TCanvas("hc2", "Difference between 3rd and 4th order", 700, 700);
       TH1D * f2n_hist = new TH1D("F2N", 
              TString::Format("#splitline{Blind Test Extracted F2N}{Random Error = %1.2f%%}", 
              100*percent_error),
              60, -0.02, 0.02); //you will have to change bin range and number of bins depending on what you are doing
       TH1D * error_hist = new TH1D("3rd vs 4th Order", 
              TString::Format("#splitline{Percentage difference between 3rd and 4th order fits:}{t' max = %1.3f. Random Error = %1.2f%%}", 
              tp_slice,100* percent_error),
              30, 0.0, 2.0);//this should most likely be fine the way it is
       for (int i = 0; i < nruns; i++) {
            float mean = 3.25021; //set to the mean value of extracted F2N ('true' value)
            f2n_hist->Fill(extracted_f2n[i]); //use if you want to look at the actual value of F2N extracted
            //f2n_hist->Fill((mean - extracted_f2n[i])/mean); //use if you want to just look at the spread in extracted F2N
            error_hist->Fill(error[i]);
       }
       hc1 -> cd();
       f2n_hist -> SetFillColor(45);
       f2n_hist -> Draw();
       f2n_hist -> GetXaxis() -> SetTitle("Extracted F2N");
       cout << f2n_hist -> GetMean() << endl;//print the extracted mean (more sig figs than plot)

       hc2 -> cd();
       error_hist -> SetFillColor(35);
       error_hist -> Draw();
       error_hist -> GetXaxis() -> SetTitle("Percentage Difference");

       
   }
   else{
        vector<Double_t> extracted_f2n, error;
        TGraph * grph1 = new TGraphErrors(size, &TP_vec[0], &FSIG1_vec[0], &XErrors[0], &YErrors[0]);
        do_fits_mi(nfits, grph2, tp_start, tp_bin, nbin, tp_max, f2n);
        do_slice(nfits, bin_slice, nbin, f2n, slice, n, extracted_f2n, error);
   }

   if (single_run_graphs){ //draw graphs corresponding to the last run of randomization or non-randomized data
       TCanvas *c3 = new TCanvas("c3", "F2N vs TP fit", 700, 700);
       c3 -> cd();
       TMultiGraph * mgf = new TMultiGraph();
       mgf-> SetTitle("Extrapolated F2N vs Max t' in Fit");//an overview plot for different orders of fit

       TGraph * o0_grph = new TGraph(nbin, tp_max, f2n[0]);
       o0_grph -> SetMarkerStyle(20);
       o0_grph -> SetMarkerColor(2);
       o0_grph -> SetTitle("0th Order Fit");

       TGraph * o1_grph = new TGraph(nbin, tp_max, f2n[1]);
       o1_grph -> SetMarkerStyle(20);
       o1_grph -> SetMarkerColor(3);
       o1_grph -> SetTitle("1st Order Fit");

       TGraph * o2_grph = new TGraph(nbin, tp_max, f2n[2]);
       o2_grph -> SetMarkerStyle(20);
       o2_grph -> SetMarkerColor(4);
       o2_grph -> SetTitle("2nd Order Fit");

       TGraph * o3_grph = new TGraph(nbin, tp_max, f2n[3]);
       o3_grph -> SetMarkerStyle(20);
       o3_grph -> SetMarkerColor(5);
       o3_grph -> SetTitle("3rd Order Fit");

       TGraph * o4_grph = new TGraph(nbin, tp_max, f2n[4]);
       o4_grph -> SetMarkerStyle(20);
       o4_grph -> SetMarkerColor(6);
       o4_grph -> SetTitle("4th Order Fit");

    /*
     * only include up to the max order you want to show. Must be less than
     * nfits - 1
     
       TGraph * o5_grph = new TGraph(nbin, tp_max, f2n[5]);
       o5_grph -> SetMarkerStyle(20);
       o5_grph -> SetMarkerColor(7);
       o5_grph -> SetTitle("5th Order Fit");

       TGraph * o6_grph = new TGraph(nbin, tp_max, f2n[6]);
       o6_grph -> SetMarkerStyle(20);
       o6_grph -> SetMarkerColor(8);
       o6_grph -> SetTitle("6th Order Fit");
    */
       //mgf -> Add(o0_grph);
       mgf -> Add(o1_grph);
       mgf -> Add(o2_grph);
       mgf -> Add(o3_grph);
       mgf -> Add(o4_grph);
       //mgf -> Add(o5_grph);
       //mgf -> Add(o6_grph);
       mgf -> Draw("APL");
       c3 -> Update();

       //if you want to keep scale consistant
       //mgf-> SetMaximum(2.6);
       //mgf -> SetMinimum(2.0);

       mgf->GetXaxis()->SetTitle("t'max (GeV^2)");
       mgf->GetXaxis()->CenterTitle();
       mgf->GetYaxis()->SetTitle("Extracted F2N");
       mgf->GetYaxis()->CenterTitle();

       c3 -> BuildLegend(0.16, 0.16, 0.35, 0.30);

       TCanvas *c4 = new TCanvas("c4", "Order Comparison", 700, 700);

       TGraph * ncmp_grph = new TGraph(nfits, n, slice); //look at slice of t' range and compare all the orders of fit
       ncmp_grph->SetMarkerStyle(20);
       ncmp_grph->SetMarkerColor(4);
       ncmp_grph->SetTitle(TString::Format("Fit Order Comparison t' max = %1.4f", tp_slice));

       ncmp_grph->GetXaxis()->SetTitle("Order of Fit");
       ncmp_grph->GetXaxis()->CenterTitle();
       ncmp_grph->GetYaxis()->SetTitle("Extracted F2N");
       ncmp_grph->GetYaxis()->CenterTitle();
       ncmp_grph -> Draw("APL");

   }

}

bool get_lines(char line[nfiles][80], FILE * out [nfiles]){
    //read the next line of output and store it in the line array
    //returns true if there was another line
    //returns false if one of the files did not have another line of data
    bool is_done = true;
    for (int i = 0; i <nfiles; i++){
        is_done = is_done && fgets(line[i], 80, out[i]);
    }
    return is_done;
}



string int_to_string(int i){
   //a hack to work around a bug with to_string
   //see http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
   ostringstream ss;
   ss << i;
   return ss.str();
}


string get_outfile_path(int index){
    // returns the path to the EVTP.OUT file
    string outfile_path = "../../christians_model/extract_f2n/EVTP"+ int_to_string(index) +".OUT";
    //cout <<  outfile_path << endl;
    return outfile_path;
}
void print_my_tree(TTree *tree){
   //helpful for debugging if you want to see what you ended up storing
   //in your tree
   Float_t tp_val = 0.0;
   Float_t fsig=0.0;
   TBranch *tp_branch = 0;
   TBranch *fsig_branch =0; 

   tree->SetBranchAddress("TP2", &tp_val, &tp_branch);
   tree ->SetBranchAddress("FSIG1", &fsig, &fsig_branch);
   Int_t entries = tree->GetEntries();
   for(int i=0; i<entries; i++){
        tp_branch->GetEntry(i);
        fsig_branch->GetEntry(i);
        cout << tp_val <<" , "<< fsig<<endl;
    }

}


int c_model_extraction_mc(bool print=false){


  //define values we will read
    Float_t TP[nfiles];
    Float_t PR2[nfiles];
    Float_t PTR[nfiles];
    Float_t FSIG[nfiles];
    Float_t UNUM[nfiles];
    Float_t F2_SPOL[nfiles];

    vector<Float_t> TP_vec;
    vector<Float_t> FSIG1_vec;
    vector<Float_t> FSIG2_vec;
    vector<Float_t> FSIG_rand;
    vector<Float_t> YErrors;
    vector<Float_t> XErrors;
    //open the  output files 
    FILE * out[nfiles]; 
    for(int i = 0; i<nfiles; i++){
        int idx = i + 1; //files indexed starting from 1
        string out_path = get_outfile_path(idx);         
         out[i] = fopen(out_path.c_str(), "r"); 
         if (print) cout << "Opened File at " << out_path << endl;
    }
        //make the tree and setup the branches
    TTree *tree = new TTree("T", "tag_data_tree");
    //set all the branches
    for(int i = 0; i<nfiles; i++){
        int idx = i + 1; //files indexed starting from 1

       tree -> Branch(TString::Format("TP%i",idx), TP+i, 
                      TString::Format("TP%i/F",idx));
       tree -> Branch(TString::Format("PR2%i",idx), PR2+i, 
                      TString::Format("PR2%i/F", idx));
       tree -> Branch(TString::Format("PTR%i", idx), PTR+i,
                      TString::Format("PTR%i/F",idx));
       tree -> Branch(TString::Format("FSIG%i", idx), FSIG+i,
                      TString::Format("FSIG%i/F", idx));
       tree -> Branch(TString::Format("UNUM%i", idx), UNUM+i, 
                      TString::Format("UNUM%i/F", idx));
       tree -> Branch(TString::Format("F2_SPOL%i",idx), F2_SPOL+i,
                      TString::Format("F2_SPOL%i/F", idx));
        }
       char line[nfiles][80];
       //start reading the files
       if (print) cout << "Reading the files... " << endl;
       char first_char;
       while (get_lines(line,out)){
           first_char = line[0][0]; //if a comment should be same for all files
           if (first_char != '#' ){ //check that line isnt a comment
               //read data and store in variables
               for(int i=0; i<nfiles; i++){
                  sscanf(line[i], "%f %f %f %f %f %f", TP+i, PR2+i, PTR+i, FSIG+i, 
                                                   UNUM+i, F2_SPOL+i);
                  //read values from string and store in arrays

                  TP[i] = TMath::Abs(TP[i]); //c_model outputs a negative t'
                  if (print) {
                     cout<<"Writing values: " <<TP[i]<<" "<<PR2[i]<<" ";
                     cout <<PTR[i]<<" "<<FSIG[i]<<" "<<UNUM[i]<<" ";
                     cout<<F2_SPOL[i]<<endl;
                  }
               }
                 
                     FSIG1_vec.push_back(FSIG[0]);
                     FSIG2_vec.push_back(FSIG[1]); //remove this line if nfiles is 1
                     TP_vec.push_back(TMath::Abs(TP[0]));
                     //cout << FSIG[0] << " " <<  FSIG[1] << endl;
                  
               tree->Fill();
          }
       }
        
        
   make_errors(FSIG1_vec, XErrors, YErrors);//make error bars
   analyze_data(TP_vec, FSIG1_vec, FSIG2_vec, XErrors, YErrors); //do fits and graph it

   //write the tree to a root file 
   char * file_name_str = "c_model_x_data.root";
   TFile * root_file = 0;
   root_file = new TFile (file_name_str, "RECREATE");
   tree->Write();

   //close the files we used
   for (int i=0; i<nfiles; i++){
      fclose(out[i]);}
   if(print) cout << "Closed Files" << endl;
   delete root_file;
   return 1;
}


