//Written by Oz Amram 6/18/2014
//
//Reads data from  3 EVTP2.OUT, EVTP3.OUT, EVTP4.OUT located in ../christians_model/
//plots the data on the same graph as well as storing it in a root file

#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <vector>

using namespace std;

string int_to_string(int i){
   //a hack to work around a bug with to_string
   //see http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
   ostringstream ss;
   ss << i;
   return ss.str();
}


string get_outfile_path(int index){
    // returns the path to the EVTP.OUT file
    // ONLY WORKS IF EVTP.OUT is in ../tag from this file
    string outfile_path = "../../christians_model/extract_f2n/EVTP"+ int_to_string(index) +".OUT";
    //cout <<  outfile_path << endl;
    return outfile_path;
}
void print_my_tree(TTree *tree){
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

void HallA_style() {
  gROOT->SetStyle("Plain");
  gStyle->SetPaperSize(TStyle::kUSLetter);
  gStyle->SetPaperSize(18,22);
  gStyle->SetOptFit(1111);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505);
  
  gStyle->SetCanvasColor(10);
  gStyle->SetPadTopMargin(.05);
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



void make_graphs(vector<Float_t> TP_vec,  vector<Float_t> FSIG1_vec, 
                 vector<Float_t> FSIG2_vec){
   HallA_style(); //make it pretty
   

   
   TCanvas *c1 = new TCanvas("c1", "Cross Sections", 700, 700);
   //make t' vs cross section graph
   TVirtualPad  *pad1 = c1->cd();
   //pad1->SetLogy(); //make Y axis on log scale
   //all the data has the same TP values so we use the same
    const int xsize = 7;
    float o2[xsize] = {4.508, 4.1147, 3.831, 3.563, 3.323, 3.109, 2.913};
    float o3[xsize] = {4.540, 4.176, 3.858, 3.589, 3.347, 3.130, 2.934};
    float o4[xsize] = {4.545, 4.178, 3.862, 3.592, 3.349, 3.133, 2.938};
    float xs[xsize] = {0.05,  0.055, 0.06,  0.065, 0.07,  0.075, 0.08};

   TGraph *o2grph = new TGraph(xsize, xs, o2);
   o2grph -> SetMarkerStyle(20);
   o2grph -> SetMarkerColor(3);
   o2grph -> SetTitle("2nd Order Fit");

   TGraph * o3grph = new TGraph(xsize, xs, o2);
   o3grph -> SetMarkerStyle(20);
   o3grph -> SetMarkerColor(4);
   o3grph -> SetTitle("3rd Order Fit");


   TGraph * o4grph = new TGraph(xsize, xs, o2);
   o4grph -> SetMarkerStyle(20);
   o4grph -> SetMarkerColor(5);
   o4grph -> SetTitle("4th Order Fit");

   TMultiGraph * xmg = new TMultiGraph();
   xmg -> Add(o2grph);
   xmg -> Add(o3grph);
   //xmg -> Add(o4grph);
   xmg -> Draw("APL");


   xmg->GetXaxis()->SetTitle("xbj");
   xmg->GetXaxis()->CenterTitle();
   xmg->GetYaxis()->SetTitle("Extracted F2N");
   xmg->GetYaxis()->CenterTitle();
   xmg-> SetTitle("Extracted F2N vs. xbj");
   //c1 ->BuildLegend();

   int size = TP_vec.size();


   TGraph * grph2 = new TGraph(size, &TP_vec[0], &FSIG2_vec[0]);
   grph2 -> SetMarkerStyle(20);
   grph2 -> SetMarkerColor(3);  
   //grph2 -> Draw("APL");

   grph2->GetXaxis()->SetTitle("t'");
   grph2->GetXaxis()->CenterTitle();
   grph2->GetYaxis()->SetTitle("Cross Section / Pole Factor");
   grph2->GetYaxis()->CenterTitle();

   //do the various fits
   float xbj = 0.08;
   const int nfits = 10;
   const float tp_start = 0.02;
   const float tp_bin = 0.0025;
   const float tp_end = 0.2;
   const int nbin = ceil((tp_end - tp_start)/tp_bin);
//CINT was giving errors that nbin wasn't const so it could be an array dimmension
   const int nbin_ = nbin; 

   float tp_slice = 0.05;
   int bin_slice = (tp_slice - tp_start)/tp_bin;

   float tp;
   int ndf_last = 0;
   Double_t tp_max[nbin_], f2n[nfits][nbin_];
   TF1 * f[nfits]; 
   TF1 * cur_f;
   for (int i = 0; i<nfits ; i++){ //set up the different order fit functions
       f[i] = new TF1(TString::Format("f%i", i), TString::Format("pol%i", i));
   }
   for (int i = 0; i<nbin; i++){
      tp = tp_start + (i * tp_bin);
      tp_max[i] = tp;
      for (int j = 0; j<nfits; j++){
          grph2 -> Fit(TString::Format("f%i", j), "QO", "", 0.0, tp);
          cur_f = f[j];
          f2n[j][i] = cur_f-> GetParameter(0);
          //cout << "Adding point: " << tp << " " << cur.back() << endl;
          
      }
   }

   TCanvas *c3 = new TCanvas("c3", "F2N vs TP fit", 700, 700);
   c3 -> cd();
   TMultiGraph * mgf = new TMultiGraph();
   mgf-> SetTitle(TString::Format("Extrapolated F2N vs Max t' in Fit (Christian's Model) Xbj = %1.4f",xbj));

   TGraph * o1_grph = new TGraph(nbin, tp_max, f2n[0]);
   o1_grph -> SetMarkerStyle(20);
   o1_grph -> SetMarkerColor(2);
   o1_grph -> SetTitle("0th Order Fit");

   TGraph * o2_grph = new TGraph(nbin, tp_max, f2n[1]);
   o2_grph -> SetMarkerStyle(20);
   o2_grph -> SetMarkerColor(3);
   o2_grph -> SetTitle("1st Order Fit");

   TGraph * o3_grph = new TGraph(nbin, tp_max, f2n[2]);
   o3_grph -> SetMarkerStyle(20);
   o3_grph -> SetMarkerColor(4);
   o3_grph -> SetTitle("2nd Order Fit");

   TGraph * o4_grph = new TGraph(nbin, tp_max, f2n[3]);
   o4_grph -> SetMarkerStyle(20);
   o4_grph -> SetMarkerColor(5);
   o4_grph -> SetTitle("3rd Order Fit");


   TGraph * o5_grph = new TGraph(nbin, tp_max, f2n[4]);
   o5_grph -> SetMarkerStyle(20);
   o5_grph -> SetMarkerColor(6);
   o5_grph -> SetTitle("4th Order Fit");


   TGraph * o6_grph = new TGraph(nbin, tp_max, f2n[5]);
   o6_grph -> SetMarkerStyle(20);
   o6_grph -> SetMarkerColor(7);
   o6_grph -> SetTitle("5th Order Fit");

   TGraph * o7_grph = new TGraph(nbin, tp_max, f2n[6]);
   o7_grph -> SetMarkerStyle(20);
   o7_grph -> SetMarkerColor(8);
   o7_grph -> SetTitle("6th Order Fit");

   mgf -> Add(o1_grph);
   mgf -> Add(o2_grph);
   mgf -> Add(o3_grph);
   mgf -> Add(o4_grph);
   mgf -> Add(o5_grph);
   //mgf -> Add(o6_grph);
   //mgf -> Add(o7_grph);
   mgf -> Draw("APL");

   mgf->GetXaxis()->SetTitle("t'max");
   mgf->GetXaxis()->CenterTitle();
   mgf->GetYaxis()->SetTitle("Extracted F2N");
   mgf->GetYaxis()->CenterTitle();

   c3 -> BuildLegend(0.2, 0.3, 0.45, 0.45);

   TCanvas *c4 = new TCanvas("c4", "Order Comparison", 700, 700);
   //choose a particular tp max and compare the F2Ns of the different orders of extractions
   Double_t comp[nfits], n[nfits];
   for (int i = 0; i<nfits; i++){
       comp[i] = f2n[i][bin_slice];
       n[i] = i;
       cout << TString::Format("Order is %i. F2N is %f", i, comp[i]) << endl;
   }

   TGraph * ncmp_grph = new TGraph(nfits, n, comp);
   ncmp_grph->SetMarkerStyle(20);
   ncmp_grph->SetMarkerColor(4);
   ncmp_grph->SetTitle(TString::Format("Fit Order Comparison t' max = %1.2f xbj = %1.4f", tp_slice, xbj));

   ncmp_grph->GetXaxis()->SetTitle("Order of Fit");
   ncmp_grph->GetXaxis()->CenterTitle();
   ncmp_grph->GetYaxis()->SetTitle("Extracted F2N");
   ncmp_grph->GetYaxis()->CenterTitle();
   ncmp_grph -> Draw("APL");



}
bool get_lines(char line[nfiles][80], FILE * out [nfiles]){
    bool is_done = true;
    for (int i = 0; i <nfiles; i++){
        is_done = is_done && fgets(line[i], 80, out[i]);
    }
    return is_done;
}


const int nfiles = 1; //global const needed to build arrays

int c_model_f2n_xvar_reader(bool print=false){

//  Add randomness to results of simulation
    bool do_rand = false;
    float percent_error = 0.01;

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
    vector<Float_t> FSIG4_vec;
    vector<Float_t> FSIG_ratio_2_3;

    //open the  output files 
    FILE * out[nfiles]; 
    for(int i = 0; i<nfiles; i++){
        int idx = i + 7; //files indexed starting from 1
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
       TRandom2 rand;
       while (get_lines(line,out)){
           first_char = line[0][0]; //if a comment should be same for all files
           if (first_char != '#' ){ //check that line isnt a comment
               //read data and store in variables
               for(int i=0; i<nfiles; i++){
                  sscanf(line[i], "%f %f %f %f %f %f", TP+i, PR2+i, PTR+i, FSIG+i, 
                                                   UNUM+i, F2_SPOL+i);

                  TP[i] = TMath::Abs(TP[i]); //c_model outputs a negative t'
                  if (print && false) {
                     cout<<"Writing values: " <<TP[i]<<" "<<PR2[i]<<" ";
                     cout <<PTR[i]<<" "<<FSIG[i]<<" "<<UNUM[i]<<" ";
                     cout<<F2_SPOL[i]<<endl;
                  }
               }
                 if (do_rand){
                     //Add add a slight gaussian error to the point instead of perfect value
                     FSIG2_vec.push_back(rand.Gaus(FSIG[0], percent_error*FSIG[1]));
                 }
                 else{
                     FSIG2_vec.push_back(FSIG[0]);
                 }
                     TP_vec.push_back(TMath::Abs(TP[0]));
                     FSIG1_vec.push_back(FSIG[0]);
                     //cout << FSIG[0] << " " <<  FSIG[1] << endl;
                  
               tree->Fill();

          }
       }
        
        
   //print_my_tree(tree);
   make_graphs(TP_vec, FSIG1_vec, FSIG2_vec);
   //write the tree to a root file 
   char * file_name_str = "c_model_x_data.root";
   TFile * root_file = 0;
   root_file = new TFile (file_name_str, "RECREATE");
   tree->Write();
   for (int i=0; i<nfiles; i++){
      fclose(out[i]);}
   if(print) cout << "Closed Files" << endl;
   delete root_file;
 
    return 1;
}


