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
    string outfile_path = "../../christians_model/var_f2n/EVTP"+ int_to_string(index) +".OUT";
    //cout <<  outfile_path << endl;
    return outfile_path;
}
void print_my_tree(TTree *tree){
   Float_t tp_val = 0.0;
   Float_t fsig=0.0;
   TBranch *tp_branch = 0;
   TBranch *fsig_branch =0; 

   tree->SetBranchAddress("TP2", &tp_val, &tp_branch);
   tree ->SetBranchAddress("FSIG2", &fsig, &fsig_branch);
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



void make_graphs(vector<Float_t> TP_vec,  vector<Float_t> FSIG2_vec, 
                 vector<Float_t> FSIG3_vec, vector<Float_t> FSIG4_vec){
   HallA_style(); //make it pretty
   TCanvas *c1 = new TCanvas("c1", "Cross Sections", 700, 700);
   //make t' vs cross section graph
   TVirtualPad  *pad1 = c1->cd();
   //pad1->SetLogy(); //make Y axis on log scale
   //all the data has the same TP values so we use the same
   int size = TP_vec.size();
   TGraph * grph1 = new TGraph(size, &TP_vec[0], &FSIG2_vec[0]);
   grph1 -> SetMarkerStyle(20);
   grph1 -> SetMarkerColor(2);
   grph1 -> SetTitle("NOMINAL");
   TGraph * grph2 = new TGraph(size, &TP_vec[0], &FSIG3_vec[0]);
   grph2 -> SetMarkerStyle(20);
   grph2 -> SetMarkerColor(3);  
   grph2 -> SetTitle("NOMINAL/SPOL");
/*   TGraph * grph3 = new TGraph(size, &TP_vec[0], &FSIG4_vec[0]);
   grph3 -> SetMarkerStyle(20);
   grph3 -> SetMarkerColor(4);  
   grph3 -> SetTitle("MOTT * SPOL");
*/

   //do the various fits
   vector<Double_t> tp_max, f2n_1o, f2n_2o, f2n_3o;
   TF1 *f1 = new TF1("f1", "pol1");
   TF2 *f2 = new TF1("f2", "pol2");
   TF3 *f3 = new TF1("f3", "pol3");
   float tp_start = 0.05;
   float tp_bin = 0.005;
   float tp_end = 0.5;
   int nbin = ceil((tp_end - tp_start)/tp_bin);
   float tp;
   int ndf_last = 0;
   int ndf;
   char * f_str;
   TF1 *f;
   vector<Double_t> * cur;
   for (int i = 0; i<nbin; i++){
      tp = tp_start + (i * tp_bin);
      tp_max.push_back(tp);
      for (int j = 0; j<3; j++){
          switch(j){
              case(0): //first order fit
                f = f1;
                f_str = "f1";
                cur = &f2n_1o;
                break;
              case(1): //2nd order fit
                f = f2;
                f_str = "f2";
                cur = &f2n_2o;
                break;
              case(2): //3rd order fit
                f = f3;
                f_str = "f3";
                cur = &f2n_3o;
                break;
           }
          grph2 -> Fit(f_str, "QO", "", 0.0, tp);
          ndf = f -> GetNDF();
          if (ndf == ndf_last) continue; //didn't include any new points so fit is redundant
          else {
              ndf_last = ndf;
              cur->push_back(f->GetParameter(0));
              //cout << "Adding point: " << tp << " " << cur.back() << endl;
          }
      }
   }

   TMultiGraph * mgr = new TMultiGraph();
   mgr -> SetTitle("Normalization: Christian's Model");
   //mgr -> Add(grph1);
   mgr -> Add(grph2);
   //mgr -> Add(grph3);
   mgr -> Draw("APL");

 
   mgr->GetXaxis()->SetTitle("t'");
   mgr->GetXaxis()->CenterTitle();
   mgr->GetYaxis()->SetTitle("Cross Section / Pole Factor");
   mgr->GetYaxis()->CenterTitle();

   c1->Update();

   TCanvas *c3 = new TCanvas("c3", "F2N vs TP fit", 700, 700);
   c3 -> cd();
   TMultiGraph * mgf = new TMultiGraph();
   mgf-> SetTitle("Extrapolated F2N vs Max t' in Fit (Christian's Model)");

   TGraph * o1_grph = new TGraph(tp_max.size(), &tp_max[0], &f2n_1o[0]);
   o1_grph -> SetMarkerStyle(20);
   o1_grph -> SetMarkerColor(2);
   o1_grph -> SetTitle("1st Order Fit");

   TGraph * o2_grph = new TGraph(tp_max.size(), &tp_max[0], &f2n_2o[0]);
   o2_grph -> SetMarkerStyle(20);
   o2_grph -> SetMarkerColor(3);
   o2_grph -> SetTitle("2nd Order Fit");

   TGraph * o3_grph = new TGraph(tp_max.size(), &tp_max[0], &f2n_3o[0]);
   o3_grph -> SetMarkerStyle(20);
   o3_grph -> SetMarkerColor(4);
   o3_grph -> SetTitle("3rd Order Fit");



   mgf -> Add(o1_grph);
   mgf -> Add(o2_grph);
   mgf -> Add(o3_grph);
   mgf -> Draw("APL");

   c3 -> BuildLegend(0.2, 0.3, 0.45, 0.45);

//   TCanvas *c2 = new TCanvas("c2", "Ratio", 700, 700);
//   c2->cd();
//   TGraph * rgrph = new TGraph(size, &TP_vec[0], &FSIG_ratios_2_3[0]);
//   rgrph -> SetMarkerSize(20);
//   rgrph -> SetMarkerColor(3);
//   rgrph -> SetTitle("SPOL");
//   rgrph -> Draw("APL");



}
bool get_lines(char line[nfiles][80], FILE * out [nfiles]){
    bool is_done = true;
    for (int i = 0; i <nfiles; i++){
        is_done = is_done && fgets(line[i], 80, out[i]);
    }
    return is_done;
}

const int nfiles = 2;

int c_model_f2n_reader(bool print=false){

  //define values we will read
    Float_t TP[nfiles];
    Float_t PR2[nfiles];
    Float_t PTR[nfiles];
    Float_t FSIG[nfiles];
    Float_t UNUM[nfiles];
    Float_t F2_SPOL[nfiles];

    vector<Float_t> TP_vec;
    vector<Float_t> FSIG2_vec;
    vector<Float_t> FSIG3_vec;
    vector<Float_t> FSIG4_vec;
    vector<Float_t> FSIG_ratio_2_3;

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

                  TP[i] = TMath::Abs(TP[i]); //c_model outputs a negative t'
                  if (print) {
                     cout<<"Writing values: " <<TP[i]<<" "<<PR2[i]<<" ";
                     cout <<PTR[i]<<" "<<FSIG[i]<<" "<<UNUM[i]<<" ";
                     cout<<F2_SPOL[i]<<endl;
                  }
               }
                 if (FSIG[0] != 0){
                     //put the values in the arrays too
                    
                     TP_vec.push_back(TMath::Abs(TP[0]));
                     FSIG2_vec.push_back(FSIG[0]);
                     FSIG3_vec.push_back(FSIG[1]);
                     FSIG4_vec.push_back(FSIG[2]);
                     //cout << FSIG[0] << " " <<  FSIG[1] << endl;
                  }
               tree->Fill();

          }
       }
        
        
   //print_my_tree(tree);
   make_graphs(TP_vec, FSIG2_vec, FSIG3_vec, FSIG4_vec);
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


