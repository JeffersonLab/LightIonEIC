//Written by Oz Amram 6/18/2014
//
//Reads data from  MISAK_DATA_1.OUT to MISAK_DATA_9.OUT
//They must be located in ../misak_model/var_x
//plots the data on the same graph as well as storing it in a root file

#include <stdio.h>
#include <math.h>
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
   string outfile_path = "../../misak_model/var_f2n/MISAK_DATA_"
                                                     +  int_to_string(index)
                                                     +   ".OUT";
    return outfile_path;
}
void print_my_tree(TTree *tree){
   Float_t tp_val = 0.0;
   Float_t fsig=0.0;
   TBranch *tp_branch = 0;
   TBranch *fsig_branch =0; 

   tree->SetBranchAddress("prt1", &tp_val, &tp_branch);
   tree ->SetBranchAddress("si1", &fsig, &fsig_branch);
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



void make_graphs(vector<Float_t> tp_vec,  vector<Float_t> si1_vec, 
                 vector<Float_t> si2_vec){
   HallA_style(); //make it pretty
   TCanvas *c1 = new TCanvas("c1", "Cross Sections", 700, 700);
   //make  vs cross section graph
   TVirtualPad  *pad1 = c1->cd();
   //pad1->SetLogy(); //make Y axis on log scale
   //all the data has the same tp values so we use the same vec
   int size = tp_vec.size();
   TGraph * grph1 = new TGraph(size, &tp_vec[0], &si1_vec[0]);
   grph1 -> SetMarkerStyle(20);
   grph1 -> SetMarkerColor(2);
   grph1 -> SetTitle("NOMINAL");
   TGraph * grph2 = new TGraph(size, &tp_vec[0], &si2_vec[0]);
   grph2 -> SetMarkerStyle(20);
   grph2 -> SetMarkerColor(3);  
   grph2 -> SetTitle("Normalized");



   grph2 -> SetTitle("Normalization: Misak's Model");
   grph2 -> Draw("APL");

 
   grph2->GetXaxis()->SetTitle("t'");
   grph2->GetXaxis()->CenterTitle();
   grph2->GetYaxis()->SetTitle("Cross Section / Pole Factor");
   grph2->GetYaxis()->CenterTitle();

   c1->Update();


}
bool get_lines(char  line[nfiles][100], FILE * out[nfiles]){
    bool not_end = true;
    for (int i =0; i<nfiles; i++){
        not_end = not_end && fgets(line[i], 100, out[i]);
    }
    return not_end;
}

const int nfiles = 2;

int m_model_f2n_reader(bool print = false){
  //define values we will read
    Float_t prt[nfiles];
    Float_t tp[nfiles];
    Float_t si[nfiles];

    //set vectors for storing values for graphing
    // prt is same for every run so only use values from DATA_1.OUT
    vector<Float_t> tp_vec;
    vector<Float_t> si1_vec;
    vector<Float_t> si2_vec;

    //open the 3 MISAK_DATA files 
    FILE * out[nfiles]; 
    for(int i = 0; i<nfiles; i++){
        int idx = i + 1; //files indexed starting from 2
        // DATA_1.OUT is the 0.9 * F2N data, DATA_2.OUT is the nominal F2N data, 
        // DATA_3.OUT is the 1.1 * F2N data
        string out_path = get_outfile_path(idx);         
         out[i] = fopen(out_path.c_str(), "r"); 
         if (print) cout << "Opened File at " << out_path << endl;
    }
        //make the tree and setup the branches
    TTree *tree = new TTree("T", "tag_data_tree");
    //set all the branches
    for(int i = 0; i<nfiles; i++){
        int idx = i + 1; //files indexed starting from 1

       tree -> Branch(TString::Format("prt%i",idx), prt+i, 
                      TString::Format("prt%i/F",idx));
       tree -> Branch(TString::Format("TP%i", idx), tp + i,
                      TString::Format("TP%i/F", idx));
       tree -> Branch(TString::Format("si%i",idx), si+i, 
                      TString::Format("si%i/F", idx));
       }
       char line[nfiles][100];
       //start reading the files
       if (print) cout << "Reading the files... " << endl;
       char first_char;
       while (get_lines(line, out)){
           first_char = line[0][0]; //if a comment should be same for all files
           if (first_char != '#' ){ //check that line isnt a comment
               //read data and store in variables
               for(int i=0; i<nfiles; i++){
                  sscanf(line[i], "%f %f %f", tp+i, prt+i, si+i);

                  if (print) {
                     cout<<"Writing values to "<< i<<": " <<tp[i]<<" "<< prt[i] << " ";
                     cout << si[i]<<" "<< endl;
                  }
               }
                     //put the values in the arrays too
                    
                     tp_vec.push_back(TMath::Abs(tp[0]));
                     si1_vec.push_back(si[0]);
                     si2_vec.push_back(si[1]);

                     tree->Fill();

          }
       }
        
        
   //if(print) print_my_tree(tree);
   make_graphs(tp_vec, si1_vec, si2_vec);
   //write the tree to a root file 
   char * file_name_str = "m_model_f2n_data.root";
   TFile * root_file = 0;
   root_file = new TFile (file_name_str, "RECREATE");
   tree->Write();
   for (int i=0; i<nfiles; i++){
      fclose(out[i]);}
   if(print) cout << "Closed Files" << endl;
   delete root_file;
 
    return 1;
}


