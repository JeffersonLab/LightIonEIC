//Written by Oz Amram 6/18/2014
//
//Reads data from  MISAK_DATA_1.OUT, MISAK_DATA_2.OUT, MISAK_DATA_3.OUT
//They must be located in ../misak_model/
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
    char cur_path[FILENAME_MAX];
    getcwd(cur_path, FILENAME_MAX);
    string str_path = string(cur_path);
    unsigned end_dir = str_path.find_last_of("/\\");
    string outfile_path = str_path.substr(0,end_dir) + "/misak_model/MISAK_DATA_"
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



void make_graphs(vector<Float_t> prt_vec,  vector<Float_t> si1_vec, 
                 vector<Float_t> si2_vec,  vector<Float_t> si3_vec,
                 vector<Float_t> si_ratio_1_2, vector<Float_t> si_ratio_3_2){
   HallA_style(); //make it pretty
   TCanvas *c1 = new TCanvas("c1", "Cross Sections", 700, 700);
   //make  vs cross section graph
   TVirtualPad  *pad1 = c1->cd();
   pad1->SetLogy(); //make Y axis on log scale
   //all the data has the same prt values so we use the same
   int size = prt_vec.size();
   TGraph * grph1 = new TGraph(size, &prt_vec[0], &si1_vec[0]);
   grph1 -> SetMarkerStyle(21);
   grph1 -> SetMarkerColor(2);
   grph1 -> SetTitle("0.9 * F2N");
   TGraph * grph2 = new TGraph(size, &prt_vec[0], &si2_vec[0]);
   grph2 -> SetMarkerStyle(21);
   grph2 -> SetMarkerColor(3);  
   grph2 -> SetTitle("Nominal F2N");
   TGraph * grph3 = new TGraph(size, &prt_vec[0], &si3_vec[0]);
   grph3 -> SetMarkerStyle(21);
   grph3 -> SetMarkerColor(4);
   grph3 -> SetTitle("1.1 * F2N");

   TMultiGraph * mgc = new TMultiGraph();

   mgc -> SetTitle("Cross Sections");
   mgc -> Add(grph1);
   mgc -> Add(grph2);
   mgc -> Add(grph3);
   mgc -> Draw("APL");

 
   mgc->GetXaxis()->SetTitle("t' (GEV/C)");
   mgc->GetXaxis()->CenterTitle();
   mgc->GetYaxis()->SetTitle("Integrated Cross Section");
   mgc->GetYaxis()->CenterTitle();

   c1->BuildLegend();
   c1->Update();

   TCanvas *c2 = new TCanvas("c2", "Ratios", 700, 700);
   c2-> cd();
   TGraph * ratio1 = new TGraph(size, &prt_vec[0], &si_ratio_1_2[0]);
   ratio1 -> SetMarkerStyle(21);
   ratio1 -> SetMarkerColor(2);
   ratio1 -> SetTitle("0.9 * F2N / Nominal");

   TGraph * ratio2 = new TGraph(size, &prt_vec[0], &si_ratio_3_2[0]);
   ratio2 -> SetMarkerStyle(21);
   ratio2 -> SetMarkerColor(3);
   ratio2 -> SetTitle("1.1 * F2N / Nominal");

   TMultiGraph * mgr = new TMultiGraph();
   mgr -> SetTitle("Ratios");
   mgr -> Add(ratio1);
   mgr -> Add(ratio2);
   mgr -> Draw("APL");
   
   mgr->GetXaxis()->SetTitle("t' (GEV/C)");
   mgr->GetXaxis()->CenterTitle();
   mgr->GetYaxis()->SetTitle("Ratio of Cross Sections");
   mgr->GetYaxis()->CenterTitle();
   mgr->GetYaxis() -> SetTitleOffset(1.8);


   c2-> BuildLegend();


}

Float_t to_tp(Float_t prt){
    //converts a prt (proton recoil transverse momentum) to tp (proton tranverse momentum
    //in nucleus)
    return prt;
}

int m_model_multi_reader(bool print=false){
  //define values we will read
    Float_t prt[3];
    Float_t tp[3];
    Float_t si[3];

    //set vectors for storing values for graphing
    // prt is same for every run so only use values from DATA_1.OUT
    vector<Float_t> tp_vec;
    vector<Float_t> si1_vec;
    vector<Float_t> si2_vec;
    vector<Float_t> si3_vec;
    vector<Float_t> si_ratio_1_2;
    vector<Float_t> si_ratio_3_2;

    //open the 3 output files 
    FILE * out[3]; 
    for(int i = 0; i<3; i++){
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
    for(int i = 0; i<3; i++){
        int idx = i + 1; //files indexed starting from 1

       tree -> Branch(TString::Format("prt%i",idx), prt+i, 
                      TString::Format("prt%i/F",idx));
       tree -> Branch(TString::Format("TP%i", idx), tp + i,
                      TString::Format("TP%i/F", idx));
       tree -> Branch(TString::Format("si%i",idx), si+i, 
                      TString::Format("si%i/F", idx));
       }
       char line[3][100];
       //start reading the files
       if (print) cout << "Reading the files... " << endl;
       char first_char;
       while (fgets(line[0], 100, out[0]) && fgets(line[1], 100, out[1]) && 
             fgets(line[2], 100, out[2])){
           first_char = line[0][0]; //if a comment should be same for all files
           if (first_char != '#' ){ //check that line isnt a comment
               //read data and store in variables
               for(int i=0; i<3; i++){
                  sscanf(line[i], "%f %f", prt+i, si+i);
                  tp[i] = to_tp(prt[i]);
//                  prt[i]+= 0.18;

                  if (print) {
                     cout<<"Writing values to "<< i<<": " <<prt[i]<<" "<< tp[i] << " ";
                     cout << si[i]<<" "<< endl;
                  }
               }
                     //put the values in the arrays too
                    
                     tp_vec.push_back(TMath::Abs(tp[0]));
                     si1_vec.push_back(si[0]);
                     si2_vec.push_back(si[1]);
                     si3_vec.push_back(si[2]);
                     si_ratio_1_2.push_back(si[0]/si[1]);
                     si_ratio_3_2.push_back(si[2]/si[1]);

                     tree->Fill();

          }
       }
        
        
   if(print) print_my_tree(tree);
   make_graphs(tp_vec, si1_vec, si2_vec, si3_vec, si_ratio_1_2, si_ratio_3_2);
   //write the tree to a root file 
   char * file_name_str = "m_model_multi_data.root";
   TFile * root_file = 0;
   root_file = new TFile (file_name_str, "RECREATE");
   tree->Write();
   for (int i=0; i<3; i++){
      fclose(out[i]);}
   if(print) cout << "Closed Files" << endl;
   delete root_file;
 
    return 1;
}


