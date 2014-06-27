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
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

string int_to_string(int i){
   //a hack to work around a bug with to_string
   //see http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
   ostringstream ss;
   ss << i;
   return ss.str();
}


string get_m_path(){
    // returns the path to the EVTP.OUT file
    // ONLY WORKS IF EVTP.OUT is in ../tag from this file
    char cur_path[FILENAME_MAX];
    getcwd(cur_path, FILENAME_MAX);
    string str_path = string(cur_path);
    cout << str_path + "/m_model_multi_data.root" << endl;
    return str_path + "/m_model_multi_data.root";
}


string get_c_path(){
    // returns the path to the EVTP.OUT file
    // ONLY WORKS IF EVTP.OUT is in ../tag from this file
    char cur_path[FILENAME_MAX];
    getcwd(cur_path, FILENAME_MAX);
    string str_path = string(cur_path);
    cout << str_path + "/c_model_multi_data.root" << endl;
    return str_path + "/c_model_multi_data.root";
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



void draw_graphs(TTree * c_tree, TTree * m_tree){
   HallA_style(); //make it pretty
   int marker_style = 20;
   TCanvas *c1 = new TCanvas("c1", "Cross Sections", 700, 700);
   //make  vs cross section graph
   TVirtualPad  *pad1 = c1->cd();
   pad1->SetLogx();
   pad1->SetLogy(); //make Y axis on log scale
   //all the data has the same prt values so we use the same

   int c_size = c_tree -> Draw("FSIG3:TP3", "PTR2 !=0", "goff");
   //the 3 index is the nominal F2N for christians_model
   TGraph * c_grph = new TGraph(c_size, c_tree->GetV2(), c_tree->GetV1());
   c_grph -> SetMarkerStyle(marker_style);
   c_grph -> SetMarkerColor(2);
   c_grph -> SetTitle("Christian's Model");
   
   int m_size = m_tree -> Draw("si2:TP2", "", "goff");
   // the 2 index is the nominal F2N for misaks model
   TGraph * m_grph = new TGraph(m_size, m_tree->GetV2(), m_tree->GetV1());
   m_grph -> SetMarkerStyle(marker_style);
   m_grph -> SetMarkerColor(3);  
   m_grph -> SetTitle("Misak's Model");

   TMultiGraph * mgc = new TMultiGraph();

   mgc -> SetTitle("Model Comparison");
   mgc -> Add(c_grph);
   mgc -> Add(m_grph);
   mgc -> Draw("APL");

 
   mgc->GetXaxis()->SetTitle("t' (GEV/C)");
   mgc->GetXaxis()->CenterTitle();
   mgc->GetYaxis()->SetTitle("Cross Section (nB/GEV**4)");
   mgc->GetYaxis()->CenterTitle();

   c1->BuildLegend();
   c1->Update();

   TCanvas * c2 = new TCanvas("c2", "Ratio", 700, 700);
   c2 -> cd();
   const Double_t * tp_vals = c_tree ->GetV2();
   const Double_t * c_vals = c_tree->GetV1();
   const Double_t * m_vals = m_tree->GetV1();
   //compute ratios
   vector<Double_t> ratio;
   for(int i = 0; i <200; i++){
        ratio.push_back(m_vals[i]/c_vals[i]);
   }
    
   TGraph * r_grph = new TGraph(m_size, m_tree->GetV2(), &ratio[0]);

   r_grph -> SetMarkerStyle(marker_style);
   r_grph -> SetMarkerColor(4);
   r_grph -> SetTitle("Model Ratios");
   r_grph ->Draw("APL");
   r_grph->GetXaxis()->SetTitle("t' (GEV/C)");
   r_grph->GetXaxis()->CenterTitle();
   r_grph->GetYaxis()->SetTitle("Ratio: Misak / Christian");
   r_grph->GetYaxis()->CenterTitle();

   c2->Update();
}

int c_and_m_comparison(bool print=false){
    string c_path = get_c_path();
    string m_path = get_m_path();
    TFile * c_file = new TFile("./c_model_multi_data.root");
    TFile * m_file = new TFile("./m_model_multi_data.root");
    TTree * c_tree = (TTree*) c_file ->Get("T");
    TTree * m_tree = (TTree*) m_file ->Get("T");
    draw_graphs(c_tree, m_tree);
  }


