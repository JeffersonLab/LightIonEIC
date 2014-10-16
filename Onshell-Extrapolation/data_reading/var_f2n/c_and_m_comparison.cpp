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
    cout << str_path + "/m_model_f2n_data.root" << endl;
    return str_path + "/m_model_f2n_data.root";
}


string get_c_path(){
    // returns the path to the EVTP.OUT file
    // ONLY WORKS IF EVTP.OUT is in ../tag from this file
    char cur_path[FILENAME_MAX];
    getcwd(cur_path, FILENAME_MAX);
    string str_path = string(cur_path);
    cout << str_path + "/c_model_f2n_data.root" << endl;
    return str_path + "/c_model_f2n_data.root";
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

const int ngraphs = 1;

void draw_graphs(Float_t tp[200], Float_t c_cs[ngraphs][200], Float_t m_cs[ngraphs][200]){
   HallA_style(); //make it pretty
   int marker_style = 20;
   TList can1, can2, mgcs, r_grphs, c_grphs, m_grphs;
   Float_t ratio[ngraphs][200];
   for (int i = 0; i<ngraphs; i++){
       int idx = i + 1;
       can1.Add(new TCanvas(TString::Format("c1%i",idx), TString::Format("Comparison%i", idx), 700, 700));
       can2.Add(new TCanvas(TString::Format("c2%i",idx), TString::Format("Ratio%i", idx), 700, 700));

       mgcs.Add(new TMultiGraph());
       c_grphs.Add(new TGraph(200, tp,c_cs[i]));
       m_grphs.Add(new TGraph(200, tp, m_cs[i]));

       //compute ratios
       for(int j = 0; j <200; j++){
            ratio[i][j] = 2*m_cs[i][j]/c_cs[i][j];
       }
        
       r_grphs.Add(new TGraph(100, tp, ratio[i]));
   }
   for (int i = 0; i<ngraphs ; i++){
       //set correct values for this graph
       int idx = i + 1;
       TCanvas * c1 = can1.At(i);
       TCanvas * c2 = can2.At(i);
       TMultiGraph * mgc = mgcs.At(i);
       TGraph * r_grph = r_grphs.At(i);
       TGraph * c_grph = c_grphs.At(i);
       TGraph * m_grph = m_grphs.At(i);

       //make  comparison graph
       c1->cd()->SetLogy();//make Y axis on log scale
       mgc -> SetTitle(TString::Format("Model Comparison: Graph #%i", idx));
       c_grph -> SetMarkerStyle(marker_style);
       c_grph -> SetMarkerColor(2);
       c_grph -> SetTitle("Christian's Model");
       
       m_grph -> SetMarkerStyle(marker_style);
       m_grph -> SetMarkerColor(3);  
       m_grph -> SetTitle("Misak's Model");


       mgc -> Add(c_grph);
       mgc -> Add(m_grph);
       mgc -> Draw("APL");

     
       mgc ->GetXaxis()->SetTitle("t' (GEV^2)");
       mgc ->GetXaxis()->CenterTitle();
       mgc ->GetYaxis()->SetTitle("Cross Section (nB/GEV**4)");
       mgc ->GetYaxis()->CenterTitle();

       c1-> BuildLegend();
       c1-> Update();


       c2->cd();
       r_grph -> SetMarkerStyle(marker_style);
       r_grph -> SetMarkerColor(4);
       r_grph -> SetTitle("Model Cross Section Ratios");
       r_grph ->Draw("APL");
       r_grph->GetXaxis()->SetTitle("t' (GEV^2)");
       r_grph->GetXaxis()->CenterTitle();
       r_grph->GetYaxis()->SetTitle("Ratio: Realistic Model / Simple Model");
       r_grph->GetYaxis()->CenterTitle();


    }
}

void get_data(TTree * c_tree, TTree * m_tree, Float_t c_cs[ngraphs][200], Float_t m_cs[ngraphs][200],
              Float_t tp[200]){
    Float_t tp_val, c_val[ngraphs], m_val[ngraphs];
    c_tree->SetBranchAddress("TP1", &tp_val);
    for(int i = 0; i <ngraphs; i ++){
        int idx = i + 1;
        c_tree->SetBranchAddress(TString::Format("FSIG%i", idx), &c_val[i]);
        m_tree->SetBranchAddress(TString::Format("si%i", idx), &m_val[i]);
    }
    int entries = c_tree->GetEntries();
    for(int j=0; j<entries; j++){
         c_tree->GetEntry(j);
         m_tree->GetEntry(j);
         tp[j] = tp_val;
         for(int i=0; i<ngraphs; i++){
            c_cs[i][j] = c_val[i];
            m_cs[i][j] = m_val[i];

      }
    }
    return;
}

int c_and_m_comparison(bool print=false){
    string c_path = get_c_path();
    string m_path = get_m_path();
    TFile * c_file = new TFile(c_path.c_str());
    TFile * m_file = new TFile(m_path.c_str());
    TTree * c_tree = c_file->Get("T");
    TTree * m_tree = m_file->Get("T");
    Float_t c_cs[ngraphs][200], m_cs[ngraphs][200], tp[200];
    get_data(c_tree, m_tree, c_cs, m_cs, tp);
    draw_graphs(tp, c_cs, m_cs);
  }


