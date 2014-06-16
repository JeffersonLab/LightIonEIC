//Written by Oz Amram 6/13/2014
//
//Reads data from ../tag/EVTP.OUT and plots it as well as storing it in a root file

#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>


using namespace std;

string get_outfile_path(int index){
    // returns the path to the EVTP.OUT file
    // ONLY WORKS IF EVTP.OUT is in ../tag from this file
    char cur_path[FILENAME_MAX];
    getcwd(cur_path, FILENAME_MAX);
    string str_path = string(cur_path);
    unsigned end_dir = str_path.find_last_of("/\\");
    string outfile_path = str_path.substr(0,end_dir) + "/tag/EVTP"+ string(index) +".OUT";
    //cout <<  outfile_path << endl;
    return outfile_path;
}
void print_my_tree(TTree *tree){
   Float_t tp_val = 0.0;
   TBranch *tp_branch = 0;

   tree->SetBranchAddress("TP", &tp_val, &tp_branch);
   Int_t entries = tree->GetEntries();
   for(int i=0; i<entries; i++){
        tp_branch->GetEntry(i);
        cout << tp_val <<endl;
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



void make_graph(const TTree *tree){
   HallA_style(); //make it pretty
   TCanvas *c1 = new TCanvas("c1", "First Graph", 700, 700);
   c1 ->Size(0,0);
   TCanvas *c2 = new TCanvas("c2", "2nd Graph", 700, 700);
   //make t' vs cross section graph
   TVirtualPad  *pad1 = c1->cd();
   pad1->SetLogy(); //make Y axis on log scale
   int size = tree -> Draw("FSIG:abs(TP)","PR2 != 0", "goff");//easiest way to turn tree data to graph
   TGraph * grph1 = new TGraph(size, tree->GetV2(), tree->GetV1());
   grph1->GetXaxis()->SetTitle("t' (GEV^2)");
   grph1->GetXaxis()->CenterTitle();
   grph1->GetYaxis()->SetTitle("Cross Section (NB/GEV^4)");
   grph1->GetYaxis()->CenterTitle();
   grph1 ->SetMarkerStyle(7);
   grph1 -> Draw("APL");
   grph1 ->SetTitle("Cross Section vs t'");

   //make t' vs F2N/POL graph
   c2 -> cd();
   int size = tree -> Draw("F2_SPOL:abs(TP)","PR2 !=0", "goff");//easiest way to turn tree data to graph
   TGraph * grph2 = new TGraph(size, tree->GetV2(), tree->GetV1());
   grph2->GetXaxis()->SetTitle("t' (GEV^2)");
   grph2->GetXaxis()->CenterTitle();
   grph2->GetYaxis()->SetTitle("Tagged F2/POLE Factor");
   grph2->GetYaxis()->CenterTitle();
   grph2 ->SetMarkerStyle(7);
   grph2 -> Draw("APL");
   grph2 ->SetTitle("F2/POLE vs t'");


   c1 -> Update();
   c2 -> Update();
}

int multi_tag_reader(bool print=false){
  //define values we will read
    Float_t TP;
    Float_t PR2;
    Float_t PTR;
    Float_t FSIG;
    Float_t UNUM;
    Float_t F2_SPOL;
    
    TFile * root_files[3];
    for(int idx = 0; idx<3; idx++){

        TFile * root_file[idx] = 0;
        TTree *trees[3] = {new TTree("T2", "tag_data2_tree"), 
                           new TTree("T3", "tag_data3_tree"),
                           new TTree("T4", "tag_data4_tree")}
                

        string out_path = get_outfile_path(idx+2); //files indexed starting from 2
        FILE * out = fopen(out_path.c_str(), "r"); //open data (have to reconvert to a char *)
        //create root file in current directory
        root_file = new TFile ("tag_data" + string(idx+2) + ".root", "RECREATE");
        //make the tree and setup the branches
        tree -> Branch("TP", &TP,"TP/F");
        tree -> Branch("PR2", &PR2,"PR2/F");
        tree -> Branch("PTR", &PTR,"PTR/F");
        tree -> Branch("FSIG", &FSIG,"FSIG/F");
        tree -> Branch("UNUM", &UNUM,"UNUM/F");
        tree -> Branch("F2_SPOL", &F2_SPOL,"F2_SPOL/F");

        if (print) cout << "Opened File at " << out_path << endl;
        char line [80];
        //start reading the file
        if (print) cout << "Reading the file... " << endl;
        char first_char;
        while (fgets(line, 80, out)){
            first_char = line[0];
            if (first_char != '#' ){ //check that line isnt a comment
                //read data and store in variables
                sscanf(line, "%f %f %f %f %f %f", &TP, &PR2, &PTR, &FSIG, &UNUM, &F2_SPOL);
                if (print) cout<<"Writing values: " <<TP<<" "<<PR2<<" ";
                if (print) cout <<PTR<<" "<<FSIG<<" "<<UNUM<<" "<<F2_SPOL<<endl;
                tree->Fill();
            }
        }
        if(print) {tree->Print();print_my_tree(tree);}
        make_graph(tree);
        tree->Write();
        fclose(out);
        if(print) cout << "Closed File" << endl;
        delete root_file;
    return 1;
}
