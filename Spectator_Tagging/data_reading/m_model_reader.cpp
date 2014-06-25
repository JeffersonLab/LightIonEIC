//Written by Oz Amram 6/20/2014
//
//Reads data from ../misak_model/MISAK_DATA.OUT and plots it as well as 
//storing it in a root file

#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>


using namespace std;

string get_outfile_path(){
    // returns the path to the EVTP.OUT file
    // ONLY WORKS IF EVTP.OUT is located at ../christians_model/EVTP.OUT
    char cur_path[FILENAME_MAX];
    getcwd(cur_path, FILENAME_MAX);
    string str_path = string(cur_path);
    unsigned end_dir = str_path.find_last_of("/\\");
    string outfile_path = str_path.substr(0,end_dir) + "/misak_model/MISAK_DATA.OUT";
    //cout <<  outfile_path << endl;
    return outfile_path;
}
void print_my_tree(TTree *tree){
   Float_t prt_val = 0.0;
   TBranch *prt_branch = 0;

   tree->SetBranchAddress("prt", &prt_val, &prt_branch);
   Int_t entries = tree->GetEntries();
   for(int i=0; i<entries; i++){
        prt_branch->GetEntry(i);
        cout << prt_val <<endl;
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
   //make t' vs cross section graph
   TVirtualPad  *pad1 = c1->cd();
   pad1->SetLogy(); //make Y axis on log scale
   int size = tree -> Draw("si:prt","", "goff");//easiest way to turn tree data to graph
   TGraph * grph1 = new TGraph(size, tree->GetV2(), tree->GetV1());
   grph1->GetXaxis()->SetTitle("Transverse Recoil Momentum");
   grph1->GetXaxis()->CenterTitle();
   grph1->GetYaxis()->SetTitle("Integrated Cross Section");
   grph1->GetYaxis()->CenterTitle();
   grph1 ->SetMarkerStyle(7);
   grph1 -> Draw("APL");
   grph1 ->SetTitle("Cross Section vs prt");

   
   c1 -> Update();
}

int m_model_reader(bool print=false){
  //define values we will read
    Float_t prt;
    Float_t si;
    Float_t x;
    Float_t f2d;
    
    TFile * root_file = 0;

    string out_path = get_outfile_path();
    FILE * out = fopen(out_path.c_str(), "r"); //open data (have to reconvert to a char *)
    //create root file in current directory
    root_file = new TFile ("m_model_data.root", "RECREATE");
    //make the tree and setup the branches
    TTree *tree = new TTree("T", "m_data_tree");
    tree -> Branch("prt", &prt,"prt/F");
    tree -> Branch("si", &si,"si/F");

    if (print) cout << "Opened File at " << out_path << endl;
    char line [100];
    //start reading the file
    if (print) cout << "Reading the file... " << endl;
    char first_char;
    while (fgets(line, 100, out)){
        first_char = line[0];
        if (first_char != '#' ){ //check that line isnt a comment
            //read data and store in variables
            sscanf(line, "%f %f", &prt, &si);
            if (print) cout<<"Writing values: " <<prt<<" "<<si<<" "<<endl;
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
