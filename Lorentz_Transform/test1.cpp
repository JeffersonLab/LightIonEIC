//Made by Oz Amram
//June 2014
// This file is an example of how to use the lorentz transformation functions provided here
// compile it with something like "g++ test1.cpp lorentz_transform.cpp -o test1.exe -I. `root-config --cflags --libs` "
// and then run with ./test1.exe 
#include <lorentz_transform.h>



int main(){
    TLorentzVector e_vec, d_vec, e_ft, d_ft, e_cm1, d_cm1, e_rot, d_rot;
    Double_t e_en = 5.0;
    Double_t e_p =  5.0;
    Double_t d_p = 100.0;
    Double_t d_en = TMath::Sqrt(TMath::Power(deut_mass, 2) + TMath::Power(d_p, 2));

    e_vec.SetPxPyPzE(0.0, 0.0, e_p, e_en);
    d_vec.SetPxPyPzE(0.0, 0.0, -d_p, d_en);

    Tran_to_fixed_target(e_vec, d_vec, e_ft, d_ft);
    Tran_to_cm(e_vec, d_vec, e_cm1, d_cm1);

//    Float_t theta = (TMath::Pi)/36.0;
//    Float_t phi = (TMath::Pi)/24.0;
//    e_rot = rotate(e_vec, theta, 0);
//    d_rot = rotate(d_vec, 0, phi);
    
    cout << "e_vec: ";
    print_vec(e_vec);
    cout << "d_vec: ";
    print_vec(d_vec);

    cout << "e_vec in FT frame: ";
    print_vec(e_ft);
    cout << "d_vec in FT frame: ";
    print_vec(d_ft);

    cout << "e_vec in cm frame: ";
    print_vec(e_cm1);
    cout << "c_vec in cm frame: ";
    print_vec(d_cm1);
    return 0;


}
