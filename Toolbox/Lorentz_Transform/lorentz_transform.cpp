//Made by Oz Amram
//June 2014
#include <lorentz_transform.h>



void Tran_to_fixed_target (TLorentzVector e_vec, TLorentzVector d_vec,
                            TLorentzVector &e_ft, TLorentzVector  &d_ft){
   //input: 
   //   e_vec - electron 4 vector in some frame
   //   d_vec - ion 4 vector in some frame
   //output:
   //   e_ft - electron 4 vector in rest frame of ion (fixed target frame)
   //   d_ft - ion 4 vector in its rest frame 
   e_ft = e_vec;
   d_ft = d_vec;
   
   d_ft.Boost(-d_vec.BoostVector()); //should be (0,0,0, mass)
   e_ft.Boost(-d_vec.BoostVector());
   return;


}
void Tran_to_cm(TLorentzVector e_vec, TLorentzVector d_vec, 
                TLorentzVector &e_cm, TLorentzVector &d_cm){
     //input: 
     //   e_vec - electron 4 vector in some frame
     //   d_vec - ion 4 vector in some frame
     //output:
     //   e_ft - electron 4 vector in center of mass frame
     //   d_ft - ion 4 vector in center of mass frame
     e_cm = e_vec;
     d_cm = d_vec;

     TLorentzVector cm_vec = e_vec + d_vec;

     e_cm.Boost(-cm_vec.BoostVector());
     d_cm.Boost(-cm_vec.BoostVector());
     return;

}

TLorentzVector rotate(TLorentzVector v, Float_t theta, Float_t phi){
    v.RotateX(theta);
    v.RotateY(phi);
    return v;

}


void print_vec(TLorentzVector v){
    cout <<"("<< v.Px() << " , " << v.Py() << " , " << v.Pz() << " , " << v.E() << ")" <<endl;
}

void print_se(TLorentzVector e_vec, TLorentzVector p_vec){
    double_t se = 2.0 * e_vec.E() * (abs(p_vec.Px()) + p_vec.E()) + deut_mass*deut_mass;
    cout << "se: " << se <<endl;
}

//int main(){
//    //example use of these functions
//    TLorentzVector e_vec, d_vec, e_ft, d_ft, e_cm1, d_cm1;
//    Double_t e_en = 5.0;
//    Double_t e_p =  5.0;
//    Double_t d_p = 100.0;
//    Double_t d_en = TMath::Sqrt(TMath::Power(deut_mass, 2) + TMath::Power(d_p, 2));
//
//    e_vec.SetPxPyPzE(e_p, 0.0, 0.0, e_en);
//    d_vec.SetPxPyPzE(-d_p, 0.0, 0.0, d_en);
//
//    Tran_to_fixed_target(e_vec, d_vec, e_ft, d_ft);
//    Tran_to_cm(e_vec, d_vec, e_cm1, d_cm1);
//
//    print_vec(e_vec);
//    print_vec(d_vec);
//    print_vec(e_ft);
//    print_vec(d_ft);
//
//    print_vec(e_cm1);
//    print_vec(d_cm1);
//    return 0;
//
//
//}
