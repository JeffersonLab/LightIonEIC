  void batch_collMC()
  {
      gROOT->ProcessLine(".L SpectatorMC_collMC.cpp+");
      gROOT->ProcessLine("mainx()");
  }
