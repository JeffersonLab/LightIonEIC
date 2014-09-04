  void batch_fixMC()
  {
      gROOT->ProcessLine(".L SpectatorMC_fixMC.cpp+");
      gROOT->ProcessLine("mainx()");
  }
