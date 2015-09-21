void batch()
{
    gROOT->ProcessLine(".L SpectatorMC_collMC.cpp+");
    gROOT->ProcessLine("mainx(0.01,0.0125892541179417,1,1.25892541179417)");
}
