{
    gROOT->Macro("Utils/CPlot.cc+");
    gROOT->Macro("Utils/MitStyleRemix.cc+");
    gROOT->Macro("Utils/CEffUser1D.cc+");
    gROOT->Macro("Utils/CEffUser2D.cc+");
    gROOT->Macro("Utils/RooGaussDoubleSidedExp.cc+");
    gROOT->Macro("Utils/RooFitter.cc+");
    // Show which process needs debugging
    gInterpreter->ProcessLine(".! ps |grep root.exe");
}
