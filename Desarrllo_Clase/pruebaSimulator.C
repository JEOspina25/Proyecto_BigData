/*Aca probaremos nuestra clase Simulator creada*/
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "Simulator.h"


{
    using namespace RooFit;

    RooRealVar Mass("Mass", "Mass", 6.05, 6.5);
    RooRealVar Mean("Mean", "Mean", 6, 6.6);
    RooRealVar Sigma("Sigma", "Sigma", 0.00001, 1);
    RooGaussian Signal("Signal", "Signal", Mass, Mean, Sigma);

    RooRealVar C("C", "C", -10, 10);
    RooExponential Backg("Backg", "Backg", Mass, C);

    RooRealVar Nsig("Nsig", "Nsig", 0, 13000);
    RooRealVar Nbkg("Nbkg", "Nbkg", 0, 13000);
    RooAddPdf Model("Model", "Model", RooArgList(Signal, Backg), RooArgList(Nsig, Nbkg));
    

    Simulator MySimulation(Model , Mass , 100);
}