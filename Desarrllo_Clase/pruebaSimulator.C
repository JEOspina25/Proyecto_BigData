/*Aca probaremos nuestra clase Simulator creada*/
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "Simulator.h"


{
    using namespace RooFit;

    // Crear una instancia de la clase Simulator
    RooRealVar x("x","x",0,10) ;
    RooRealVar mean1("mean1","mean of gaussian 1",2) ;
    RooRealVar mean2("mean2","mean of gaussian 2",3) ;
    RooRealVar sigma("sigma","width of gaussians",1) ;
    RooGaussian gauss1("gauss1","gaussian PDF",x,mean1,sigma) ;
    RooGaussian gauss2("gauss2","gaussian PDF",x,mean2,sigma) ;

    // Build Argus background PDF
    RooRealVar argpar("argpar","argus shape parameter",-1.0) ;
    RooRealVar cutoff("cutoff","argus cutoff",9.0) ;
    RooArgusBG argus("argus1","Argus PDF",x,cutoff,argpar) ;

    // Add the components

    RooRealVar g1frac("g1frac","fraction of gauss1",0.5) ;
    RooRealVar g2frac("g2frac","fraction of gauss2",0.1) ;
    RooAddPdf model("sum","g1+g2+a",RooArgList(gauss1,gauss2,argus),RooArgList(g1frac,g2frac)) ;
    Simulator simulator(model, x);

    // Generar datos sintéticos
    simulator.Generate();

    // ... código posterior ...
}