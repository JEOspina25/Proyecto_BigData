#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"
#include "TObject.h"
#include "Simulator.h"

//using namespace RooFit;

ClassImp(Simulator)

/*Descripcion*/

Simulator::Simulator() : TObject(){}
Simulator::Simulator( RooAbsPdf & Model , RooRealVar const& Obs): TObject(){

    this->DataSet = Model.generate(Obs , 10000);
    this->FitResult = Model.fitTo(*DataSet , RooFit::Extended(true) , RooFit::Save(true));
};




