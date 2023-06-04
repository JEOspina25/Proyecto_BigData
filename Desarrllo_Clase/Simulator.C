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
Simulator::Simulator( RooAbsPdf& Model_ , RooRealVar & Obs_ , int const& nbin_): TObject(){
    nbin = nbin_;
    this->Obs = &Obs_;
    this->Model = &Model_;
    this->DataSet = Model->generate(*Obs , 10000); // Default constructor
    this->FitResult = Model->fitTo(*DataSet , RooFit::Extended(true) , RooFit::Save(true));
}

TCanvas* Simulator::MainPlot(RooAbsPdf *Model ,RooDataSet *DataSet){
    int x = 800;
    int y = 600;

    TCanvas *c = new TCanvas("c", "c", x, y);
    c->Divide(1, 2 , 0 ,0);

    c->cd(1); gPad->SetRightMargin(0.01);
    RooPlot* frame = Obs->frame();
    DataSet->plotOn(frame);
    Model->plotOn(frame);
    frame->Draw("A");
    frame->GetYaxis()->CenterTitle(); 
    frame->GetYaxis()->SetTitleSize(0.07); 
    frame->GetYaxis()->SetTitleOffset(0.5);
    frame->GetYaxis()->SetLabelSize(0.045);

    c->cd(2);gPad->SetRightMargin(0.01); gPad->SetBottomMargin(0.3);
    RooHist* pullHist = frame->pullHist();
    RooPlot* pullFrame = Obs->frame();
    pullFrame->addPlotable(pullHist, "P");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitle("(Data - Fit) /#sigma");
    pullFrame->GetYaxis()->CenterTitle(); 
    pullFrame->GetYaxis()->SetTitleSize(0.07); 
    pullFrame->GetYaxis()->SetTitleOffset(0.5);
    pullFrame->GetYaxis()->SetLabelSize(0.045);
    pullFrame->GetXaxis()->CenterTitle(); 
    pullFrame->GetXaxis()->SetTitleSize(0.07); 
    pullFrame->GetXaxis()->SetTitleOffset(0.9);
    pullFrame->GetXaxis()->SetLabelSize(0.045);

    TLine* zeroLine = new TLine(6.05, 0, 6.5, 0);
    zeroLine->SetLineStyle(2);
    pullFrame->Draw("A");
    zeroLine->Draw("same");
    c->Update();

    return c;
}






