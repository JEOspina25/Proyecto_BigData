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

using namespace std;
using namespace RooFit;

// Funcion para guardar los parametros del fit

// Funcion para crear figura
TCanvas* CreateCanvas(TString cname, RooFitResult* result, RooDataSet* data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel, RooAddPdf sumgau, RooPolynomial bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean) 
{   
    // Número de bines
    Double_t nbin = ((supM-infM)/0.005)+1;

    int H = 500;
    int W = 600;

    // Creando Canvas
    TCanvas *c1 = new TCanvas(cname,cname,500,50,W,H);
    c1->cd() ;  
    c1->SetLeftMargin(0.11);
    c1->SetRightMargin(0.01);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.1);

    // Creando frame
    RooPlot* Mframe = M.frame(infM,supM,nbin);
    // Se dibujan los datos y el ajuste
    data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
    MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
    
    // Dibujar el frame
    MassModel.plotOn(Mframe,Components(bkg1),LineColor(kBlue),LineWidth(2),LineStyle(kDashed),Name("bkg")); 
    data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
    MassModel.plotOn(Mframe);
    Mframe->SetTitle(""); 
    
    // Ajustes del frame
    Mframe->SetYTitle("Events / 4 MeV"); 
    Mframe->SetLabelSize(0.05,"XY");
    Mframe->SetTitleSize(0.05,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();
    Mframe->GetYaxis()->SetNdivisions(505,1);
    Mframe->GetXaxis()->SetNdivisions(505,1);
    Mframe->GetXaxis()->SetTickLength(0.01);    
    Mframe->GetXaxis()->SetDecimals(1); 
    Mframe->SetTitleOffset(0.85,"X");
    Mframe->SetTitleOffset(1.1,"Y");
    // Mframe->SetMinimum(0.5); 
    Mframe->Draw();
    
    // Legend 1
    TLegend *Legend1 = new TLegend(0.18,0.18,0.38,0.48); 
    Legend1->SetTextSize(0.04);
    Legend1->SetFillColor(0);
    Legend1->SetBorderSize(0);
    Legend1->SetFillStyle(0);
    Legend1->AddEntry("", "29 nb^{-1}(13 Tev)","");
    Legend1->AddEntry(Mframe->findObject("Data")," Data","ep"); 
    Legend1->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
    Legend1->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
    Legend1->Draw();

    TLegend *Legend2 = new TLegend(0.65,0.8,0.85,0.9);
    Legend2->AddEntry("", "16 < p_{T} < 24 Gev, |#eta| < 2.1","");
    Legend2->SetBorderSize(0);
    Legend2->SetTextSize(0.04);
    Legend2->SetFillStyle(0);
    Legend2->SetMargin(0.1);
    Legend2->Draw();

    c1->Modified();
    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
    c1->Update();
    return c1; 
  
}

int DoubleGaussiansMC(){

    // Valores Límete para las variables
    Double_t Mmin = 1.8; 
    Double_t Mmax = 1.975;

    // Variables a usar
    RooRealVar M("M", "Mass (K#pi#pi) (GeV)", Mmin, Mmax);

    // ---- MassModel ----

    // -Parámetros Señal-
    RooRealVar mean("mean"," Mass mean",1.875,1.85,1.9,"GeV");

    // Gausiana 1
    RooRealVar width("width"," Mass width",0.02,0.001,0.025,"GeV"); // Sigma1
    RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

    // Gausiana 2
    RooRealVar width2("width2"," Mass width2 ",0.025,0.001,0.05,"GeV"); // Sigma2
    RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

    // -Parámetros Background-
    RooRealVar c0("c","c",0.0,1000.0);
    RooRealVar c1("c","c",0.0,1.0);
    RooPolynomial Bkg("Bkg","Exp. Background",M,RooArgList(c0,c1));

    // Cantidad de datos por cada componente
    RooRealVar Ns("Ns","Ns",0.,500);
    RooRealVar Nb("Nb","Nb",1500.,2000);   
    RooRealVar fs("fs","fs",0.8,0.,1.);

    // Suma de las dos gausianas
    RooAddPdf Sumgaus("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));

    // Modelo de masa (2 Gausianas + Exponencial)
    RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sumgaus,Bkg),RooArgList(Ns,Nb));

    // Generación
    RooDataSet* Data_M = MassModel.generate(M,40000);

    // ---- Fitting ----
    RooFitResult* ResultFit = MassModel.fitTo(*Data_M,Extended(),Minos(kFALSE),Save(kTRUE));
 
    // // Hacer Gráfica
    Double_t supM = Mmax;
    Double_t infM = Mmin;
    TCanvas* Canvas1 = CreateCanvas("Canvas_MasaK", ResultFit, Data_M, M, supM, infM, MassModel, Sumgaus, Bkg, Ns, Nb, width, width2, fs, mean);  
    Canvas1->Print("../plots/Plot_Mass_K_2pi.png");

    Int_t bins = 38;
    Int_t datos = 800;

    return 0;
}