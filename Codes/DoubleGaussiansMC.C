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
TCanvas* CreateCanvasNomPull(TString cname, RooFitResult* result, RooDataSet* data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel, RooAddPdf sumgau, RooPolynomial bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean) 
{   
    // Número de bines
    Double_t nbin = ((supM-infM)/0.005)+1;

    int H = 500;
    int W = 600;

    // Creando Canvas
    TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
    c1->cd() ;  
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.1);
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
    Mframe->SetYTitle("Events / 2.5 MeV"); 
    Mframe->SetLabelSize(0.07,"XY");
    Mframe->SetTitleSize(0.08,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();
    Mframe->GetYaxis()->SetNdivisions(505,1);
    Mframe->GetXaxis()->SetNdivisions(505,1);
    Mframe->GetXaxis()->SetTickLength(0.0);    
    Mframe->GetXaxis()->SetDecimals(1); 
    Mframe->SetTitleOffset(0.8,"X");
    Mframe->SetTitleOffset(0.6,"Y");
    Mframe->SetMinimum(0.5); 
    Mframe->Draw();
    
    // Legend 1
    TLegend *leg = new TLegend(0.08,0.58,0.28,0.88); 
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
    leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
    leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
    leg->Draw();
    
    // Valores y errores en valor medio y ancho de la señal
    Double_t Mpsi = mean.getVal()*1000.0;
    Double_t MpsiE = mean.getError()*1000.0;

    Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
    Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
    
    // Legend Parámetros
    TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
    legpar->SetTextSize(0.04);
    legpar->SetTextFont(42);
    legpar->SetFillColor(0);
    legpar->SetBorderSize(0);
    legpar->SetFillStyle(0);
    legpar->AddEntry("",Form("M(B_{c}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
    legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
    legpar->AddEntry("",Form("N_{B_{c}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
    legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
    legpar->Draw();

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
    RooRealVar mean("mean"," Mass mean",1.855,1.85,1.9,"GeV");

    // Gausiana 1
    RooRealVar width("width"," Mass width",0.02,0.001,0.05,"GeV"); // Sigma1
    RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

    // Gausiana 2
    RooRealVar width2("width2"," Mass width2 ",0.025,0.001,0.1,"GeV"); // Sigma2
    RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

    // -Parámetros Background-
    RooRealVar c("c","c",0.0,1000.0);
    RooPolynomial Bkg("Bkg","Exp. Background",M,c);

    // Cantidad de datos por cada componente
    RooRealVar Ns("Ns","Ns",0.,2000);
    RooRealVar Nb("Nb","Nb",0.,2000);   
    RooRealVar fs("fs","fs",0.8,0.,1.);

    // Suma de las dos gausianas
    RooAddPdf Sumgaus("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));

    // Modelo de masa (2 Gausianas + Exponencial)
    RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sumgaus,Bkg),RooArgList(Ns,Nb));

    // Generación
    RooDataSet* Data_M = Sumgaus.generate(M,10000);

    // ---- Fitting ----
    RooFitResult* ResultFit = MassModel.fitTo(*Data_M,Extended(),Minos(kFALSE),Save(kTRUE));
 
    // // Hacer Gráfica
    Double_t supM = Mmax;
    Double_t infM = Mmin;
    TCanvas* Canvas1 = CreateCanvasNomPull("Canvas_MasaK", ResultFit, Data_M, M, supM, infM, MassModel, Sumgaus, Bkg, Ns, Nb, width, width2, fs, mean);  
    Canvas1->Print("../plots/MassK.png");

    Int_t bins = 38;
    Int_t datos = 800;

    // TCanvas *c1 = new TCanvas("c1","",600,500);

    // RooPlot *xframe2 = M.frame();
    // Data_M->plotOn(xframe2,Name("Data_M"), MarkerSize(1 ),MarkerStyle(8), DrawOption("P"),Binning(bins),DataError(RooAbsData::SumW2),XErrorSize(0));
    // //MassModel.plotOn(xframe2,LineColor(kBlue),LineWidth(2),Name("fit"));
    // MassModel.plotOn(xframe2,Components(Sumgaus),LineColor(kBlue+1),LineWidth(3),Name("signal")); 
    // MassModel.plotOn(xframe2,Components(Bkg),LineColor(kBlue),LineWidth(4), LineStyle(kDashed) ,Name("Bkg")); 


    // xframe2->GetXaxis()->SetNdivisions(6);
    // xframe2->SetMinimum(-1); 
    // xframe2->GetYaxis()->SetRangeUser(0, 250);   
    // xframe2->SetTitle(" ");
    // xframe2->SetYTitle("Events/0.4 Mev"); 
    // xframe2->SetXTitle("#Delta M(Gev)");
    // xframe2->Draw();

    // TLatex *tex2 = new TLatex(0.12,0.836,"CMS");
    // tex2->SetNDC();
    // tex2->SetTextFont(60);
    // tex2->SetTextSize(0.05); 
    // tex2->SetLineWidth(2);
    // tex2->Draw("L");

    // auto legend = new TLegend(1,1.5,.62,.38);
    // legend->AddEntry("", "29 nb^{-1}(13 Tev)","");
    // legend->SetBorderSize(0);
    // legend->SetTextSize(0.05);
    // legend->SetFillStyle(0);
    // legend->SetMargin(0.1);
    // legend->Draw();

    // auto legend2 = new TLegend(0.5,0.8,.4,0.9);
    // legend2->AddEntry("", "16 < p_{T} < 24 Gev, |#eta| < 2.1","");
    // legend2->SetBorderSize(0);
    // legend2->SetTextSize(0.04);
    // legend2->SetFillStyle(0);
    // legend2->SetMargin(0.1);
    // legend2->Draw();




    // //Leyenda: objetos de la figura
    // TLegend *leg = new TLegend(0.4,0.49,0.83,0.71); 
    // leg->SetTextSize(0.04);
    // leg->SetFillColor(0);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);
    // leg->AddEntry(xframe2->findObject("data")," Data","ep"); 
    // leg->AddEntry(xframe2->findObject("signal")," Fit ","l");
    // leg->AddEntry(xframe2->findObject("bkg1"),"combinatorial background","l");
    // leg->Draw();
    // //Se almacena el lianzo en la carpeta plots
    // c1->Update();

    // c1->Draw();
    // c1->Print("plots/MassK.png");
    return 0;
}