//credit: https://root.cern/doc/master/rf101__basics_8C.html

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;

void plot_Proyecto()
{
Int_t bins = 48;
Int_t datos = 900;
// S e t u p   m o d e l
// ---------------------

// Declare variables x,mean,sigma with associated name, title, initial value and allowed range
RooRealVar x("x", "x", 0.14, 0.16);
RooRealVar mean("mean", "mean of gaussian", 0.1455, 0.14, 0.16);
RooRealVar sigma("sigma", "width of gaussian", 0.0012, 0, 0.2);

RooRealVar c("c","c",0);
RooPolynomial bkg1("bkg1","Background",x,RooArgSet(c),0);

RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

//Pesos de Background y señal
RooRealVar Ns("Ns","Ns",0.,500);
RooRealVar Nb("Nb","Nb",0.,500); 
//Modelo para la masa
RooAddPdf MassModel("MassModel","MassModel",RooArgList(gauss,bkg1),RooArgList(Ns,Nb));
// Build gaussian pdf in terms of x,mean and sigma

// G e n e r a t e   e v e n t s
// -----------------------------

// Generate a dataset of 1000 events in x from gauss
RooDataSet *data = gauss.generate(x, datos,Binning(bins));


// F i t   m o d e l   t o   d a t a
// -----------------------------

// Fit pdf to data
MassModel.fitTo(*data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
bkg1.Print("v"); 


//mean.Print();
//sigma.Print();

// Draw frame on a canvas
// -----------------------------
TCanvas *c1 = new TCanvas("c1","",600,500);

RooPlot *xframe2 = x.frame();
data->plotOn(xframe2, MarkerSize(0.8), DrawOption("P"),Binning(bins),DataError(RooAbsData::SumW2),XErrorSize(0));
//MassModel.plotOn(xframe2,LineColor(kBlue),LineWidth(2),Name("fit"));
MassModel.plotOn(xframe2,Components(gauss),LineColor(kBlue),LineWidth(2),Name("signal")); 
MassModel.plotOn(xframe2,Components(bkg1),LineColor(kBlue),LineWidth(4), LineStyle(kDashed) ,Name("bkg1")); 

xframe2->GetXaxis()->SetNdivisions(6);
xframe2->SetMinimum(-1); 
xframe2->GetYaxis()->SetRangeUser(0, 250);   
xframe2->Draw();

TLatex *tex2 = new TLatex(0.2,0.926,"CMS");
tex2->SetNDC();
tex2->SetTextFont(61);
tex2->SetTextSize(0.05); 
tex2->SetLineWidth(2);
tex2->Draw();

//Leyenda: objetos de la figura
TLegend *leg = new TLegend(0.4,0.49,0.83,0.71); 
leg->SetTextSize(0.04);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(xframe2->findObject("Data")," Data","ep"); 
leg->AddEntry(xframe2->findObject("signal")," Fit result","l");
leg->AddEntry(xframe2->findObject("bkg1"),"Background.","l");
leg->Draw();
//Se almacena el lianzo en la carpeta plots
c1->Draw();
c1->Print("plots/Plot_DeltaM.png");
}
