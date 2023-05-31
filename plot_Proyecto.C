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
Int_t bins = 40;
Int_t datos = 900;
// S e t u p   m o d e l
// ---------------------

// Declare variables x,mean,sigma with associated name, title, initial value and allowed range
RooRealVar x("x", "x", 0.14, 0.16);
RooRealVar mean("mean", "mean of gaussian", 0.1455, 0.14, 0.16);
RooRealVar sigma("sigma", "width of gaussian", 0.0012, 0, 0.2);

RooRealVar c("c","c",0,1   );
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
MassModel.fitTo(*data);

//mean.Print();
//sigma.Print();

// Draw frame on a canvas
// -----------------------------
TCanvas *c1 = new TCanvas("c1","",600,800);

RooPlot *xframe2 = x.frame();
data->plotOn(xframe2, MarkerSize(0.8), DrawOption("P"),Binning(bins),DataError(RooAbsData::SumW2),XErrorSize(0));
MassModel.plotOn(xframe2);
bkg1.plotOn(xframe2,LineColor(kGreen),LineWidth(2));

xframe2->GetXaxis()->SetNdivisions(6);
xframe2->Draw();

TLatex *tex2 = new TLatex(0.2,0.926,"CMS");
tex2->SetNDC();
tex2->SetTextFont(61);
tex2->SetTextSize(0.05); 
tex2->SetLineWidth(2);
tex2->Draw();
//Se almacena el lianzo en la carpeta plots
c1->Draw();
c1->Print("plots/Datos_fit_pull.png");
}
