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

using namespace RooFit;
using namespace std;

class Simulator{
    /*Esta clase recibe como input un modelo pdf con sus respectivos argumentos y genera a partir de ello
    1.) datos sinteticos
    2.) un fit a esos datos sinteticos
    3.) un pull y montecarlo de dichos datos
    
    todo esto solo con findes educativos de formacion*/

    private: /*Atributos*/

    RooAbsPdf Model; // Modelo que pasa el usuario
    RooRealVar Obs; //Observale asociado al modelo

    RooDataSet Sintetic_Data; //Datos sinteticos 
    int Num_Data = 10000; //numero de datos generados sinteticamente
    int H = 800 , W = 600; //ancho del canvas
    RooFitResult FitResult; //resultado del fit

    public:

    Simulator(RooAbsPdf _Model , RooRealVar _Obs){  //Constructor

    this -> Model = _Model;
    this -> Obs = _Obs;

    }

    /*Metodos-----------------------------------------------------------------------------------------------------*/
    Generate(){
    
    RooDataSet Sintetic_Data = Model.generate(Obs,Num_Data);
    }

    ~Simulator(){}//Destructor
};





