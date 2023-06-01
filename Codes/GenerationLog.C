#include<iostream>

{
    using namespace RooFit;
    
    // Observable:
   RooRealVar Delta_Mass("Delta_Mass","m_{ES} (GeV)",.14,.16);

// Parameters:
   RooRealVar sigmean("sigmean","B^{#pm} mass", .146  , .140 , .155 );
   RooRealVar sigwidth("sigwidth","B^{#pm} width", 01, 1, 10);

// Build a Gaussian PDF:
   RooGaussian signalModel("signal","signal PDF",Delta_Mass,sigmean,sigwidth);

// Build  background PDF:
   
    
    RooRealVar a("a","a",100, 0 ,10000) ;
    RooRealVar b("b","b", 100 , 0, 10000) ;
    RooGenericPdf background("background","background", "a*log(b*x)", RooArgSet(Delta_Mass,a,b)) ;

    /*
    RooRealVar m0("m0", "Media", .5, 0.00001, 1);  // Ajusta la media según la escala en la que se encuentran tus datos
    RooRealVar k("k", "Desviación estándar", 1.5, 0.00001, 10); 
    RooLognormal background("background", "Fondo Log-Normal", Delta_Mass, m0, k);*/
    

// Construct a signal and background PDF:.q
   RooRealVar nsig("nsig","#signal events",200,0.,10000);
   RooRealVar nbkg("nbkg","#background events",800,0.,10000);
   RooAddPdf model("model","g+a",RooArgList(signalModel,background),RooArgList(nsig,nbkg));

// The PDF is used to generate an un-binned toy data set, then the PDF is fit to that data set using an un-binned maximum likelihood fit.
// Then the data are visualized with the PDF overlaid.

// Generate a toy MC sample from composite PDF:
   RooDataSet *data = model.generate(Delta_Mass, 10000);

// Perform extended ML fit of composite PDF to toy data:
   model.fitTo(*data , Extended(true));

// Plot toy data and composite PDF overlaid:
   RooPlot* mesframe = Delta_Mass.frame();
   data->plotOn(mesframe);
   model.plotOn(mesframe);
   model.plotOn(mesframe, Components(background), LineStyle(ELineStyle::kDashed));

   mesframe->Draw();
}



