/**
 * avalanche.cc
 * General program flow based on example code from the Garfield++ website.
 *
 * Demonstrates electron avalanche and induced signal readout with
 * 2D finite-element visualization in Garfield++ with a LEM.  LEM 
 * parameters are from: 
 * C. Shalem et. al. Nucl. Instr. Meth. A, 558, 475 (2006).
 *
*/
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <string>

#include "MediumMagboltz.hh"
#include "ComponentElmer.hh"
#include "ComponentBase.hh"
#include "Sensor.hh"
#include "ViewField.hh"
#include "Plotting.hh"
#include "ViewFEMesh.hh"
#include "ViewSignal.hh"
#include "ViewMedium.hh"
#include "GarfieldConstants.hh"
#include "Random.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"

#include <TApplication.h>
#include "TFile.h"
#include "TNtuple.h"
#include "TTUBE.h"
#include "TNode.h"

using namespace std;
using namespace Garfield;

void CheckElectricField(ComponentElmer *elm); // visualize よう 関数。

double visualXlw=-0.06 * 1;
double visualXup= 0.06 * 1;
double visualZlw=-0.10 * 1;
double visualZup= 0.10 * 1;
int points = 200;



//unit [cm]
const double gkRhole = 0.0085;
const double gkThickness = 0.06851;


int main(int argc, char *argv[]=NULL) 
{
  TApplication app("app", &argc, argv);
  cerr << "input parameters = " <<  argc << endl;
  for (int i=0; i<argc; i++) {
    cerr << i << " parameter = " << argv[i] << endl; 
  }
  
  // ----- initial settings
  if ( !argv[1] ) argv[1]=0;
  int runNum   = atoi( argv[1] ); 
  int nRuns    = atoi( argv[2] ); 
  string gemdir = "./gemcell";
  string gasdir = "./gas";
  //-------------------------//

  nRuns = 1; 
  const double pitch= 1.00;  // LEM pitch in cm
//  const double bx = 0;
//  const double by = 0;
//  const double bz = 3.5;

  // Variables describing signal binning.
  double tEnd = 600.0;
  int nsBins = 500;

  // ----- define the medium.
  MediumMagboltz* gas = new MediumMagboltz();
  //  gas->SetTemperature(293.);                // Set the temperature (K)
  //  gas->SetPressure(760.);                   // Set the pressure (Torr)
  //  gas->SetComposition("Ar", 95.0, "CF4", 3.0, "IC4H10", 2.0); // t2k
  //  gas->SetComposition("Ar", 90., "CH4", 10.); // Specify the gas mixture (Ar/CH4 90:10)
  //  gas->SetComposition("ar", 90., "co2", 10.); // Specify the gas mixture (Ar/CO2 90:10)
  //  gas->LoadGasFile(gasdir+"/Ar90CH410/ar_90_ch4_10.gas");
  gas->LoadGasFile(gasdir+"/ar_90_ch4_10.gas");
  //  gas->LoadGasFile(gasdir+"/ar_95_cf4_3_ic4h10_2.gas");
  gas->LoadIonMobility(gasdir+"/IonMobility_Ar+_Ar.txt");
  gas->EnablePenningTransfer (0.25, 0,"ar");
  gas->EnableDrift();                         // Allow for drifting in this medium

  gas->Initialise(true);
  // ----- import an Elmer-created LEM and the weighting field for the readout electrode.
  cout << "--------------------" << endl;
  ComponentElmer * elm = new ComponentElmer(
					    gemdir+"/mesh.header",
					    gemdir+"/mesh.elements",
					    gemdir+"/mesh.nodes",
					    gemdir+"/dielectrics.dat",
					    gemdir+"/gemcell.result",
					    "cm");

  elm->SetMedium(0,gas);
  elm->EnableMirrorPeriodicityX();
  elm->EnableMirrorPeriodicityY();
  elm->PrintMaterials();
  elm->PrintRange();

  //  CheckElectricField(elm);

  // ----- Set up a sensor object.
   Sensor* sensor = new Sensor();
   sensor->AddComponent(elm);
   double defvol = 10*pitch;	
   double lwZvol = -0.4000; // 1mm: volume definition on Z
   sensor->SetArea(-defvol, -defvol, lwZvol,
		   defvol,  defvol, defvol);
   sensor->AddElectrode(elm,"wtlel");
   sensor->SetTimeWindow(0.,tEnd/nsBins,nsBins);
   
   // ----- Create an avalanche object
   AvalancheMicroscopic* aval = new AvalancheMicroscopic();
   aval->SetSensor(sensor);
   aval->EnableSignalCalculation();

   aval->EnableAvalancheSizeLimit(3);

//   elm->SetMagneticField(bx,by,bz);
//   aval->EnableMagneticField();

   AvalancheMC* drift = new AvalancheMC();
   drift->SetSensor(sensor);
   drift->SetDistanceSteps(0.0005); // 0.001 10um
   drift->UnsetTimeWindow();


   // data 保存用 ntuple for each event
   TNtuple *tn0 = new TNtuple("tn0","tn0","nElecsAbs:nElecsEff:nIons");
   // data 保存用 ntuple for each electron
   TNtuple *tn1 = new TNtuple("tn1","tn1","evnt:eleX:eleY:eleZ:xi:yi:zi"); 

   TString nRun;
   nRun += runNum;   
   cerr << "\n\n----- iteration nruns = " << nRuns << endl;

   for (int i=0; i<nRuns; i++) {
     TCanvas *cFE;
     cFE = new TCanvas("cFE","cFE", 500,500);
     cFE->Draw();

     ViewFEMesh * vFE = new ViewFEMesh();
     ViewDrift* vD = new ViewDrift();
     vD->SetArea(-0.6,-0.5,-0.5,
		 0.6, 0.5, 0.5);
     aval ->EnablePlotting(vD);
     drift->EnablePlotting(vD);

     vFE->SetCanvas(cFE);
     vFE->SetComponent(elm);
     vFE->SetArea(visualXlw, visualZlw, -0.1,
		  visualXup, visualZup, +0.1);
     vFE->SetFillMesh(-1);
     vFE->SetColor(1,kYellow+2);
     vFE->SetColor(2,kGreen+3);
     vFE->SetColor(3,kGreen+3);
     vFE->SetViewDrift(vD);
     vFE->EnableAxes();
     vFE->SetXaxisTitle("x (cm)");
     vFE->SetYaxisTitle("z (cm)");    
     vFE->SetPlane(0,-1,0,0,0,0);


//     double ri = (1.0E-4*300) * sqrt(RndmUniform() );//[mm]
//     double thetai = TwoPi * RndmUniform();
//     double xi = ri * cos(thetai); // random
//     double yi = ri * sin(thetai); // random
     //     double zi = 0.0100*3;
     double xi = 0;
     double yi = 0;
     double zi = 0.10;
     double e0 = 0.1;


     //initial position
     double x0 = xi;
     double y0 = yi;
     double z0 = zi;

     double xe1, ye1, ze1, te1, e1;
     double xe2, ye2, ze2, te2, e2; 
     int status;

     cerr << "----- start avalanche" << endl;
     aval -> AvalancheElectron(x0, y0, z0, 0., e0, 0., 0., 0);     
     // initialize
     int nElecsAbs=0;
     int nElecsEff=0;
     int nIons=0;
     aval ->GetAvalancheSize(nElecsAbs, nIons);
     cerr << "size elec:" << nElecsAbs << " ion:" << nIons << endl;
     int nElecEndPoints = aval->GetNumberOfElectronEndpoints();
     std::cerr << "... avalanche complete with " << nElecEndPoints << " electron tracks." << std::endl;



     for (int j = nElecEndPoints; j--;) {
       aval->GetElectronEndpoint(j,xe1,ye1,ze1,te1,e1,
				 xe2,ye2,ze2,te2,e2, status);

       tn1->Fill(i,xe1,ye1,ze1,xi,yi,zi);

       if ( ze2<-0.3 ) nElecsEff++; // count electrons which reach to sufficienty lower reagion (effective gain) 
       
     }
     tn0->Fill(nElecsAbs,nElecsEff,nIons); // 絶対ゲイン: どんだけ、増幅したか、
     std::cerr << "... nElecs effective = " << nElecsEff << std::endl;
     //     vD->Plot();
     vFE->Plot();
     //     vField->Plot();

   }//end nRuns loop

   TFile *tf = new TFile("./output/root" + nRun + ".root","recreate");
   tf->cd();
   tn0->Write();
   tn1->Write();

   tf->Close();	

   // Extract the calculated signal.
   double bscale = tEnd/nsBins;  // time per bin
   double sum = 0.;              // to keep a running sum of the integrated signal
   TH1F * hS = new TH1F("hh","hh",nsBins,0,tEnd);               // total signal
   TH1F * hInt = new TH1F("hInt","hInt",nsBins,0,tEnd);         // integrated signal
   // Fill the histograms with the signals.
   //  Note that the signals will be in C/(ns*binWidth), and we will divide by e to give a signal in e/(ns*binWidth).
   //  The total signal is then the integral over all bins multiplied by the bin width in ns.
   for(int i = 0; i < nsBins; i++) {
     double wt = sensor->GetSignal("wtlel",i) / ElementaryCharge;
     sum += wt;
     hS->Fill(i*bscale,wt);
     hInt->Fill(i*bscale,sum);
   }

   // ----- Plot fields.
   std::cout << "#FINISHED: " << std::endl;

   app.Run(kTRUE);

   return 0;
}

void CheckElectricField(ComponentElmer *elm)
{
  gStyle->SetOptStat(0);
  TH2F *h2Ex = new TH2F("h2Ex","h2Ex",points,visualXlw,visualXup,points, visualZlw,visualZup);
  TH2F *h2Ey = new TH2F("h2Ey","h2Ey",points,visualXlw,visualXup,points, visualZlw,visualZup);
  TH2F *h2Ez = new TH2F("h2Ez","h2Ez",points,visualXlw,visualXup,points, visualZlw,visualZup);
  TH2F *h2E  = new TH2F("h2E" ,"h2E" ,points,visualXlw,visualXup,points, visualZlw,visualZup);
  TH2F *h2V  = new TH2F("h2V" ,"h2V" ,points,visualXlw,visualXup,points, visualZlw,visualZup);
  
  TGraph *gE0 = new TGraph();
  TGraph *gE1 = new TGraph();
  
  Medium* medium;
  int status;
  //   double ex, ey, ez, v, absE;
  double ex, ey, ez, v;
  for (int i=0; i<points; i++) { // x
    //    cerr << i << endl;
    for (int j=0; j<points; j++) { // z
      double x = h2Ez->GetXaxis()->GetBinCenter(i+1);
      double z = h2Ez->GetYaxis()->GetBinCenter(j+1);
      //double y = h2Ez->GetYaxis()->GetBinCenter(j+1);
      //         elm->ElectricField(x, y, 0.0, ex, ey, ez, v, medium, status);
      elm->ElectricField(x, 0.0, z, ex, ey, ez, v, medium, status);
      
      h2Ex->SetBinContent(i+1,j+1,ex );
      h2Ey->SetBinContent(i+1,j+1,ey );
      h2Ez->SetBinContent(i+1,j+1,ez );
      h2E->SetBinContent(i+1,j+1, sqrt(ex*ex+ey*ey+ez*ez));
      h2V->SetBinContent(i+1,j+1, v);
      if (i==points/2   ) gE0->SetPoint(gE0->GetN(),z,sqrt(ex*ex+ey*ey+ez*ez));
      if (i==points/2+20) gE1->SetPoint(gE1->GetN(),z,sqrt(ex*ex+ey*ey+ez*ez));
      
    }
  }
  TCanvas * cE = new TCanvas("cE","cE",1250,900);	
  cE->Divide(3,2);
  cE->cd(1);	
  gE0->Draw("ap*"); // unit [V]
  gE1->Draw(" p*"); // unit [V]
  gE1->SetMarkerColor(4);
  cE->cd(2);
  h2V->SetTitle(";x axis;y axis");
  h2V->Draw("colz");
  cE->cd(3);
  h2E->SetTitle(";x axis;y axis");
  h2E->Draw("colz");
  
  cE->cd(4); h2Ex->Draw("colz"); // Efield x component	
  cE->cd(5); h2Ey->Draw("colz"); // Efield y component
  cE->cd(6); h2Ez->Draw("colz"); // Efield z component	
  
  cE->Print("output/field.pdf");

#if 0
  // ----- Set up the object for FE mesh visualization: takes a bit longer time

  TCanvas * cFE = new TCanvas("cFE","cFE",700,700);
  ViewFEMesh * vFE = new ViewFEMesh();
  vFE->SetCanvas(cFE);
  vFE->SetComponent(elm);
  // note: here the x-y axes correspond to projection 
  //       chosen and z-axis is irrelevant for 2D projections
  vFE->SetArea(visualXlw, visualZlw, -0.1,
	       visualXup, visualZup, +0.1);
  vFE->SetPlane(0,-1,0,0,0,0);
  vFE->SetFillMesh(-1);
  vFE->SetColor(1,kYellow+2);
  vFE->SetColor(2,kGreen+3);
  vFE->SetColor(3,kGreen+3);
  vFE->EnableAxes(); // comment this to disable creation of independent axes when contours are plotted
  vFE->Plot();
  cFE->Print("output/cFE.pdf");
#endif


};

