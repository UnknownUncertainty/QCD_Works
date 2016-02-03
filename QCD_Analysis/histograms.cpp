#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <time.h>
#include <stdio.h>



#include <TSystem.h>
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TChain.h>
#include <TBranch.h>
#include <TEventList.h>
#include "TClonesArray.h"
#include <TFile.h>
#include <TH1F.h>
#include <TPad.h>
#include <TF1.h>
#include <TH2F.h>

#include "TLorentzVector.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TAttMarker.h"
#include <TAttMarker.h>
//#include "TClone.h"
//#include "TROOT.h"
#include <TAttFill.h>
#include "TStyle.h"
#include "TColor.h"
#include "TF2.h"
#include "TExec.h"
#include "TCanvas.h"
#include <TAxis.h>
#include "TColor.h"
#include <TMath.h>
#include <TProfile.h>
#include "TStyle.h"
#include <TString.h>
#include "TLatex.h"
#include "TLegend.h"
#include "TVectorD.h"
//#include "TProfile2D.h."
#include "TLine.h"
#include <algorithm>
#include <time.h>
//#include <vector>
#include "TGraph.h"
#include <cmath>
#include <string>
#include <vector>
//#include "/Users/nsonmez/root/macros/tdrStyle.C"

using namespace std;

int main()
{	//gROOT->ProcessLine(".L /Users/nsonmez/root/macros/tdrStyle.C");
	//gStyle->SetOptStat("ne");
	
	//setTDRStyle();
	
	//______________________________________________________________
	// read the root file


	TFile * mainroot = new TFile("~/Desktop/Programs/qcd_delphes/delphes_analysis/result.root","READ");
	

	int PT_BIN_ARRAY[]={0,50,100,150,250,500,750,1000,1500,2000,3000,4000,5000,6500};	
	string Pt_Reso_Cut_Name[] = {"No_Cut","Eta_Cut","R_Cut","R_Eta_Cut"};
	string JET_NAME[] = {"Leading_Jet","Second_Jet","Third_Jet","Fourth_Jet"};

//	char PT_Reso_canvas_name[100];
//	char PT_Reso_canvas_title[100];
//	char PT_Reso_canvas_print[1000];
//	char PT_Reso_Get_Hist[100];
	
	char get_hist_name[100];
	char canvas_name[100];
	char canvas_title[1000];
	char canvas_print[1000];
////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------Rapidity Resolution for All Jets ---------------------------------//
////////////////////////////////////////////////////////////////////////////////////////////////
/*	
	//char fit_projection_name[100];
	//----Dizin İçerisindeki Histogramları Okuma İşlemi Yapılmakta--------//////////////////

	TH2F *hist_DeltaY_vs_GenY[14];
	TCanvas *Canvas_Rap_Reso_All_Y_[14];
	TProfile *hist_Rapi_Prof[14];
	TH1D * projRapiY[14];
	TF1 *fitProjectionRapi[14];
        
	//TGraph *gr = new TGraph(n,x,y);
	double Rapidity_Reso_All_Jets_Sigma[13];
	double Rapidity_Reso_All_Jets_Pt_Ort[13];
	double Rapidity_Reso_All_Jets_Mean[13];
	int Rapidity_Reso_All_Jets_n=13;

	//double Rapidity_Reso_All_Jets_Pt_Ort[13];

	//double sigma[13];
	//double Pt_value[13];
	

	for(int i=0; i<13; i++)
	{
		
	    sprintf(get_hist_name, "Rap_Resolution_TH2F/DeltaY_vs_GenY_%ipt%i_y%i" ,PT_BIN_ARRAY[i] ,PT_BIN_ARRAY[i+1] ,i);
	    sprintf(canvas_name,  "Rapidity_Reso_%i< PT <%i",PT_BIN_ARRAY[i] ,PT_BIN_ARRAY[i+1]);
	    sprintf(canvas_title, "Rapidity_Reso_%i< PT <%i",PT_BIN_ARRAY[i] ,PT_BIN_ARRAY[i+1]);

	    hist_DeltaY_vs_GenY[i]=(TH2F*)mainroot->Get(get_hist_name); 
	    Canvas_Rap_Reso_All_Y_[i] = new TCanvas(canvas_name,canvas_title,0,0,800,1200);

            Canvas_Rap_Reso_All_Y_[i]->SetFillColor(41);
            Canvas_Rap_Reso_All_Y_[i]->Divide(1,3);
            Canvas_Rap_Reso_All_Y_[i]->cd(1);
	
	    hist_DeltaY_vs_GenY[i]->GetXaxis()->SetTitle("Ygen");
            hist_DeltaY_vs_GenY[i]->GetYaxis()->SetTitle("Ygen-Ycalo");
            hist_DeltaY_vs_GenY[i]->Draw("cont0");

	    Canvas_Rap_Reso_All_Y_[i]->cd(2);

            hist_Rapi_Prof[i] = hist_DeltaY_vs_GenY[i]->ProfileX();
	    //rapidity_profile->SetTitle("Rapidity_Res_60_100");
            hist_Rapi_Prof[i]->SetMinimum(-0.2);
            hist_Rapi_Prof[i]->SetMaximum(0.2);
            hist_Rapi_Prof[i]->SetTitle("Profile");
            hist_Rapi_Prof[i]->GetXaxis()->SetTitle("Ygen");
            hist_Rapi_Prof[i]->GetYaxis()->SetTitle("Ygen-Ycalo");
            hist_Rapi_Prof[i]->Draw("E1");
	
		//TLine *line1;
		//TLine *line2;
	//(c1->GetUxmin()
	    Canvas_Rap_Reso_All_Y_[i]->Update();
	    TLine *line1 = new TLine(2.5,-0.2,2.5,0.2);
            line1->SetLineColor(kRed);
            line1->Draw("SAME");

            TLine *line2 = new TLine(-2.5,-0.2,-2.5,0.2);
            line2->SetLineColor(kRed);
            line2->Draw("SAME");
  	 
	    Canvas_Rap_Reso_All_Y_[i]->cd(3);

//        gStyle->SetOptStat(1111);
	
    	   projRapiY[i] = hist_DeltaY_vs_GenY[i]->ProjectionY();
	   projRapiY[i]->SetTitle("Projection Y is Fitted to Gaussian");

           fitProjectionRapi[i] = new TF1("fitProjectionRapi[i]","gaus",-0.2,0.2);
           projRapiY[i]->Fit("fitProjectionRapi[i]","WW","",-0.2,0.2);
           projRapiY[i]->GetYaxis()->SetTitle("Events");
           projRapiY[i]->Draw();
 	   
	   sprintf(canvas_print, "/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/Rapidity_Reso_%i_PT_%i",PT_BIN_ARRAY[i] ,PT_BIN_ARRAY[i+1]);
           Canvas_Rap_Reso_All_Y_[i]->Print(canvas_print);



	   
	   Rapidity_Reso_All_Jets_Sigma[i] = fitProjectionRapi[i]->GetParameter(2);
           Rapidity_Reso_All_Jets_Mean[i] = fitProjectionRapi[i]->GetParameter(1);
	   Rapidity_Reso_All_Jets_Pt_Ort[i] = ( PT_BIN_ARRAY[i] + PT_BIN_ARRAY[i+1] )/2;
	   	 
	}

cout << " *************************Rapidity Grafikleri Bitti************************* " << endl;
cout << " *************************Rapidity Grafikleri Bitti************************* " << endl;
cout << " *************************Rapidity Grafikleri Bitti************************* " << endl;
*/

//TCanvas *canvas1v1 = new TCanvas("canvas1v1","title",0,0,600,500);
/*
TCanvas *c1 = new TCanvas("canvas_name","canvas_title",0,0,800,1200);
TGraph *gr = new TGraph(Rapidity_Reso_All_Jets_n,Rapidity_Reso_All_Jets_Pt_Ort,Rapidity_Reso_All_Jets_Sigma);
gr->Draw("AC*");
c1->Print("/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/deneme.eps");

TCanvas *c2 = new TCanvas("canvas_name","canvas_title",0,0,800,1200);
TGraph *gr2 = new TGraph(Rapidity_Reso_All_Jets_n,Rapidity_Reso_All_Jets_Pt_Ort,Rapidity_Reso_All_Jets_Mean);
gr2->Draw("AC*");
c2->Print("/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/deneme2.eps");
*/
cout << " *************************TGraph Grafikleri Bitti************************* " << endl;
cout << " *************************TGraph Grafikleri Bitti************************* " << endl;
cout << " *************************TGraph Grafikleri Bitti************************* " << endl;

  //////////////////////////
  //R32 for genjet/////////
  /////////////////////////
	
 
	string Ratio_32_canvas_name[]  = {"3_Jets_2_Jets_Gen","3_Jets_2_Jets_Calo"};
	string Ratio_32_canvas_title[] = {"3_Jets_2_Jets_Gen","3_Jets_2_Jets_Calo"};
	string Ratio_32_X_Axis[] = {"Gen HT (GeV)","Calo HT (GeV)"};
	string Ratio_32_Y_Axis[] = {"Ratio_{32} for GenJet","Ratio_{32} for CaloJet"};

	TCanvas *Ratio_R32_canvas[2];
	//TCanvas *Ratio_R32_canvas[2];
	TH1F *Ratio_GenJet_R32_Clone[2]; 
	//TCanvas *Ratio_CaloJet_R32 = new TCanvas("3jets2jetscalo","3jets and 2 jets calo",0,0,800,1200);
	TH1F *Ht_R2[2];
	TH1F *Ht_R3[2];
	Ht_R2[0] = (TH1F*)mainroot->Get("Ht_hist/Ht_GenJet_R2");
  	Ht_R3[0] = (TH1F*)mainroot->Get("Ht_hist/Ht_GenJet_R3");
	Ht_R2[1] = (TH1F*)mainroot->Get("Ht_hist/Ht_CaloJet_R2");
  	Ht_R3[1] = (TH1F*)mainroot->Get("Ht_hist/Ht_CaloJet_R3");
            //Canvas_Rap_Reso_All_Y_[i]->Divide(1,3);
	for(int i=0; i<2; i++)
	{
		sprintf(canvas_name,"%s",Ratio_32_canvas_name[i].c_str());
		sprintf(canvas_title,"%s",Ratio_32_canvas_title[i].c_str());
		Ratio_R32_canvas[i] = new TCanvas(canvas_name,canvas_title,0,0,1200,1500);

		Ratio_R32_canvas[i]->Divide(1,2);
		Ratio_R32_canvas[i]->cd(1);
		gPad->SetGrid();
		gPad->SetLogy();
		gPad->SetLogx();
	
		Ht_R2[i]->SetMarkerStyle(23);
		Ht_R2[i]->SetLineColor(kRed);
		Ht_R2[i]->SetMarkerColor(kRed);
		Ht_R2[i]->GetXaxis()->SetRangeUser(70,2000);
      	  //Ht_R2[i]->SetMinimum(1);
		Ht_R2[i]->SetTitle("gen HT dist. jets #geq 3 and #geq 2");
		Ht_R2[i]->GetXaxis()->SetTitle(Ratio_32_X_Axis[i].c_str());
		Ht_R2[i]->Draw("E1");

		Ht_R3[i]->SetMarkerStyle(24);
		Ht_R3[i]->SetMarkerColor(kBlue);
		Ht_R3[i]->SetLineColor(kBlue);
		Ht_R3[i]->Draw("SAME E1");

       		Ratio_R32_canvas[i]->cd(2);
       		gPad->SetGrid();

       		Ratio_GenJet_R32_Clone[i]  = (TH1F*)Ht_R3[i]->Clone(Ratio_32_Y_Axis[i].c_str());
        	Ratio_GenJet_R32_Clone[i]->SetLineColor(kBlue);
        	Ratio_GenJet_R32_Clone[i]->SetMinimum(0);
        	Ratio_GenJet_R32_Clone[i]->SetMaximum(1);
        	Ratio_GenJet_R32_Clone[i]->GetXaxis()->SetRangeUser(0,2000);
        	Ratio_GenJet_R32_Clone[i]->Divide(Ht_R2[i]);
        	Ratio_GenJet_R32_Clone[i]->SetMarkerStyle(21);
        	Ratio_GenJet_R32_Clone[i]->SetTitle("Ratio of Njets>3 over Njets>2");
        	Ratio_GenJet_R32_Clone[i]->GetXaxis()->SetTitle(Ratio_32_X_Axis[i].c_str());
        	Ratio_GenJet_R32_Clone[i]->GetYaxis()->SetTitle(Ratio_32_Y_Axis[i].c_str());
        	Ratio_GenJet_R32_Clone[i]->Draw();

	gPad->BuildLegend();
	sprintf(canvas_print,"/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/%s",Ratio_32_canvas_name[i].c_str());
        Ratio_R32_canvas[i]->Print(canvas_print);
	}
	
	TH1F *R32_GenJet_over_CaloJet_Clone; 
	TCanvas *R32_GenJet_over_CaloJet_canvas = new TCanvas("R32_GenJet_over_CaloJet","R32_GenJet_over_CaloJet",0,0,800,1200);

	R32_GenJet_over_CaloJet_canvas->Divide(1,2);
	R32_GenJet_over_CaloJet_canvas->cd(1);
	gPad->SetGrid();

	Ratio_GenJet_R32_Clone[0]->SetLineColor(kBlue);
	Ratio_GenJet_R32_Clone[0]->Draw();
	Ratio_GenJet_R32_Clone[0]->SetTitle("R32 for GenJet and CaloJet");

	Ratio_GenJet_R32_Clone[1]->SetLineColor(kRed);
	Ratio_GenJet_R32_Clone[1]->Draw("SAME");
	Ratio_GenJet_R32_Clone[1]->SetTitle("R32 for GenJet and CaloJet");
	//gPad->BuildLegend();
	Ratio_GenJet_R32_Clone[0]->SetTitle("R32 for GenJet and CaloJet");
	R32_GenJet_over_CaloJet_canvas->cd(2);
	gPad->SetGrid();

	R32_GenJet_over_CaloJet_Clone = (TH1F*)Ratio_GenJet_R32_Clone[1]->Clone("Ratio CaloJet(R32) over GenJet(R32)");
	R32_GenJet_over_CaloJet_Clone->SetLineColor(kCyan);
	R32_GenJet_over_CaloJet_Clone->Divide(Ratio_GenJet_R32_Clone[0]);
	R32_GenJet_over_CaloJet_Clone->GetXaxis()->SetTitle("HT(GeV)");
	R32_GenJet_over_CaloJet_Clone->GetYaxis()->SetTitle("CaloJet(R32) over GenJet(R32)");
	R32_GenJet_over_CaloJet_Clone->SetMinimum(0);
        R32_GenJet_over_CaloJet_Clone->SetMaximum(2);
	R32_GenJet_over_CaloJet_Clone->GetXaxis()->SetRangeUser(0,2000);
	R32_GenJet_over_CaloJet_Clone->SetTitle("CaloJet(R32) over GenJet(R32)");
	R32_GenJet_over_CaloJet_Clone->Draw();
	R32_GenJet_over_CaloJet_canvas->Print("/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/R32_Ratio_GenJet_over_CaloJet");

//-----------------------------------------------------------------------------------------
//PT Resolution ------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
/*	
	TH2F     *PT_Reso_hist[4][4];
	TCanvas  *PT_Reso_canvas[4][4];
	TProfile *PT_Reso_Prof[4][4];
	TProfile *PT_Reso_ProfY[4][4];
	//TH1D     *PT_Reso_Projection[4][4];
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{

		 //sprintf(hist_name,  "PT_Reso_%s_%s" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());

	   	 sprintf(canvas_name,  "PT_Reso_%s_%s" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());
	   	 sprintf(canvas_title, "PT_Reso_%s_%s" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());
		 sprintf(get_hist_name,  "PT_Resolution_TH2F/PT_Reso_%s_%s" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());

		 PT_Reso_hist[i][j]   = (TH2F*)mainroot->Get(get_hist_name);
	   	 PT_Reso_canvas[i][j] = new TCanvas(canvas_name, canvas_title, 0,0,800,1200);


		// PT_Reso_canvas[i][j]->SetFillColor(41);
         	 PT_Reso_canvas[i][j]->Divide(1,3);
                 PT_Reso_canvas[i][j]->cd(1);

	  	 PT_Reso_hist[i][j]->SetTitle("Pt Resolution of the Leading Jet");
            	 PT_Reso_hist[i][j]->GetXaxis()->SetTitle("PT of the Leading jet");
            	 PT_Reso_hist[i][j]->GetYaxis()->SetTitle("Resolution");
            	 PT_Reso_hist[i][j]->Draw("COLZ");

     	   	 PT_Reso_canvas[i][j]->cd(2);

	         PT_Reso_Prof[i][j] = PT_Reso_hist[i][j]->ProfileX();
        	 PT_Reso_Prof[i][j]->SetMinimum(0);
            	 PT_Reso_Prof[i][j]->SetMaximum(0.5);
		 if(j >= 2)PT_Reso_Prof[i][j]->GetXaxis()->SetRangeUser(0,400);
	         else PT_Reso_Prof[i][j]->GetXaxis()->SetRangeUser(0,1000);				 
       		 PT_Reso_Prof[i][j]->SetTitle("Pt Resolution of the Leading Jet (Profile)");
            	 PT_Reso_Prof[i][j]->GetXaxis()->SetTitle("PT of the Leading jet");
            	 PT_Reso_Prof[i][j]->GetYaxis()->SetTitle("Resolution");
            	 PT_Reso_Prof[i][j]->Draw("E1");
//][ hist
		 PT_Reso_canvas[i][j]->cd(3);
//PT_Reso_ProfY
		 PT_Reso_ProfY[i][j] = PT_Reso_hist[i][j]->ProfileY();
    	  	// PT_Reso_Projection[i][j]->SetMinimum(0);
            	//PT_Reso_Projection[i][j]->SetMaximum(0.5);
		// if(j >= 2)PT_Reso_Projection[i][j]->GetXaxis()->SetRangeUser(0,400);
	        //else PT_Reso_Projection[i][j]->GetXaxis()->SetRangeUser(0,1000);				 
       		 PT_Reso_ProfY[i][j]->SetTitle("Pt Resolution of the Leading Jet (Profile)");
            	 PT_Reso_ProfY[i][j]->GetXaxis()->SetTitle("PT of the Leading jet");
            	 PT_Reso_ProfY[i][j]->GetYaxis()->SetTitle("Resolution");
            	 PT_Reso_ProfY[i][j]->Draw();
PT_Reso_Projection[i][j] = PT_Reso_hist[i][j]->ProfileY();
    	  	// PT_Reso_Projection[i][j]->SetMinimum(0);
            	//PT_Reso_Projection[i][j]->SetMaximum(0.5);
		// if(j >= 2)PT_Reso_Projection[i][j]->GetXaxis()->SetRangeUser(0,400);
	        //else PT_Reso_Projection[i][j]->GetXaxis()->SetRangeUser(0,1000);				 
       		 PT_Reso_Projection[i][j]->SetTitle("Pt Resolution of the Leading Jet (Profile)");
            	 PT_Reso_Projection[i][j]->GetXaxis()->SetTitle("PT of the Leading jet");
            	 PT_Reso_Projection[i][j]->GetYaxis()->SetTitle("Resolution");
            	 PT_Reso_Projection[i][j]->Draw();
	
	        sprintf(canvas_print,"/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/PT_Reso_%s_%s",Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());
    	        PT_Reso_canvas[i][j]->Print(canvas_print);	    	

		}
	}

*/	TH2F     *PT_Reso_hist[4][4];
	TH1 *PT_Reso_ProjectionY;
	//TH1 *PT_Reso_ProjectionY;
	TH1D *PT_Reso_ProfileY;
	int nbins=0;
	//double mean =0;
	//double mean2 =0;
	//double bincontent = 0;	
	//double sigma=0;
	TCanvas *PT_Reso_RMS;
	TCanvas *PT_Reso_Mean;
	TGraph *gr_Mean;
	TGraph *gr_RMS;
	vector<double> mean;
	vector<double> rms;
	vector<double> pt_mean_rms;
	double rms_control=0;
	double mean_control=0;

	double rms_control_prof=0;
	double mean_control_prof=0;

	int bin_counter=0;

	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{

		 //sprintf(hist_name,  "PT_Reso_%s_%s" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());

	   	sprintf(canvas_name,  "PT_Reso_%s_%s_mean_rms" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());
	        sprintf(canvas_title, "PT_Reso_%s_%s_mean_rms" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());

		 sprintf(get_hist_name,  "PT_Resolution_TH2F/PT_Reso_%s_%s" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());
		 PT_Reso_hist[i][j]   = (TH2F*)mainroot->Get(get_hist_name);

		 
					nbins = PT_Reso_hist[i][j]->GetYaxis()->GetNbins();
					mean_control=0;
					rms_control=0;
					mean.clear();
					rms.clear();
					bin_counter++;

					for (int k=1; k <= nbins; k++)
					{	
						PT_Reso_ProjectionY = PT_Reso_hist[i][j]->ProjectionY("",k,k);
						PT_Reso_ProfileY    = PT_Reso_hist[i][j]->ProfileY();

						mean_control  = PT_Reso_ProjectionY->GetMean();
						rms_control   = PT_Reso_ProjectionY->GetRMS();

						//mean_control_prof  = PT_Reso_ProfileY->GetMean();
						//rms_control_prof   = PT_Reso_ProfileY->GetRMS();
						
						//printf(" Bins %d, Mean Prof = %g, RMS Prof = %g\n",k,mean_control_prof,rms_control_prof);
						if(mean_control !=0 && rms_control !=0)
						{
						bin_counter++;
						mean.push_back((double )PT_Reso_ProjectionY->GetMean()); 
						rms.push_back((double) PT_Reso_ProjectionY->GetRMS());  						
						pt_mean_rms.push_back((double)k);

						//printf(" Bins %d, Mean = %g, Sigma = %g\n",k,mean[k],rms[k]);
						
						}
      						delete PT_Reso_ProjectionY; 
						delete PT_Reso_ProfileY;
						
					}
					PT_Reso_RMS = new TCanvas(canvas_name,canvas_title,0,0,800,1200);
					gr_RMS  = new TGraph(bin_counter,&(pt_mean_rms[0]),&(mean[0]));
					gr_RMS->Draw("AC*");
				     	sprintf(canvas_print,"/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/PT_Reso_%s_%s_rms",Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());
					PT_Reso_RMS->Print(canvas_print);

					PT_Reso_Mean = new TCanvas(canvas_name,canvas_title,0,0,800,1200);
					gr_Mean = new TGraph(bin_counter,&(pt_mean_rms[0]),&(rms[0]));
					gr_Mean->Draw("AC*");	
					sprintf(canvas_print,"/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/PT_Reso_%s_%s_mean",Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());
					PT_Reso_Mean->Print(canvas_print);

					//delete PT_Reso_RMS;
					//delete PT_Reso_Mean;
					//delete gr_Mean;
					//delete gr_RMS;

		}
	}

//TProfile *profile = new TProfile((const char*)Pplot,"Profile plot",nbins,-5,5,"S"); 
//for(j= 0, j < max; j++){ 
//    chain->GetEntry(j); 
 //   profile->Fill(w->at(1),f->at(1)); 
//} 
//TH1F *plot = new TH1F((const char*)Ddistribution,"Distribution",nbins,-5,5); 
//for(Int_t binN = 1, binN <= nbins, binN++){ 
 //   Double_t rms = profile->GetBinError(binN); 
 //   plot->Fill(plot->GetBinCenter(binN), rms); 
//}
	 
cout << " *************************PT Resolution Grafikleri Bitti************************* " << endl;
cout << " *************************PT Resolution Grafikleri Bitti************************* " << endl;
cout << " *************************PT Resolution Grafikleri Bitti************************* " << endl;


//-----------------------------------------------------------------------------------------
//HT Resolution ------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

	int jet_num2[] = {1,2,3,4};
	TH2F     *HT_Reso_hist[4][4];
	TCanvas  *HT_Reso_canvas[4][4];
	TProfile *HT_Reso_Prof[4][4];

	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{

		 //sprintf(hist_name,  "PT_Reso_%s_%s" ,Pt_Reso_Cut_Name[i].c_str(),JET_NAME[j].c_str());

	   	 sprintf(canvas_name,  "HT_Reso_Jet_Number_%i_%s" ,jet_num2[i],Pt_Reso_Cut_Name[j].c_str());
	   	 sprintf(canvas_title, "HT_Reso_Jet_Number_%i_%s" ,jet_num2[i],Pt_Reso_Cut_Name[j].c_str());
		 sprintf(get_hist_name,  "HT_Resolution_TH2F/HT_Reso_Jet_Number_%i_%s" ,jet_num2[i],Pt_Reso_Cut_Name[j].c_str());

		 HT_Reso_hist[i][j]   = (TH2F*)mainroot->Get(get_hist_name);
	   	 HT_Reso_canvas[i][j] = new TCanvas(canvas_name, canvas_title, 0,0,800,1200);


		 //HT_Reso_canvas[i][j]->SetFillColor(41);
         	 HT_Reso_canvas[i][j]->Divide(1,2);
                 HT_Reso_canvas[i][j]->cd(1);

	  	 HT_Reso_hist[i][j]->SetTitle("Pt Resolution of the Leading Jet");
            	 HT_Reso_hist[i][j]->GetXaxis()->SetTitle("PT of the Leading jet");
            	 HT_Reso_hist[i][j]->GetYaxis()->SetTitle("Resolution");
            	 HT_Reso_hist[i][j]->Draw("COLZ");

     	   	 HT_Reso_canvas[i][j]->cd(2);


	         HT_Reso_Prof[i][j] = HT_Reso_hist[i][j]->ProfileX();
		 HT_Reso_Prof[i][j]->GetXaxis()->SetRangeUser(0,2000);
        	 HT_Reso_Prof[i][j]->SetMinimum(0);
            	 HT_Reso_Prof[i][j]->SetMaximum(0.2);
            	 HT_Reso_Prof[i][j]->SetTitle("Ht Resolution of the Leading Jet (Profile)");
            	 HT_Reso_Prof[i][j]->GetXaxis()->SetTitle("HT of the Leading jet");
                 HT_Reso_Prof[i][j]->GetYaxis()->SetTitle("Resolution");
		// gPad->SetGrid();
		 //gPad->SetLogy();
		 //gPad->SetLogx();
            	 HT_Reso_Prof[i][j]->Draw("E1");


	        sprintf(canvas_print,"/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis/Analysis/Analysis_Histograms/HT_Reso_Jet_Number_%i_%s" ,jet_num2[i],Pt_Reso_Cut_Name[j].c_str());
    	        HT_Reso_canvas[i][j]->Print(canvas_print);	    	

		}
	}

}
 

