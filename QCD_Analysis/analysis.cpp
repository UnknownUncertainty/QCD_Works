#include <TSystem.h>
#include <TChain.h>
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TEventList.h>
#include "TClonesArray.h"
#include <TH1D.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <TMath.h>
#include <TProfile.h>
#include "TStyle.h"
#include <time.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "Utilities.h"
#include "TVectorD.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "DelphesClasses.h"
#include "skimming_events.h"

#include <algorithm>
#include <time.h>
#include "TLatex.h"
#include "TLegend.h"

//#include "Math/GenVector/VectorUtil.h"



using namespace std;

const int PT_BIN_ARRAY[]={0,50,100,150,250,500,750,1000,1500,2000,3000,4000,5000,6500};
const int JET_NUMBER = 0;
//char JET_NAME[] = {'a','b','c','d'};


//JET_NAME[1].push_back("Second_Jet");
//JET_NAME[2].push_back("Third_Jet");
//JET_NAME[3].push_back("Fourth_Jet");
//"Leading_Jet","Second_Jet","Third_Jet","Fourth_Jet"};

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

int main(int argc, char*argv[])
{
	//----------------------------------------------------------------------------------------
	/////////////////////////////////////////////  Read Input Parameters  ////////////////////
	//----------------------------------------------------------------------------------------

	if (argc < 2) {printf("******Error in input parameters\n");return 1;}

	CommandLine c1;
	c1.parse(argc,argv);
	gROOT->ProcessLine("#include >vector<");

    string InputFileName	= c1.getValue<string>	("InputFileName");
    string OutputFileName	= c1.getValue<string>	("OutputFileName");
    //string OutputFileTag	= c1.getValue<string>	("OutputFileTag");
    //string JetAlgo		= c1.getValue<string>	("JetAlgo");
    //vector<string> HLTbit	= c1.getVector<string> 	("HLTbit","");
    int NENTRIES			= c1.getValue<int>		("NEntries");
    //bool DATA				= c1.getValue<bool>		("DATA");
   // bool SaveData			= c1.getValue<bool>		("SaveData");
    //double MH2			= c1.getValue<double>	("MH2");
    //double cross			= c1.getValue<double>	("CrossSection");
    //double efficiency		= c1.getValue<double>	("Efficiency");

    //double MUONPT_CUT=23.;
    //double ELECPT_CUT=33.;
    //double MET_CUT=40.;

	// input parametrelerini ekrana yazdiriyor
	if (!c1.check()) return 0;
	c1.print(); // Printing the options

    string outputfile = OutputFileName + ".root";
	// output dosyasini olustur.
	TFile *outf = new TFile(outputfile.c_str(),"RECREATE");


	cout << "________________________________________________________________\n";
	cout << "\n";
	cout << "time start  " << endl;
	gSystem->Exec("date '+%H:%M:%S'");


	//----------------------------------------------------------------------------------------
	///////////////////////////////////////////////  Histogram Output  ///////////////////////
	//----------------------------------------------------------------------------------------
	// Book histograms


	//----------------------------------------------------------------------------------------
	//output dosyasinda jets adinda bir klasor olustur
	//jetdir adinda olusan dizine git ve asagidaki histogramlari, ilgili degiskenleri uyarinca tanimla
	//----------------------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////
//////////Eta'sı 2.5 ' den küçük Jetlerin Pt Kesimine Göre Grafikleri//
///////////////////////////////////////////////////////////////////////

/*	char jet_ptcut_etacut_hist_title[100];
	char jet_ptcut_etacut_hist_name[100];
 	for(int b=0; b<4;b++)
	{
		for (int a =0; a < 4; a++)
		{
		sprintf(jet_ptcut_etacut_hist_title,"Ptcut-%i%i Jet Eta < 2.5",a,b);
		sprintf(jet_ptcut_etacut_hist_name ,"Jet_ptcut_etacut_%i%i",a,b);
		jet_ptcut_etacut_[a][b]         = new TH1F(jet_ptcut_etacut_hist_name ,   jet_ptcut_etacut_hist_title , 100, 0, 350);
		}
	}

*/
//vector<char> hist_name(100);
//vector<char> hist_title(1000);
/*for(int i=0; i<4; i++){
cout << JET_NAME[i] << endl;
}
*/
char hist_name[100];
char hist_title[1000];
string JET_NAME[] = {"Leading_Jet","Second_Jet","Third_Jet","Fourth_Jet"};
string Cut_Name[] = {"No_Cut","Eta_Cut","R_Cut","R_Eta_Cut"};

///////////////////////////////////////////////////////////////////////////////////////////
//////////GenJet and CaloJet Ht Distribution //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	TDirectory *Ht_hist= outf->mkdir("Ht_hist");
	Ht_hist->cd();

	TH1F  *hist_genjetht  = new TH1F("Ht_GenJet_All"     , "Ht_GenJet_All ",   500.0, 0.0, 10000.0);
	TH1F  *hist_calojetht = new TH1F("Ht_CaloJet_All"    , "Ht_CaloJet_All",   500.0, 0.0, 10000.0);

	TH1F  *hist_genjet_R2   = new TH1F("Ht_GenJet_R2"    , "Ht_GenJet_R2 ",   500.0, 0.0, 10000.0);
	TH1F  *hist_genjet_R3   = new TH1F("Ht_GenJet_R3"    , "Ht_GenJet_R3 ",   500.0, 0.0, 10000.0);
	TH1F  *hist_calojet_R2  = new TH1F("Ht_CaloJet_R2"   , "Ht_GenJet_R2 ",   500.0, 0.0, 10000.0);
	TH1F  *hist_calojet_R3  = new TH1F("Ht_CaloJet_R3"   , "Ht_GenJet_R3 ",   500.0, 0.0, 10000.0);
	
	//hist_genjet_R3

//////////////////////////////////////////////////////////////////////////////////////////
//--------------------------Rapidity Resolution for PT Intervals------------------------//
//////////////////////////////////////////////////////////////////////////////////////////
	
	//------------For Each Eta Intervals Between 0 - 2.5, PT Intervals(0-6000)------//

	//Satırlar PT Aralıklarını ve Sütunlar ise Rapidity Aralıklarını Göstermekedir.--/

	TDirectory *Rap_Resolution = outf->mkdir("Rap_Resolution");
	Rap_Resolution->cd();
 	TH1F *hist_deltaY[14][5];

	for (int i=0; i<12; i++)
	{
	  for (int j=0; j<5; j++)
	  {
	    sprintf(hist_name, "DeltaY_%ipt%i_y%i" ,PT_BIN_ARRAY[i] ,PT_BIN_ARRAY[i+1] ,j);
	    sprintf(hist_title,"%i<Pt<%i, Rapidity resolution for %.1f < y < %.1f",PT_BIN_ARRAY[i],PT_BIN_ARRAY[i+1],j*0.5,(j+1)*0.5);
	    hist_deltaY[i][j]=new TH1F(hist_name, hist_title, 100, -0.2, +0.2 );
	  }
	}

	//------------For All Eta Values, PT Intervals(0-6000)---------------------------//
	TDirectory *Rap_Reso_TH2F = outf->mkdir("Rap_Resolution_TH2F");
	Rap_Reso_TH2F->cd();

	TH2F *hist_DeltaY_vs_GenY[14];
	for(int i=0; i<13; i++)
	{
	    sprintf(hist_name, "DeltaY_vs_GenY_%ipt%i_y%i" ,PT_BIN_ARRAY[i] ,PT_BIN_ARRAY[i+1] ,i);
	    sprintf(hist_title,"%i<Pt<%i, Rapidity Reso all Y" ,PT_BIN_ARRAY[i],PT_BIN_ARRAY[i+1]);
	    hist_DeltaY_vs_GenY[i] = new TH2F(hist_name, hist_title, 150, -5.0, 5.0 , 150.0, -0.5, 0.5 );
	}

//////////////////////////////////////////////////////////////////////////////////////////
//--------------------------PT Resolution for Each Jet ---------------------------------//
//////////////////////////////////////////////////////////////////////////////////////////
	//hist_name.clear();
	//hist_title.clear();
	
	TDirectory *PT_Resolution_TH2F= outf->mkdir("PT_Resolution_TH2F");
	PT_Resolution_TH2F->cd();

	//Satırlar Yapılan Cutları ve Sütunlar İse Kaçıncı Jet Olduğunu Göstermektedir.-//
	TH2F *hist_PT_Reso[4][4];
	for(int i=0; i<4; i++)
	{	
		for(int j=0; j<4; j++)
		{
	   	sprintf(hist_name,  "PT_Reso_%s_%s" ,Cut_Name[i].c_str(),JET_NAME[j].c_str());
	   	sprintf(hist_title, "%s_DeltaPT_over_GenPT_vs_GenPT_for_%s" ,Cut_Name[i].c_str(),JET_NAME[j].c_str());
                hist_PT_Reso[i][j] = new TH2F(hist_name, hist_title, 250.0, 0.0, 6500.0 ,250.0, -0.5, 0.5);	    	    	
		}	
	}
//////////////////////////////////////////////////////////////////////////////////////////
//--------------------------HT Resolution for Each Jet ---------------------------------//
//////////////////////////////////////////////////////////////////////////////////////////
	//hist_name.clear();
	//hist_title.clear();
	int jet_num2[] = {1,2,3,4};
	TDirectory *HT_Resolution_TH2F= outf->mkdir("HT_Resolution_TH2F");
	HT_Resolution_TH2F->cd();

	//Satırlar Jet Sayısını Belirtir ve Sütunlar İse Yapılan Cutları Göstermektedir.-//
	TH2F *hist_HT_Reso[4][4];
	for(int i=0; i<4; i++)
	{	
		for(int j=0; j<4; j++)
		{
	   	sprintf(hist_name,  "HT_Reso_Jet_Number_%i_%s" ,jet_num2[i],Cut_Name[j].c_str());
	   	sprintf(hist_title, "HT_Resolution_for_Jet_Number >= %i and %s PT_Cut 50 GeV" ,jet_num2[i],Cut_Name[j].c_str());
                hist_HT_Reso[i][j] = new TH2F(hist_name, hist_title, 250.0, 0.0, 6500.0 ,250.0, -0.5, 0.5);	    	    	
		}	
	}





	//----------------------------------------------------------------------------------------
	////////////////////////////////////////  My Data STructure  /////////////////////////////
	//----------------------------------------------------------------------------------------

	TH1F::SetDefaultSumw2(true);

	gSystem->Load("libExRootAnalysis");
	gSystem->Load("libDelphes");


	// Create chain of root trees
	TChain chain("Delphes");
	//chain.Add(InputFileName.c_str());
	char filename[1000];
	FILE *input;
	input = fopen(InputFileName.c_str(),"r");
	if (input != NULL)
	{
		// lets read each line and get the filename from it
		while (fscanf(input,"%s\n",filename) != EOF)
		{
			printf("%s\n",filename);
			chain.Add(filename);
		}
	}


	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Get pointers to branches used in this analysis
	TClonesArray *branchJet    = treeReader->UseBranch("Jet");
 	TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
 	//TClonesArray *branchSHT    = treeReader->UseBranch("ScalarHT");


	//----------------------------------------------------------------------------------------
	/////////////////////////////////////  LOOP  Over the EVENTS  ////////////////////////////
	//----------------------------------------------------------------------------------------

	cout << "Set Branch Addresses" << endl;

	int decade = 0;
	unsigned entries = 0;

	if (NENTRIES==-1) entries=numberOfEntries;
	else entries=NENTRIES;

	cout << "Reading TREE: " << numberOfEntries << " events available, \n";
	cout << "\t\t" << NENTRIES << " of them will be analyzed." << endl;



//Jet *genjet;
//Jet *calojet;
vector<Jet *> mycalojet;
vector<Jet *> mygenjet;


//TLorentzVector jet_2[2];
//TLorentzVector jet_3[2];
//TLorentzVector jet_4[2];
//vector<Jet *> m2Jets;
//vector<Jet *> m3Jets;
//vector<Jet *> m4Jets;


//GenJet *genjet_2[25];
//GenJet *genjet_3[25];

	// Loop over all events
	for(Long64_t entry = 0; entry < entries; ++entry)
	{

		double progress = 10.0*entry/(1.0*entries);
		int k = TMath::FloorNint(progress);
		if (k > decade) {   cout << 10*k << "%\t"; gSystem->Exec("date '+%H:%M:%S'"); cout << endl;	}
		decade = k;

		// ROOT dosyalari icindeki her bir veri dallarini okuyor.
		treeReader->ReadEntry(entry);


//#################################################################################################################################################
//#################################################################################################################################################


	///////////////////////////////////////////////////////////////////////////////
	//////////////////////////////Jet Grafikleri///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////

	//CaloJet İçin Ht Hesaplama ve -----------------------------------------------
	double calo_jet_ht_sum = 0.0;
	int jet_counter=0;
	mycalojet.clear();

	if(branchJet->GetEntries() > 0) // Bütün Jetler İçin
	{
		
		for (int i=0; i< branchJet->GetEntries(); i++)
		{
		
		mycalojet.push_back((Jet*) branchJet->At(i));
		calo_jet_ht_sum = calo_jet_ht_sum + mycalojet[i]->PT;

			if(mycalojet[i]->PT > 50 && abs( mycalojet[i]->Eta ) <= 2.5)
			{	
				jet_counter++;				
			}
					
		}
		hist_calojetht->Fill(calo_jet_ht_sum);
		if(jet_counter >=2)hist_calojet_R2->Fill(calo_jet_ht_sum);
		if(jet_counter >=3)hist_calojet_R3->Fill(calo_jet_ht_sum);

	}
	

	//GenJet İçin Ht Hesaplama----------------------------------------------------
	double gen_jet_ht_sum = 0.0;
	jet_counter=0;

	mygenjet.clear();

	if(branchGenJet->GetEntries() > 0) // Bütün Jetler İçin
	{		
		for (int i=0; i< branchGenJet->GetEntries(); i++)
		{		
		mygenjet.push_back((Jet*) branchGenJet->At(i));
		gen_jet_ht_sum = gen_jet_ht_sum + mygenjet[i]->PT;

			if(mygenjet[i]->PT > 50 && abs( mygenjet[i]->Eta ) <= 2.5)
			{	
				jet_counter++;
			}
		}
		hist_genjetht->Fill(gen_jet_ht_sum);
		if(jet_counter >=2)hist_genjet_R2->Fill(gen_jet_ht_sum);
		if(jet_counter >=3)hist_genjet_R3->Fill(gen_jet_ht_sum);
		
	}	


	//---------------------------------------------------------------------
	//Çözünürlük Grafikleri------------------------------------------------
	//---------------------------------------------------------------------


	//-----------Eta Çözünürlüğü----------------------------------------------------------------------------------------------
	mycalojet.clear();
	mygenjet.clear();
	double delta_Phi = 0.0;
	double delta_Eta = 0.0;
	double delta_R = 0.0;
	//double pt_reso_etacut = 0.0;
	double pt_reso = 0.0;
	//int counterjet = 0;

	jet_counter=0;
	calo_jet_ht_sum = 0.0;
	gen_jet_ht_sum  = 0.0;
	double calo_jet_ht_sum_Eta = 0.0;
	double gen_jet_ht_sum_Eta  = 0.0;
	double calo_jet_ht_sum_delta_R = 0.0;
	double gen_jet_ht_sum_delta_R  = 0.0;
	double calo_jet_ht_sum_Eta_delta_R = 0.0;
	double gen_jet_ht_sum_Eta_delta_R  = 0.0;

	double ht_reso = 0.0;
	double ht_reso_Eta = 0.0;
	double ht_reso_delta_R = 0.0;
	double ht_reso_Eta_delta_R = 0.0;
	




	if(branchGenJet->GetEntries() > JET_NUMBER && branchJet->GetEntries() > JET_NUMBER && branchGenJet->GetEntries() == branchJet->GetEntries() )
	{

		for (int i=0; i< branchGenJet->GetEntries(); i++)
		{
		mycalojet.push_back((Jet*) branchJet->At(i));
		mygenjet.push_back((Jet*) branchGenJet->At(i));

		jet_counter++;

		  for (int j=0; j <12; j++)
	          {
	            if(mygenjet[i]-> PT > PT_BIN_ARRAY[j] && mygenjet[i]-> PT < PT_BIN_ARRAY[j+1] && mycalojet[i]-> PT > PT_BIN_ARRAY[j] && mycalojet[i]-> PT < PT_BIN_ARRAY[j+1])
	            {
	              if (abs(mygenjet[i]->Eta)<0.5) hist_deltaY[j][0]->Fill(mygenjet[i]->Eta - mycalojet[i]->Eta);
	              if (abs(mygenjet[i]->Eta)>0.5 && abs(mygenjet[i]->Eta)<1.0) hist_deltaY[j][1]->Fill(mygenjet[i]->Eta - mycalojet[i]->Eta);
	              if (abs(mygenjet[i]->Eta)>1.0 && abs(mygenjet[i]->Eta)<1.5) hist_deltaY[j][2]->Fill(mygenjet[i]->Eta - mycalojet[i]->Eta);
	              if (abs(mygenjet[i]->Eta)>1.5 && abs(mygenjet[i]->Eta)<2.0) hist_deltaY[j][3]->Fill(mygenjet[i]->Eta - mycalojet[i]->Eta);
	              if (abs(mygenjet[i]->Eta)>2.0 && abs(mygenjet[i]->Eta)<2.5) hist_deltaY[j][4]->Fill(mygenjet[i]->Eta - mycalojet[i]->Eta);

		      hist_DeltaY_vs_GenY[j]->Fill( mygenjet[i]->Eta, mygenjet[i]->Eta - mycalojet[i]->Eta);
	            }
	          }
	
			//-----PT Resolution Histogramları İÇin---------------------------------//
		
			delta_Phi = abs(mygenjet[i]->Phi - mycalojet[i]->Phi);
			delta_Eta = abs(mygenjet[i]->Eta - mycalojet[i]->Eta);
			delta_R = sqrt( pow(delta_Phi,2) + pow(delta_Eta,2));
			pt_reso = ( mygenjet[i]->PT - mycalojet[i]->PT ) / (mygenjet[i]->PT);
			

			if(i<=3)     hist_PT_Reso[0][i]->Fill(mygenjet[i]->PT, pt_reso);
			if(i<=3 && abs(mycalojet[i]->Eta) < 2.5 && abs(mygenjet[i]->Eta) < 2.5 ) hist_PT_Reso[1][i]->Fill(mygenjet[i]->PT, pt_reso);
			if(delta_R < 0.25 && i<=3)    hist_PT_Reso[2][i]->Fill(mygenjet[i]->PT, pt_reso);
			if(delta_R < 0.25 && i<=3 && abs(mycalojet[i]->Eta) < 2.5 && abs(mygenjet[i]->Eta) < 2.5) hist_PT_Reso[3][i]->Fill(mygenjet[i]->PT, pt_reso);

			//-----HT Hesapla For Each Cut PT_CUT 50 GeV --------------------------------//
			
			if(mygenjet[i]->PT > 50) gen_jet_ht_sum  = gen_jet_ht_sum  + mygenjet[i]->PT;
			if(mycalojet[i]->PT > 50)calo_jet_ht_sum = calo_jet_ht_sum + mycalojet[i]->PT;

			if(mygenjet[i]->PT > 50  && abs(mygenjet[i]->Eta)  < 2.5)  gen_jet_ht_sum_Eta   = gen_jet_ht_sum_Eta  + mygenjet[i]->PT;
			if(mycalojet[i]->PT > 50 && abs(mycalojet[i]->Eta) < 2.5)  calo_jet_ht_sum_Eta  = calo_jet_ht_sum_Eta + mycalojet[i]->PT;

			if(mygenjet[i]->PT > 50 && delta_R < 0.25 ) gen_jet_ht_sum_delta_R  = gen_jet_ht_sum_delta_R  + mygenjet[i]->PT;
			if(mycalojet[i]->PT > 50 && delta_R < 0.25 ) calo_jet_ht_sum_delta_R = calo_jet_ht_sum_delta_R + mycalojet[i]->PT;

			if(mygenjet[i]->PT > 50 && delta_R < 0.25 && abs(mygenjet[i]->Eta) < 2.5) gen_jet_ht_sum_Eta_delta_R  = gen_jet_ht_sum_Eta_delta_R  + mygenjet[i]->PT;
			if(mycalojet[i]->PT > 50 && delta_R < 0.25 && abs(mycalojet[i]->Eta) < 2.5)calo_jet_ht_sum_Eta_delta_R = calo_jet_ht_sum_Eta_delta_R + mycalojet[i]->PT;



		
	 	}


			//-----HT Histogramlar ve Resolution Hesapla For Each Cut --------------------------------//
			ht_reso             = abs((gen_jet_ht_sum             - calo_jet_ht_sum))             / (gen_jet_ht_sum);
			ht_reso_Eta         = abs((gen_jet_ht_sum_Eta         - calo_jet_ht_sum_Eta ))        / (gen_jet_ht_sum_Eta );
			ht_reso_delta_R     = abs((gen_jet_ht_sum_delta_R     - calo_jet_ht_sum_delta_R)     / (gen_jet_ht_sum_delta_R));
			ht_reso_Eta_delta_R = abs((gen_jet_ht_sum_Eta_delta_R - calo_jet_ht_sum_Eta_delta_R) / (gen_jet_ht_sum_Eta_delta_R));
				for(int i=0; i<4; i++)
				{	
					if(jet_counter >= (i+1))     hist_HT_Reso[i][0]->Fill(gen_jet_ht_sum, ht_reso);
					if(jet_counter >= (i+1))     hist_HT_Reso[i][1]->Fill(gen_jet_ht_sum_Eta, ht_reso_Eta);	
					if(jet_counter >= (i+1))     hist_HT_Reso[i][2]->Fill(gen_jet_ht_sum_delta_R, ht_reso_delta_R);			
					if(jet_counter >= (i+1))     hist_HT_Reso[i][3]->Fill(gen_jet_ht_sum_Eta_delta_R, ht_reso_Eta_delta_R);
				}
				
	}



	

}





	//end of event loop

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////
	//---------------------------------------------------------------------------------------------------------
	// end of event loop

	//myfile.close();
	outf->Write();
	outf->Close();

	//end of main loop
}
