#include "makePlots.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TImage.h"

//Constructor
makePlots::makePlots(TChain* inchain):Chain1(inchain) 
{
	readmap();
	cout << "Constructor of makePlot ... \n\n" << endl;
}
//Destructor
makePlots::~makePlots()
{
	delete c; // delete canvas first before deleting app
	delete app;
	cout << "\n\n";
	cout << "Destructor of makePlot ... " << endl;
}

///
/// ==================== Init ==================== ///
///
void makePlots::Init( string pedfile, string gainfile, string noisyfile ){
	cout << "----------Init start----------" << endl;
	// read init files and initialize
	yamlReader();
	GainFactorReader( gainfile );
	noisyChannelReader( noisyfile );
	read_P_and_N( pedfile );
	//Crosstalk(injCh);
	P.root_logon();

	// initialize root branch
	Chain1->SetBranchAddress("event",&event);
	Chain1->SetBranchAddress("chip",&chip);
	Chain1->SetBranchAddress("roll",&roll);
	Chain1->SetBranchAddress("dacinj",&dacinj);
	Chain1->SetBranchAddress("timesamp",&timesamp);
	Chain1->SetBranchAddress("hg",&hg);
	Chain1->SetBranchAddress("lg",&lg);
	Chain1->SetBranchAddress("tot_fast",&tot_fast);
	Chain1->SetBranchAddress("tot_slow",&tot_slow);
	Chain1->SetBranchAddress("toa_rise",&toa_rise);
	Chain1->SetBranchAddress("toa_fall",&toa_fall);
	
	// initialize global parameters 
	hg_sig = new double*[NCH];
	lg_sig = new double*[NCH];
	for(int i = 0; i < NCH; i++){
		hg_sig[i] = new double[NSCA];
		lg_sig[i] = new double[NSCA];
	}

	app = new TApplication("app",0,0);
	c = new TCanvas();
	cout << "----------Init complete----------" << endl << endl;
}
///
/// ==================== PlotProducer ==================== ///
///
void makePlots::PlotProducer(){

	char title[200];

	/// Set Output Root File
	int start = input_fileName.find_last_of("/");
	int end   = input_fileName.find(".root");
	string outf = input_fileName.substr(start+1,end-start-1);

	sprintf(title,"root_plot/plot_%s.root",outf.c_str());
	TFile *outfile = new TFile(title,"recreate");
	cout << "output file = " <<  title << endl << endl;
  
	/// Define Parameters 
	int TotalEntries = Chain1->GetEntries();
	int Nevents = TotalEntries/NCHIP;
	cout << "Total Events = " << Nevents << endl;
	int MaxTS = 2; //choose this time sample to be the peak
	int AverageEvents = 0;

	double *XTalkCoupling_Average    = new double[NCHANNEL];
	double *dac_ctrl                 = new double[Nevents];
	double *hg_NoisyChannel          = new double[Nevents];
	double **hg_allCh       = new double*[NCHANNEL];
	double **lg_allCh       = new double*[NCHANNEL];
	double **tot_allCh      = new double*[NCHANNEL];
	double **mip_allCh      = new double*[NCHANNEL];
	double **XTalkCoupling  = new double*[NCHANNEL];
	double **hgFitMean      = new double*[NCHANNEL];
	double **lgFitMean      = new double*[NCHANNEL];
	double **hgFitSigma     = new double*[NCHANNEL];
	double **lgFitSigma     = new double*[NCHANNEL];
	float **hgMean          = new float*[NCHANNEL];
	float **lgMean          = new float*[NCHANNEL];
	float **hgSigma         = new float*[NCHANNEL];
	float **lgSigma         = new float*[NCHANNEL];
	for(int i = 0; i < NCHANNEL; i++){
		hg_allCh[i]      = new double[Nevents];
		lg_allCh[i]      = new double[Nevents];
		tot_allCh[i]     = new double[Nevents];
		mip_allCh[i]     = new double[Nevents];
		XTalkCoupling[i] = new double[Nevents];
		hgFitMean[i]     = new double[NSCA];
		lgFitMean[i]     = new double[NSCA];
		hgFitSigma[i]    = new double[NSCA];
		lgFitSigma[i]    = new double[NSCA];
		hgMean[i]        = new float[NSCA];
		lgMean[i]        = new float[NSCA];
		hgSigma[i]       = new float[NSCA];
		lgSigma[i]       = new float[NSCA];
	}
	double **mip_Ring_1Chip = new double*[NRings];
	for(int i = 0; i < NRings; i++){
		mip_Ring_1Chip[i] = new double[Nevents];
	}
	double **hg_SubPed = new double*[NSCA];
	double **lg_SubPed = new double*[NSCA];
	for(int i = 0; i < NSCA; i++){
		hg_SubPed[i] = new double[NCH];
		lg_SubPed[i] = new double[NCH];
	}
	double **hg_SubPedCM = new double*[NCH];
	double **lg_SubPedCM = new double*[NCH];
	for(int i = 0; i < NCH; i++){
		hg_SubPedCM[i] = new double[NSCA];
		lg_SubPedCM[i] = new double[NSCA];
	}
	double ***mip_Ring_4Chip           = new double**[NRings];
	double ***XTalkCoupling_Ring_4Chip = new double**[NRings];
	for(int i = 0; i < NRings; i++){
		mip_Ring_4Chip[i]           = new double*[NCHIP];
		XTalkCoupling_Ring_4Chip[i] = new double*[NCHIP];
		for(int j = 0; j < NCHIP; j++){
			mip_Ring_4Chip[i][j]           = new double[Nevents];
			XTalkCoupling_Ring_4Chip[i][j] = new double[Nevents];
		}
	}

	/// Declare TDirectories
	sprintf(title,"injCh%d",injCh);
	TDirectory *cdinjCh = outfile->mkdir(title);
	sprintf(title,"allCh");
	TDirectory *cdallCh = cdinjCh->mkdir(title);
	sprintf(title,"Pedestal");
	TDirectory *cdPedestal = cdinjCh->mkdir(title);
	sprintf(title,"hgPedestal");
	TDirectory *cdhgPedestal = cdPedestal->mkdir(title);
	sprintf(title,"lgPedestal");
	TDirectory *cdlgPedestal = cdPedestal->mkdir(title);
	sprintf(title,"hgNoise");
	TDirectory *cdhgNoise = cdPedestal->mkdir(title);
	sprintf(title,"lgNoise");
	TDirectory *cdlgNoise = cdPedestal->mkdir(title);
  
	cdinjCh->cd();
  
	/// Declare Histograms
	TH1D *h_hgPedestal[NSCA][NCHANNEL];
	TH1D *h_lgPedestal[NSCA][NCHANNEL];

	/// Initialize Parameters
	for(int ich = 0; ich < NCHANNEL; ich++) { XTalkCoupling_Average[ich] = 0; }

	for(int sca = 0; sca < NSCA; ++sca) {
		for(int ichannel = 0; ichannel < NCHANNEL; ichannel++) {
			char h_title[50];
			cdPedestal->cd();
			sprintf(h_title,"h_hgPedestal_Ch%d_SCA%d", ichannel, sca);
			h_hgPedestal[sca][ichannel] = new TH1D(h_title,h_title,200,-500,500);
			sprintf(h_title,"h_lgPedestal_Ch%d_SCA%d", ichannel, sca);
			h_lgPedestal[sca][ichannel] = new TH1D(h_title,h_title,200,-500,500);
			hgMean [ichannel][sca] = 0;
			lgMean [ichannel][sca] = 0;
			hgSigma[ichannel][sca] = 0;
			lgSigma[ichannel][sca] = 0;
		}
	}

  
	/// --------------- Start of Loop --------------- ///
	for(int entry = 0; entry < TotalEntries ; ++entry) {
    
		if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }
		Chain1 -> GetEntry(entry);
		dac_ctrl[event] = dacinj;

		// Timesample 
		int TS[NSCA];
		int TS0_sca, MaxTS_sca;
		for(int sca = 0 ; sca < NSCA ; sca++) {
			TS[sca] = timesamp[sca];
			if (timesamp[sca] == 0) { TS0_sca = sca ; }
			if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
		}

		/// Passing hg, lg to self define array
		for (int ich = 0; ich < NCH; ich++){
			for (int sca = 0; sca < NSCA; sca++){
				hg_sig[ich][sca] = hg[sca][ich];
				lg_sig[ich][sca] = lg[sca][ich];
			}
		}
		if(subPed_flag) { Pedestal_CM_Subtractor( chip ); } // Ped & CM Subtraction
		
		for (int ich = 0; ich < NCH; ich++){
			for (int sca = 0; sca < NSCA; sca++){
				//cout << event << " " << ich << " " << sca << " " <<hg_sig[ich][sca] << " " << lg_sig[ich][sca] << endl;
			}
		}

		/// Pedestal histograms
		for (int ich = 0; ich < NCH; ich++){
			for (int sca = 0; sca < NSCA; sca++){
				int ichannel = chip*64 + ich;
				//cout << chip << " " << ich << " " << sca << " " << hg_sig[ich][sca] << endl;
				h_hgPedestal[sca][ichannel]->Fill( hg_sig[ich][sca] );
				h_lgPedestal[sca][ichannel]->Fill( lg_sig[ich][sca] );
				hgMean [ichannel][sca] += hg_sig[ich][sca];
				lgMean [ichannel][sca] += lg_sig[ich][sca];
				hgSigma[ichannel][sca] += ( hg_sig[ich][sca] * hg_sig[ich][sca] );
				lgSigma[ichannel][sca] += ( lg_sig[ich][sca] * lg_sig[ich][sca] );
			}
		}
	
		/// Injection & Cross Talk Analysis
		for(int ich = 0; ich < NCH; ich++){
			int channel      = ich + chip*64;
			double hg = hg_sig[ich][MaxTS_sca];
			double lg = lg_sig[ich][MaxTS_sca];
			double tot = tot_slow[ich];
			
			hg_allCh[channel][event]  = hg;
			lg_allCh[channel][event]  = lg;
			tot_allCh[channel][event] = tot;

			/// mip conversion
			double energy_mip = mipConverter( hg, lg, tot, channel);
			mip_allCh[channel][event] = energy_mip;
		}

		if ( chip == 3 ) {
			for(int ichannel = 0; ichannel < NCHANNEL; ichannel++){
				/// Injection XTalk calculation
				int ichip = ichannel / NCH;
				int inj_channel;
				if ( maskCh_flag )
					inj_channel = ( injChip * NCH ) + injCh;
				else
					inj_channel = ( ichip * NCH ) + injCh;

				XTalkCoupling[ichannel][event] = mip_allCh[ichannel][event] / mip_allCh[inj_channel][event];
				// cout << " event = " << event << " channel = " << ichannel << " energy = " << mip_allCh[ichannel][event] << " Xtalk = " << XTalkCoupling[ichannel][event] << endl;
				if( event>50 && event<=700 ){
					XTalkCoupling_Average[ichannel] += XTalkCoupling[ichannel][event];
					AverageEvents++;
				}
				// Calulate ring Energy
				int iring;
				iring = ringPositionFinder( inj_channel, ichannel );
				if( iring > -1 ) {
					if ( maskCh_flag )
						mip_Ring_1Chip[iring][event] += mip_allCh[ichannel][event];
					else
						mip_Ring_4Chip[iring][ichip][event] += mip_allCh[ichannel][event];
				}
			}

			for(int ichip = 0; ichip < NCHIP; ichip++){
				for(int iring = 1; iring < NRings; iring++) {
					XTalkCoupling_Ring_4Chip[iring][ichip][event] = mip_Ring_4Chip[iring][ichip][event] / mip_Ring_4Chip[0][ichip][event];
				}
			}
		}
	
	}

	/// --------------- End of Loop --------------- ///

	for(int ichannel = 0; ichannel < NCHANNEL; ichannel++){
		XTalkCoupling_Average[ichannel] /= (AverageEvents/NCHANNEL);
		for (int sca = 0; sca < NSCA; sca++){
			if (ichannel%2 == 1) continue;
			h_hgPedestal[sca][ichannel]->Fit("gaus","Q");
			hgFitMean [ichannel][sca] = h_hgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(1);
			hgFitSigma[ichannel][sca] = h_hgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(2);
			//if ( hgFitMean[ichannel][sca] > 200 ) 
			//cout << " ichannel " << ichannel << " sca " << sca << " Mean " << hgFitMean[ichannel][sca] << " Sigma " << hgFitSigma[ichannel][sca] << endl;
			h_lgPedestal[sca][ichannel]->Fit("gaus","Q");
			lgFitMean [ichannel][sca] = h_lgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(1);
			lgFitSigma[ichannel][sca] = h_lgPedestal[sca][ichannel]->GetFunction("gaus")->GetParameter(2);
	  
			hgMean [ichannel][sca] /= Nevents;
			lgMean [ichannel][sca] /= Nevents;
		}
	}
	for(int ichannel = 0; ichannel < NCHANNEL; ichannel++){
		for(int sca = 0; sca < NSCA; sca++){
			hgSigma[ichannel][sca] /= Nevents;
			lgSigma[ichannel][sca] /= Nevents;
			hgSigma[ichannel][sca] -= ( hgMean[ichannel][sca] * hgMean[ichannel][sca] );
			lgSigma[ichannel][sca] -= ( lgMean[ichannel][sca] * lgMean[ichannel][sca] );
			hgSigma[ichannel][sca] = sqrt( hgSigma[ichannel][sca] );
			lgSigma[ichannel][sca] = sqrt( hgSigma[ichannel][sca] );
		}
	}

	/// Plots!!!!!
	cdinjCh->cd();
	for(int ichip = 0; ichip < NCHIP; ichip++){
		int inj_channel = (ichip*64) + injCh;
	
		TGraph* ginjCh_hg  = new TGraph( Nevents, dac_ctrl, hg_allCh[inj_channel] );
		sprintf(title,"hg");
		ginjCh_hg->SetTitle(title);
		ginjCh_hg->SetName(title);
		ginjCh_hg->SetMarkerColor(P.Color(0));
		TGraph* ginjCh_lg  = new TGraph( Nevents, dac_ctrl, lg_allCh[inj_channel] );
		sprintf(title,"lg");
		ginjCh_lg->SetTitle(title);
		ginjCh_lg->SetName(title);
		ginjCh_lg->SetMarkerColor(P.Color(1));
		TGraph* ginjCh_tot = new TGraph( Nevents, dac_ctrl, tot_allCh[inj_channel] );
		sprintf(title,"tot");
		ginjCh_tot->SetTitle(title);
		ginjCh_tot->SetName(title);
		ginjCh_tot->SetMarkerColor(P.Color(2));
		TGraph* ginjCh_mip = new TGraph( Nevents, dac_ctrl, mip_allCh[inj_channel] );
		sprintf(title,"mip_InjCh%d_chip%d", injCh, ichip);
		ginjCh_mip->SetTitle(title);
		ginjCh_mip->SetName(title);
		ginjCh_mip->Write();
	
		TMultiGraph *multig_InjCh_hltot = new TMultiGraph();
		multig_InjCh_hltot->Add(ginjCh_hg);
		multig_InjCh_hltot->Add(ginjCh_lg);
		multig_InjCh_hltot->Add(ginjCh_tot);
		sprintf(title,"hglgtot_InjCh%d_chip%d", injCh, ichip);
		multig_InjCh_hltot->SetTitle(title);
		multig_InjCh_hltot->SetName(title);
		multig_InjCh_hltot->Write();


		/// Xtalk vs dac_ctrl
		TMultiGraph *multig_XTalkCoupling_ring = new TMultiGraph();
		for(int iring = 1; iring < NRings; iring++){
			TGraph* gXTalkCoupling = new TGraph(Nevents, mip_allCh[inj_channel], XTalkCoupling_Ring_4Chip[iring][ichip] );
			sprintf(title,"ring %d", iring);
			gXTalkCoupling->SetTitle(title);
			gXTalkCoupling->SetName(title);
			gXTalkCoupling->SetMarkerColor(P.Color(iring-1));
			gXTalkCoupling->SetLineWidth(0);
			gXTalkCoupling->SetFillColor(0);
			multig_XTalkCoupling_ring->Add(gXTalkCoupling);
		}
		sprintf(title,"XtalkCoupling_InjCh%d_chip%d", injCh, ichip);
		multig_XTalkCoupling_ring->SetTitle(title);
		multig_XTalkCoupling_ring->SetName(title);
		multig_XTalkCoupling_ring->Draw("AP");
		c->Update();
		multig_XTalkCoupling_ring->GetYaxis()->SetRangeUser(-0.01,0.5);
		multig_XTalkCoupling_ring->Write();
	}


  
	/// 2D Average Xtalk 
	int NNoisy = 8;
	int NoisyChannel[8] = {248,186,214,120,126,42,254,190};

	TH2Poly *poly = new TH2Poly;
	InitTH2Poly(*poly);
	poly->SetMinimum(-0.02);
	for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
		float X, Y;
		int forCH = ichannel / 2;
		bool NoisyBool = false;
		X = CHmap[forCH].first;
		Y = CHmap[forCH].second;
		if(ichannel%64 == injCh){
			poly->Fill(X,Y,-2);
		}
		else {
			if(!NoisyBool){
				poly->Fill(X,Y,XTalkCoupling_Average[ichannel]);
				//poly->Fill(X,Y,forCH);
			}
		}
	}
	sprintf(title,"XtalkCoupling_Poly_InjCh%d", injCh);
	poly->SetTitle(title);
	poly->SetName(title);
	poly->Write();

  
	/// 2D Pedestal after CM
	cdPedestal->cd();
	TH2Poly *polyhgPed[NSCA];
	TH2Poly *polylgPed[NSCA];
	TH2Poly *polyhgErr[NSCA];
	TH2Poly *polylgErr[NSCA];
	for (int sca = 0; sca < NSCA; sca++){
		cdhgPedestal->cd();
		polyhgPed[sca] = new TH2Poly();
		InitTH2Poly(*polyhgPed[sca]);
		sprintf(title,"hgPedestal_SCA%d",sca);
		polyhgPed[sca]->SetTitle(title);
		polyhgPed[sca]->SetName(title);
		polyhgPed[sca]->SetMarkerSize(1);

		cdhgNoise->cd();
		polyhgErr[sca] = new TH2Poly();
		InitTH2Poly(*polyhgErr[sca]);
		sprintf(title,"hgSigma_SCA%d",sca);
		polyhgErr[sca]->SetTitle(title);
		polyhgErr[sca]->SetName(title);
		polyhgErr[sca]->SetMarkerSize(1);
	
		cdlgPedestal->cd();
		polylgPed[sca] = new TH2Poly();
		InitTH2Poly(*polylgPed[sca]);
		sprintf(title,"lgPedestal_SCA%d",sca);
		polylgPed[sca]->SetTitle(title);
		polylgPed[sca]->SetName(title);
		polylgPed[sca]->SetMarkerSize(1);

		cdlgNoise->cd();
		polylgErr[sca] = new TH2Poly();
		InitTH2Poly(*polylgErr[sca]);
		sprintf(title,"lgSigma_SCA%d",sca);
		polylgErr[sca]->SetTitle(title);
		polylgErr[sca]->SetName(title);
		polylgErr[sca]->SetMarkerSize(1);
	}

	for(int sca = 0; sca < NSCA; sca++){
		for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
			float X, Y;
			int forCH = ichannel / 2;
			X = CHmap[forCH].first;
			Y = CHmap[forCH].second;
			polyhgPed[sca]->Fill( X, Y, (int)hgFitMean [ichannel][sca]);
			polyhgErr[sca]->Fill( X, Y, (int)hgFitSigma[ichannel][sca]);
			polylgPed[sca]->Fill( X, Y, (int)lgFitMean [ichannel][sca]);
			polylgErr[sca]->Fill( X, Y, (int)lgFitSigma[ichannel][sca]);
			//polyhgPed[sca]->Fill( X, Y, hgMean [ichannel][sca]);
			//polyhgErr[sca]->Fill( X, Y, hgSigma[ichannel][sca]);
			//polylgPed[sca]->Fill( X, Y, lgMean [ichannel][sca]);
			//polylgErr[sca]->Fill( X, Y, lgSigma[ichannel][sca]);
		}
		cdhgPedestal->cd();
		polyhgPed[sca]->Write();
		cdlgPedestal->cd();
		polylgPed[sca]->Write();
		cdhgNoise->cd();
		polyhgErr[sca]->Write();
		cdlgNoise->cd();
		polylgErr[sca]->Write();
	}
  
	/// 1D Average Xtalk
	// Unconnect channels

	double labelCnct[NCHANNEL/2], labelUnCnct[NCHANNEL/2];
	double XTalkCoupling_Cnct[NCHANNEL/2], XTalkCoupling_UnCnct[NCHANNEL/2];
	for( int ichannel = 0; ichannel < NCHANNEL; ichannel++ ) {
		if ( ichannel%2 == 1 ) {
			labelCnct[ichannel/2] = ichannel;
			XTalkCoupling_UnCnct[ichannel/2] = XTalkCoupling_Average[ichannel];
		}
		else {
			labelUnCnct[ichannel/2] = ichannel;
			XTalkCoupling_Cnct[ichannel/2] = XTalkCoupling_Average[ichannel];
		}
	}
	TGraph* gXTalkCoupling_Cnct   = new TGraph(NCHANNEL/2, labelCnct , XTalkCoupling_Cnct);
	sprintf(title,"InjCh%d_XtalkCoupling_Connected_Channel", injCh);
	gXTalkCoupling_Cnct->SetTitle(title);
	gXTalkCoupling_Cnct->SetName(title);
	gXTalkCoupling_Cnct->Write();
	TGraph* gXTalkCoupling_UnCnct = new TGraph(NCHANNEL/2, labelUnCnct , XTalkCoupling_UnCnct);
	sprintf(title,"InjCh%d_XtalkCoupling_UnConnected_Channel", injCh);
	gXTalkCoupling_UnCnct->SetTitle(title);
	gXTalkCoupling_UnCnct->SetName(title);
	gXTalkCoupling_UnCnct->Write();


	cdallCh->cd();
	for(int ichannel = 0; ichannel < NCHANNEL; ichannel++){
		int ichip = ichannel / 64;
		int inj_channel = (ichip*64) + injCh;
	
		TGraph* ginjCh_hg  = new TGraph( Nevents, dac_ctrl, hg_allCh[ichannel] );
		sprintf(title,"hg");
		ginjCh_hg->SetTitle(title);
		ginjCh_hg->SetName(title);
		ginjCh_hg->SetMarkerColor(P.Color(0));
		TGraph* ginjCh_lg  = new TGraph( Nevents, dac_ctrl, lg_allCh[ichannel] );
		sprintf(title,"lg");
		ginjCh_lg->SetTitle(title);
		ginjCh_lg->SetName(title);
		ginjCh_lg->SetMarkerColor(P.Color(1));
		TGraph* ginjCh_tot = new TGraph( Nevents, dac_ctrl, tot_allCh[ichannel] );
		sprintf(title,"tot");
		ginjCh_tot->SetTitle(title);
		ginjCh_tot->SetName(title);
		ginjCh_tot->SetMarkerColor(P.Color(2));
		TGraph* ginjCh_mip = new TGraph( Nevents, dac_ctrl, mip_allCh[ichannel] );
		sprintf(title,"mip_Ch%d", ichannel);
		ginjCh_mip->SetTitle(title);
		ginjCh_mip->SetName(title);
		ginjCh_mip->Write();
	
		TMultiGraph *multig_InjCh_hltot = new TMultiGraph();
		multig_InjCh_hltot->Add(ginjCh_hg);
		multig_InjCh_hltot->Add(ginjCh_lg);
		multig_InjCh_hltot->Add(ginjCh_tot);
		sprintf(title,"hglgtot_Ch%d", ichannel);
		multig_InjCh_hltot->SetTitle(title);
		multig_InjCh_hltot->SetName(title);
		multig_InjCh_hltot->Write();

		TGraph* gXTalkCoupling = new TGraph(Nevents, mip_allCh[inj_channel], XTalkCoupling[ichannel] );
		sprintf(title,"xtalk_Ch%d", ichannel);
		gXTalkCoupling->SetTitle(title);
		gXTalkCoupling->SetName(title);
		gXTalkCoupling->Write();
	}


	outfile->Write();
	outfile->Close();


	/// deallocate 
	for (int i = 0; i < NRings; i++){
		for (int j = 0; j < NCHIP; j++) {
			delete[] mip_Ring_4Chip[i][j];
			delete[] XTalkCoupling_Ring_4Chip[i][j];
		}
		delete[] mip_Ring_4Chip[i];
		delete[] XTalkCoupling_Ring_4Chip[i];
	}
	delete[] mip_Ring_4Chip;
	delete[] XTalkCoupling_Ring_4Chip;
	
	for (int i = 0; i < NCHANNEL; i++){
		delete[] hg_allCh[i];       
		delete[] lg_allCh[i];       
		delete[] tot_allCh[i];      
		delete[] mip_allCh[i];      
		delete[] XTalkCoupling[i];  
	}
	delete[] hg_allCh;       
	delete[] lg_allCh;       
	delete[] tot_allCh;      
	delete[] mip_allCh;      
	delete[] XTalkCoupling;  
	delete[] XTalkCoupling_Average;   
	delete[] dac_ctrl;                
	delete[] hg_NoisyChannel;         

}
///
/// ==================== cosmicAnalyzer ==================== ///
///
void makePlots::cosmicAnalyzer(){
	
	char title[200];
	/// Set Output Root File
	int start = input_fileName.find_last_of("/");
	int end   = input_fileName.find(".root");
	string outf = input_fileName.substr(start+1,end-start-1);

	sprintf(title,"cosmicAnalysis/plot_%s.root",outf.c_str());
	TFile *outfile = new TFile(title,"recreate");
	cout << "output file = " << title << endl;
  
	/// Declare Parameters
	int TotalEntries = Chain1->GetEntries();
	int Nevents = TotalEntries/NCHIP;
	cout << "Total Events = " << Nevents << endl;
	int MaxTS = 2;              //choose this time sample to be the peak
	int mipCount = 0;

	/// Define Histograms
	TH1D *h_mipAllCh = new TH1D("h_mipAllCh","",50,0,400);
  
  
	//==================== Loop over the events ====================
   
	for(int entry = 0; entry < TotalEntries ; ++entry){
		if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }
		Chain1 -> GetEntry(entry);

		/// Timesample
		int TS[NSCA];
		int TS0_sca, MaxTS_sca;
		for(int sca = 0 ; sca < NSCA ; sca++) {
			TS[sca] = timesamp[sca];
			if (timesamp[sca] == 0) { TS0_sca = sca ; }
			if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
		}

		/// hg, lg signal passed to self defined parameters 
		for (int ich = 0; ich < NCH; ich++){
			for (int sca = 0; sca < NSCA; sca++){
				hg_sig[ich][sca] = hg[sca][ich];
				lg_sig[ich][sca] = lg[sca][ich];
			}
		}
		if(subPed_flag) { Pedestal_CM_Subtractor( chip ); }  // Ped & CM Subtraction


		/// Calculate hit per chip
		int hit = 0;
		for(int ich = 0; ich < NCH; ich++){
			if ( mipSigCheck(hg_sig[ich], TS ) ) { hit++; }
		}

		/// Fill histogram
		for (int ich = 0; ich < NCH; ich+=2) {
			if ( ich + chip*NCH == 44 ) continue;
			//if ( ich != 18 ) continue;
			if ( mipSigCheck(hg_sig[ich], TS ) && hit == 1 ) {
				h_mipAllCh->Fill( hg_sig[ich][MaxTS_sca] );
				pulsePlotter( hg_sig[ich], TS , event, chip, ich, -1, -1);
				mipCount++;
			}
		}
	}

	//... ==================== End of Loop ==================== ...

	cout << endl << "totalEvent# = " << Nevents << " signal# = " << mipCount << endl;
	cout << "efficiency = " << (float)mipCount / Nevents << endl;
 
	outfile->Write();
	outfile->Close();

}


///
/// ==================== Pulse_display ==================== ///
///
void makePlots::Pulse_display( int displayChannel, int pulseDisplay_type, int lowerR, int upperR, int startEv ){

	int Nevents = Chain1->GetEntries();
	cout << "Total Events = " << Nevents << endl;

	/// --------------- Loop Over Events --------------- ///
	for(int ev = 0; ev < Nevents ; ++ev){
		Chain1 -> GetEntry(ev);
		if( event < startEv ) continue;
		int TS[NSCA];
		for(int i = 0 ; i < NSCA ; ++i)  TS[i] = timesamp[i];

		/// Passing hg lg signal to self defined array
		for (int ich = 0; ich < NCH; ich++){
			for (int sca = 0; sca < NSCA; sca++){
				hg_sig[ich][sca] = hg[sca][ich];
				lg_sig[ich][sca] = lg[sca][ich];
			}
		}
		if(subPed_flag) { Pedestal_CM_Subtractor( chip ); } // Ped & CM Subtraction

		/// --------------- Pulse Plots --------------- ///
		
		/// Loop Over every Channel
		if ( pulseDisplay_type == 0 && displayChannel == -1) { 
			for( int ich =0; ich < 64; ich+=2){
				pulsePlotter( hg_sig[ich], TS, event, chip, ich, lowerR, upperR);
			}
		}
		/// injCh display
		else if ( pulseDisplay_type == 1 ) { 
			pulsePlotter( lg_sig[injCh], TS, event, chip, injCh, lowerR, upperR);
		}
		/// find signal and display 
		else if ( pulseDisplay_type == 2 ) {  
			for( int ich =0; ich < 64; ich+=2){
				if ( mipSigCheck( hg_sig[ich], TS ) ) {
					pulsePlotter( hg_sig[ich], TS, event, chip, ich, lowerR, upperR );
				}
			}
		}
		/// selected channel display
		else {  
			int ichip = displayChannel / 64;
			int ich   = displayChannel % 64;
			if ( chip != ichip ) continue;
			pulsePlotter( hg_sig[ich], TS, event, ichip, ich, lowerR, upperR);
		}
	}

}
///
/// ==================== mipConverter ==================== ///
///
double makePlots::mipConverter( double hg_SubPed, double lg_SubPed, double tot , int channel){
  
	double mip;
	int ichip = channel / NCH;
	int ich   = channel % NCH;
  
	if( hg_SubPed < HGTP[ichip][ich]){
		mip = hg_SubPed * ADC2MIP;
	}
	else{
		//if( lg_sig < LGTP[channel]){
		if( lg_SubPed < LGTP_default){
			mip = ( lg_SubPed * LG2HG_Conversion[ichip][ich] * ADC2MIP);
		}
		else{
			mip = ( (tot - TOTOffSet[ichip][ich]) * TOT2LG_Conversion[ichip][ich] * LG2HG_Conversion[ichip][ich] * ADC2MIP);
		}
	}
	return mip;
}

///
/// ==================== ringPositionFinder ==================== ///
///
int makePlots::ringPositionFinder( int inj_channel, int ichannel){

	int formatCh = ichannel / 2;
	int formatInjCh = inj_channel / 2;
	double X     = CHmap[formatCh].first;
	double Y     = CHmap[formatCh].second;
	double inj_X = CHmap[formatInjCh].first;
	double inj_Y = CHmap[formatInjCh].second;

	int ring = -1;
  
	if ( ichannel % 2 == 0 ){ 
		double dx = X - inj_X;
		double dy = Y - inj_Y;
		double dR = sqrt(dx*dx + dy*dy);
	
		if ( formatCh == formatInjCh ){ ring = 0; }
		else {
			ring = 0;
			while(true) {
				ring++;
				if( ring == 5 ) {
					ring = -1;
					break; 
				}
				if( dR < 1.12455*ring*1.2 ) break;
			}
		}
	}

	return ring;
}

///
/// ==================== yamlReader ==================== ///
///
void makePlots::yamlReader(){

	int start = input_fileName.find_last_of("/");
	int end = input_fileName.find("_pedestal.root");
	string f_substr = input_fileName.substr(start+1,end-start-1);
	string rootFileName(f_substr);
	end = input_fileName.find("ana_output");
	f_substr = input_fileName.substr(0,end-1);
	string yamlPath(f_substr);
	char yamlFileName[100];
	sprintf(yamlFileName,"%s/yaml/%s.yaml",yamlPath.c_str(), rootFileName.c_str());
  
	string searchstr;
	string line;
	ifstream yamlFile(yamlFileName);
	if(!yamlFile.is_open()){
		cout << "Did not find injection file " << yamlFileName
			 << ".\n Take this run as pedestal.(Inj_dac = 0)" << endl;
	}
	if(yamlFile.is_open()){
		cout << "yamlFile = " << yamlFileName << endl;
		while( true ) {
			if ( yamlFile.eof() ) break;
			getline (yamlFile, line);
	  
			if ( line.find("channelIds:") != -1 ){
				string tmp;
				//yamlFile >> tmp >> searchstr;
				start = line.find("[");
				end = line.find("]");
				searchstr = line.substr(start+1,end-start+1);
				injCh = atoi(searchstr.c_str());
				cout << "InjCh = " << injCh << endl;
			}
			if ( line.find("acquisitionType") != -1 ){
				start = line.find(":");
				searchstr = line.erase(0, start+2);
				acquisitionType = searchstr;
				cout << "acquisitionType = " << acquisitionType << endl;
			}
			else if ( maskCh_flag == true && line.find("channelIdsToMask:") != -1 ) {
				getline(yamlFile, line);
				int count = 3; 
				while ( line.find("[]") == -1) {
					getline(yamlFile, line);
					count--;
					if ( count < 0 ) break;
				}
				if ( count > -1 ) {
					injChip = count;
					cout << "1 channel injection = " << injChip << endl;
				}
			}
		}
	}

	cout << endl;
}

///
/// ==================== GainFactorReader ==================== ///
///
void makePlots::GainFactorReader( string gainfile ){

	string TB_GainFactors("src_txtfile/TPro_fittingoutput.txt");
  
	ifstream GainFile(gainfile);
	string line;
	char tmp[50];
	int ichip, ich;
	if(!GainFile.is_open()){
		cout << "Did not find GainFactor file " << gainfile
			 << ".\n Take this run's GainFactor as default & calculate the Gainfactors " << endl;
		for(ichip = 0; ichip < NCHIP; ichip++){
			for(ich = 0; ich < NCH; ich++){
				HGTP[ichip][ich] = 1500;
				LG2HG_Conversion[ichip][ich] = 8.5;
				LGTP[ichip][ich] = 900;
				TOT2LG_Conversion[ichip][ich] = 3.8;
				TOTOffSet[ichip][ich] = 180;
			}
		}
		//Gain_factor_producer();
	}
  
	if(GainFile.is_open()){
		cout << "gainFile = " << gainfile << endl;
		getline(GainFile,line);
		while(!GainFile.eof()){
			GainFile >> tmp >> tmp >> ichip >> ich >> tmp;
			GainFile >> HGTP[ichip][ich] >> LG2HG_Conversion[ichip][ich] >> LGTP[ichip][ich] >> TOT2LG_Conversion[ichip][ich] >> TOTOffSet[ichip][ich];
			getline(GainFile,line);
		}
	}
	cout << endl;
}
///
/// ==================== noisyChannelReader ==================== ///
///
void makePlots::noisyChannelReader( string noisyFileName ) {
	ifstream noisyFile( noisyFileName );
	string line;
	int ichip, ich;
	if(!noisyFile.is_open()){
		cout << "Did not find noisyChannelFile " << noisyFileName
			 << ".\n not setting noisyChannels " << endl;
	}
	else{
		cout << "noisyFile = " << noisyFileName << endl;
		while(!noisyFile.eof()){
			noisyFile >> ichip >> ich;
			int ichannel = ichip*NCH + ich;
			noisyChannel.push_back(ichannel);
		}
		noisyChannel.pop_back();
	}
	cout << "noisy channels = ";
	for(int i = 0; i < noisyChannel.size(); i++){
		cout << noisyChannel.at(i) << " " ;
	}
	cout << endl << endl;
}


///
/// ==================== Pedestal_CM_Subtractor ==================== ///
///
void makePlots::Pedestal_CM_Subtractor( int chip ){
	
  	for (int ich = 0; ich < NCH; ich++){
		for (int sca = 0; sca < NSCA; sca++){
			hg_sig[ich][sca] -= avg_HG[chip][ich][sca];  // Pedestal Subtraction
			lg_sig[ich][sca] -= avg_LG[chip][ich][sca];
		}
	}

	double *hgCM_sca, *lgCM_sca;
	hgCM_sca = CMCalculator_v2( hg_sig, chip ); // Calculate CM for each sca
	lgCM_sca = CMCalculator_v2( lg_sig, chip );
	

	for (int ich = 0; ich < NCH; ich++){
		for (int sca = 0; sca < NSCA; sca++){
			hg_sig[ich][sca] -= hgCM_sca[sca]; // CM subtraction
			lg_sig[ich][sca] -= lgCM_sca[sca];
			//cout << hg_sig[ich][sca] << endl;
		}
	}

}

///
/// ==================== CMCalculator_v2 ==================== ///
///
double* makePlots::CMCalculator_v2 ( double **sig_subPed, int chip ) {
	// Calculate CM for each TS
	static double meanChipPedestal[NSCA];

	int scaCount[NSCA];
	for (int sca = 0; sca < NSCA; sca++) {
		meanChipPedestal[sca] = 0;
		scaCount[sca] = 0;
	}

	for (int ich = 0; ich < NCH; ich+=2) {
		int ichannel = ich + chip*NCH;
		bool noisy_flag = false;
		
		// Determine if the channel is a noisy channel
		for(int i = 0; i < noisyChannel.size(); i++){
			if ( ichannel == noisyChannel.at(i) ) { noisy_flag = true; }
		}
		if (noisy_flag) continue;

		// Calculate mean pedestal for each sca 
		for (int sca = 0; sca < NSCA; sca++) {
			meanChipPedestal[sca] += sig_subPed[ich][sca];
			scaCount[sca]++;
		}
	}

	for (int sca = 0; sca < NSCA; sca++) {
		meanChipPedestal[sca] /= scaCount[sca];
	}
	return meanChipPedestal;
}

///
/// ==================== CMCalculator ==================== ///
///
double makePlots::CMCalculator( double **sig_subPed, int *TS ) {
	// Common mode calculation, per event, per chip
	// Assume common mode is identical for all time samples
	// Pick 5-9 TS for CM
	int scaCount = 0;
	double meanChipPedestal = 0;
  
	for(int ich = 0; ich < NCH; ich+=2){
		for( int timesample = 0; timesample <= 9; timesample++){
			int sca = 0;
			while(true){
				if ( TS[sca] == timesample ) break;
				else sca++;
			}
			meanChipPedestal += sig_subPed[sca][ich];
			scaCount++;
		}
	}
	meanChipPedestal /= scaCount;
	return meanChipPedestal;
}

///
/// ==================== mipSegCheck ==================== ///
///
bool makePlots::mipSigCheck( double *sig, int *TS ) {
	int p_noisy_cut  = 2000;
	int n_noisy_cut  = -1000;
	double noSignal_cut = 30;
	bool sig_flag = false;
	double sig_ts[NSCA];
  
	for( int sca = 0; sca < NSCA; sca++) {
		sig_ts[ TS[sca] ] = sig[sca];
		if ( TS[sca] > 1 && TS[sca] < 5 && sig[sca] < p_noisy_cut && sig[sca] > noSignal_cut ) sig_flag = true;
	}
	for( int sca = 0; sca < NSCA; sca++) {
		if ( sig[sca] < n_noisy_cut ) sig_flag = false;
	}
	if ( sig_ts[0] > sig_ts[2] || sig_ts[0] > sig_ts[3] ) sig_flag = false;
	if ( sig_ts[0] < -200 || sig_ts[0] > 150 ) sig_flag = false;
	if ( sig_ts[11] < -200 || sig_ts[11] > 150 ) sig_flag = false;
	if ( sig_ts[2] < 15 ) sig_flag = false;
  
	return sig_flag;
}

///
/// ==================== pulsePlotter ==================== ///
///
void makePlots::pulsePlotter( double *sig, int *TS, int ev, int ichip, int ich, int lowerR, int upperR ) {
	
	// This function plots the input 13 timesamples and show it on the screen 
	double double_TS[NSCA];
	for(int i = 0; i < NSCA; i++) double_TS[i] = (double)TS[i];
	TGraph *gr = new TGraph(NSCA, double_TS, sig );
	char plot_title[50];
	gr->SetMarkerColor(ichip+1);
	gr->SetMarkerStyle(22);
	gr->SetMarkerSize(1.2);
	if ( upperR != -1 ) gr->GetYaxis()->SetRangeUser(lowerR,upperR);
	gr->Draw("AP");
	sprintf(plot_title,"evt %d chip %d channel%d", ev, ichip, ich);
	gr->SetTitle(plot_title);
	gr->GetXaxis()->SetTitle("TS");
	gr->GetYaxis()->SetTitle("ADC");
	c->Update();
	gPad->WaitPrimitive();
  
}

///
/// ==================== read_P_and_N ==================== ///
///
void makePlots::read_P_and_N(string ped_file){

	int end = input_fileName.find("ana_output");
	string pedPath = input_fileName.substr(0,end-1);
	char pedFileName[100];
	sprintf(pedFileName,"%s/pedestal",pedPath.c_str());

	char HG_name[100],LG_name[100];
	sprintf(HG_name,"%sHG.txt",pedFileName);
	sprintf(LG_name,"%sLG.txt",pedFileName);
	ifstream inHG(HG_name);
	ifstream inLG(LG_name);
	if( !inHG.is_open() || !inLG.is_open()){
		cout << "File not found! Neither " << HG_name << " or " << LG_name
			 << " exist!" << endl;
		return;}
	else{
		cout << "pedFile = " << endl;
		cout << "1. " << HG_name << "\n" << "2. "<< LG_name << endl;
		string line;
		int ichip,ich;
		int testeof;

	
		while(true){
			inHG >> testeof;
			if( inHG.eof() ) break;
			else{
				inHG >> ichip; 
				inHG >> ich;
				for(int sca = 0; sca < NSCA ; ++sca){
					inHG >> avg_HG[ichip][ich][sca];
					inHG >> sigma_HG[ichip][ich][sca];
				}
			}
		}

		while(true){
			inLG >> testeof;
			if( inLG.eof() ) break;
			else{
				inLG >> ichip;
				inLG >> ich;
				for(int sca = 0; sca < NSCA ; ++sca){
					inLG >> avg_LG[ichip][ich][sca];
					inLG >> sigma_LG[ichip][ich][sca];
				}
			}
		}

    
		cout << "Reading pedestal file done!" << endl << endl;
		inHG.close();
		inLG.close();
	}  
}

///
/// ==================== readmap ==================== ///
///
void makePlots::readmap(){
	ifstream file("./src_txtfile/CH_map.txt");
	string line;
	int ichip,ich,itype,iformatCH;
	double iposx, iposy;
	while(true){
		getline(file,line);
		if( file.eof() ) break;
		file >> ichip >> ich >> iposx >> iposy >> itype;
		iformatCH = ichip*32 + ich/2;
		CHmap[iformatCH] = make_pair(iposx,iposy);}
	file.close();
	//Since there is no such pad, assign a unreasonable value
	CHmap[2*32+60/2] = make_pair(1000.,1000.);

}

///
/// ==================== InitTH2Poly ==================== ///
///
void makePlots::InitTH2Poly(TH2Poly& poly)
{
	int MAXVERTICES = 6;
	double HexX[MAXVERTICES];
	double HexY[MAXVERTICES];
	int iu,iv,CellXYsize;
	ifstream file("src_txtfile/poly_frame.txt");
	string line;
  
	for(int header = 0; header < 4; ++header )     getline(file,line);
  
	while(true){
		getline(file,line);
		if( file.eof() ) break;
		file >> iu >> iv >> CellXYsize;    
		for(int i = 0; i < CellXYsize ; ++i){
			getline(file,line);
			file >> HexX[i] >> HexY[i];
		}
		poly.AddBin(CellXYsize, HexX, HexY);
	}
	file.close();
}



///
/// ==================== Gain_factor_producer ==================== ///
///
void makePlots::Gain_factor_producer(){
  
	//-------------------- Define Parameters --------------------

	int TotalEntries = Chain1->GetEntries();
	int Nevents = TotalEntries/NCHIP;
	int MaxTS = 2;

	int hg_injCh[NCHIP][Nevents], lg_injCh[NCHIP][Nevents], tot_injCh[NCHIP][Nevents];
	int dac_ctrl[Nevents];
	int HGLGfitmax[NCHIP], TOTLGfitmax[NCHIP];
	bool HGTP_flag[NCHIP], LGTP_flag[NCHIP];

	//==================== Set Output File ====================

  


	//==================== Initialize ====================

	for(int ichip = 0; ichip < NCHIP; ichip++){
		HGTP_flag[ichip] = false;
		LGTP_flag[ichip] = false;
	}


	//==================== Loop Over Events ====================

	for(int ev = 0; ev < TotalEntries; ev++){
		if(ev%1000 == 0){ cout << "Now Processing = " << ev << endl;}
		Chain1->GetEntry(ev);
		dac_ctrl[event] = dacinj;
    
		int TS0_sca, MaxTS_sca;
		for(int sca = 0 ; sca < NSCA ; sca++) {
			if (timesamp[sca] == 0) { TS0_sca = sca ; }
			if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
		}
    
    
		hg_injCh[chip][event] = ( hg[MaxTS_sca][injCh] - hg[TS0_sca][injCh] );
		lg_injCh[chip][event] = ( lg[MaxTS_sca][injCh] - lg[TS0_sca][injCh] );
		tot_injCh[chip][event] = tot_slow[injCh];


		if(hg_injCh[chip][event] > HGTP[chip][injCh] && HGTP_flag[chip] == false){
			HGTP_flag[chip] = true;
			HGLGfitmax[chip] = lg_injCh[chip][event];
		}
		if(lg_injCh[chip][event] > LGTP[chip][injCh] && LGTP_flag[chip] == false){
			LGTP_flag[chip] = true;
			TOTLGfitmax[chip] = tot_injCh[chip][event];
		}
	}

  
	//==================== End Loop ====================

	//...

	//==================== Plots ====================

	char pltTit[100];
	string Xtit, Ytit, Opt;
	int MkSty, MkClr, LClr, fitmin, fitmax;
	float MkSize, LWid;
	bool Stat, Wait, SavePlot;
	TCanvas *c1 = new TCanvas();

	TGraph** gh = new TGraph*[NCHIP];
	TGraph** gl = new TGraph*[NCHIP];
	TGraph** gtot = new TGraph*[NCHIP];
	TGraph** LG2HG = new TGraph*[NCHIP];
	TGraph** TOT2LG = new TGraph*[NCHIP];

	for(int ichip = 0; ichip < NCHIP; ichip++){
		gh[ichip] = new TGraph(Nevents,dac_ctrl,hg_injCh[ichip]);
		gl[ichip] = new TGraph(Nevents,dac_ctrl,lg_injCh[ichip]);
		gtot[ichip] = new TGraph(Nevents,dac_ctrl,tot_injCh[ichip]);
		LG2HG[ichip] = new TGraph(Nevents,lg_injCh[ichip],hg_injCh[ichip]);
		TOT2LG[ichip] = new TGraph(Nevents,tot_injCh[ichip],lg_injCh[ichip]);
    
		LG2HG[ichip]->Fit("pol1","","",fitmin = 0,fitmax = HGLGfitmax[ichip]);
		//TOT2LG[ichip]->Fit("pol1","","",fitmin = TOTOffSet,fitmax = TOTLGfitmax[ichip]);
		TOT2LG[ichip]->Fit("pol1","","",fitmin = 200,fitmax = 300);
    
		TF1* Linear_fit_LG2HG = LG2HG[ichip]->GetFunction("pol1");
		TF1* Linear_fit_TOT2LG = TOT2LG[ichip]->GetFunction("pol1");
		LG2HG_Conversion[ichip][injCh] = Linear_fit_LG2HG->GetParameter(1);
		TOT2LG_Conversion[ichip][injCh] = Linear_fit_TOT2LG->GetParameter(1);
		TOTOffSet[ichip][injCh] = -Linear_fit_TOT2LG->GetParameter(0)/Linear_fit_TOT2LG->GetParameter(1);
		cout << LG2HG_Conversion[ichip][injCh] << " " <<  TOT2LG_Conversion[ichip][injCh] << " " << TOTOffSet[ichip][injCh] << endl;
		/*
		  TOT2LG[ichip]->Draw("AP");
		  c1->Update();
		  gPad->WaitPrimitive();
		*/
		sprintf(pltTit,"HG_Chip%d",ichip);
		P.GStd(*gh[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);

		sprintf(pltTit,"LG_Chip%d",ichip);
		P.GStd(*gl[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);

		sprintf(pltTit,"TOT_Chip%d",ichip);
		P.GStd(*gtot[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);
    
		sprintf(pltTit,"LG2HG_Chip%d",ichip);
		P.GStd(*LG2HG[ichip], pltTit, Xtit = "LG", Ytit = "HG", Opt = "AP", Wait = 0, SavePlot = 1);
    
		sprintf(pltTit,"TOT2LG_Chip%d",ichip);
		P.GStd(*TOT2LG[ichip], pltTit, Xtit = "TOT", Ytit = "LG", Opt = "AP", Wait = 0, SavePlot = 1);
	}
}

///
/// ==================== Crosstalk ==================== ///
///
void makePlots::Crosstalk(Int_t CH){

	TCanvas* c1 = new TCanvas();
	for(int ichip = 0; ichip < NCHIP; ichip++){
    
		int cross_num = 6;
		int formatInj_Ch = CH/2+ichip*32;
		float Xdist = 0.974452;
		float Ydist = 0.5626;
		float cross_posx[cross_num], cross_posy[cross_num];
		float X = CHmap[formatInj_Ch].first;
		float Y = CHmap[formatInj_Ch].second;
  
		cross_posx[0] = X - Xdist;
		cross_posy[0] = Y + Ydist;
		cross_posx[1] = X;
		cross_posy[1] = Y + 2*Ydist;
		cross_posx[2] = X + Xdist;
		cross_posy[2] = Y + Ydist;
		cross_posx[3] = X + Xdist;
		cross_posy[3] = Y - Ydist;
		cross_posx[4] = X;
		cross_posy[4] = Y - 2*Ydist;
		cross_posx[5] = X - Xdist;
		cross_posy[5] = Y - Ydist;
		
		cout << "Chip" << ichip << " FirstRing Channels = " ;
		for(int icross = 0; icross < 6; icross++){
			//cout << cross_posx[i] << " " << cross_posy[i] << endl;
			int ch = 0;
			bool good_channel = true;
			while(abs(CHmap[ch].first-cross_posx[icross]) > 1e-4 || abs(CHmap[ch].second-cross_posy[icross]) > 1e-4){
				if( ch > 127) {
					good_channel = false;
					break;
				}
				ch++;	    
				//cout << ch << " "<<abs(CHmap[ch].first - cross_posx[i])<< " " << abs(CHmap[ch].second - cross_posy[i])<< endl;
				/*
				  forCH = chip*32 + cross_ch[i]/2;
				  X = CHmap[forCH].first;
				  Y = CHmap[forCH].second;	      
				  poly->Fill(X,Y,cross_ch[i]);
				*/
      
			}
			cross_ch_FirstRing[ichip][icross] = ch * 2;
			cout << cross_ch_FirstRing[ichip][icross] << " " ;
		}
		cout << endl;
	}
	cout << endl;
	TH2Poly *poly = new TH2Poly;
	InitTH2Poly(*poly);
	/*  
		for(int i=0; i<6; i++){
		poly->Draw("colztext0");
		c1->Update();
		//gPad->WaitPrimitive();
		}
	*/
	delete c1;
  
}

///
/// ==================== Cut ==================== ///
///
int makePlots::Cut(Long64_t entry, Long64_t sigma)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	if(entry<10*sigma){
		return -1;
	}
	else return 1;
}



/*
  double pulseShape_fcn_v2(double t, double tmax, double amp, double amp0 = 0., double tau = 22., int n_ord = 3){

  const double ampl_norm = 1.608;
  const double alpha = 0.931;

  if( t>tmax-_trise )
  return (amp*ampl_norm * (1 - ((t-(tmax-_trise))/tau)/(n_ord+1)) * std::pow((t-(tmax-_trise))/tau, n_ord) * std::exp(-alpha*(t-(tmax-_trise))/tau)) + amp0;
  else return 0;

  }
*/


/*
///
///==================== LEDAnalyzer ====================///
///
void makePlots::LEDAnalyzer(){
char title[200];

/// Set Output Root File
int start = input_fileName.find_last_of("/");
int end   = input_fileName.find(".root");
string outf = input_fileName.substr(start+1,end-start-1);

sprintf(title,"cosmicAnalysis/plot_%s.root",outf.c_str());
TFile *outfile = new TFile(title,"recreate");
cout << "output file = " << title << endl;
  
/// Declare Parameters
int TotalEntries = Chain1->GetEntries();
int Nevents = TotalEntries/NCHIP;
cout << "Total Events = " << Nevents << endl;
int MaxTS = 2;              //choose this time sample to be the peak
int mipCount = 0;

double *dac_ctrl   = new double[Nevents];
double **hg_allCh       = new double*[NCHANNEL];
double **lg_allCh       = new double*[NCHANNEL];
double **tot_allCh      = new double*[NCHANNEL];
double **mip_allCh      = new double*[NCHANNEL];
for(int i = 0; i < NCHANNEL; i++){
hg_allCh[i]      = new double[Nevents];
lg_allCh[i]      = new double[Nevents];
tot_allCh[i]     = new double[Nevents];
mip_allCh[i]     = new double[Nevents];
}
double **hg_SubPed = new double*[NSCA];
double **lg_SubPed = new double*[NSCA];
for(int i = 0; i < NSCA; i++){
hg_SubPed[i] = new double[NCH];
lg_SubPed[i] = new double[NCH];
}
double **hg_sig = new double*[NCH];
double **lg_sig = new double*[NCH];
for(int i = 0; i < NCH; i++){
hg_sig[i] = new double[NSCA];
lg_sig[i] = new double[NSCA];
}
  
/// Declare directories
  
/// Define Histograms
TH1D *h_mipAllCh = new TH1D("h_mipAllCh","",50,0,400);
  
/// Initialize
  
/// --------------- Start of Loop -------------- ///
   
for(int entry = 0; entry < TotalEntries ; ++entry){
    
if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }
Chain1 -> GetEntry(entry);
dac_ctrl[event] = dacinj;

int TS[NSCA];
int TS0_sca, MaxTS_sca;
for(int sca = 0 ; sca < NSCA ; sca++) {
TS[sca] = timesamp[sca];
if (timesamp[sca] == 0) { TS0_sca = sca ; }
if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
}

for(int ich = 0; ich < NCH; ich++ ){
for(int sca = 0; sca < NSCA; sca++ ){
hg_SubPed[sca][ich] = hg[sca][ich] - avg_HG[chip][ich][sca]; // Pedestal Subtraction
lg_SubPed[sca][ich] = lg[sca][ich] - avg_LG[chip][ich][sca];
}
}
double hgCM = CMCalculator( hg_SubPed, TS ); // Calculate CM for the chip
double lgCM = CMCalculator( lg_SubPed, TS );
		
int hit = 0;
for(int ich = 0; ich < NCH; ich++){
for (int sca = 0; sca < NSCA; sca++){
if ( subPed_flag ){
hg_sig[ich][sca] = hg_SubPed[sca][ich] - hgCM; // CM subtraction 
lg_sig[ich][sca] = lg_SubPed[sca][ich] - lgCM;
}
else {
hg_sig[ich][sca] = hg[sca][ich];
lg_sig[ich][sca] = lg[sca][ich];
}
}
if ( mipSigCheck(hg_sig[ich], TS ) ) hit++;
}

for (int ich = 0; ich < NCH; ich+=2) {
if ( ich + chip*NCH == 44 ) continue;
if ( mipSigCheck(hg_sig[ich], TS ) && hit < 2) {
h_mipAllCh->Fill( hg_sig[ich][MaxTS_sca] );
pulsePlotter( hg_sig[ich], TS , event, chip, ich, -1, -1);
mipCount++;
}
}

// mip conversion
//double energy_mip = mipConverter( hg_sig, lg_sig, tot, channel);
//mip_allCh[channel][event] = energy_mip;
}

/// --------------- End of Loop --------------- ///

cout << endl << "totalEvent# = " << Nevents << " signal# = " << mipCount << endl;
cout << "efficiency = " << (float)mipCount / Nevents << endl;
 
// Plots!!!!!

outfile->Write();
outfile->Close();


// deallocate 
for (int i = 0; i < NCHANNEL; i++){
delete[] hg_allCh[i];       
delete[] lg_allCh[i];       
delete[] tot_allCh[i];      
delete[] mip_allCh[i];      
}
delete[] hg_allCh;       
delete[] lg_allCh;       
delete[] tot_allCh;      
delete[] mip_allCh;      
delete[] dac_ctrl;                

  
}
*/





