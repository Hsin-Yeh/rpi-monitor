// Arthor: Hsin-Yeh Wu
// Email : thankyouyou06@gmail.com
//
// This class is the main class for analyzing the rpi data 


#ifndef makePlots_h
#define makePlots_h

#include "TChain.h"
#include "TH2Poly.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "PlotSetting.h"
#include <string>
#include <utility> //std::pair
#include <map>     //std::map

using namespace std;

const int NCHIP = 4;
const int NCH = 64;
const int NSCA = 13;
const int NformatCH = 128;
const int NCHANNEL = 256;
const int NRings = 5;

class makePlots{
 public:

	makePlots (TChain* inchain);
	~makePlots();

	//public function
	void Init( string pedfile, string gainfile, string noisyfile );
	void PlotProducer();
	void cosmicAnalyzer();
	void Pulse_display( int displayChannel = -1 , int acq_type = 0, int lowerR = -1, int upperR = -1 , int startEv = 0 );
  
	//public parameter
	string input_fileName;
	bool subPed_flag;
	bool maskCh_flag;

  
	//PlotSetting class 
	PlotSetting P;

  
 private:

	// private function
	double  mipConverter( double hg_SubPed, double lg_SubPed, double tot , int channel);
	int     ringPositionFinder( int inj_channel, int channel);
	double  CMCalculator( double **sig_SubPed, int *TS );
	double* CMCalculator_v2(double **sig_SubPed, int chip );
	void    Pedestal_CM_Subtractor( int chip );
	bool    mipSigCheck( double *sig, int *TS );
	void    pulsePlotter( double *sig, int *TS, int ev, int ichip, int ich, int lowerR, int upperR );
	void    Crosstalk(Int_t ch);
	void    InitTH2Poly(TH2Poly& poly);  //Give frame to TH2Poly
	void    Gain_factor_producer();
	int     Cut(Long64_t entry, Long64_t sigma);

	// src_txtfile reader 
	void    yamlReader();
	void    GainFactorReader( string gainfile );
	void    noisyChannelReader( string noisyFileName );
	void    read_P_and_N(string ped_file);
	void    readmap();

	//private parameter
	TApplication   *app;
	TCanvas        *c;
	TTree          *Chain1;
	int            cross_ch_FirstRing[NCHIP][6];
	double         **hg_sig;
	double         **lg_sig;
	
	//pedestal parameter
	float          avg_HG[NCHIP][NCH][NSCA];
	float          sigma_HG[NCHIP][NCH][NSCA];
	float          avg_LG[NCHIP][NCH][NSCA];
	float          sigma_LG[NCHIP][NCH][NSCA];

	//yaml parameter
	int injCh;
	string acquisitionType;
	string ModuleNumber;
	int injChip;

	//gainFactor parameter
	double LG2HG_Conversion[NCHIP][NCH];
	double TOT2LG_Conversion[NCHIP][NCH];
	double HGTP[NCHIP][NCH];
	double LGTP[NCHIP][NCH];
	double TOTOffSet[NCHIP][NCH];
	double ADC2MIP = 0.0227;
	double LGTP_default = 900;

	//noisy parameter
	vector<int> noisyChannel;


  
	///////////////////////////////
	// Declaration of leaf types //
	///////////////////////////////
  
	//one event consists of 4 entries, every entry consists of one chip
	Int_t         event;          // event number, 
	Int_t         chip;           // chip number 
	Int_t         roll;           // rollposition
	Int_t         dacinj;         // injection dac number 
	Int_t         timesamp[13];   // Time sample [sca]
	Int_t         hg[13][64];     // high gain [ sca ] [ CHANNEL ]
	Int_t         lg[13][64];     // low  gain [ sca ] [ CHANNEL ]
	Int_t         tot_fast[64];   // not used 
	Int_t         tot_slow[64];   // usually use this to be the number for tot 
	Int_t         toa_rise[64];   // Timing 
	Int_t         toa_fall[64];   // Timing 
  
	// map < key = chip*32+ch/2 , pair <x, y> > 
	map<int,pair < double,double > > CHmap;
};

#endif
