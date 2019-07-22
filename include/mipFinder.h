
#ifndef makePlots_h
#define makePlots_h

#include "makeplots.h"

using namespace std;

class makePlots{
 public:

  makePlots (TChain* inchain);
  ~makePlots();

  //public function
  void Init( string pedfile, string gainfile );
  void PlotProducer();
  void Inj_Pulse_display( int displayChannel = -1 , bool inj_flag = false );

  //public parameter
  string input_fileName;
  
 private:
  
  ///////////////////////////////
  // Declaration of leaf types //
  ///////////////////////////////
  
  Int_t         event;
  Int_t         chip;
  Int_t         roll;
  Int_t         dacinj;
  Int_t         timesamp[13];
  Int_t         hg[13][64];
  Int_t         lg[13][64];
  Int_t         tot_fast[64];
  Int_t         tot_slow[64];
  Int_t         toa_rise[64];
  Int_t         toa_fall[64];
  
  // map < key = chip*32+ch/2 , pair <x, y> > 
  map<int,pair < double,double > > CHmap;
};

#endif

