// Author : Hsin-Yeh Wu
// E-mail : thankyouyou06@gmail.com
//
// This class is for general plot settings.
// If you have a kind of plotsetting frequently used, you could set a function here for it 

#ifndef PlotSetting_h
#define PlotSetting_h

#include "TROOT.h"
#include "TH2Poly.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <string>
using namespace std;

class PlotSetting{
 public:

	PlotSetting ();
	~PlotSetting ();

	char plotfolder_path[200];

	// public function
	void root_logon();  // A general setting for the frame 
	int Color(int c);   // A set of colors I like to use ( not root default colors ) 
	
	// Frequently used plot settings 
	void HStd(TH1& h, char* plot_title, string Xtitle, string Ytitle, bool Wait, bool SavePlot); 
	void GStd(TGraph& g, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot); 
	void G(TGraph& g, char* plot_title, string Xtitle, string Ytitle,int MarkerStyle,int MarkerColor,
		   int MarkerSize, int LineColor, int LineWidth, string Option,bool OptStat, bool Wait, bool SavePlot);
	void Multi(TMultiGraph& g, TLegend& legend, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot);
	void MultiAdd(TMultiGraph& multig, TGraph& g, TLegend& legend, char* plot_title, int MarkerStyle,int MarkerColor, int MarkerSize);
	void Poly(TH2Poly& poly, char* plot_title, string Xtitle, string Ytitle, string Option,bool OptStat, bool Wait, bool SavePlot);


 private:
  
};


#endif
