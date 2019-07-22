#include "PlotSetting.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TImage.h"
#include "TMultiGraph.h"
#include "TLatex.h"

//Constructo
PlotSetting::PlotSetting()
{
}
//Destructor
PlotSetting::~PlotSetting()
{
}

void PlotSetting::HStd(TH1& h, char* plot_title, string Xtitle, string Ytitle, bool Wait, bool SavePlot) 
{
	if(!Wait){ gROOT->SetBatch(kTRUE); }
	else { gROOT->SetBatch(kFALSE); }
	TCanvas* c = new TCanvas();
	h.SetTitle(plot_title);
	h.GetXaxis()->SetTitle(Xtitle.c_str());
	h.GetYaxis()->SetTitle(Ytitle.c_str());
	h.GetYaxis()->SetTitleOffset(1.3);
	h.Draw();
	c->Update();
	if(Wait){ gPad->WaitPrimitive(); }
	if(SavePlot){
		char title[200];
		sprintf(title,"%s/%s.png",plotfolder_path,plot_title);
		TImage *img = TImage::Create();
		img->FromPad(c);
		img->WriteImage(title);
		cout << title << " has been saved" << endl;
	}
	delete c;  
}

void PlotSetting::GStd(TGraph& g, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot) 
{
	if(!Wait){ gROOT->SetBatch(kTRUE); }
	else { gROOT->SetBatch(kFALSE); }
	TCanvas* c = new TCanvas();
	g.SetTitle(plot_title);
	g.GetXaxis()->SetTitle(Xtitle.c_str());
	g.GetYaxis()->SetTitle(Ytitle.c_str());
	g.GetYaxis()->SetTitleOffset(1.3);
	g.Draw(Option.c_str());
	c->Update();
	if(Wait){ gPad->WaitPrimitive(); }
	if(SavePlot){
		char title[200];
		sprintf(title,"%s/%s.png",plotfolder_path,plot_title);
		TImage *img = TImage::Create();
		img->FromPad(c);
		img->WriteImage(title);
		cout << title << " has been saved" << endl;
	}
	delete c;  
}

void PlotSetting::G(TGraph& g, char* plot_title, string Xtitle, string Ytitle,int MarkerStyle,int MarkerColor,
					int MarkerSize, int LineColor, int LineWidth, string Option,bool OptStat, bool Wait, bool SavePlot) 
{
	if(!Wait){ gROOT->SetBatch(kTRUE); }
	else { gROOT->SetBatch(kFALSE); }
	TCanvas* c = new TCanvas();
	gStyle->SetOptStat(OptStat);
	g.SetTitle(plot_title);
	g.GetXaxis()->SetTitle(Xtitle.c_str());
	g.GetYaxis()->SetTitle(Ytitle.c_str());
	g.GetYaxis()->SetTitleOffset(1.2);
	g.SetMarkerStyle(MarkerStyle);
	g.SetMarkerColor(MarkerColor);
	g.SetLineColor(LineColor);
	g.SetLineWidth(LineWidth);
	g.SetMarkerSize(MarkerSize);
	g.Draw(Option.c_str());
	c->Update();
	if(Wait){ gPad->WaitPrimitive(); }
	if(SavePlot){
		char title[200];
		sprintf(title,"%s/%s.png",plotfolder_path,plot_title);
		TImage *img = TImage::Create();
		img->FromPad(c);
		img->WriteImage(title);
		cout << title << " has been saved" << endl;
	}
	delete c;  
}


void PlotSetting::Multi(TMultiGraph& g, TLegend& legend, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot) 
{
	if(!Wait){ gROOT->SetBatch(kTRUE); }
	else { gROOT->SetBatch(kFALSE); }
	TCanvas* c = new TCanvas();
	g.Draw(Option.c_str());
	g.SetTitle(plot_title);
	g.GetXaxis()->SetTitle(Xtitle.c_str());
	g.GetYaxis()->SetTitle(Ytitle.c_str());
	g.GetYaxis()->SetTitleOffset(1.3);
	legend.SetBorderSize(0);
	legend.Draw();
	c->Update();
	if(Wait){ gPad->WaitPrimitive(); }
  
	if(SavePlot){
		char title[200];
		sprintf(title,"%s/%s.png",plotfolder_path,plot_title);
		TImage *img = TImage::Create();
		img->FromPad(c);
		img->WriteImage(title);
		cout << title << " has been saved" << endl;
	}
	delete c;
}

void PlotSetting::MultiAdd(TMultiGraph& multig, TGraph& g, TLegend& legend, char* plot_title, int MarkerStyle,int MarkerColor, int MarkerSize)
{
	g.SetTitle(plot_title);
	g.SetMarkerStyle(MarkerStyle);
	g.SetMarkerColor( Color(MarkerColor) );
	g.SetMarkerSize(MarkerSize);
	g.SetLineColor( Color(MarkerColor) );
	g.SetLineWidth(4);
	multig.Add(&g);
	TLegendEntry* lentry = legend.AddEntry(&g, plot_title,"l");
	//lentry->SetTextColor(MarkerColor);
}

void PlotSetting::Poly(TH2Poly& poly, char* plot_title, string Xtitle, string Ytitle, string Option,bool OptStat, bool Wait, bool SavePlot) 
{
	if(!Wait){ gROOT->SetBatch(kTRUE); }
	else { gROOT->SetBatch(kFALSE); }
	TCanvas* c = new TCanvas();
	TLatex latex;
	latex.SetTextSize(0.05);
	latex.SetTextAlign(13);  //align at top
	gStyle->SetOptStat(OptStat);
	gStyle->SetPalette(55);
	poly.SetTitle(plot_title);
	poly.GetXaxis()->SetTitle(Xtitle.c_str());
	poly.GetYaxis()->SetTitle(Ytitle.c_str());
	poly.GetYaxis()->SetTitleOffset(1.);
	poly.Draw(Option.c_str());
	//  latex.DrawLatex(5.2,7.8,"<E / EInj>_{EInj=200}^{EInj=1000}");
	latex.DrawLatex(6,7.5,"<E / EInj>");
	c->Update();
     
	if(Wait){ gPad->WaitPrimitive(); }

	if(SavePlot){
		char title[200];
		sprintf(title,"%s/%s.png",plotfolder_path,plot_title);
		TImage *img = TImage::Create();
		img->FromPad(c);
		img->WriteImage(title);
		cout << title << " has been saved" << endl;
	}
	delete c;
}

void PlotSetting::root_logon(){

	cout << endl << "Welcome to the ATLAS rootlogon.C" << endl;

	//
	// based on a style file from BaBar
	//

	//..BABAR style from RooLogon.C in workdir
	TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");

	// use plain black on white colors
	Int_t icol=0;
	atlasStyle->SetFrameBorderMode(icol);
	atlasStyle->SetCanvasBorderMode(icol);
	atlasStyle->SetPadBorderMode(icol);
	atlasStyle->SetPadColor(icol);
	atlasStyle->SetCanvasColor(icol);
	atlasStyle->SetStatColor(icol);
	//atlasStyle->SetFillColor(icol);

	// set the paper & margin sizes
	atlasStyle->SetPaperSize(20,26);
	atlasStyle->SetPadTopMargin(0.05);
	atlasStyle->SetPadTopMargin(0.11);
	atlasStyle->SetPadRightMargin(0.13);
	atlasStyle->SetPadBottomMargin(0.15);
	atlasStyle->SetPadLeftMargin(0.12);

	// use large fonts
	//Int_t font=72;
	Int_t font=42;
	Double_t tsize=0.045;
	atlasStyle->SetTextFont(font);



	atlasStyle->SetTextSize(tsize);
	atlasStyle->SetLabelFont(font,"x");
	atlasStyle->SetTitleFont(font,"x");
	atlasStyle->SetLabelFont(font,"y");
	atlasStyle->SetTitleFont(font,"y");
	atlasStyle->SetLabelFont(font,"z");
	atlasStyle->SetTitleFont(font,"z");

	atlasStyle->SetLabelSize(tsize,"x");
	atlasStyle->SetTitleSize(tsize,"x");
	atlasStyle->SetLabelSize(tsize,"y");
	atlasStyle->SetTitleSize(tsize,"y");
	atlasStyle->SetLabelSize(tsize,"z");
	atlasStyle->SetTitleSize(tsize,"z");

	atlasStyle->SetLabelOffset(tsize/4,"x");
	atlasStyle->SetLabelOffset(tsize/4,"y");
  


	//use bold lines and markers
	atlasStyle->SetMarkerStyle(20);
	atlasStyle->SetMarkerSize(0.5);
	atlasStyle->SetHistLineWidth(2.);
	atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

	//
	atlasStyle->SetLegendFillColor(0);
	atlasStyle->SetLegendFont(62);
	//atlasStyle->SetLegendTextSize(0.03);

	//get rid of X error bars and y error bar caps
	//atlasStyle->SetErrorX(0.001);

	//do not display any of the standard histogram decorations
	//  atlasStyle->SetOptTitle(0);
	//atlasStyle->SetOptStat(1111);
	atlasStyle->SetOptStat(0);
	//atlasStyle->SetOptFit(1111);
	atlasStyle->SetOptFit(0);

	// put tick marks on top and RHS of plots
	atlasStyle->SetPadTickX(1);
	atlasStyle->SetPadTickY(1);

	//  gROOT->SetStyle("Plain");

	//gStyle->SetPadTickX(1);
	//gStyle->SetPadTickY(1);

	gROOT->SetStyle("ATLAS");
	gROOT->ForceStyle();
}

int PlotSetting::Color(int c)
{
	if(c == 0){ return 633;}
	else if(c == 1){ return 600;}
	else if(c == 2){ return 800;}
	else if(c == 3){ return 419;}
	else if(c == 4){ return 880;}
	else if(c == 5){ return 803;}
	else{ return 0; }
}

