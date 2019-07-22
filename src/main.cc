#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <time.h>

#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <cstdio>

#include "TGraph.h"
#include "TH2Poly.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TStyle.h"

using namespace std;
void readFile(string &fileName);
std::fstream& GotoLine(std::fstream& file, unsigned int num);
int decode_raw();
int format_channels();
int roll_position(); //Output is the location of TS 0
void InitTH2Poly(TH2Poly& poly);  //Give frame to TH2Poly
void readmap();


#define RAWSIZE 30784
#define totSCA 15  // 13SCAs + TOT + TOA
#define NSCA 13    // Real SCAs
#define NCHANNEL 256

unsigned char raw[RAWSIZE/2]; 
unsigned int ev[4][1924];
int dati[4][128][totSCA];
int global_TS[NSCA];
int dati_sum[4][128][NSCA];
double dati_sumsq[4][128][NSCA];
int mem_counter[4][128][NSCA];
map<int,pair < double,double > > CHmap;
TCanvas *c1;
//TApplication  *app;

float avg_hg[4][64];
float avg_lg[4][64];
int evt_counter;


int printData(){
	for(int chip = 0; chip < 4; chip++){
		for(int ch = 0; ch < 128; ch++){
			printf("chip %d ch %d : ", chip, 63-(ch%64));
			for(int sca = 0; sca < 15; sca++){
				printf("%d ", dati[chip][ch][sca]);
			}
			printf("\n");
		}
	}
	return(0);
}


int countTotalEvent(fstream& file){
	int totalEvent = 0;
	unsigned char tmp;
	while(!file.eof()){
		for(int i = 0 ; i < RAWSIZE/2; ++i){
			tmp = file.get();
		}
		tmp = file.get();
		tmp = file.get();
		totalEvent++;
	}
	file.clear();
	file.seekg(0,ios::beg);
	return totalEvent;
}


std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
	return file;
}

void pulsePlotter(){
	// This function plots the input 13 timesamples and show it on the screen
	int isca[13];
	for(int sca = 0; sca < 13; sca++){
		isca[sca] = sca;
	}

	for(int ichannel = 0; ichannel < 256; ichannel+=2){
		int ichip = ichannel / 64;
		int ich = ichannel % 64;
		if(ich != 38 ) continue;
		//int TS[NSCA];
		//for(int i = 0; i < NSCA; i++) TS[i] = global_TS[i];
		
		TGraph *gr = new TGraph(13, global_TS, dati[ichip][(63-ich)]);
		char plot_title[50];
		gr->SetMarkerColor(ichip+1);
		gr->SetMarkerStyle(22);
		gr->SetMarkerSize(1.2);
		gr->Draw("AP");
		sprintf(plot_title,"chip %d channel%d", ichip, ich);
		gr->SetTitle(plot_title);
		gr->GetXaxis()->SetTitle("TS");
		gr->GetYaxis()->SetTitle("ADC");
		c1->Update();
		gPad->WaitPrimitive();
	}
}


void drawPlots(int eventcount){

	char title[200];
	int maxTS;
	TH2Poly *poly = new TH2Poly;
	InitTH2Poly(*poly);
	double double_TS[NSCA];
	for(int sca = 0; sca < NSCA; sca++){
		if (global_TS[sca] == 2 )
			maxTS = sca;
	}
		
	for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
		int ichip = ichannel / 64;
		int ich = ichannel % 64;
		float X, Y;
		int forCH = ichannel / 2;
		bool NoisyBool = false;
		X = CHmap[forCH].first;
		Y = CHmap[forCH].second;
		poly->SetMinimum(-0.02);
		poly->SetMaximum(2000);
		poly->Fill(X,Y,avg_hg[ichip][ich]);
		//if(ich == 6 && dati[ichip][(63-ich)][maxTS] < 100)
		//printData();
	}
	sprintf(title,"test%d", eventcount);
	poly->SetTitle(title);
	poly->SetName(title);
	poly->Draw("colztext");
	c1->Update();
	//gPad->WaitPrimitive();
}
 


int decode_raw(){
    int i, j, k, m;
	bool m_compressedData = true;
    unsigned char x, y;
	unsigned int t;
	unsigned int bith, bit11, bit10, bit9, bit8, bit7, bit6, bit5, bit4, bit3, bit2, bit1, bit0;
	for( i = 0; i < 1924; i = i+1){
		for (k = 0; k < 4; k = k + 1){
			ev[k][i] = 0;
		}
   	}

	if( !m_compressedData ){
		for(int  i = 0; i < 1924; i++){
			for (int j = 0; j < 16; j++){
				x = raw[i*16 + j];
				x = x&0xf;
				for (int sk = 0; sk < 4; sk++)
					ev[sk][i] = ev[sk][i] | (uint16_t) ( ((x>>sk) & 1) << (15 - j));
			}
		}
	}
	else{
		for(int  i = 0; i < 1924; i++){
			for (int j = 0; j < 8; j++){
				x = raw[i*8 + j];
				//cout << (i+1)*(j+1) << " " << (int)x << endl;
				y = (x>>4)&0xf;
				x = x&0xf;
				for (int sk = 0; sk < 4; sk++){
					ev[sk][i] = ev[sk][i] | (uint16_t) ( ((x>>sk) & 1) << (14 - j*2));
					ev[sk][i] = ev[sk][i] | (uint16_t) ( ((y>>sk) & 1) << (15 - j*2));
				}
			}
		}
	}

    
    /*****************************************************/
    /*    Gray to binary conversion                      */
    /*****************************************************/
   	for(k = 0; k < 4 ; k = k +1 ){
   		for(i = 0; i < 1920; i = i + 1){
   			bith = ev[k][i] & 0x8000;
   			t = ev[k][i] & 0x7fff;
   			bit11 = (t >> 11) & 1;
        	bit10 = bit11 ^ ((t >>10) &1);
        	bit9 = bit10 ^ ((t >>9) &1);
        	bit8 = bit9 ^ ((t >>8) &1);
        	bit7 = bit8 ^ ((t >>7) &1);
        	bit6 = bit7 ^ ((t >>6) &1);
        	bit5 = bit6 ^ ((t >>5) &1);
        	bit4 = bit5 ^ ((t >>4) &1);
        	bit3 = bit4 ^ ((t >>3) &1);
        	bit2 = bit3 ^ ((t >>2) &1);
        	bit1 = bit2 ^ ((t >>1) &1);
        	bit0 = bit1 ^ ((t >>0) &1);
        	ev[k][i] =  bith | ((bit11 << 11) + (bit10 << 10) + (bit9 << 9) + (bit8 << 8) + (bit7 << 7) + (bit6 << 6) + (bit5 << 5) + (bit4 << 4) + (bit3  << 3) + (bit2 << 2) + (bit1  << 1) + bit0);
        }
    }
	return 0;
}


int format_channels(){
	for(int chip = 0; chip < 4; chip = chip +1 ){
        for(int ch = 0; ch < 128; ch = ch +1 ){
            for(int sca = 0 ; sca < 15 ; sca = sca +1){
                dati[chip][ch][sca] = ev[chip][sca*128+ch] & 0x0FFF;
            }
        }
    }
    return(0);
}


int roll_position(){
	unsigned int roll_check;  //Just to check if 4 chip has same rollmask
	int chip,first,sec,rollpos;
	bool skip_evt = false;
	for (chip =0; chip < 4; chip = chip +1 ){
		unsigned int roll;
		roll = ev[chip][1920] & 0x1FFF;
		//printf("roll pos = %x\n",roll);
		if(chip == 0) roll_check = roll;
		else if( roll_check != roll ){
			cout << "Problematic event!( No. " << evt_counter <<") Chip " << chip
				 << "has different rollMask! Skip event!" << endl;
			skip_evt = true;
		}

		unsigned char bits[13];
		first = -1; // first is actually second XD
		sec   = -1;
		for(int bit = 0; bit < 13 ; ++bit){
			bits[bit] = (roll >> bit) & 1;}
		for(int bit = 0 ; bit < 13; ++bit) {if((int)bits[bit] == 1) first = bit;}
		for(int bit = 0 ; bit < 13; ++bit) {
			if((int)bits[bit] == 1 && first != bit) sec = bit;}
		if(first == 12 && sec == 11) rollpos = 0;
		else if(first == 12 && sec == 0)  rollpos = 1;
		else rollpos = first+1;
		
		//for(int bit = 0; bit < 13 ; ++bit)      cout << (int)bits[bit] << " ";
		//cout << first << " , " << sec << ", rollpos = " << rollpos << endl;
		//getchar();
	}
	if(skip_evt){ return (-1); }

	for(int sca = 0; sca < NSCA; sca++){
		global_TS[sca] = (12 - sca + abs(12 - rollpos - 10)) % 13;
		//cout << sca << " " << global_TS[sca] << " ";
	}

	return rollpos;
}


int average_data(){
	int maxTS;
	for(int sca = 0; sca < NSCA; sca++){
		if (global_TS[sca] == 2 )
			maxTS = sca;
	}

	for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
		int ichip = ichannel / 64;
		int ich = ichannel % 64;
		avg_hg[ichip][ich] += dati[ichip][127-ich][maxTS];
		avg_lg[ichip][ich] += dati[ichip][63 -ich][maxTS];
	}
	return(0);
}

/*
int writeData(){
	ofstream outf("HGTS2");
	
	
}
*/


void InitTH2Poly(TH2Poly& poly)
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

void readmap(){
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



int main(int argc, char** argv){
	string s1;
	string filename;
	string line;
	string arg_string;
	int avgEventNumber = 5;
	vector<string> arg_list;
  
	for(int i = 0 ; i < argc ; ++i){
		arg_string = argv[i];
		arg_list.push_back(arg_string);
	}
	TApplication theApp("App",&argc, argv);
	c1 = new TCanvas();
	gStyle->SetOptStat(false);
	readmap();

	clock_t t;
	
	while(true){
		t = clock();
		unsigned char tmp;
		fstream myFile(arg_list[1]);
		//string s = "scp -q pi@192.168.50.196:~/rpi-daq/test.txt ./";
		//system(s.c_str());
		int totalEvent = countTotalEvent(myFile);
		cout << totalEvent << endl;
		if(totalEvent < avgEventNumber+1) continue;

		for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
			int ichip = ichannel / 64;
			int ich = ichannel % 64;
			avg_hg[ichip][ich] = 0;
			avg_lg[ichip][ich] = 0;
		}
		
		int event = 0;
		while( event < totalEvent-avgEventNumber-2){
			for(int i = 0 ; i < RAWSIZE/2; ++i){
				raw[i] = myFile.get();
			}
			tmp = myFile.get();
			tmp = myFile.get();
			event++;
		}

		event = 0;
		while( event < avgEventNumber ){
			for(int i = 0 ; i < RAWSIZE/2; ++i){
				raw[i] = myFile.get();
			}
			tmp = myFile.get();
			tmp = myFile.get();
			
			decode_raw();
			format_channels();
			int rollpos = roll_position();
			average_data();
			//cout << rollpos << endl;

			//printData();
			//pulsePlotter();
			//myFile.close();
			event++;
		}
		gSystem->ProcessEvents();
		for(int ichannel = 0; ichannel < NCHANNEL; ichannel+=2){
			int ichip = ichannel / 64;
			int ich = ichannel % 64;
			avg_hg[ichip][ich] /= avgEventNumber;
			avg_lg[ichip][ich] /= avgEventNumber;
		}
		drawPlots(totalEvent);
		int dt = clock() - t;
		cout << (float)dt/CLOCKS_PER_SEC << endl;
		//gPad->WaitPrimitive();
	}
	//myFile >> s1;
	return 0;
}
