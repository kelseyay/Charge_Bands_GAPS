//Hopefully this script will provide a primary beta from energy depositions in Not That Much Time.
//This will start simple with beta from at least three tracker hits with energy deposition > 0.4

//WORK IN PROGRESS!! Use on grond data only
//Apply TOF cuts (either by slimming the data or by implementing them here, or maybe both as a failsafe! I like the idea of both because if you run on a slimmed dataset, it doesn't take much time)

// How to use:
// ./GndTres


using namespace std;

#include "KYtools.C"
#include <boost/math/special_functions/lambert_w.hpp>
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
#include "TStyle.h"
#include <Math/RootFinder.h>
#include <Math/Functor.h>

#include "GGeometry.hh"
using namespace ROOT::Math;
using namespace std;
//Headers included in KYtools are included here, so this is nice.

//FIXME: does this work on mac?
#include <sys/stat.h>

void histplot1f_fit(float& mean, float& std, string ctitle, TH1F* h1, string title, string xtitle, string ytitle, string savename){
    TCanvas * c1 = new TCanvas(ctitle.c_str(), ctitle.c_str(), 200, 10, 900, 900);
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.1);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.1);

    h1->SetTitle(title.c_str());
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->SetLineWidth(2);

    h1->Fit("gaus");
    TF1 *fitFunc = h1->GetFunction("gaus");
    mean = fitFunc->GetParameter(1);
    std = fitFunc->GetParameter(2);
    //cout << "Mean: " << fitFunc->GetParameter(1) << endl;
    //cout << "Sigma?: " << fitFunc->GetParameter(2) << endl;
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy(0);
    gStyle->SetTitleW(0.9);
    gStyle->SetOptFit();
    h1->Draw();
    //h1->Draw("SAME L");

    char histname[400];
    sprintf(histname, "%s.png",savename.c_str());
    c1->SaveAs(histname);
}

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<string>("end_name", "name at the end of pngs","","e");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.2,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1.0,"u");
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->ParseCommandLine(argc, argv);
parser->Parse();

double betacut = parser->GetOption<double>("beta_low");
if(betacut <= 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
string end_name = parser->GetOption<string>("end_name");
cout << "reco_path " << reco_path << endl;
cout << "out path " << out_path << endl;
//cout << "out path last string " << out_path[out_path.length()-1] << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

int MainLoopScaleFactor = parser->GetOption<int>("MainloopScale");

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

string txtname = "MCEdep_Beta.txt";

//Prepare textile for saving values
std::ofstream myfile;
myfile.open(out_path + txtname);
myfile << TString::Format( "Filename : %s", reco_path.c_str() )  << endl;
myfile << TString::Format( "Beta High : %f", betahigh) << endl;
myfile << TString::Format( "Beta Low : %f", betacut) << endl;
myfile.close();

//Prepare reconstructed event
CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare MC event
CEventMc* MCEvent = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
TreeMC->SetBranchAddress("Mc", &MCEvent); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);

//auto tkr_edep_beta_z1 = new TF1("tkr_edep_beta_z1","OneZ_TKR(x)",0.01,0.95);

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0.1; //No low Tof cut right now
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
int NHitsmin = 2;

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8

int TRG = parser->GetOption<int>("TRG");

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

const int nBins = 1;

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

TH1F*HBeta = new TH1F ("BetaRec","Ground Brec",70, 1.0-0.5, 1.0+0.5);

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 1000; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
		    cout << "Event number " << i << endl;
	}

	//Cuts are implemented in this chunk:
	//Currently there is space for the 1PP/3PP hit but no calculation
	if(Event->GetNTracks() == 1){  //First select the single track event
	    bool Outer_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
		bool Inner_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
		bool TKRflag = 0;
		int Layer_Hits_Tracker[7] = {}; //Seven layers

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if( pt != nullptr && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) && pt->GetChi2()/pt->GetNdof() < 3.2 /*&& -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh*/ && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(Event->GetPrimaryBeta()) > betacut && fabs(Event->GetPrimaryBeta()) <  betahigh ){
			//-----------EVENT LEVEL CUT APPLIED
			    double Brec = Event->GetPrimaryBeta();
				cout << "Event is " << i << " beta reco is " << Brec << endl;
				HBeta->Fill(Brec);


			}

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i


float mu_rec;
float sigma_rec;
histplot1f_fit(mu_rec,sigma_rec,"Hrec", HBeta, "Beta Rec", "Beta Rec", "Nentries", out_path + "BRecGround"+ end_name);

myfile.open(out_path + txtname,std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}
