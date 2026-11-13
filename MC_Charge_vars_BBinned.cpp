//Hopefully this script will provide a primary beta from energy depositions in Not That Much Time.
//This will start simple with beta from at least three tracker hits with energy deposition > 0.4

//WORK IN PROGRESS!!

// How to use:
// ./BinZvars -i /home/kelsey/simulations/simdat/proton/v3.0.0/slimmed_TL1/MC_alpha_TL1 -o Alpha_v300_TL1/More_Test/ -u 1.2


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

float Ztof = 5.574;
float Atof = 10.3;
float Ltof = 0.635;
float rtof = 1.032;
float iontof = (0.000016 * pow(Ztof,0.9));

float z = 1;
float Ztkr = 14;
float Atkr = 28;
float rtkr = 2.33; //Density g/cm^2
float Ltkr = 0.25; //Mmk so for the MC, 0.95 beta (limit of solver) corresponds to 0.787 MeV

float iontkr = (0.000016 * pow(Ztkr,0.9));
float C_1 = 0.3071/2; //MeV/ g/cm^2 #2*pi*constants not 4*pi*constants for MPV
float me = 0.511; //mass of electron * c^2

double ZOne_TKR(double x)
{
    return C_1*pow(z,2)*(rtkr*Ltkr*Ztkr/Atkr)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/iontkr ) + log(C_1*pow(z,2)*(rtkr*Ltkr*Ztkr/Atkr)*(1/pow(x,2))/iontkr)  + 0.2 - pow(x,2) );
}

double ZOne_TOF(double x)
{
    return C_1*pow(z,2)*(rtof*Ltof*Ztof/Atof)*(1/pow(x,2))* ( log( (1.022 * pow(x,2)/(1 - pow(x,2)) )/iontof ) + log(C_1*pow(z,2)*(rtof*Ltof*Ztof/Atof)*(1/pow(x,2))/iontof)  + 0.2 - pow(x,2) );
}

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

float solvermin_tkr = ZOne_TKR(0.95); //0.787;
float solvermin_tof = ZOne_TOF(0.95);

//Try figuring out why this worked: https://root.cern/doc/v620/classTF1.html
TF1 *Z1_tkr_Solve = new TF1("Z1_tkr_Solve", [](double *x, double *p){ return ZOne_TKR(x[0]); }, 0.01, 0.95, 0);
TF1 *Z1_tof_Solve = new TF1("Z1_tof_Solve", [](double *x, double *p){ return ZOne_TOF(x[0]); }, 0.01, 0.95, 0);

//Example:
//TGraph * g = new TGraph(npointx, xvec, yvec);
//TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*g->Eval(x[0]); }, xmin, xmax, 1);

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<string>("end_name", "name at the end of pngs","","e");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.2,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1.0,"u");
parser->AddCommandLineOption<double>("tkr_factor", "tkr factor",1,"k");
parser->AddCommandLineOption<double>("tof_factor", "tof factor",1,"f");
parser->AddCommandLineOption<bool>("TKR", "tkr use",1,"s");
parser->AddCommandLineOption<bool>("TF", "tof use",1,"t");
parser->AddCommandLineOption<int>("TRG", "Which trigger?",0,"r");
parser->AddCommandLineOption<bool>("print", "print out info?",0,"n");
parser->AddCommandLineOption<int>("MainloopScale", "Main loop scale factor",1,"m");
parser->AddCommandLineOption<bool>("MC_Weighting", "Monte Carlo Weighting on or off?",0,"w");
parser->ParseCommandLine(argc, argv);
parser->Parse();

bool MC_Weight = parser->GetOption<bool>("MC_Weighting");
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
myfile << TString::Format( "Charge : %f", z) << endl;
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

bool TKR = parser->GetOption<bool>("TKR");
bool TF = parser->GetOption<bool>("TF");
int TRG = parser->GetOption<int>("TRG");
int print = parser->GetOption<bool>("print");

double tk = parser->GetOption<double>("tkr_factor");
double tf = parser->GetOption<double>("tof_factor");

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

const int NBbins = 4;

TH1F * HBetaProxy[NBbins];
TH1F * HTOF_Prim[NBbins];
TH1F * HTKR_TrM[NBbins];

const float Bedges[] = {0.5, 0.61, 0.73, 0.9, 1.1}; //Edges of the histogram

for(int i = 0; i < NBbins; i++){
    HBetaProxy[i] = new TH1F (Form("HBetaProxy%i",i),Form("HBetaProxy %f",Bedges[i]),70, 0.1,1.2);
    HTOF_Prim[i] = new TH1F (Form("HTOF_Prim%i",i),Form("HTOF_Prim %f",Bedges[i]),80, 0,30);
    HTKR_TrM[i] = new TH1F (Form("HTKR_TrM%i",i),Form("HTKR_TrM %f",Bedges[i]), 70, 0.5,4);
}

int bcounts = 0;
TH2D * HRecB_vs_GenB = new TH2D("HRecB_vs_GenB","Rec_Beta vs Gen_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1.1);
TH2D * HProxB_vs_GenB = new TH2D("HProxB_vs_GenB","Prox_B vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1, 50, 0.1, 1.1);
TH2D * HRecB_vs_ProxB = new TH2D("HRecB_vs_ProxB","Rec_Beta vs Prox_B",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1.1);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 1000; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

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
		if( pt != nullptr && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) /*&& pt->GetChi2()/pt->GetNdof() < 3.2 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh*/ && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) > betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
			//-----------EVENT LEVEL CUT APPLIED

			//cout << "Event: " << i << endl;
			//First iteration over event for flags and EnergyDepositionMip vector filling
			double Bgen = MCEvent->GetPrimaryBeta();
			//double gamma = sqrt( 1/(1 - pow(Bgen,2)) );
			//cout << "Bgen " << Bgen << endl;

			for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){

                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                //cout << "VolmeId = " << VolumeId << endl;
                if(volspec(VolumeId,0,2) == 20 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                    int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                    Layer_Hits_Tracker[layer]++;
                    TKRflag = 1;
                    //cout << "TKR Hit" << endl;
                }
                //cout << "Event is " << i << endl;
                if(volspec(VolumeId,0,2) == 10 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){ Outer_TOF_flag[volspec(VolumeId,2,1)] = 1; }
                if(volspec(VolumeId,0,2) == 11 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){ Inner_TOF_flag[volspec(VolumeId,2,1)] = 1; }
			}

			//Do we want to only run this on certain tracks? Yeah probably. Can remove the TOF flags. Just run on whatever lol.
			if( Event->GetTrack(0)->GetEnergyDeposition().size() >= NHitsmin ) {
			    //Note, runnning on the slimmed dataset already has the required cuts.
			//if( (Inner_TOF_flag[1] || TKRflag) && (Outer_TOF_flag[0] && Inner_TOF_flag[0]) ) { //Guarantees 3 hit requirement, if the energy deposition is lower than the solver minimum, don't use. Still try to calculate.
			//CBE_bot and TKR required. For this study I'm doing right this second no UMB and CBEtop
			    //cout << "Event " << i << " is a 2-hit event! Try to calculate BP! " << endl;
			//Sides only below, exclude events with UMB and CBEtop
			//if( !(Outer_TOF_flag[0] && Inner_TOF_flag[0]) && ( (Outer_TOF_flag[2] && Inner_TOF_flag[2]) || (Outer_TOF_flag[3] || Inner_TOF_flag[3]) || (Outer_TOF_flag[4] && Inner_TOF_flag[4]) || (Outer_TOF_flag[5] && Inner_TOF_flag[5]) ) && (Inner_TOF_flag[1] || TKRflag) ) { //Guarantees 3 hit requirement, if the energy deposition is lower than the solver minimum, don't use. Still try to calculate.

			//All things below
			//if( ((Outer_TOF_flag[0] && Inner_TOF_flag[0]) || (Outer_TOF_flag[2] && Inner_TOF_flag[2]) || (Outer_TOF_flag[3] || Inner_TOF_flag[3]) || (Outer_TOF_flag[4] && Inner_TOF_flag[4]) || (Outer_TOF_flag[5] && Inner_TOF_flag[5]) ) && (Inner_TOF_flag[1] || TKRflag) ) { //Guarantees 3 hit requirement, if the energy deposition is lower than the solver minimum, don't use. Still try to calculate.
			    if(print) cout << endl << "Event is " << i << endl;
				if(print) cout << "Zenith angle is " << Event->GetPrimaryMomentumDirection().CosTheta() << endl;
			    if(print) cout << "Bgen is " << Bgen << endl;
				if(print) cout << "Breco is " << Event->GetPrimaryBeta() << endl;

				vector<double> TOF_prim; //Want a vector for the TOF_prim variable
			    vector<double> Beta_Proxy;
                vector<double> Tr_Mean_TKR;

                vector<double> Combine_Edeps;
                vector<int> Combine_layer;
                vector<int> Combine_mod;
                vector<int> Combine_detector;
                vector<int> Combine_strip;

                double Brec = Event->GetPrimaryBeta();

				//cout << "Event is " << i << endl;

				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
				    unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
					float Edep = Event->GetTrack(0)->GetEnergyDeposition(isig);
					float Ang_Edep = 0;

					if( (volspec(VolumeId,0,3) == 100 || volspec(VolumeId,0,3) == 110 ) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){
					    //Add UMB and CBEtop hits to the TOF_prim vector
						TOF_prim.push_back(Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta()));
					}

					if(TF && GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){
					    if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){ //It's a flat paddle.
					        Ang_Edep = Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta());
							if(Ang_Edep > solvermin_tof){
							    double root = Z1_tof_Solve->GetX(Ang_Edep);
								if(print)cout << "Hit is " << isig << endl;
								if(print)cout << "Energy Deposition at Volid " << VolumeId << " is " << Edep << " angle corrected is " << Ang_Edep << endl;
                                if(print)cout << "Calculated Beta at Volid = "  << root << endl;
                                Beta_Proxy.push_back(root);
							}
                        }else{ //It's a vertical paddle
                            //The calculation of the correction factor depends on the orientation of the paddle
                            float xy_path_corr = 1;
                            float x_vec = Event->GetTrack(0)->GetMomentumDirection()[0][0];
                            float y_vec = Event->GetTrack(0)->GetMomentumDirection()[0][1];
                            if(volspec(VolumeId,2,1) == 2 ||volspec(VolumeId,2,1) == 3){ //+/- X paddles
                                if(print)cout << "xvec, yvec = " << x_vec << " " << y_vec << endl;
                                if(print)cout << "cos(phi) = " << fabs(x_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) ) << endl;
                                xy_path_corr = fabs(x_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) );
                            }
                            if(volspec(VolumeId,2,1) == 4 ||volspec(VolumeId,2,1) == 5){ //+/- Y paddles
                                if(print)cout << "xvec, yvec = " << x_vec << " " << y_vec << endl;
                                if(print)cout << "cos(phi) = " << fabs(y_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) ) << endl;
                                xy_path_corr = fabs(y_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) );
                            }

                            Ang_Edep = Edep*fabs( sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2)))*xy_path_corr;
                            if(Ang_Edep > solvermin_tof){
							    double root = Z1_tof_Solve->GetX(Ang_Edep);
								if(print)cout << "Hit is " << isig << endl;
								if(print)cout << "Energy Deposition at Volid " << VolumeId << " is " << Edep << " angle corrected is " << Ang_Edep << endl;
                                if(print)cout << "Calculated Beta = "  << root << endl;
                                Beta_Proxy.push_back(root);
							}
                        }
					}

					if(TKR && volspec(VolumeId,0,2) == 20 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
						int layer = GGeometryObject::GetTrackerLayer(VolumeId);
						int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);
						int sdmod = GGeometryObject::GetLayerModule(VolumeId);
						int sddet = GGeometryObject::GetModuleDetector(VolumeId);
						Ang_Edep = Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta());

						if(Layer_Hits_Tracker[layer] == 1){
						    if(Ang_Edep > solvermin_tkr){
						        double root = Z1_tkr_Solve->GetX(Ang_Edep);
						        Beta_Proxy.push_back(root);
							}
							Tr_Mean_TKR.push_back(sqrt(Ang_Edep));
						}

						if(Layer_Hits_Tracker[layer] == 2){
						    //cout << "Event " << i << " has a double hit in the tracker! lmds = " << layer << sdmod << sddet << sdstrip << endl;
							bool new_com = 1; //If the combined hit is not found in the vecors, then it must be added. Assume it is new unless it's found!
							for(int com = 0; com < Combine_Edeps.size();com++){ //Check the vectors that track if hits are combined
							    if(layer == Combine_layer.at(com) && sdmod == Combine_mod.at(com) && sddet == Combine_detector.at(com) ){
									//cout << "Hewwo, I am probably the friend of that hit!" << endl;
									//cout << "With our powers combined " << Combine_Edeps.at(com) << " + " << Ang_Edep << endl;
									Tr_Mean_TKR.push_back(sqrt(Ang_Edep + Combine_Edeps.at(com)));
					                if(Ang_Edep + Combine_Edeps.at(com) > solvermin_tkr){
					                    double root = Z1_tkr_Solve->GetX(Ang_Edep + Combine_Edeps.at(com) );
					                    Beta_Proxy.push_back(root);
									}
								}
							}
							if(new_com){
							    Combine_Edeps.push_back(Ang_Edep);
								Combine_layer.push_back(layer);
								Combine_mod.push_back(sdmod);
								Combine_detector.push_back(sddet);
								Combine_strip.push_back(sdstrip);
							}

 						}

						/*
				        if(sdstrip != 0 && sdstrip != 7 && Layer_Hits_Tracker[layer] < 2){
							if(Ang_Edep > solvermin_tkr){
							    double root = Z1_tkr_Solve->GetX(Ang_Edep);
								if(print)cout << "Hit is " << isig << endl;
								if(print)cout << "Energy Deposition at Volid " << VolumeId << " is " << Edep << " angle corrected is " << Ang_Edep << endl;
                                if(print)cout << "Calculated Beta at Volid = "  << root << endl;
                                Beta_Proxy.push_back(root);
							}
						}*/
				    }
				} //Closed bracket for iteration over isig with TOF cuts

                double TOF_prim_Ave = 0;
                double CtrTOF_prim_Ave = 0;
                if(TOF_prim.size() < 2){
                    TOF_prim_Ave = 0;
                }else{
                    for(unsigned int isig = 0; isig < double(TOF_prim.size()); isig++){
              		TOF_prim_Ave += TOF_prim.at(isig);
              		CtrTOF_prim_Ave++;
               	}
               	if(CtrTOF_prim_Ave == 0){ TOF_prim_Ave = 0;
               	}else{ TOF_prim_Ave /= CtrTOF_prim_Ave;}
                }

                for(int i = 0; i < NBbins+1; i++){
                    if(Brec > Bedges[i] && Brec < Bedges[i+1]){
                        HTOF_Prim[i]->Fill(TOF_prim_Ave/Ltof);
                    }
                }


				std::sort(Tr_Mean_TKR.begin(), Tr_Mean_TKR.begin()+Tr_Mean_TKR.size());
                double TrM_TKR = 0;
                double CtrTrM_TKR = 0;
                if(Tr_Mean_TKR.size() < NHitsmin){
                    TrM_TKR = 0;
                }else{
                    for(unsigned int isig = 0; isig < double(Tr_Mean_TKR.size())/2; isig++){
              		TrM_TKR += Tr_Mean_TKR.at(isig);
              		CtrTrM_TKR++;
               	}
               	if(CtrTrM_TKR == 0){ TrM_TKR = 0;
               	}else{ TrM_TKR /= CtrTrM_TKR;}
                }

                for(int i = 0; i < NBbins+1; i++){
                    if(Brec > Bedges[i] && Brec < Bedges[i+1]){
                        HTKR_TrM[i]->Fill(TrM_TKR);
                    }
                }

				if(Beta_Proxy.size() >= NHitsmin){
				    std::sort(Beta_Proxy.begin(), Beta_Proxy.end(),std::greater<>());

                        for(unsigned int isig = 0; isig < double(Beta_Proxy.size()); isig++){
                            if(print)cout << "Sorted Beta Proxy " << isig << " is " << Beta_Proxy.at(isig) << endl;
                        }

                        double TrBP = 0;
                        double CtrTrBP = 0;
                        //for(unsigned int isig = 0; isig < floor(double(Beta_Proxy.size())/2); isig++){
                        for(unsigned int isig = 0; isig < double(Beta_Proxy.size())/2; isig++){
                            TrBP += Beta_Proxy.at(isig);
                            CtrTrBP++;
                        }if(CtrTrBP == 0){ TrBP = 0;
                        }else{ TrBP /= CtrTrBP;}

                        //To floor or not to floor!!!
                        //cout << "Beta Prox floor(BP/2) = " << TrBP << endl;
                        if(print)cout << "Beta Prox BP/2 not floor = " << TrBP << endl;
                        if(TrBP > 0){
                            bcounts++;

                            for(int i = 0; i < NBbins+1; i++){
                                if(Brec > Bedges[i] && Brec < Bedges[i+1]){
                                    HBetaProxy[i]->Fill(TrBP);
                                }
                            }

                        }

				}


			} //Closed bracket for if statement for TOF cuts

			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut (beta, cos, pt exists)

	} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i

for(int i = 0; i < NBbins; i++){
    histplot1f(Form("HTKR_TrM%i",i), HTKR_TrM[i], "TKR TrM " + roundstr_d(Bedges[i],2) + "-" + roundstr_d(Bedges[i+1],2), "TrM", "Nentries", out_path + "Bedge_TKR_TrM"+ roundstr_d(Bedges[i],2) + "-" + roundstr_d(Bedges[i+1],2)+ end_name);
    histplot1f(Form("HTOF_PrimM%i",i), HTOF_Prim[i], "TOF Primary " + roundstr_d(Bedges[i],2) + "-" + roundstr_d(Bedges[i+1],2), "TOF Primary dE/dx", "Nentries", out_path + "Bedge_TOF_prim"+ roundstr_d(Bedges[i],2) + "-" + roundstr_d(Bedges[i+1],2)+ end_name);
    histplot1f(Form("HBeta_ProxM%i",i), HBetaProxy[i], "Beta Proxy " + roundstr_d(Bedges[i],2) + "-" + roundstr_d(Bedges[i+1],2), "Bprox", "Nentries", out_path + "Bedge_Bprox"+ roundstr_d(Bedges[i],2) + "-" + roundstr_d(Bedges[i+1],2)+ end_name);
    HTKR_TrM[i]->SaveAs((out_path + to_string(i) + end_name + "HTKR_TrM.root").c_str());
    HTOF_Prim[i]->SaveAs((out_path + to_string(i) + end_name + "HTOF_Prim.root").c_str());
    HBetaProxy[i]->SaveAs((out_path + to_string(i) + end_name + "HBeta_Prox.root").c_str());
}

myfile.open(out_path + txtname,std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile.close();

//Histogram section
//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}
