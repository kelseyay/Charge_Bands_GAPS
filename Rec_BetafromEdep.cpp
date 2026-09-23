//Hopefully this script will provide a primary beta from energy depositions in Not That Much Time.
//This will start simple with beta from at least three tracker hits with energy deposition > 0.4
//Alright so I named this sometihng stupid. We have Edep vs beta and this is Edep FROM beta...


// How to use:
// /home/kelsey/GAPS_TrackerEnergyDep/build/RecEdepBeta -i /home/kelsey/simulations/simdat/proton/v3.0.0/triggerlevel1/ -l 0.2 -u 1.5
// /home/kelsey/GAPS_TrackerEnergyDep/build/RecEdepBeta -i /home/kelsey/simulations/simdat/flight/251221/FPSI/ -l 0.2 -u 1.5 -r 2

// How to use with factors:
// /home/kelsey/GAPS_TrackerEnergyDep/build/RecEdepBeta -i /home/kelsey/simulations/simdat/flight/251221/FPSI/ -l 0.2 -u 1.2 -r 2 -f 1.4 -k 1.15 -o ZCuts_Flight_Data/

using namespace std;

#include "KYtools.C"
#include <boost/math/special_functions/lambert_w.hpp>
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
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
float Ltkr = 0.22; //For instrument tracker, set Ltkr = 0.22

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

float solvermin_tkr = ZOne_TKR(0.95); //0.787;
float solvermin_tof = ZOne_TOF(0.95);

//Try figuring out why this worked: https://root.cern/doc/v620/classTF1.html
TF1 *Z1_tkr_Solve = new TF1("Z1_tkr_Solve", [](double *x, double *p){ return ZOne_TKR(x[0]); }, 0.01, 0.95, 0);
TF1 *Z1_tof_Solve = new TF1("Z1_tof_Solve", [](double *x, double *p){ return ZOne_TOF(x[0]); }, 0.01, 0.95, 0);

bool charge_cut_z1(double Proxy_Beta, double Reco_Beta){
    return true;
}

/*
bool charge_cut_z1(double Proxy_Beta, double Reco_Beta){
    if(Reco_Beta > 0.75 && Proxy_Beta > 0.66){ return 1; }else if(Reco_Beta < 0.75 && Proxy_Beta >  0.83*Reco_Beta + 0.0328){return 1;}else{return 0;}
}*/

bool charge_cut_z2(double Proxy_Beta, double Reco_Beta){
    if(Reco_Beta > 0.75 && Proxy_Beta < 0.55){ return 1; }else if(Reco_Beta < 0.75 && Proxy_Beta <  0.78*Reco_Beta - 0.033){return 1;}else{return 0;}
}

bool std_Rec_event_cuts(const CEventRec* Event){
    bool PassAllCuts = false;
    if(Event->GetNTracks() == 1){  //First select the single track event
    	CTrackRec* pt = Event->GetPrimaryTrack();
    	uint pt_index = 0;
   	    for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

    	//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
    	if( pt != nullptr && pt->GetChi2()/pt->GetNdof() < 3.2 && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0){
            PassAllCuts = true;
    		//-----------EVENT LEVEL CUT APPLIED
        } //End selecton cuts

    } //Single Track cut
     return PassAllCuts;
}

//Example:
//TGraph * g = new TGraph(npointx, xvec, yvec);
//TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*g->Eval(x[0]); }, xmin, xmax, 1);

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<string>("end_name", "name at the end of pngs","","e");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.8,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1,"u");
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


string end_name = parser->GetOption<string>("end_name");
bool MC_Weight = parser->GetOption<bool>("MC_Weighting");
double betacut = parser->GetOption<double>("beta_low");
if(betacut <= 0 || betacut >=1){ betacut = 0.8; cout << "Error with low beta choice. Setting Beta low to 0.8" << endl; }
cout << "beta cut = " << betacut << endl;

double betahigh = parser->GetOption<double>("beta_high");
if(betahigh <= 0 || betahigh >=2 || betahigh < betacut){ betahigh = 1; cout << "Error with high/low beta choice! Setting Beta upper to 1" << endl; }
cout << "beta high = " << betahigh << endl;

string reco_path = parser->GetOption<string>("in_path");
string out_path = parser->GetOption<string>("out_file");
cout << "reco_path " << reco_path << endl;
cout << "out path " << out_path << endl;
//cout << "out path last string " << out_path[out_path.length()-1] << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

int MainLoopScaleFactor = parser->GetOption<int>("MainloopScale");

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

string txtname = "Rec_Beta_from_Edep.txt";

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

//auto tkr_edep_beta_z1 = new TF1("tkr_edep_beta_z1","OneZ_TKR(x)",0.01,0.95);

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0.1; //No low Tof cut right now
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
int NHitsmin = 3;

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

//For this plot, only want Beta Proxy and Beta reconstruction
TH1F * HBetaProxy = new TH1F("HBetaProxy","HBetaProxy",40, betacut, betahigh);
TH1F * HBetaRec = new TH1F("HBetaRec","HBetaRec",40, betacut, betahigh);
TH2D * HRecB_vs_ProxB = new TH2D("HRecB_vs_ProxB","Rec_Beta vs Prox_B",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1);

//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 10000; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
		    cout << "Event number " << i << endl;
	}

    if(LG_std_event_selection_cuts(Event)){
        //cout << "This event passes LG cuts! " << i << endl;

        if(std_Rec_event_cuts(Event) &&  fabs(Event->GetPrimaryBeta()) > betacut && fabs(Event->GetPrimaryBeta()) <  betahigh && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) )){
            //cout << "Event passes the standard cuts! " << endl;

            bool Outer_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
     	    bool Inner_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
  	        bool TKRflag = 0;
  	        int Layer_Hits_Tracker[7] = {}; //Seven layers

    		for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
                if(volspec(VolumeId,0,2) == 20 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                    int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                    Layer_Hits_Tracker[layer]++;
                    TKRflag = 1;
                }
                if(volspec(VolumeId,0,2) == 10 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){ Outer_TOF_flag[volspec(VolumeId,2,1)] = 1; }
                if(volspec(VolumeId,0,2) == 11 && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){ Inner_TOF_flag[volspec(VolumeId,2,1)] = 1; }
    		}

            if( (Outer_TOF_flag[0] && Inner_TOF_flag[0])  && Inner_TOF_flag[1] && !Inner_TOF_flag[2] && !Inner_TOF_flag[3] && !Inner_TOF_flag[4] && !Inner_TOF_flag[5] && !Outer_TOF_flag[2] && !Outer_TOF_flag[3] && !Outer_TOF_flag[4] && !Outer_TOF_flag[5]  ){ //If UMB and CBE top and CBE bot and NOTHING ELSE
                    //cout << "Event also passes the TOF cuts!! " << endl;
                    if(print) cout << "This event also passes HG cuts! " << endl;
                    if(print) cout << endl << "Event is " << i << endl;
    				if(print) cout << "Breco is " << Event->GetPrimaryBeta() << endl;
    			    vector<double> Beta_Proxy;

    				for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
    				    unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event
    					float Edep = Event->GetTrack(0)->GetEnergyDeposition(isig);
    					float Ang_Edep = 0;

    					if(TF && GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){
    					    if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){ //It's a flat paddle.
    					        Ang_Edep = tf*Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta());
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

                                Ang_Edep = tf*Edep*fabs( sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2)))*xy_path_corr;
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
    				        if(sdstrip != 0 && sdstrip != 7 && Layer_Hits_Tracker[layer] < 2){
    			                Ang_Edep = tk*Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta());
    							if(Ang_Edep > solvermin_tkr){
    							    double root = Z1_tkr_Solve->GetX(Ang_Edep);
    								if(print)cout << "Hit is " << isig << endl;
    								if(print)cout << "Energy Deposition is " << Edep << " angle corrected is " << Ang_Edep << endl;
                                    if(print)cout << "Calculated Beta at Volid " << VolumeId << " = "  << root << endl;
                                    Beta_Proxy.push_back(root);
    							}
    						}
    				    }
    				} //Closed bracket for iteration over isig with TOF cuts

    				if(Beta_Proxy.size() > 1){
    				    std::sort(Beta_Proxy.begin(), Beta_Proxy.end(),std::greater<>());

                            //for(unsigned int isig = 0; isig < double(Beta_Proxy.size()); isig++){
                            //    if(print)cout << "Sorted Beta Proxy " << isig << " is " << Beta_Proxy.at(isig) << endl;
                            //}

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
                            if(TrBP > 0 && charge_cut_z1(TrBP,Event->GetPrimaryBeta()) ){
                                HBetaProxy->Fill(TrBP);
                                HBetaRec->Fill(Event->GetPrimaryBeta());
                                HRecB_vs_ProxB->Fill(Event->GetPrimaryBeta(),TrBP);
                            }

    				} //Closed bracket Beta Proxy calculation

                } //TOF Cuts

    			} //Closed bracket for if statement for standard event selection cuts

    			//-----------EVENT LEVEL CUTS END

    }//Closed bracket for LG cuts

}  //Closed bracket for iteration through tree events, move on to the next event i

HBetaRec->SetMaximum(HBetaRec->GetEntries());
HBetaProxy->SetMaximum(HBetaProxy->GetEntries());

//HBetaProxy->Scale( 1./HBetaProxy->Integral(),"WIDTH");

histplot1f("c1",HBetaProxy,"Proxy Beta","Proxy Beta","NEntries", out_path + "Rec_BetaProxy"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR) + "TF_Factor" + roundstr_d(tf,2) + "TK_Factor" + roundstr_d(tk,2) + end_name );
histplot1f("c2",HBetaRec,"Reconstructed Beta","Reconstructed Beta","NEntries", out_path + "Rec_BetaRec " + "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2) + "TOF" + to_string(TF) + "TKR" + to_string(TKR) + "TF_Factor" + roundstr_d(tf,2) + "TK_Factor" + roundstr_d(tk,2) + end_name );
histplot2d("c2_5",HRecB_vs_ProxB,"Prox_B versus Rec_B","Reconstructed Beta", "Proxy Beta","NEntries", out_path + "RecBProxB" + "B" +  roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR) + "TF_Factor" + roundstr_d(tf,2) + "TK_Factor" + roundstr_d(tk,2) + end_name );


HBetaRec->SaveAs((out_path + "HRecData.root").c_str());
HBetaProxy->SaveAs((out_path + "HProxyData.root").c_str());


myfile.open(out_path + txtname,std::ios::app);
myfile << TString::Format( "Trigger : %i", TRG) << endl;
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile << "TOF Factor " << tf << endl;
myfile << "TKR Factor " << tk << endl;
myfile.close();


//Histogram section
//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}
