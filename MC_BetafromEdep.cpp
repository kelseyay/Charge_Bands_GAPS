//Hopefully this script will provide a primary beta from energy depositions in Not That Much Time.
//This will start simple with beta from at least three tracker hits with energy deposition > 0.4
// HEY MESSING WITH THIS! SINGLE TRACK CRITERIA IS BEING FUTZED ABOUT

// How to use:
// /home/kelsey/GAPS_TrackerEnergyDep/build/MCEdepBeta -i /home/kelsey/simulations/simdat/proton/v3.0.0/triggerlevel1/ -l 0.4 -u 0.8


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

TH1D * HBetaProxy = Plotting.DefineTH1D("HBetaProxy",40, betacut, betahigh, "Beta Proxy", "entries", 10, 100000);
TH1D * HBetaRec = Plotting.DefineTH1D("HBetaRec",40, betacut, betahigh, "Beta Reconstructed", "entries", 10, 100000);
TH1D * HBetaGen = Plotting.DefineTH1D("HBetaGen",40, betacut, betahigh, "Generated Beta", "entries", 10, 100000);
TH2D * HRecB_vs_GenB = new TH2D("HRecB_vs_GenB","Rec_Beta vs Gen_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1);
TH2D * HProxB_vs_GenB = new TH2D("HProxB_vs_GenB","Prox_B vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1, 50, 0.1, 1);
TH2D * HRecB_vs_ProxB = new TH2D("HRecB_vs_ProxB","Rec_Beta vs Prox_B",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1);

TH1F * HBetaGen_Weight_noTOFcuts = new TH1F("HBetaGen_Weight_noTOFcuts","HBetaGen_Weight_noTOFcuts",40, betacut, betahigh);
TH1F * HBetaGen_Weight = new TH1F("HBetaGen_Weight","HBetaGen_Weight",40, betacut, betahigh);
TH1F * HBetaProxy_Weight = new TH1F("HBetaProxy_Weight","HBetaProx_Weight",40, betacut, betahigh);
TH1F * HBetaRec_Weight = new TH1F("HBetaRec_Weight","HBetaRec_Weight",40, betacut, betahigh);
float bcounts_Weight = 0;

//TH2D * HRecB_vs_GenB_Weight = new TH2D("HRecB_vs_GenB_Weight","Rec_Beta * Tr_Mean vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 1);
//TH2D * HProxB_vs_GenB_Weight = new TH2D("HProxB_vs_GenB_Weight","Prox_B vs Gen_Beta",50, betacut - 0.1, betahigh + 0.1, 50, 0.1, 1);


//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;

//MC Weighting things below:
TChain*TreeSimulationParameter = new TChain("SimulationParameterTree");
TreeSimulationParameter->Add(FilenameRoot);
GSimulationParameter * Parameter = new GSimulationParameter;
TreeSimulationParameter->SetBranchAddress("SimulationParameter", &Parameter);
TreeSimulationParameter->GetEntry(0);

double FluxScaleFactor = 71.1552/32.058750*0.25;

int BetaBins = 25;
double StartingPlaneAcceptance = 1;
TH1D* HPrimaryBeta = nullptr;
double BinWidthFactor = 1;
std::vector<double> PrimaryBetaLowHigh;

vector<pair<double, double> > CosZenithCut;
CosZenithCut.push_back(make_pair(-0.75, -1));
CosZenithCut.push_back(make_pair(-0.5, -0.75));
CosZenithCut.push_back(make_pair(-0.25, -0.5));
CosZenithCut.push_back(make_pair(0, -0.25));

vector<TGraph*> GMuonTotalFluxUnscaled;
GMuonTotalFluxUnscaled = MC_weighting_framework_muon(GMuonTotalFluxUnscaled);

vector<TGraph*> GProtonTotalFluxUnscaled;
GProtonTotalFluxUnscaled = MC_weighting_framework_proton(GProtonTotalFluxUnscaled);

vector<TGraph*> GAlphaTotalFluxUnscaled;
GAlphaTotalFluxUnscaled = MC_weighting_framework_alpha(GAlphaTotalFluxUnscaled);

CAnalysisManager AnalysisManagerRec;
AnalysisManagerRec.SetGSimulationParameterTChain(TreeSimulationParameter);
StartingPlaneAcceptance = AnalysisManagerRec.GetStartingPlaneAcceptance();
PrimaryBetaLowHigh = AnalysisManagerRec.GetPrimaryBetaLowHigh();
HPrimaryBeta = AnalysisManagerRec.GetHPrimaryBeta();
BinWidthFactor = (PrimaryBetaLowHigh.at(1)-PrimaryBetaLowHigh.at(0))/double(BetaBins) / HPrimaryBeta->GetBinWidth(1);

//End MC weighting things

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 10; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);
    TreeMC->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
		    cout << "Event number " << i << endl;
	}

	//Cuts are implemented in this chunk:
	//if(Event->GetNTracks() == 1){  //First select the single track event
		bool Umbflag = 0;
		bool CBEtopflag = 0;
		bool CBEbotflag = 0;
		bool TKRflag = 0;
		int Layer_Hits_Tracker[7] = {}; //Seven layers

		CTrackRec* pt = Event->GetPrimaryTrack();
		uint pt_index = 0;
      	        for( ; pt_index < Event->GetNTracks(); pt_index++) if( Event->GetTrack(pt_index)->IsPrimary() ) break;

		//Note downwards beta enforced by Event->GetPrimaryBeta() (should be positive) multiplied by Event->GetPrimaryMomentumDirection()[2] (z trajectory of particle)
		if( pt != nullptr && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) ) && pt->GetChi2()/pt->GetNdof() < 3.2 && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) > -coslow && -fabs(Event->GetPrimaryMomentumDirection().CosTheta()) < -coshigh && Event->GetPrimaryBeta()*Event->GetPrimaryMomentumDirection()[2] < 0 && fabs(MCEvent->GetPrimaryBeta()) > betacut && fabs(MCEvent->GetPrimaryBeta()) <  betahigh ){
			//-----------EVENT LEVEL CUT APPLIED

			//cout << "Event: " << i << endl;
			//First iteration over event for flags and EnergyDepositionMip vector filling
			double Bgen = MCEvent->GetPrimaryBeta();
			//double gamma = sqrt( 1/(1 - pow(Bgen,2)) );
			//cout << "Bgen " << Bgen << endl;

			//cout << "Particle spec is " << MCEvent->GetTrack(0)->GetPdg() << endl;
            double RateScale = 0;
            if( MC_Weight ){
                //cout << "Particle species track 0 " << MCEvent->GetTrack(0)->GetPdg() << endl;
                if(MCEvent->GetTrack(0)->GetPdg() == 13) {RateScale = RateScale_muon(TreeSimulationParameter, MainLoopScaleFactor, StartingPlaneAcceptance, BinWidthFactor, HPrimaryBeta, CosZenithCut, GMuonTotalFluxUnscaled, FluxScaleFactor, Event); }
                if(MCEvent->GetTrack(0)->GetPdg() == 2212) { FluxScaleFactor = 1;  RateScale = RateScale_TOA(TreeSimulationParameter, MainLoopScaleFactor, StartingPlaneAcceptance, BinWidthFactor, HPrimaryBeta, CosZenithCut, GProtonTotalFluxUnscaled, FluxScaleFactor, Event);}
                if(MCEvent->GetTrack(0)->GetPdg() == 1000020040) {FluxScaleFactor = 1;  RateScale = RateScale_TOA(TreeSimulationParameter, MainLoopScaleFactor, StartingPlaneAcceptance, BinWidthFactor, HPrimaryBeta, CosZenithCut, GAlphaTotalFluxUnscaled, FluxScaleFactor, Event);}
                //cout << "BGen = " << Bgen << " Cos(theta) = " << Event->GetPrimaryMomentumDirectionGenerated().CosTheta() << " RateScale = " << RateScale << endl;
            }
            HBetaGen_Weight_noTOFcuts->Fill(MCEvent->GetPrimaryBeta(),RateScale);

			for(uint isig=0; isig<Event->GetPrimaryTrack()->GetEnergyDeposition().size(); isig++){
                unsigned int VolumeId  = Event->GetPrimaryTrack()->GetVolumeId(isig); //Check the VolumeId of the event
                if(volspec(VolumeId,0,2) == 20 && Event->GetPrimaryTrack()->GetEnergyDeposition(isig) > TrackerCut){
                    int layer = GGeometryObject::GetTrackerLayer(VolumeId);
                    Layer_Hits_Tracker[layer]++;
                    TKRflag = 1;
                }
                if(volspec(VolumeId,0,3) == 100 && Event->GetPrimaryTrack()->GetEnergyDeposition(isig) > TofCutLow){ Umbflag = 1; } // cout << "UMB hit!" <<endl ;
                if(volspec(VolumeId,0,3) == 110 && Event->GetPrimaryTrack()->GetEnergyDeposition(isig) > TofCutLow) { CBEtopflag = 1; }// cout << "CBE top hit!" << endl;
                if(volspec(VolumeId,0,3) == 111 && Event->GetPrimaryTrack()->GetEnergyDeposition(isig) > TofCutLow) { CBEbotflag = 1; }// cout << "CBE bot hit!" << endl;
			}

			//Do we want to only run this on certain tracks? Yeah probably. Can remove the TOF flags. Just run on whatever lol.
			if(Umbflag && CBEtopflag && (CBEbotflag || TKRflag)){ //Guarantees 3 hit requirement, if the energy deposition is lower than the solver minimum, don't use. Still try to calculate.
			    if(print) cout << endl << "Event is " << i << endl;
				if(print) cout << "Zenith angle is " << Event->GetPrimaryMomentumDirection().CosTheta() << endl;
			    if(print) cout << "Bgen is " << Bgen << endl;
				if(print) cout << "Breco is " << Event->GetPrimaryBeta() << endl;
			    vector<double> Beta_Proxy;

				for(uint isig=0; isig<Event->GetPrimaryTrack()->GetEnergyDeposition().size(); isig++){
				    unsigned int VolumeId  = Event->GetPrimaryTrack()->GetVolumeId(isig); //Check the VolumeId of the event
					float Edep = Event->GetPrimaryTrack()->GetEnergyDeposition(isig);
					float Ang_Edep = 0;

					if(TF && GGeometryObject::IsTofVolume(VolumeId) && Event->GetPrimaryTrack()->GetEnergyDeposition(isig) > TofCutLow){
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
                            float x_vec = Event->GetPrimaryTrack()->GetMomentumDirection()[0][0];
                            float y_vec = Event->GetPrimaryTrack()->GetMomentumDirection()[0][1];
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

					if(TKR && volspec(VolumeId,0,2) == 20 && Event->GetPrimaryTrack()->GetEnergyDeposition(isig) > TrackerCut){
						int layer = GGeometryObject::GetTrackerLayer(VolumeId);
						int sdstrip = GGeometryObject::GetDetectorStrip(VolumeId);
				        if(sdstrip != 0 && sdstrip != 7 && Layer_Hits_Tracker[layer] < 2){
			                Ang_Edep = Edep*fabs(Event->GetPrimaryMomentumDirection().CosTheta());
							if(Ang_Edep > solvermin_tkr){
							    double root = Z1_tkr_Solve->GetX(Ang_Edep);
								if(print)cout << "Hit is " << isig << endl;
								if(print)cout << "Energy Deposition at Volid " << VolumeId << " is " << Edep << " angle corrected is " << Ang_Edep << endl;
                                if(print)cout << "Calculated Beta at Volid = "  << root << endl;
                                Beta_Proxy.push_back(root);
							}
						}
				    }
				} //Closed bracket for iteration over isig with TOF cuts

				if(Beta_Proxy.size() > 2){
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
                        if(TrBP > 0){
                            HBetaProxy->Fill(TrBP);
                            HBetaGen->Fill(Event->GetPrimaryBetaGenerated());
                            HBetaRec->Fill(Event->GetPrimaryBeta());
                            HRecB_vs_GenB->Fill(Event->GetPrimaryBetaGenerated(),Event->GetPrimaryBeta());
                            HProxB_vs_GenB->Fill(Event->GetPrimaryBetaGenerated(),TrBP);
                            HRecB_vs_ProxB->Fill(Event->GetPrimaryBeta(),TrBP);
                            /*if(Event->GetNTracks() > 1){
                                cout << "Event " << i << " Multi-track event with good BP " << endl;
                                cout << "Beta proxy = " << TrBP << endl;
                                cout << " Brec = " << Event->GetPrimaryBeta() << endl;
                                cout << "Bgen = " << Event->GetPrimaryBetaGenerated() << endl;
                            }*/

                            //Weighted MC filling below

                            if(MC_Weight){
                                bcounts_Weight = bcounts_Weight+RateScale;
                                HBetaGen_Weight->Fill(Event->GetPrimaryBetaGenerated(),RateScale);
                                HBetaProxy_Weight->Fill(TrBP,RateScale);
                                HBetaRec_Weight->Fill(Event->GetPrimaryBeta(),RateScale);
                            }

                        }

				}


			} //Closed bracket for if statement for TOF cuts

			//-----------EVENT LEVEL CUTS END


		} //Closed bracket for event level cut (beta, cos, pt exists)

	//} //Closed bracket for single track cut

}  //Closed bracket for iteration through tree events, move on to the next event i

histplot2d("c1",HRecB_vs_GenB,"Rec_B versus Gen_B","Generated Beta", "Reconstructed Beta","NEntries", out_path + "GenBRec" + "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2) + "TOF" + to_string(TF) + "TKR" + to_string(TKR) );
histplot2d("c2",HProxB_vs_GenB,"Prox_B versus Gen_B","Generated Beta", "Proxy Beta","NEntries", out_path + "GenBProxB" + "B" +  roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );
histplot2d("c2_5",HRecB_vs_ProxB,"Prox_B versus Rec_B","Reconstructed Beta", "Proxy Beta","NEntries", out_path + "RecBProxB" + "B" +  roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );


HBetaRec->SetMaximum(HBetaRec->GetEntries());
HBetaGen->SetMaximum(HBetaGen->GetEntries());
HBetaProxy->SetMaximum(HBetaProxy->GetEntries());
histplot1d("c3",HBetaProxy,"Proxy Beta","Proxy Beta","NEntries", out_path + "BetaProxy"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );
histplot1d("c4",HBetaGen,"Generated Beta","Generated Beta","NEntries", out_path + "BetaGen"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );
histplot1d("c5",HBetaRec,"Reconstructed Beta","Reconstructed Beta","NEntries", out_path + "BetaRec " + "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );

//HBetaProxy_Weight->Scale( 1./HBetaProxy_Weight->Integral(),"WIDTH");
HBetaProxy_Weight->SetMinimum(1);
HBetaProxy_Weight->SetMaximum(bcounts_Weight);

HBetaGen_Weight->SetMinimum(1);
HBetaGen_Weight->SetMaximum(bcounts_Weight);

HBetaRec_Weight->SetMaximum(bcounts_Weight);
HBetaRec_Weight->SetMinimum(1);

HBetaGen_Weight_noTOFcuts->SetMaximum(bcounts_Weight*100);
HBetaGen_Weight_noTOFcuts->SetMinimum(1);

histplot1f("c6",HBetaGen_Weight,"Generated Beta Weighted","Generated Beta Weighted","NEntries", out_path + "Weighted_BetaGen"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );
histplot1f("c7",HBetaProxy_Weight,"Proxy Beta Weighted","Proxy Beta Weighted","NEntries", out_path + "Weighted_BetaProx"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );
histplot1f("c8",HBetaRec_Weight,"Reconstructed Beta Weighted","Reconstructed Beta Weighted","NEntries", out_path + "Weighted_BetaRec"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );
histplot1f("c9",HBetaGen_Weight_noTOFcuts,"Generated Beta Weighted (No TOF Cuts)","Generated Beta Weighted","NEntries", out_path + "Weighted_BetaGen_NoTOFCuts"+ "B" + roundstr_d(betacut,2) + "-" + roundstr_d(betahigh,2)+ "TOF" + to_string(TF) + "TKR" + to_string(TKR)  );


HBetaGen_Weight->SaveAs((out_path + "HGenWeight.root").c_str());
HBetaProxy_Weight->SaveAs((out_path + "HProxyWeight.root").c_str());
HBetaRec_Weight->SaveAs((out_path + "HRecWeight.root").c_str());

myfile.open(out_path + txtname,std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile.close();


//Histogram section
//--------------------------------------

cout << endl << "I am done" << endl;

return 1;

}



/*
//Theoretically save on compute time if put weighting stuff at last possible place.
double AcceptanceScale;
if (TreeSimulationParameter != nullptr){
    //acceptance scaling factor based on beta of the primary
    AcceptanceScale = MainLoopScaleFactor*StartingPlaneAcceptance/(BinWidthFactor*HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(Event->GetPrimaryBetaGenerated())));

    if(HPrimaryBeta->GetBinContent(HPrimaryBeta->FindBin(Event->GetPrimaryBetaGenerated())) == 0) AcceptanceScale = 0;
    //Find the bin associated with the generated beta in the vector, see if it doesn't exist?
}else AcceptanceScale = 1; //Set it = to 1 if there's a nullptr? That's a surprise...
//cout << "AcceptanceScale = " << AcceptanceScale << endl;

int AngularRegion = -1;
for(unsigned int a = 0; a < CosZenithCut.size(); a++) if(Event->GetPrimaryMomentumDirectionGenerated().CosTheta() < CosZenithCut.at(a).first && Event->GetPrimaryMomentumDirectionGenerated().CosTheta() > CosZenithCut.at(a).second) AngularRegion = a;
//This is just checking which "bin" the generated cos(theta) is in
if(AngularRegion < 0) continue; //Don't bother if angular region wasn't found
RateScale = FluxScaleFactor*AcceptanceScale*GMuonTotalFluxUnscaled.at(AngularRegion)->Eval(Event->GetPrimaryBetaGenerated());
//cout << "Beta? = " << Event->GetPrimaryBetaGenerated() << endl;
//cout << "Angle? = " << Event->GetPrimaryMomentumDirectionGenerated().CosTheta() << endl;
cout << "RateScale?? = " << RateScale << endl;
cout << "RateScale from KYTools = " << RateScale_muon(TreeSimulationParameter, MainLoopScaleFactor, StartingPlaneAcceptance, BinWidthFactor, HPrimaryBeta, CosZenithCut, GMuonTotalFluxUnscaled, FluxScaleFactor, Event) << endl;
*/
