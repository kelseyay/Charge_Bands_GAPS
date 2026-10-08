//How to use: ./DZedep3 -i /data1/nextcloud/cra_data/data/2025/production/v26.01/reconstructed/flight/251226/starlink251226_17 -s 0 -t 1 -f 1 -k 1 -r 2 -o test/
//Also use: ./DZedep3 -i /home/kelsey/simulations/simdat/flight/251226/26.01/starlink251226_1 -s 0 -t 1 -f 1 -k 1 -r 2 -l 0.9 -u 0.99 -o test/

using namespace std;

#include "KYtools.C"
#include <boost/math/special_functions/lambert_w.hpp>
//Headers included in KYtools are included here, so this is nice.

//FIXME: does this work on mac?
#include <sys/stat.h>

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<string>("out_file", "name of output root file", "", "o");
parser->AddCommandLineOption<double>("beta_low", "low Beta Cut",0.2,"l");
parser->AddCommandLineOption<double>("beta_high", "upper Beta Cut",1.0,"u");
parser->AddCommandLineOption<double>("tkr_factor", "tkr factor",1,"k");
parser->AddCommandLineOption<double>("tof_factor", "tof factor",1,"f");
parser->AddCommandLineOption<bool>("TKR", "tkr use",1,"s");
parser->AddCommandLineOption<bool>("TF", "tof use",1,"t");
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
cout << "reco_path " << reco_path << endl;
cout << "out path " << out_path << endl;
//cout << "out path last string " << out_path[out_path.length()-1] << endl;

if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash!" << endl; out_path = out_path + '/'; }

int MainLoopScaleFactor = parser->GetOption<int>("MainloopScale");

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;

//Prepare textile for saving values
std::ofstream myfile;
string txtname = "RecChargeZedep3.txt";
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

//int MainLoopScaleFactor = 1; //Set this number to scale the step size. Larger means runs faster and fewer events
double TofCutLow = 0; //No low Tof cut right now
double TrackerCut = 0.4; //Threshold for an energy deposition to be considered a hit
int NHitsmin = 1;

double xlow = 0.1; //Low range for histogram MeV
double xhigh = 3; //High range for histogram MeV

double coshigh = 0.54; //0.995; //0.92 //0.54 is the highest angle that can hit UMB, CBEtop, CBEbot
double coslow = 1; //0.62 //0.8
double fitlow = 0.5;
double fithigh = 3.5;
double acutlow = 1.5;
const Int_t NBins = 50;

//Full tracker histogram range
double mpvmin = 0.66;
double mpvmax = 0.75;

bool TKR = parser->GetOption<bool>("TKR");
bool TF = parser->GetOption<bool>("TF");
int TRG = parser->GetOption<int>("TRG");

double tk = parser->GetOption<double>("tkr_factor");
double tf = parser->GetOption<double>("tof_factor");

int passctr = 0;
int strkctr = 0;

int alphactr = 0;
int pctr = 0;
int muctr = 0;
int yesid_actr = 0;
int yesid_pctr = 0;

float Ztof = 5.574;
float Atof = 10.3;
float Ltof = 0.635;
float rtof = 1.032;
float Gtofnew = tf; //1/0.9; //(1/0.71);
//cout << "sqrt(Gtofnew) = " << sqrt(tf) << endl;

float Ztkr = 14;
float Atkr = 28.3;
float Ltkr = 0.22; //Should = 0.25 for the MC simulation!
float rtkr = 2.33;
float Gtkrnew = tk; //1; // (1/0.83);
//cout << "sqrt(Gtkrnew) = " << sqrt(tk) << endl;

//Prepare cuts:
map<int, unsigned int> TofIndexVolumeIdMap;

//For Plotting purposes
ca::GPlottingTools Plotting;
char text[400]; //This variable is used later to name the plots

TH1D * HTruncatedMeanEnergyDepositionMip = Plotting.DefineTH1D("HTruncatedMeanEnergyDepositionMip",100, 0, 10, "sqrt(sqrt(truncated mean E))nergy deposition downgoing MIP [MIP]", "entries", 0.5, 1e4);
TH1D * HChargeMip = Plotting.DefineTH1D("HChargeMip",200, 0.1, 4, "particle charge for downgoing MIP", "entries", 0.5, 1e4);
TH1D * HTofMult = Plotting.DefineTH1D("HTofMult",200, 0.1, 2, "Multiplicative Factor TOF", "entries", 0.5, 1e4);
TH1D * HTkrMult = Plotting.DefineTH1D("HTkrMult",200, 0.1, 2, "Multiplicative Factor TKR", "entries", 0.5, 1e4);

//TH2D * HGenB_vs_GenZ = new TH2D("HGenB_vs_GenZ","Gen_Beta * Gen_Z vs Gen_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 );
TH2D * HRecB_vs_CalcZ = new TH2D("HRecB_vs_CalcZ","Calculated Z vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 4);

TH2D * HRecB_vs_TrM_TOF = new TH2D("HRecB_vs_TrM_TOF","sqrt(TrM_TOF) Z vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 4);
TH2D * HRecB_vs_TrM_TKR = new TH2D("HRecB_vs_TrM_TKR","sqrt(TrM_TKR) Z vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 4);

TH2D * HRecB_vs_dedx_TOF = new TH2D("HRecB_vs_dedx_TOF","dE/dx Z vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 20);
TH2D * HRecB_vs_dedx_TKR = new TH2D("HRecB_vs_dedx_TKR","dE/dx Z vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.1 , 20);



//TH2D * HRecB_vs_Cal = new TH2D("HRecB_vs_RecBTrunM","Rec_Beta * Tr_Mean vs Rec_Beta",50,betacut - 0.1, betahigh + 0.1, 50, 0.5*betacut -0.1 , 1.5*betahigh*2 + 0.1 );

TH1D * hedep; //Full tracker energy deposition plot
hedep = new TH1D ("h0", ("Edep l Beta " + to_string(betacut) + " - " + to_string(betahigh) ).c_str(), NBins, xlow,xhigh);


//How many entries:
cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
cout << "Beta high " << betahigh << endl;
cout << "Beta cut " << betacut << endl;
cout << "MPV range on 2D Full Tracker will be " << mpvmin << " - " << mpvmax << endl;

//Now we can go over the loop
TreeRec->GetEntry(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;

//Using i to loop over every event in the tree
//for(unsigned int i = 0; i < 10; i+=MainLoopScaleFactor){
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor){
    TreeRec->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(MainLoopScaleFactor*10))) == 0){
	    cout << "Event number " << i << endl;
	}

    vector<int> lg_hits = {};
    bool TofTriggerSingleTrack = false;
    LG_std_event_selection_cuts(Event,lg_hits,TofTriggerSingleTrack);

    if(TofTriggerSingleTrack){

        if(std_Rec_event_cuts(Event) &&  fabs(Event->GetPrimaryBeta()) > betacut && fabs(Event->GetPrimaryBeta()) <  betahigh && ( (TRG == 0) || ((int)Event->GetTriggerSources().at(0) == TRG) )){

            int Outer_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
            int Inner_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
            int TKRflag = 0;
            int Layer_Hits_Tracker[7] = {}; //Seven layers

            vector<int> tof_hg_hits = {};
            tof_flags(Event, tof_hg_hits, Outer_TOF_flag, Inner_TOF_flag, Layer_Hits_Tracker, TKRflag );

            if( std::is_permutation(lg_hits.begin(), lg_hits.end(), tof_hg_hits.begin(), tof_hg_hits.end())  ){ //If UMB and CBE top and CBE bot and NOTHING ELSE

                vector<double> Zevent;
                vector<double> Tr_Mean_TOF;
                vector<double> Tr_Mean_TKR;
                double Brec = Event->GetPrimaryBeta();

                for(uint isig=0; isig<Event->GetTrack(0)->GetEnergyDeposition().size(); isig++){
                    unsigned int VolumeId  = Event->GetTrack(0)->GetVolumeId(isig); //Check the VolumeId of the event

                    if(TF && GGeometryObject::IsTofVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TofCutLow){
                        //A hit in the COR or CBE_sides needs to be multiplied by sin(theta) instead of cos(theta)
                        double ion = 0.000016*pow(Ztof,0.9);
                        double eps_noz = 0.15355*(Ztof/Atof)*(rtof*Ltof)*(1/pow(Brec,2));
                        double const_terms = 0.2 - pow(Brec,2) + log( 1.022*pow(Brec,2)/( ion*(1 - pow(Brec,2))) )  ;

                        if(volspec(VolumeId,2,1) == 0 || volspec(VolumeId,2,1) == 1){ //It's a flat paddle
                            //This isn't actually the truncated mean. It's the energy deposition, the truncating happens in the Zevent phase
                            double Ang_Edep = Gtofnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()) ;
                            Zevent.push_back( sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Ang_Edep/ion)*exp(const_terms) ) - const_terms  ) ) );
                            Tr_Mean_TOF.push_back(sqrt(Ang_Edep));
                        } else{ //It's a vertical paddle
                            double Ang_Edep = Gtofnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs( sqrt(1-pow(Event->GetPrimaryMomentumDirection().CosTheta(),2)) ) ;
                            float xy_path_corr = 1;
                            float x_vec = Event->GetTrack(0)->GetMomentumDirection()[0][0];
                            float y_vec = Event->GetTrack(0)->GetMomentumDirection()[0][1];

                            if(volspec(VolumeId,2,1) == 2 ||volspec(VolumeId,2,1) == 3){ //+/- X paddles
                                xy_path_corr = fabs(x_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) );
                            }
                            if(volspec(VolumeId,2,1) == 4 ||volspec(VolumeId,2,1) == 5){ //+/- Y paddles
                                xy_path_corr = fabs(y_vec) / sqrt( pow(x_vec,2) + pow(y_vec,2) );
                            }

                            Ang_Edep = xy_path_corr*Ang_Edep ;
                            Zevent.push_back( sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Ang_Edep/ion)*exp(const_terms) ) - const_terms  ) ) );
                            Tr_Mean_TOF.push_back(sqrt(Ang_Edep));
                        }
                    }

                    if(TKR && GGeometryObject::IsTrackerVolume(VolumeId) && Event->GetTrack(0)->GetEnergyDeposition(isig) > TrackerCut){
                        //A hit in the COR or CBE_sides needs to be multiplied by sin(theta) instead of cos(theta)
                        double ion = 0.000016*pow(Ztkr,0.9);
                        double eps_noz = 0.15355*(Ztkr/Atkr)*(rtkr*Ltkr)*(1/pow(Brec,2));
                        double const_terms = 0.2 - pow(Brec,2) + log( 1.022*pow(Brec,2)/( ion*(1 - pow(Brec,2))) )  ;
                        double Ang_Edep = Gtkrnew*Event->GetTrack(0)->GetEnergyDeposition(isig)*fabs(Event->GetPrimaryMomentumDirection().CosTheta()) ;
                        Zevent.push_back( sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Ang_Edep/ion)*exp(const_terms) ) - const_terms  ) ) );
                        Tr_Mean_TKR.push_back(sqrt(Ang_Edep));
                        //cout << "z from MPV form Tr Mean = " << sqrt( (ion/eps_noz)*exp( boost::math::lambert_w0( (Ang_Edep/ion)*exp(const_terms) ) - const_terms  ) ) << endl;

                    }

                } //Iterate over events for ZEdep

                std::sort(Zevent.begin(), Zevent.begin()+Zevent.size());
                double TrZ = 0;
                double CtrTrZ = 0;

                if(Zevent.size() < NHitsmin){
                    TrZ = 0;
                }else{
                    for(unsigned int isig = 0; isig < double(Zevent.size())/2; isig++){
              		TrZ += Zevent.at(isig);
              		CtrTrZ++;
               	}
               	if(CtrTrZ == 0){ TrZ = 0;
               	}else{ TrZ /= CtrTrZ;}
                }
                HRecB_vs_CalcZ->Fill(Event->GetPrimaryBeta(),TrZ);
                std::sort(Tr_Mean_TOF.begin(), Tr_Mean_TOF.begin()+Tr_Mean_TOF.size());
                double TrM_TOF = 0;
                double CtrTrM_TOF = 0;


                //dE/dx
                //cout << "Event is " << i << " beta is " << Event->GetPrimaryBeta() << " dedx TOF " << endl;
                for(unsigned int isig = 0; isig < double(Tr_Mean_TOF.size()); isig++){
                    //cout << "Hit " << isig << " energy deposition is " << pow(Tr_Mean_TOF.at(isig),2) << endl;
                    HRecB_vs_dedx_TOF->Fill(Event->GetPrimaryBeta(),(pow(Tr_Mean_TOF.at(isig),2)/Ltof));
                }

                for(unsigned int isig = 0; isig < double(Tr_Mean_TKR.size()); isig++){
                    //cout << "Hit " << isig << " energy deposition is " << pow(Tr_Mean_TOF.at(isig),2) << endl;
                    HRecB_vs_dedx_TKR->Fill(Event->GetPrimaryBeta(),(pow(Tr_Mean_TKR.at(isig),2)/Ltkr));
                }

                //Calculated truncated mean
                if(Tr_Mean_TOF.size() < NHitsmin){
                    TrM_TOF = 0;
                }else{
                    for(unsigned int isig = 0; isig < double(Tr_Mean_TOF.size())/2; isig++){
              		TrM_TOF += Tr_Mean_TOF.at(isig);
              		CtrTrM_TOF++;
               	}
               	if(CtrTrM_TOF == 0){ TrM_TOF = 0;
               	}else{ TrM_TOF /= CtrTrM_TOF;}
                }
                HRecB_vs_TrM_TOF->Fill(Event->GetPrimaryBeta(),TrM_TOF);

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
                HRecB_vs_TrM_TKR->Fill(Event->GetPrimaryBeta(),TrM_TKR);

            } //Closed bracket for TOF cuts, using rec data this time.

        } //Closed bracket for standard cuts, using rec data this time.

    } //Closed bracket for LG TOF cuts


}  //Closed bracket for iteration through tree events, move on to the next event i


myfile.open(out_path + txtname,std::ios::app);
myfile << "Total Events/Mainscale Factor " << TreeRec->GetEntries()/MainLoopScaleFactor << endl;
myfile.close();


//Histogram section
//--------------------------------------

string name = "";
if(TF && TKR) name = "Both";
if(TF && !TKR) name = "TOF";
if(!TF && TKR) name = "TKR";

HRecB_vs_CalcZ->SaveAs( (out_path + "HRecB_vs_CalcZ.root").c_str() );
histplot1d("c1", HChargeMip, (name + ": Charge Distribution Rec Beta " + roundstr_d(betacut,2) + " - " + roundstr_d(betahigh,2)  + " TK=" + roundstr_d(tk,2) + " TF=" + roundstr_d(tf,2)).c_str(),"Charge","NEvents", out_path + name + "_Zedep3");
histplot2d("c2", HRecB_vs_CalcZ, "Z_calc versus Rec_B","Reconstructed Beta","Z_calc","NEntries", out_path + name + "Rec2D_Zedep3");
histplot1d("c5", hedep, "Tracker Energy Deposition for "+to_string(betacut)+" - "+to_string(betahigh),"Angle Corrected Energy Deposition","NEvents", out_path + name + "Hedep_Zedep3");

histplot2d("c6", HRecB_vs_TrM_TOF, "Data TOF: Sqrt(TrM_TOF) versus Rec_B","Reconstructed Beta","Sqrt(TrM_TOF)","NEntries", out_path + "BothRecBrecTrM_TOF2D");
histplot2d("c7", HRecB_vs_TrM_TKR, "Data TKR: Sqrt(TrM_TKR) versus Rec_B","Reconstructed Beta","Sqrt(TrM_TKR)","NEntries", out_path + "BothRecBrecTrM_TKR2D");
histplot2d("c7", HRecB_vs_dedx_TOF, "Data TOF: dE/dx versus Rec_B","Reconstructed Beta","dE/dx [MeV/cm]","NEntries", out_path + "BothRecBrecdedx_TOF2D");
histplot2d("c8", HRecB_vs_dedx_TKR, "Data TKR: dE/dx versus Rec_B","Reconstructed Beta","dE/dx [MeV/cm]","NEntries", out_path + "BothRecBrecdedx_TKR2D");


cout << "Hedep Max Bin Center = " << hedep->GetBinCenter(hedep->GetMaximumBin()) << endl;


//histplot2d("c6", HTrunM_vs_RecB, "sqrt(Ang_Edep) versus Rec_B","Reconstructed Beta","sqrt(truncated mean E)","NEntries", out_path + "BothRecBTrunM");
//cout << "Max bin is " << HChargeMip->GetMaximumBin() << /*" factor should be " << 1/(HChargeMip->GetMaximumBin()) << */ endl;
//cout << "Max bin NEntries is " << HChargeMip->GetBinContent(HChargeMip->GetMaximumBin()) << endl;
cout << "TOF: " << TF << " TKR: " << TKR << endl;
cout << "Max bin Center is  " << HChargeMip->GetBinCenter(HChargeMip->GetMaximumBin()) << " factor should be " << pow(1/HChargeMip->GetBinCenter(HChargeMip->GetMaximumBin()),2) << endl;

cout << endl << "I am done" << endl;

return 1;

}
