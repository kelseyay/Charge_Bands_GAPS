//Use this code to pare down an entire large dataset to a small, small dataset!!

using namespace std;

#include "KYtools.C"

#include "GDataEvent.hh"
#include "GGeometryTools.hh"
#include "GDataPoint.hh"
#include "GDataTrack.hh"
#include "GDataVertex.hh"

using namespace Crane::Analysis;
namespace ca = Crane::Analysis;
namespace cl = Crane::Common;
//using Crane::Calibration;

bool is_selected(const CEventRec* Event){

    vector<int> lg_hits = {};
    bool TofTriggerSingleTrack = false;
    bool Selected_Slim = false;
    LG_std_event_selection_cuts(Event,lg_hits,TofTriggerSingleTrack);

    if(TofTriggerSingleTrack){
        //cout << "Event ID is " << Event->GetEventNumber() << endl;
        //cout << "This event passes LG cuts! " << endl;

        if(std_Rec_event_cuts(Event) &&  fabs(Event->GetPrimaryBeta()) > 0.2 && fabs(Event->GetPrimaryBeta()) <  1.2 && (((int)Event->GetTriggerSources().at(0) == 2) )){
            //cout << "Event passes the standard cuts! " << endl;

            int Outer_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
       	    int Inner_TOF_flag[7] = {}; //TOF top = 0, bot = 1, +X = 2, -X = 3, +Y = 4, -Y = 5. 1 or 3 PPs is 6
            int TKRflag = 0;
            int Layer_Hits_Tracker[7] = {}; //Seven layers
            vector<int> tof_hg_hits = {};

            tof_flags(Event, tof_hg_hits, Outer_TOF_flag, Inner_TOF_flag, Layer_Hits_Tracker, TKRflag );

            /*
            cout << "Passed HG hits: " << endl;
            for (int element : tof_hg_hits) {
                std::cout << element << " ";
            }

            cout << endl << "How about LG hits: " << endl;
            for (int element : lg_hits) {
                std::cout << element << " ";
            }

            cout << endl;*/

            if (std::is_permutation(lg_hits.begin(), lg_hits.end(), tof_hg_hits.begin(), tof_hg_hits.end())) {
                //cout << "WOW the LG hits and HG hits all match!! " << endl;
                Selected_Slim = true;
            }

        } // Closed bracket standard event selection cuts

    } //Closed bracket TOF trigger cuts

    return Selected_Slim;
}

int main(int argc, char *argv[]){

GOptionParser* parser = GOptionParser::GetInstance();
parser->AddProgramDescription("Minimal Reproducable Example for Extracing Data from Reco Data");
parser->AddCommandLineOption<string>("in_path", "path to instrument data files", "./*", "i");
parser->AddCommandLineOption<bool>("save", "save the special events to a root file?",1,"s");
parser->AddCommandLineOption<string>("out_file", "name of output path", "", "o");
parser->AddCommandLineOption<string>("sts_root_name", "name of output root file", "sts_root_test.root", "n");
parser->ParseCommandLine(argc, argv);
parser->Parse();

bool SAVE = parser->GetOption<bool>("save");

string out_path = parser->GetOption<string>("out_file");
cout << "out path: " << out_path << endl;
if(out_path != "" && out_path[out_path.length()-1] != '/' ){ cout <<  "out path no slash! Adding! " << endl; out_path = out_path + '/'; }

string reco_path = parser->GetOption<string>("in_path");

string sts_file_name = parser->GetOption<string>("sts_root_name");
if(sts_file_name.compare(sts_file_name.length()-5,sts_file_name.length(),".root") != 0){ cout << "NO .root at the end of the root name! Adding!" << endl; sts_file_name = sts_file_name + ".root"; }

cout << reco_path << endl;
if(reco_path.compare(reco_path.length()-5,reco_path.length(),".root") == 0){ cout << ".root at the end of the reco path! Deleting!" << endl; reco_path = reco_path.substr(0,reco_path.length()-5); }

char FilenameRoot[400];
sprintf(FilenameRoot,"%s*.root",reco_path.c_str());
cout << FilenameRoot << endl;


int ev_track = 0;

TFile *f;
TFile *f_source;
TDirectoryFile *GOptions_copy;

CEventRec* Event = new CEventRec(); //New reconstructed event
TChain * TreeRec = new TChain("TreeRec"); //New TreeRec Tchain object (this is new to me)
TreeRec->SetBranchAddress("Rec", &Event); //Set the branch address using Event (defined above)
TreeRec->Add(FilenameRoot);

//Prepare FPSI Reconstruction variable for locating vertex:
Crane::Reconstruction::TrackFit::GDataEvent * reco_data_event_ = new Crane::Reconstruction::TrackFit::GDataEvent();
TChain * TreeGReco = new TChain("TreeGReco");
TreeGReco->SetBranchAddress("FindPrimaryStarIterative", &reco_data_event_); //Set the branch address using Event (defined above)
TreeGReco->Add(FilenameRoot);

CEventMc* MCEvent = new CEventMc(); //New reconstructed event
TChain * TreeMC = new TChain("TreeMc"); //New TreeMC Tchain object (this is new to me)
TreeMC->SetBranchAddress("Mc", &MCEvent); //Set the branch address using Event (defined above)
TreeMC->Add(FilenameRoot);

/*
TChain* events = new TChain("TreeRec");
CEventRec* reco_event = new CEventRec;
events->SetBranchAddress("Rec", &reco_event);
events->Add(FilenameRoot);*/

TTree *Copy_GRecoTree = new TTree("TreeGReco", "GReco Tree");
TTree *Copy_RecTree = new TTree("TreeRec", "Rec Tree");
//TTree *Copy_MCTree = new TTree("TreeMc", "MC Tree");
Copy_GRecoTree = TreeGReco->CloneTree(0);
Copy_RecTree = TreeRec->CloneTree(0);
//Copy_MCTree = TreeMC->CloneTree(0);

cout << "Total Number of events / Mainscale Factor = " << TreeRec->GetEntries() << endl;


for(uint i=0; i<TreeRec->GetEntries(); i++){
    TreeRec->GetEntry(i);
    //TreeGReco->GetEntry(i);
    //TreeMC->GetEntry(i);

    if( ((int)i % (int)ceil(TreeRec->GetEntries()/(10))) == 0){
		    cout << "Event number " << i << endl;
	}

        if(is_selected(Event)){
            /*
            cout << "Event " << i << " is unbelieveably rad " << endl;
            cout << "Event ID? " << Event->GetEventId() << endl;
            cout << "Event Number? " << Event->GetEventNumber() << endl;
            cout << "Saved Event Number = " << ev_track << endl;
            ev_track++;*/

            if(SAVE){
                Copy_RecTree->Fill();
			}

        }

}

if(SAVE){
    //Make a directory to save the root file in
    string outdir = out_path + "Slim_Trim_Skim_Search";
    char SaveDir[600];
    sprintf(SaveDir, "mkdir %s", outdir.c_str());
    int success = system(SaveDir);
    if (success == 0){std::cout << "Directory " << SaveDir <<" created!" << std::endl;};

    string full_title = out_path + "Slim_Trim_Skim_Search/" + sts_file_name;

    char SaveRootFile[600];
    f = new TFile(full_title.c_str(), "RECREATE");

    TreeRec->GetEntry(0);
    cout << "source root file " << TreeRec->GetCurrentFile()->GetName() << endl;

    TObject* geo_tree = TreeRec->GetCurrentFile()->Get("GGeometry");
    geo_tree->Write("GGeometry"); //This does work!
    //Copy_GRecoTree->Write();
    Copy_RecTree->Write();
    //Copy_MCTree->Write();
    f->Close();
}

cout << "I am done" << endl;

}
