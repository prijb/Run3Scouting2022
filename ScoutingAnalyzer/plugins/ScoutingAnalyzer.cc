// -*- C++ -*-
//
// Package:    Scouting2022Analyzer/ScoutingAnalyzer
// Class:      ScoutingAnalyzer
//
/**\class ScoutingAnalyzer ScoutingAnalyzer.cc Scouting2022Analyzer/ScoutingAnalyzer/plugins/ScoutingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Prijith Babu Pradeep
//         Created:  Mon, 19 Jun 2023 00:18:32 GMT
//
//

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

// dataformats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/fillCovariance.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/Common/interface/AssociationMap.h"

// dataformats (PAT stuff?)
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// dataformats (scouting specific?)
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"

// dataformats (trigger)
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// what?? (gen??)
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// tracking tools 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// CMSSW stuff
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace std;

using reco::TrackCollection;

class ScoutingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingAnalyzer(const edm::ParameterSet&);
  ~ScoutingAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // Stuff that involves inheritance from WatchRuns and WatchLuminosityBlocks
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void clearVars();

  // ----------member data ---------------------------
  // The darn tokens (Will be initialised later)
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> pVtxToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> sVtxToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> pfToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron>> electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonsToken;

  const edm::EDGetTokenT<double> rhoToken;
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;  //HLT results token

  // Mapping trigger paths
  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  // L1 stuff
  bool doL1;
  triggerExpression::Data triggerCache_;

  // Some malarkey involving triggers
  unsigned char trig;
  edm::InputTag algInputTag_;
  edm::InputTag extInputTag_;
  edm::EDGetToken algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string> l1Seeds_;
  std::vector<bool> l1Result_;


  // This thing is also an ntuplizer, so now initialise the variables to fill a tree
  // Primary vertex details
  UInt_t n_pVtx;
  vector<Float16_t> pVtx_x;
  vector<Float16_t> pVtx_y;
  vector<Float16_t> pVtx_z;
  vector<Float16_t> pVtx_xError;
  vector<Float16_t> pVtx_yError;
  vector<Float16_t> pVtx_zError;
  vector<Int_t> pVtx_trksize;
  vector<Float16_t> pVtx_chi2;
  vector<Int_t> pVtx_ndof;
  vector<Bool_t> pVtx_isvalidvtx;

  //Secondary vertex details
  UInt_t n_sVtx;
  vector<Float16_t> sVtx_x;
  vector<Float16_t> sVtx_y;
  vector<Float16_t> sVtx_z;
  vector<Float16_t> sVtx_dxy;
  vector<Float16_t> sVtx_dxySig;
  vector<Float16_t> sVtx_xError;
  vector<Float16_t> sVtx_yError;
  vector<Float16_t> sVtx_zError;
  vector<Int_t> sVtx_trksize;
  vector<Float16_t> sVtx_chi2;
  vector<Int_t> sVtx_ndof;
  vector<Bool_t> sVtx_isvalidvtx;

  vector<Int_t> sVtx_nMuon;
  vector<Float16_t> sVtx_mass;

  /*
  //PF candidates
  UInt_t n_pf;
  vector<Float16_t> PFParticle_pt;
  vector<Float16_t> PFParticle_eta;
  vector<Float16_t> PFParticle_phi;
  vector<Int_t> PFParticle_pdgId;
  vector<Int_t> PFParticle_vertex;
  vector<Float16_t> PFParticle_normchi2;
  vector<Float16_t> PFParticle_dz;
  vector<Float16_t> PFParticle_dxy;
  vector<Float16_t> PFParticle_dzsig;
  vector<Float16_t> PFParticle_dxysig;
  vector<UInt_t>	PFParticle_lostInnerHits;
  vector<UInt_t>	PFParticle_quality;
  vector<Float16_t> PFParticle_trkpt;
  vector<Float16_t> PFParticle_trketa;
  vector<Float16_t> PFParticle_trkphi;
  vector<Bool_t> PFParticle_relativetrkvars;
  */

  //Electron colletion details (what is max_ele for??)
  
  const static int 	max_ele = 1000;
  UInt_t n_ele;
  vector<Float16_t> Electron_pt;
  vector<Float16_t> Electron_eta;
  vector<Float16_t> Electron_phi;
  vector<Float16_t> Electron_m;
  vector<Int_t> Electron_charge;
  vector<Float16_t> Electron_detain;
  vector<Float16_t> Electron_dphiin;
  vector<Float16_t> Electron_sigmaietaieta;
  vector<Float16_t> Electron_hoe;
  vector<Float16_t> Electron_ooemoop;
  vector<Int_t>	Electron_missinghits;
  vector<Float16_t> Electron_ecaliso;
  vector<Float16_t> Electron_hcaliso;
  vector<Float16_t> Electron_tkiso;
  vector<Float16_t> Electron_r9;
  vector<Float16_t> Electron_smin;
  vector<Float16_t> Electron_smaj;
  vector<UInt_t> Electron_seedid;
  vector<Bool_t> Electron_rechitzerosuppression;
  

  //Muon collection details (added vertex id)
  const static int 	max_mu = 1000;
  UInt_t n_mu;
  vector<Float16_t> Muon_pt;
  vector<Float16_t> Muon_eta;
  vector<Float16_t> Muon_phi;
  vector<Float16_t> Muon_m;
  vector<Float16_t> Muon_ecaliso;
  vector<Float16_t> Muon_hcaliso;
  vector<Float16_t> Muon_trkiso;
  vector<Float16_t> Muon_chi2;
  vector<Float16_t> Muon_ndof;
  vector<Float16_t> Muon_charge;
  vector<Float16_t> Muon_dxy;
  vector<Float16_t> Muon_dz;
  vector<Float16_t> Muon_dxyerror;
  vector<Float16_t> Muon_dzerror;
  vector<Float16_t> Muon_nvalidmuon_hits;
  vector<Float16_t> Muon_nvalidpixelhits;
  
  vector<Float16_t> Muon_nmatchedstations;
  vector<Float16_t> Muon_type;
  vector<Float16_t> Muon_nvalidstriphits;
  vector<Float16_t> Muon_trkqoverp;
  vector<Float16_t> Muon_trklambda;
  vector<Float16_t> Muon_trkpt;
  vector<Float16_t> Muon_trkphi;
  vector<Float16_t> Muon_trketa;
  vector<Float16_t> Muon_trkqoverperror;
  vector<Float16_t> Muon_trklambdaerror;
  vector<Float16_t> Muon_trkpterror;
  vector<Float16_t> Muon_trkphierror;
  vector<Float16_t> Muon_trketaerror;
  vector<Float16_t> Muon_trkdszerror;
  vector<Float16_t> Muon_trkdsz;

  //Vertex id
  vector<std::vector<int>> Muon_vtxIndx;

  //Rho info
  UInt_t n_rhoval;
  vector<Float16_t> rho;

  //The actual tree declared. We're gonna make it in the constructor dw
  TTree* tree;

  //Run and lumisection
  Int_t run;
  Int_t event;
  Int_t lumSec;




};

//
// constructors and destructor (a lotta initialisations of all the tokens)
// Initialisation is of the form
// tokenName(consumes<type>(iConfig.getParameter<edm::InputTag>(string you pass in python config)))
ScoutingAnalyzer::ScoutingAnalyzer(const edm::ParameterSet& iConfig):
  pVtxToken(consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("primaryVtx"))),
  sVtxToken(consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("secondaryVtx"))),
  //pfToken(consumes<std::vector<Run3ScoutingParticle>>(iConfig.getParameter<edm::InputTag>("pfparticles"))),
  electronsToken(consumes<std::vector<Run3ScoutingElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonsToken(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  //Get trigger related tokens
  triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerresults")),
  triggerResultsToken(consumes<edm::TriggerResults>(triggerResultsTag)),
  //What sorcery is this?
  doL1(iConfig.existsAs<bool>("doL1")?iConfig.getParameter<bool>("doL1"):false)
  {

  //now do what ever initialization is needed
  // What does usesResource mean??
  usesResource("TFileService");

  // L1 stuff?
  if (doL1){
    algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
    extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
    algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
    l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(
    iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }

  // Access the TFileService (help...)
  edm::Service<TFileService> fs;

  // Make the tree using the file service (huh....)
  tree = fs->make<TTree>("tree", "tree");

  // These three branches contain the event number, and what run and lumisection they belong to?
  tree->Branch("lumSec", &lumSec, "lumSec/i" );
  tree->Branch("run", &run, "run/i" );
  tree->Branch("event", &event, "event/i" );

  // The interesting stuff

  // Triggers
  tree->Branch("trig", &trig, "trig/b");
  tree->Branch("l1Result", "std::vector<bool>" ,&l1Result_, 32000, 0);

  // Primary vertex info (why does only the first branch have a /i?)
  tree->Branch("n_pVtx", &n_pVtx, "n_pVtx/i");
  tree->Branch("pVtx_x", &pVtx_x);
  tree->Branch("pVtx_y", &pVtx_y);
  tree->Branch("pVtx_z", &pVtx_z);
  tree->Branch("pVtx_xError", &pVtx_xError);
  tree->Branch("pVtx_yError", &pVtx_yError);
  tree->Branch("pVtx_zError", &pVtx_zError);
  tree->Branch("pVtx_trksize", &pVtx_trksize);
  tree->Branch("pVtx_chi2", &pVtx_chi2);
  tree->Branch("pVtx_ndof", &pVtx_ndof);
  tree->Branch("pVtx_isvalidvtx", &pVtx_isvalidvtx);

  // Secondary vertex info 
  tree->Branch("n_sVtx", &n_sVtx, "n_sVtx/i");
  tree->Branch("sVtx_x", &sVtx_x);
  tree->Branch("sVtx_y", &sVtx_y);
  tree->Branch("sVtx_z", &sVtx_z);
  tree->Branch("sVtx_xError", &sVtx_xError);
  tree->Branch("sVtx_yError", &sVtx_yError);
  tree->Branch("sVtx_zError", &sVtx_zError);
  tree->Branch("sVtx_dxy", &sVtx_dxy);
  tree->Branch("sVtx_dxySig", &sVtx_dxySig);
  tree->Branch("sVtx_trksize", &sVtx_trksize);
  tree->Branch("sVtx_chi2", &sVtx_chi2);
  tree->Branch("sVtx_ndof", &sVtx_ndof);
  tree->Branch("sVtx_isvalidvtx", &sVtx_isvalidvtx);

  tree->Branch("sVtx_nMuon", &sVtx_nMuon);
  tree->Branch("sVtx_mass", &sVtx_mass);

  /*
  // Pf particles
  tree->Branch("n_pf", &n_pf, "n_pf/i");
  tree->Branch("PFParticle_pt", &PFParticle_pt);
  tree->Branch("PFParticle_eta", &PFParticle_eta);
  tree->Branch("PFParticle_phi", &PFParticle_phi);
  tree->Branch("PFParticle_pdgId", &PFParticle_pdgId);
  tree->Branch("PFParticle_vertex", &PFParticle_vertex);
  tree->Branch("PFParticle_normchi2", &PFParticle_normchi2);
  tree->Branch("PFParticle_dz", &PFParticle_dz);
  tree->Branch("PFParticle_dxy", &PFParticle_dxy);
  tree->Branch("PFParticle_dzsig", &PFParticle_dzsig);
  tree->Branch("PFParticle_dxysig", &PFParticle_dxysig);
  tree->Branch("PFParticle_lostInnerHits", &PFParticle_lostInnerHits);
  tree->Branch("PFParticle_quality", &PFParticle_quality);
  tree->Branch("PFParticle_trkpt", &PFParticle_trkpt);
  tree->Branch("PFParticle_trketa", &PFParticle_trketa);
  tree->Branch("PFParticle_trkphi", &PFParticle_trkphi);
  tree->Branch("PFParticle_relativetrkvars", &PFParticle_relativetrkvars);
  */



  // Electrons
  // Note: Some shenanigans involving electron variables
  
  tree->Branch("n_ele", &n_ele, "n_ele/i");
  tree->Branch("Electron_pt", &Electron_pt);
  tree->Branch("Electron_eta", &Electron_eta);
  tree->Branch("Electron_phi", &Electron_phi);
  tree->Branch("Electron_m", &Electron_m);
  tree->Branch("Electron_charge", &Electron_charge);
  tree->Branch("Electron_detain", &Electron_detain);
  tree->Branch("Electron_dphiin", &Electron_dphiin);
  tree->Branch("Electron_sigmaietaieta", &Electron_sigmaietaieta);
  tree->Branch("Electron_hoe", &Electron_hoe);
  tree->Branch("Electron_ooemoop", &Electron_ooemoop);
  tree->Branch("Electron_missinghits", &Electron_missinghits);
  tree->Branch("Electron_ecaliso", &Electron_ecaliso);
  tree->Branch("Electron_hcaliso", &Electron_hcaliso);
  tree->Branch("Electron_tkiso", &Electron_tkiso);
  tree->Branch("Electron_r9", &Electron_r9);
  tree->Branch("Electron_smin", &Electron_smaj);
  tree->Branch("Electron_smaj", &Electron_smin);
  tree->Branch("Electron_seedid", &Electron_seedid);
  tree->Branch("Electron_rechitzerosuppression", &Electron_rechitzerosuppression);
  
  // Muons
  tree->Branch("n_mu",&n_mu,"n_mu/i");
  tree->Branch("Muon_pt", &Muon_pt);
  tree->Branch("Muon_eta", &Muon_eta);
  tree->Branch("Muon_phi", &Muon_phi);
  tree->Branch("Muon_m", &Muon_m);
  tree->Branch("Muon_ecaliso", &Muon_ecaliso);
  tree->Branch("Muon_hcaliso", &Muon_hcaliso);
  tree->Branch("Muon_trkiso", &Muon_trkiso);
  tree->Branch("Muon_chi2", &Muon_chi2);
  tree->Branch("Muon_ndof", &Muon_ndof);
  tree->Branch("Muon_charge", &Muon_charge);
  tree->Branch("Muon_dxy", &Muon_dxy);
  tree->Branch("Muon_dz", &Muon_dz);
  tree->Branch("Muon_dxyerror", &Muon_dxyerror);
  tree->Branch("Muon_dzerror", &Muon_dzerror);
  tree->Branch("Muon_nvalidmuon_hits", &Muon_nvalidmuon_hits);
  tree->Branch("Muon_validpixelhits", &Muon_nvalidpixelhits );
  
  tree->Branch("Muon_nmatchedstations", &Muon_nmatchedstations);
  tree->Branch("Muon_type", &Muon_type);
  tree->Branch("Muon_nvalidstriphits", &Muon_nvalidstriphits);
  tree->Branch("Muon_trkqoverp", &Muon_trkqoverp);
  tree->Branch("Muon_trklambda", &Muon_trklambda);
  tree->Branch("Muon_trkpt", &Muon_trkpt);
  tree->Branch("Muon_trkphi", &Muon_trkphi);
  tree->Branch("Muon_trketa", &Muon_trketa);
  tree->Branch("Muon_trkqoverperror", &Muon_trkqoverperror);
  tree->Branch("Muon_trklambdaerror", &Muon_trklambdaerror);
  tree->Branch("Muon_trkpterror", &Muon_trkpterror);
  tree->Branch("Muon_trkphierror", &Muon_trkphierror);
  tree->Branch("Muon_trketaerror", &Muon_trketaerror);
  tree->Branch("Muon_trkdzerror", &Muon_trkdszerror);
  tree->Branch("Muon_trkdz", &Muon_trkdsz);

  tree->Branch("Muon_vtxIndx", &Muon_vtxIndx);

  // Rho
  tree->Branch("n_rhoval", &n_rhoval, "n_rhoval/i");
  tree->Branch("rho", &rho);

}

// Destructor!
ScoutingAnalyzer::~ScoutingAnalyzer() {
}

//
// member functions
//

// Fill stuff

// ------------ method called for each event  ------------
void ScoutingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  using namespace std;
  using namespace reco;

  //Token's already initalized, so now we use this token to fill a newly initialized handle
  //This is the FWLite stuff you're more familiar with but FWLite gets by label rather than token

  //Get trigger handles
  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken(triggerResultsToken, triggerResultsHandle);
  bool triggerResultsValid = triggerResultsHandle.isValid();

  Handle<vector<Run3ScoutingVertex>> pVtxHandle;
  iEvent.getByToken(pVtxToken, pVtxHandle);
  bool pVtxValid = pVtxHandle.isValid();

  Handle<vector<Run3ScoutingVertex>> sVtxHandle;
  iEvent.getByToken(sVtxToken, sVtxHandle);
  bool sVtxValid = sVtxHandle.isValid();

  //Handle<vector<Run3ScoutingParticle>> pfHandle;
  //iEvent.getByToken(pfToken, pfHandle);
  //bool pfValid = pfHandle.isValid();

  
  Handle<vector<Run3ScoutingElectron>> electronsHandle;
  iEvent.getByToken(electronsToken, electronsHandle);
  bool electronValid = electronsHandle.isValid();
  

  Handle<vector<Run3ScoutingMuon>> muonsHandle;
  iEvent.getByToken(muonsToken, muonsHandle);
  bool muonValid = muonsHandle.isValid();

  Handle<double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  bool rhoValid = rhoHandle.isValid();

  //No clue what this is but it directly feeds into the branch without handles
  run = iEvent.eventAuxiliary().run();
  event = iEvent.eventAuxiliary().event();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();
  
  //Now fill collections!

  //HLT. Why is trig summed over all events this way without resetting?
  trig = 0;
  //Check how triggers are fired?
  if(triggerResultsValid){
    for (size_t i = 0; i < triggerPathsVector.size(); i++){
      //??
      if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
      //If the 
      if (i == 0  && triggerResultsHandle->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   1; // DST_PixelTracking
      if (i == 1  && triggerResultsHandle->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   2; // DST_HLTMu
    }
  }
  
  bool mu_or_bit = 0;

  //L1 
  if (doL1) {
   //Get l1GtUtils from the event, event setup and the token
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);

  
  //Runs over all the L1 algorithms and checks if event passes or not
  /*
    for( int r = 0; r<280; r++){
	    string name ("empty");
	    bool algoName_ = false;
	    algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
	    cout << "getAlgNameFromBit = " << algoName_  << endl;
	    cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
	  }
  */

    //I think this fills the data with the L1 decisions?
    for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
      bool l1htbit = 0;	
      l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
      //cout<<string(l1Seeds_[iseed])<<"  "<<l1htbit<<endl;
      //cout << "L1 seed name: " << string(l1Seeds_[iseed]) << endl;
      //cout << "L1 decision: " <<  l1htbit << endl;
      //cout << "" << endl;
      l1Result_.push_back( l1htbit );
      mu_or_bit += l1htbit;
      }
  }


  //Fill pvtx where you count and add that to event as well
  //Check if handle is valid
  //* dereferences handle and gets collection. Then pass each element by reference?
  n_pVtx = 0;
  if(pVtxValid){
    for (auto &pVtx: *pVtxHandle){
      //Use the pointer to pVtx 
      auto *pVtx_iter = &pVtx;
      pVtx_x.push_back(pVtx_iter->x());
      pVtx_y.push_back(pVtx_iter->y());
      pVtx_z.push_back(pVtx_iter->z());
      pVtx_xError.push_back(pVtx_iter->xError());
      pVtx_yError.push_back(pVtx_iter->yError());
      pVtx_zError.push_back(pVtx_iter->zError());
      pVtx_trksize.push_back(pVtx_iter->tracksSize());
      pVtx_chi2.push_back(pVtx_iter->chi2());
      pVtx_ndof.push_back(pVtx_iter->ndof());
      pVtx_isvalidvtx.push_back(pVtx_iter->isValidVtx());
      n_pVtx++;
    }
  }

  //Fill svtx
  n_sVtx = 0;

  if(sVtxValid){
    for (auto &sVtx: *sVtxHandle){
      auto *sVtx_iter = &sVtx;
      sVtx_x.push_back(sVtx_iter->x());
      sVtx_y.push_back(sVtx_iter->y());
      sVtx_z.push_back(sVtx_iter->z());
      //Calculate dxy
      Float16_t dx = sVtx_iter->x() - pVtx_x.at(0);
      Float16_t dy = sVtx_iter->y() - pVtx_y.at(0);
      Float16_t dxy = TMath::Sqrt(dx*dx + dy*dy);
      sVtx_dxy.push_back(dxy);
      sVtx_xError.push_back(sVtx_iter->xError());
      sVtx_yError.push_back(sVtx_iter->yError());
      sVtx_zError.push_back(sVtx_iter->zError());
      //Calculate dxysig
      Float16_t dxerr = TMath::Sqrt(std::pow(sVtx_iter->xError(), 2) + std::pow(pVtx_xError.at(0), 2));
      Float16_t dyerr = TMath::Sqrt(std::pow(sVtx_iter->yError(), 2) + std::pow(pVtx_yError.at(0), 2));
      Float16_t dxyerr = TMath::Sqrt(std::pow(dx*dxerr, 2) + std::pow(dy*dyerr, 2))/dxy;
      Float16_t dxysig = dxy/dxyerr;
      sVtx_dxySig.push_back(dxysig);
      sVtx_trksize.push_back(sVtx_iter->tracksSize());
      sVtx_chi2.push_back(sVtx_iter->chi2());
      sVtx_ndof.push_back(sVtx_iter->ndof());
      sVtx_isvalidvtx.push_back(sVtx_iter->isValidVtx());

      //Fill with mass
      TLorentzVector sVtx_Sum;
      UInt_t n_sVtx_mu = 0;


      if (muonValid){
        for (auto &muo: *muonsHandle){
          auto *muons_iter = &muo;
          const auto& Muon_vtxIndx_inter = muons_iter->vtxIndx();

          if(std::find(Muon_vtxIndx_inter.begin(), Muon_vtxIndx_inter.end(), n_sVtx) != Muon_vtxIndx_inter.end()){

            TLorentzVector muon_fourvec;
            muon_fourvec.SetPtEtaPhiM(muons_iter->pt(), muons_iter->eta(), muons_iter->phi(), muons_iter->m());
            sVtx_Sum += muon_fourvec;
            n_sVtx_mu += 1;
          }
        }
      }

      sVtx_nMuon.push_back(n_sVtx_mu);
      sVtx_mass.push_back(sVtx_Sum.M());
      n_sVtx++;
    }
  }
  /*
  n_pf = 0;
  if(pfValid) {
    for (auto &pfcand : *pfHandle) {
      auto *pfcand_iter = &pfcand;
      PFParticle_pt.push_back(pfcand_iter->pt());
      PFParticle_eta.push_back(pfcand_iter->eta());
      PFParticle_phi.push_back(pfcand_iter->phi());
      PFParticle_pdgId.push_back(pfcand_iter->pdgId());
      PFParticle_vertex.push_back(pfcand_iter->vertex());
      PFParticle_normchi2.push_back(pfcand_iter->normchi2());
      PFParticle_dz.push_back(pfcand_iter->dz());
      PFParticle_dxy.push_back(pfcand_iter->dxy());
      PFParticle_dzsig.push_back(pfcand_iter->dzsig());
      PFParticle_dxysig.push_back(pfcand_iter->dxysig());
      PFParticle_lostInnerHits.push_back(pfcand_iter->lostInnerHits());
      PFParticle_quality.push_back(pfcand_iter->quality());
      PFParticle_trkpt.push_back(pfcand_iter->trk_pt());
      PFParticle_trketa.push_back(pfcand_iter->trk_eta());
      PFParticle_trkphi.push_back(pfcand_iter->trk_phi());
      PFParticle_relativetrkvars.push_back(pfcand_iter->relative_trk_vars());
      n_pf++;
    } 
  } 
  */
  
  //Fill electrons
  n_ele = 0;
  if(electronValid) {
    for (auto &ele : *electronsHandle) {
      auto *electrons_iter = &ele;
      Electron_pt.push_back(electrons_iter->pt());
      Electron_eta.push_back(electrons_iter->eta());
      Electron_phi.push_back(electrons_iter->phi());	
      Electron_m.push_back(electrons_iter->m());
      Electron_charge.push_back(electrons_iter->charge());
      Electron_detain.push_back(electrons_iter->dEtaIn());
      Electron_dphiin.push_back(electrons_iter->dPhiIn());
      Electron_sigmaietaieta.push_back(electrons_iter->sigmaIetaIeta());
      Electron_hoe.push_back(electrons_iter->hOverE());	
      Electron_ooemoop.push_back(electrons_iter->ooEMOop());
      Electron_missinghits.push_back(electrons_iter->missingHits());
      Electron_ecaliso.push_back(electrons_iter->ecalIso());
      Electron_hcaliso.push_back(electrons_iter->hcalIso());
      Electron_tkiso.push_back(electrons_iter->trackIso());
      Electron_r9.push_back(electrons_iter->r9());
      Electron_smin.push_back(electrons_iter->sMin());
      Electron_smaj.push_back(electrons_iter->sMaj());
      Electron_seedid.push_back(electrons_iter->seedId());
      Electron_rechitzerosuppression.push_back(electrons_iter->rechitZeroSuppression());
      n_ele++;
    } 
  } 
  
  //Fill muons
  n_mu = 0;
  if(muonValid) {
    for (auto &muo : *muonsHandle) {
      auto *muons_iter = &muo;
      Muon_pt.push_back(muons_iter->pt());
      Muon_eta.push_back(muons_iter->eta());
      Muon_phi.push_back(muons_iter->phi());	
      Muon_m.push_back(muons_iter->m());
      Muon_ecaliso.push_back(muons_iter->ecalIso());
      Muon_hcaliso.push_back(muons_iter->hcalIso());
      Muon_trkiso.push_back(muons_iter->trackIso());
      Muon_chi2.push_back(muons_iter->normalizedChi2());
      Muon_ndof.push_back(muons_iter->trk_ndof());
      Muon_charge.push_back(muons_iter->charge());
      Muon_dxy.push_back(muons_iter->trk_dxy());
      Muon_dxyerror.push_back(muons_iter->trk_dxyError());
      Muon_dz.push_back(muons_iter->trk_dz());
      Muon_dzerror.push_back(muons_iter->trk_dzError());
      Muon_nvalidmuon_hits.push_back(muons_iter->nValidRecoMuonHits());
      Muon_nvalidpixelhits.push_back(muons_iter->nValidPixelHits());
      Muon_nvalidstriphits.push_back(muons_iter->nValidStripHits());
      Muon_nmatchedstations.push_back(muons_iter->nRecoMuonMatchedStations());

      Muon_type.push_back(muons_iter->type());
      Muon_trkqoverp.push_back(muons_iter->trk_qoverp());
      Muon_trklambda.push_back(muons_iter->trk_lambda());
      Muon_trkpt.push_back(muons_iter->trk_pt());
      Muon_trkphi.push_back(muons_iter->trk_phi());
      Muon_trketa.push_back(muons_iter->trk_eta());
      Muon_trkqoverperror.push_back(muons_iter->trk_dxyError());
      Muon_trklambdaerror.push_back(muons_iter->trk_dzError());
      Muon_trkpterror.push_back(muons_iter->trk_qoverpError());
      Muon_trkphierror.push_back(muons_iter->trk_lambdaError());
      Muon_trketaerror.push_back(muons_iter->trk_phiError());
      Muon_trkdsz.push_back(muons_iter->trk_dsz());
      Muon_trkdszerror.push_back(muons_iter->trk_dszError());

      Muon_vtxIndx.push_back(muons_iter->vtxIndx());
      n_mu++;
    } 
  } 

  //Fill rho
  n_rhoval = 0;
  if(rhoValid) {
    rho.push_back(*rhoHandle);
    n_rhoval++;
  }

  //To study events that mass muon OR
  /*
  if(mu_or_bit == true){
  cout << "Event passes Muon OR" << endl;
  cout << "Number of muons in event: " << n_mu << endl;
  }
  */
  //Finally fill the tree! (Clear vars does the cleaning)
  tree->Fill();
  clearVars();

}

// ------------ clear variables

void ScoutingAnalyzer::clearVars() {

  l1Result_.clear();

  pVtx_x.clear();
  pVtx_y.clear();
  pVtx_z.clear();
  pVtx_xError.clear();
  pVtx_yError.clear();
  pVtx_zError.clear();
  pVtx_trksize.clear();
  pVtx_chi2.clear();
  pVtx_ndof.clear();
  pVtx_isvalidvtx.clear();
  
  sVtx_x.clear();
  sVtx_y.clear();
  sVtx_z.clear();
  sVtx_xError.clear();
  sVtx_yError.clear();
  sVtx_zError.clear();
  sVtx_dxy.clear();
  sVtx_dxySig.clear();
  sVtx_trksize.clear();
  sVtx_chi2.clear();
  sVtx_ndof.clear();
  sVtx_isvalidvtx.clear();
  sVtx_nMuon.clear();
  sVtx_mass.clear();

  /*
  PFParticle_pt.clear();
  PFParticle_eta.clear();
  PFParticle_phi.clear();
  PFParticle_pdgId.clear();
  PFParticle_vertex.clear();
  PFParticle_normchi2.clear();
  PFParticle_dz.clear();
  PFParticle_dxy.clear();
  PFParticle_dzsig.clear();
  PFParticle_dxysig.clear();
  PFParticle_lostInnerHits.clear();
  PFParticle_quality.clear();
  PFParticle_trkpt.clear();
  PFParticle_trketa.clear();
  PFParticle_trkphi.clear();
  PFParticle_relativetrkvars.clear();
  */
  
  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  Electron_m.clear();
  Electron_charge.clear();
  Electron_detain.clear();
  Electron_dphiin.clear();
  Electron_sigmaietaieta.clear();
  Electron_hoe.clear();
  Electron_ooemoop.clear();
  Electron_missinghits.clear();
  Electron_ecaliso.clear();
  Electron_hcaliso.clear();
  Electron_tkiso.clear();
  Electron_r9.clear();
  Electron_smin.clear();
  Electron_smaj.clear();
  Electron_seedid.clear();
  Electron_rechitzerosuppression.clear();
  

  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_m.clear();
  Muon_ecaliso.clear();
  Muon_hcaliso.clear();
  Muon_trkiso.clear();
  Muon_chi2.clear();
  Muon_ndof.clear();
  Muon_charge.clear();
  Muon_dxy.clear();
  Muon_dxyerror.clear();
  Muon_dz.clear();
  Muon_dzerror.clear();
  Muon_nvalidmuon_hits.clear();
  Muon_nvalidpixelhits.clear();
  Muon_nmatchedstations.clear();
  Muon_type.clear();
  Muon_nvalidstriphits.clear();
  Muon_trkqoverp.clear();
  Muon_trklambda.clear();
  Muon_trkpt.clear();
  Muon_trkphi.clear();
  Muon_trketa.clear();
  Muon_trkqoverperror.clear();
  Muon_trklambdaerror.clear();
  Muon_trkpterror.clear();
  Muon_trkphierror.clear();
  Muon_trketaerror.clear();
  Muon_trkdszerror.clear();
  Muon_trkdsz.clear();

  Muon_vtxIndx.clear();

  rho.clear();

}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingAnalyzer::endJob() {
  // please remove this method if not needed
}

void ScoutingAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  // HLT paths

  //For 2022
  triggerPathsVector.push_back("DST_Run3_PFScoutingPixelTracking_v*");
  triggerPathsVector.push_back("DST_HLTMuon_Run3_PFScoutingPixelTracking_v*");
 
  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
      std::string pathName = hltConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }
}

void ScoutingAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingAnalyzer);
