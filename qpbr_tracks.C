#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEvePointSet.h"
#include "TEveBoxSet.h"
#include "TEveLine.h"
#include "TEveScene.h"
#include "TEveEventManager.h"
#include "TEveWindow.h"

#include "TEveTrackPropagator.h"
#include "TTree.h"
#include "TEveTrack.h"
#include "TEveViewer.h"
#include "TGLViewer.h"
#include "TEveBrowser.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGTab.h"
#include "TGedEditor.h"
#include "TEveGedEditor.h"
#include "TGLSAViewer.h"
#include "TEveTrans.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TGLWidget.h"
#include "TEveLegoEventHandler.h"
#include "TGLEmbeddedViewer.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoTrack.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TTimer.h"
#include "TPad.h"
#include "TROOT.h"
#include "RtypesCore.h"
#include "TMath.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/Cylindrical3D.h"
#include "TEveCaloLegoOverlay.h"
#include "TEveCalo.h"

//==============================================================================
// Forward declarations
//------------------------------------------------------------------------------
TEveCalo3D*   MakeCalo3D(TEveCaloData* data, TEveWindowSlot* slot);
TEveCaloLego* MakeCaloLego(TEveCaloData* data, TEveWindowSlot* slot);
void          MakeViewerScene(TEveWindowSlot* slot, TEveViewer*& v, TEveScene*& s);

//==============================================================================
// Globals
//==============================================================================
bool       gdebug=kTRUE;
//bool       debug2=kFALSE;
bool       debug2=kTRUE ;

Double_t kCenter[3]  = {0., 0., 0.};                                      // GL camera center (geometry)
Double_t kfov       = 30.;
Double_t kdolly     = 0.;
Double_t khRotate   = TMath::Pi()/30.;
Double_t kvRotate   = TMath::Pi()/8.;

typedef std::list<TEveElement*>     LList_t;
typedef LList_t::iterator           LList_i;

//==============================================================================
// Constants.
//------------------------------------------------------------------------------
const Double_t kX_max  = 300;       // Top Volume Box dimensions
const Double_t kY_max  = 300;
const Double_t kZ_max  = 300;

const Double_t kR_min = 240;        // child Tube dimensions
const Double_t kR_max = 250;
const Double_t kZ_d   = 300;

const TLorentzVector vrtxo    = TLorentzVector(0.,  0.,     0.,    0.);   // collision point
const TEveVectorD    origin   = TEveVectorD(0., 0., 0.);

const Int_t kSleep = 1500;                                                 // time loop delay in [ms]

//------------------------------------------------------------------------------

void pbr_create_geo()
{
// define a new geometry object "EZB_sim_geo"
    if (gGeoManager) {
       delete gGeoManager;
    }                      
// new call below initializes gGeoManager global pointer
    new TGeoManager("EZB_sim_geo", "simple detector geometry-cube enclosure+cylinder");

//-- materials definitions
     TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
     TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
//-- media definitions
     TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
     TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
     matVacuum->SetTransparency(80);
     matAl->SetTransparency(70);

// [0] this is the top volume cube
     TGeoVolume *top = gGeoManager->MakeBox("TOP",Vacuum,kX_max,kY_max,kZ_max);
     gGeoManager->SetTopVisible();
     gGeoManager->SetTopVolume(top);
     top->SetFillColor(3);

// [1] this is the cylinder
     TGeoVolume *vcylinder= gGeoManager->MakeTube("Cylinder",Al, kR_min, kR_max, kZ_d);  // Al <-> nullptr <- placeholder for TGeoMedium*
                 vcylinder->SetLineColor(kGreen);
     TGeoShape  *scylinder= vcylinder->GetShape();

     top->AddNode(vcylinder,1);

     gGeoManager->CloseGeometry();   

}

void qpbr_tracks()
{

     Double_t rfP0    = 1.0;  // set incoming beam particle momentum
     Double_t rfE0    = 10.0; // set incoming beam particle energy 

     Double_t tauc    = 1.0;  // collision time (arbitrary)

//_ geometry rendering
     pbr_create_geo();
     TEveManager::Create();   // this initializes gEve global pointer
     TEveGeoTopNode* top = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode()); // "/TOP_1/Cylinder_1"
     top->SetVisOption(0);
     top->SetVisLevel(6);
     gEve->AddGlobalElement(top);
//_ end geometry rendering

//_ track generation
     auto propag = new TEveTrackPropagator();
     propag->SetStepper(TEveTrackPropagator::kRungeKutta);
     propag->SetMaxR(kR_min);                              // these set maximum extents of track propagation
     propag->SetMaxZ(kZ_d);
     //propag->SetRnrDecay(kFALSE);

     auto tracks  = new TEveTrackList("All Tracks");
     auto beams   = new TEveTrackList("Beams");
     auto ctracks = new TEveTrackList("Products");

     auto pmc = new TEvePathMarkD(TEvePathMarkD::kDecay, origin, tauc ); // collision vertex mark
                                                                         // at the origin
     TEveVector4D endpoint;
     Double_t     Rho;
     Double_t     Phi = 0.;
     Double_t     rfE;
     Int_t        itrk = 0;
     TTree *tep = new TTree("tep", "TTree itrk,x,y,z,t,E");              // endpoint storage for hits on Cylinder_1
     tep->Branch("itrk",&itrk,"itrk/I");
     tep->Branch("x",&endpoint.fX,"x/D");
     tep->Branch("y",&endpoint.fY,"y/D");
     tep->Branch("z",&endpoint.fZ,"z/D");
     tep->Branch("t",&endpoint.fT,"t/D");
     tep->Branch("E",&rfE,"E/D");

     ////////////////////////////////
     // 2x primary particle beams  //
     ////////////////////////////////
     const Int_t Nptracks = 2;             
     auto p1 = new TParticle();
     auto p2 = new TParticle();
     // primary vertex
     p1->SetProductionVertex(0.,0., kZ_max   ,0.);
     p2->SetProductionVertex(0.,0.,-kZ_max   ,0.);
     p1->SetMomentum(0.,0.,-rfP0,rfE0);
     p2->SetMomentum(0.,0.,+rfP0,rfE0);

     auto track = new TEveTrack(p1,0,propag);
     track->AddPathMark(*pmc);
     track->MakeTrack();
//_generate track time (track 0)
     auto trackp = propag->GetLastPoints();     // const std::vector<TEveVector4D>&
     endpoint  = trackp.back();
     endpoint += TEveVector4D(0.0,0.0,0.0, tauc);
     if(gdebug){
     endpoint.Dump(); cout<< " Track 0" << endl;
     }
//_end generate track time
     track->SetElementName("Beam 0");
     track->SetMainColor(kYellow);
     track->PrintPathMarks();
     tracks->AddElement(track);
     TEveTrack *track0 = track;     

          track = new TEveTrack(p2,1,propag);
     track->AddPathMark(*pmc);
     track->MakeTrack();
//_generate track time (track 1)
          trackp = propag->GetLastPoints();
     endpoint  = trackp.back();
     endpoint += TEveVector4D(0.0,0.0,0.0, tauc);
     if(gdebug){
     endpoint.Dump(); cout<< " Track 1" << endl;
     }
//_end generate track time
     track->SetElementName("Beam 1");
     track->SetMainColor(kYellow);
     tracks->AddElement(track);
     TEveTrack *track1 = track;     


     /////////////////////////////////////////////////////////////////
     // 4x secondary particle products (arbitrary Energy+Momentum)  //
     /////////////////////////////////////////////////////////////////
     auto rnd = gRandom;

     const Int_t Nstracks = 4;
     for (Int_t i =Nptracks; i<(Nstracks+Nptracks); i++) {

       auto p = new TParticle();
       rfE = rnd->Uniform(rfE0);
       p->SetProductionVertex(vrtxo);
       p->SetMomentum(4.*pow(-1,i),0.,(1.+i)*pow(-1,i),rfE);

       auto track = new TEveTrack(p,i,propag);
       track->MakeTrack();
//_establish parent-child linkage
       track0->AddElement(track);
       track1->AddElement(track);
        track->AddParent(track0);
        track->AddParent(track1);
//_end establish parent-child linkage
//_generate track time (track i)
       trackp = propag->GetLastPoints();
       endpoint  = trackp.back();
       endpoint += TEveVector4D(0.0,0.0,0.0, tauc+i);  // detection time (arbitrary)
       endpoint.Dump();
//_end generate track time
       track->SetMainColor(kBlue);
       track->SetName(Form("p %d",i-Nptracks+1));      // children track labeling starts at 1
       track->SetLineColor((pow(-1,i) > 0 ? kBlue : kRed));

       track->SetRnrPoints(kTRUE);
       track->SetMarkerStyle(4);
       track->PrintPathMarks();

       ctracks->AddElement(track);

//_detect Cylinder_1 hits

       Rho = TMath::Sqrt(endpoint.fX*endpoint.fX+endpoint.fY*endpoint.fY);
       if(Rho == kR_min) {
          // we have a hit, store x,y,z,t,E
          itrk = i;
          cout << "track " << itrk << " matches" << endl;      
          tep->Fill();    // store to tree
       }

//_end detect Cylinder_1 hits

     }

     gEve->GetEventScene()->AddElement(tracks);
     gEve->GetDefaultViewer()->SetElementName("GL-Full");

//_ end track generation

     // setup camera view
     auto vGL  = gEve->GetDefaultGLViewer();
     //vGL->CurrentCamera().SetExternalCenter(kTRUE);
     vGL->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
     vGL->CurrentCamera().Configure(kfov*1.5,kdolly,kCenter,khRotate,kvRotate);
     vGL->RequestDraw(TGLRnrCtx::kLODHigh);

//_create lego plots (on TCanvas)

     if(gdebug){
     tep->Print();
     tep->Show(1);
     tep->Scan("z:E");
     }

     const Double_t k1eps = 1.0 + 0.01;         // small tolerance to ensure values at the min/max boundaries are counted

     TH2D *h2  = new TH2D("ZPhiE","Z,Phi vs E",10,-(kZ_d*k1eps),(kZ_d*k1eps),10,-(TMath::Pi()*k1eps),(TMath::Pi()*k1eps));
     ROOT::Math::Cylindrical3D< Double_t > p3D; // Cartesian to Cylindrical coordinate transform
     Int_t nentries = (Int_t)tep->GetEntries();

     for (Int_t i=0; i<nentries; i++) {

        tep->GetEntry(i);
        p3D.SetXYZ(endpoint.fX,endpoint.fY,endpoint.fZ);
        Rho = p3D.Rho();
        Phi = p3D.Phi();
        if(gdebug) cout<< "Rho Z Phi rfE "<<Rho<< " "<< p3D.Z()<<" " << Phi<<" "<<rfE<<endl;
        h2->Fill(p3D.Z(),Phi,rfE);

     }

//_end create lego plots on TCanvas

//_create two-tabbed TEve display
     // Create the second tab 
     TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
     gEve->GetBrowser()->GetTabRight()->SetTab(1);

//_create two vertical slots in tab #2
     // Set new window for 2 vertical sub-slots
     TEveWindowPack* packH = slot->MakePack();
     packH->SetElementName("Composite");
     packH->SetShowTitleBar(kFALSE);
     slot = packH->NewSlot();
     TEveWindowPack* pack0 = slot->MakePack();
     pack0->SetShowTitleBar(kFALSE);
     auto  slotTop     = pack0->NewSlot();
     auto  slotsBottom = pack0->NewSlot();

//_create two horizontal slots in bottom tab #2
     TEveWindowPack* pack1 = slotsBottom->MakePack();
     pack1->SetHorizontal();
     pack1->SetShowTitleBar(kFALSE);
     auto  slotBottLeft    = pack1->NewSlot();
     auto  slotBottRight   = pack1->NewSlot();

//create top track plot in Eve 
     auto  v  = new TEveViewer("GL-Hits");
     //auto  sb = new  TEveScene("Track-hits");
     v->SpawnGLEmbeddedViewer(gEve->GetEditor());
     vGL = v->GetGLViewer();
     vGL->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
     vGL->CurrentCamera().Configure(1.1*kfov,kdolly,kCenter,khRotate,kvRotate);
     vGL->RequestDraw(TGLRnrCtx::kLODHigh);
     slotTop->ReplaceWindow(v);

     gEve->GetViewers()->AddElement(v);

     v->AddScene(gEve->GetGlobalScene());    // This is the Geometry
     auto sb = gEve->SpawnNewScene("Track-Hits");
      v->AddScene(sb);   
     gSystem->ProcessEvents();

     TEveTrackList *cbeams  = new TEveTrackList("Collider Beams");
     TEveTrack *jtrack0 = new TEveTrack(*track0);    // shallow copies, no children
     TEveTrack *jtrack1 = new TEveTrack(*track1);
     jtrack0->SetElementName("b0");
     jtrack1->SetElementName("b1");
     jtrack0->SetRnrPoints(kTRUE);
     jtrack1->SetRnrPoints(kTRUE);
     jtrack0->SetMarkerStyle(4);
     jtrack1->SetMarkerStyle(4);
     jtrack0->MakeTrack();
     jtrack1->MakeTrack();
     cbeams->AddElement(jtrack0);
     cbeams->AddElement(jtrack1);
     sb->AddElement(cbeams);
     gSystem->ProcessEvents();

//_end create top track plot in Eve 

//create bottom lego plot in Eve 
     auto data = new TEveCaloDataHist();
     data->AddHistogram((TH2F *)h2);
     data->RefSliceInfo(0).Setup("ECAL", 0.000, kRed );
     data->GetEtaBins()->SetTitleFont(100);
     data->GetEtaBins()->SetTitle("Z ");
     data->GetPhiBins()->SetTitleFont(100);
     data->GetPhiBins()->SetTitle("Phi ");
     data->IncDenyDestroy();

     gEve->AddToListTree(data, kFALSE);

     auto lego   = MakeCaloLego(data, slotBottLeft);
     gSystem->ProcessEvents();

//_end create bottom lego plot in Eve 

//_create lego/track animation (on TCanvas)

     slotBottRight->StartEmbedding();
     auto c2 = new TCanvas("Animation","Animation E Distribution");
     slotBottRight->StopEmbedding();
     gSystem->ProcessEvents();

     std::string s;
     TPaveText *pt = new TPaveText(-0.95,1.09,-0.33,0.92);
     TH2D *h2a = new TH2D("ZPhiEa","Z,Phi vs E",10,-(kZ_d*k1eps),(kZ_d*k1eps),10,-(TMath::Pi()*k1eps),(TMath::Pi()*k1eps));

     //Sort tep TTree in order of endpoint.fT (time)
     Int_t *index = new Int_t[nentries];
     tep->Draw("t","","goff");                           // this extracts
     TMath::Sort(nentries,tep->GetV1(),index, kFALSE);   // and sorts the t-values

     // (pseudo) time-loop
     for (Int_t i=0;i<nentries;i++) {

        tep->GetEntry(index[i]);
        p3D.SetXYZ(endpoint.fX,endpoint.fY,endpoint.fZ);
        h2a->Fill(p3D.Z(),p3D.Phi(),rfE);
        if(gdebug)cout << "i= " << i << "  index ="<< index[i] << " "<< p3D.Z() << " " << p3D.Rho() << " " << p3D.Phi() << " " << rfE << " " << endpoint.fT << endl;
        h2a->Draw("lego2 0");
        pt->Clear();
        s = "t = " + std::to_string(endpoint.fT);
        pt->AddText(s.c_str());
        pt->Draw();
        c2->Update();

//_start track-by-track addition top track plot

//_start FindTrack
     if(!debug2){
       track = ctracks->FindTrackByLabel(itrk);
     }else{
       LList_i iel;
       track=0;
       for (LList_i iel=ctracks->BeginChildren(); iel!=ctracks->EndChildren(); ++iel)
         if (((TEveTrack*)(*iel))->GetLabel() == itrk) {
           
            track=(TEveTrack*)(*iel);
            break;
         }
     }
//_end FindTrack

     if(track){
        track->SetRnrSelf(kTRUE);
        cout<< "++ itrk="<<itrk<<endl;
        sb->AddElement(track);          // update scene top track plot
        sb->Changed();
        sb->Repaint(kTRUE);
        //gEve->GetScenes()->RepaintChangedScenes(kFALSE);
        //gEve->GetScenes()->RepaintAllScenes(kTRUE);
        //gSystem->ProcessEvents();
     }

//_end track-by-track addition top track plot

        gSystem->ProcessEvents();
        gSystem->Sleep(kSleep); // in milliseconds
     if(gdebug)cout<< "t-loop "<< i<<endl;
     }

//_end create lego/track animation (on TCanvas)

//_end create two horizontal slots in bottom tab #2
//_end create two vertical slots in tab #2
//_end create two-tabbed TEve display

     gEve->GetBrowser()->GetTabRight()->SetTab(2);
}

//______________________________________________________________________________
//__________                                                           _________
//__________    Following taken from tutorials/eve/calorimeters.C      _________
//______________________________________________________________________________

//______________________________________________________________________________

void MakeViewerScene(TEveWindowSlot* slot, TEveViewer*& v, TEveScene*& s)
{
   // Create a scene and a viewer in the given slot.

   v = new TEveViewer("Viewer lego");
   v->SpawnGLViewer(gEve->GetEditor());
   slot->ReplaceWindow(v);
   gEve->GetViewers()->AddElement(v);
   s = gEve->SpawnNewScene("Scene lego");
   v->AddScene(s);
}

//______________________________________________________________________________
TEveCaloLego* MakeCaloLego(TEveCaloData* data, TEveWindowSlot* slot)
{
   // Eta-phi lego view.

   TEveViewer* v;
   TEveScene* s;
   if (slot) {
      MakeViewerScene(slot, v, s);
   } else {
      v = gEve->GetDefaultViewer();
      s = gEve->GetEventScene();
   }
   v->SetElementName("Final E Distribution");
   //v->SetElementName("Viewer - Lego");
   //s->SetElementName("Scene - Lego");

   auto lego = new TEveCaloLego(data);
   s->AddElement(lego);

   // By the default lego extends is (1x1x1). Resize it to put in 'natural'
   // coordinates, so that y extend in 2*Pi and set height of lego two times
   //  smaller than y extend to have better view in 3D perspective.

   lego->InitMainTrans();                                 // Set rotation and scaling for lego plot
   lego->RefMainTrans().SetRotByAngles( 6.8 , 0. , 0. );
   lego->RefMainTrans().SetScale( 0.1, 10.0 , 7.0 );

   auto glv = v->GetGLViewer();

   // draws scales and axis on borders of window
   TEveCaloLegoOverlay* overlay = new TEveCaloLegoOverlay();
   glv->AddOverlayElement(overlay);
   overlay->SetCaloLego(lego);

   // set event handler to move from perspective to orthographic view.
   //glv->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   //glv->SetCurrentCamera(TGLViewer::kCameraOrthoXnOZ);
   //glv->SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
   glv->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
   Double_t kfov       = 0.85; 
   Double_t kdolly     = 0.;
   Double_t kcenter[3] = {2.,0., 2.};                       // GL camera center (lego plot)
   Double_t khrotate   = TMath::Pi();
   Double_t kvrotate   = TMath::Pi();
   khrotate*= 0.33;   // ~1/3*Pi
   kvrotate*= 1.50;   //  3/2*Pi
   glv->CurrentCamera().Configure( kfov, kdolly, kcenter, khrotate, kvrotate);
   glv->RequestDraw(TGLRnrCtx::kLODHigh);

   glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, lego));
   gEve->AddToListTree(lego, kTRUE);
   gEve->FullRedraw3D(kFALSE, kFALSE);

   return lego;
}
