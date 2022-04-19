#include <iostream>
#define newMClowpu_cxx
#include "newMClowpu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include <set>
#include "TRandom3.h"
#include <vector>
#include <string>
#include <map>
#include <utility>
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "settings.h"
#include <TF1.h>

#include <time.h>


using namespace std;

double binmcweight(int i){
//=================================== "5-10",   "10-15",  "15-7000", "15-7000", "15-7000", "15-7000",
//vector<double> sigmas=      { 30590000000, 3578000000, 1375000000, 1375000000, 1375000000, 1375000000 }
//vector<double>   noevents=    { 3299019, 7416745, 19864163, 19864163, 19864163, 19864163}
// Sum mWeights from NANOAOD    { 45800000, 6072000, 5812983*0.01647 + 3816768*0.01642 + 5465981*0.01649 + 4768431*0.01647 = 327081.245830 }
// Lumi sumweight/sigma         {45800000/30590000000, 6072000/3578000000, 327081.245830/1375000000 }
vector<double> binnedmcweight = { 0.0014972213, 0.0016970375, 0.00023787727, 0.00023787727, 0.00023787727, 0.00023787727 };
  //  vector<double> binnedmcweight = {9272.4534, 482.42187, 69.220133, 69.220133, 69.220133, 69.220133 }; //eski degerler
return binnedmcweight[i];
}

double DPhi(double phi1, double phi2){
    double dphi_ = fabs(phi1-phi2);
    return (dphi_ <= TMath::Pi())? dphi_ : TMath::TwoPi() -dphi_;
}

double oplus(double a, double b) {
  return sqrt(a*a + b*b);
}


void newMClowpu::Loop()
{
    time_t start,end; double dif;  time (&start); // time zimbirtisi
   
    if (fChain == 0) return; //MakeClass
    Long64_t nentries = fChain->GetEntriesFast();//MakeClass
    Long64_t nbytes = 0, nb = 0; //MakeClass
    
    cout<<"nentries"<<nentries<<endl;
    //nentries = 800000;
 
    vector<TF1*>fit_crossfit ;
     for ( int k = 0; k<9; ++k)
     {
         TF1 *fiteta=(TF1*)crossfit->Get(Form("frs_%d",k));      fit_crossfit.push_back(fiteta);
     }
    
    TH2D *hot_coldzonemap = (TH2D*)hot_cold->Get("h2hot_hot_cold");
    
    //====================== root dosyamı ac
    TFile myFile("11Nisan_MCs_V5_mpf.root", "RECREATE");
    
    //CondFormat'ın icinden cektigimiz text dosyaları,degistirdiklerim var!!///
    JetCorrectorParameters *pfchs_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_V5_MC_L1FastJet_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_V5_MC_L2Relative_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_V5_MC_L3Absolute_AK4PFchs.txt");
    //JetCorrectorParameters *pfchs_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_V4_MC_L2L3Residual_AK4PFchs.txt");
    
    vector<JetCorrectorParameters> vParam_pfchs;
    vParam_pfchs.push_back(*pfchs_l1);
    vParam_pfchs.push_back(*pfchs_l2);
    vParam_pfchs.push_back(*pfchs_l3);
    //vParam_pfchs.push_back(*pfchs_l2l3res);
    FactorizedJetCorrector *pfchs_jec = new FactorizedJetCorrector(vParam_pfchs);
    
    //static const double etamax = 1.3;
    static const int netabins = 9;
    //static const double etabins[netabins+1] = {0,0.5};
    //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
    //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,3.7,4.2,4.7};
    static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7,1.3};
    //static const double etabins[netabins+1] = {0,1.3,3.2,4.7};
    const double x[9][65]=
    {
       
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-0.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3637, 5220, 5492}, // Eta_0.5-1.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2941, 3832}, // Eta_1.0-1.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2500, 2640}, // Eta_1.5-2.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684}, // Eta_2.0-2.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_2.5-3.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_3.0-3.2
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032},// Eta_3.2-4.7
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-1.3
        /*{10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_4.0-4.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032} // Eta_4.5-5.0*/
         
    };
    const int nx[9] = {64,63,58,54,48,40,40,40,64};
    //const int nx[9] = {64};
    auto *y = jp::dptrangevseta;  //settings.h 'dan bilgileri alacagin zaman ac;)
    const int ny[9] = {128,126,116,108,96,80,80,80,128};       //settings.h 'dan bilgileri alacagin zaman ac;)
    
     TH1::SetDefaultSumw2();
     TH2::SetDefaultSumw2();
        
        //std::vector<TH1F*> vectptt;
        //std::vector<TH1F*> vectgen;
        //======================tag and probe ========================
        vector<TProfile*> vpchf; vector<TProfile*> vpnhf; vector<TProfile*> vpnef; vector<TProfile*> vpcef; vector<TProfile*> vpmuf;
        vector<TProfile*> vppuf; vector<TProfile*> vphhf; vector<TProfile*> vphef;
        
        vector<TProfile*> vpchftp; vector<TProfile*> vpnhftp; vector<TProfile*> vpneftp; vector<TProfile*> vpceftp; vector<TProfile*> vpmuftp;
        vector<TProfile*> vppuftp; vector<TProfile*> vphhftp; vector<TProfile*> vpheftp;
        //butun eventlerın histogramlarını buraya tanımlıyoruz//
        vector<TH2D*>vmpf; vector<TH2D*>vmpfx;
        TH1D *hpt(nullptr), *hpt_1(nullptr), *hpt_2(nullptr), *hpt_3(nullptr),*hpt_N(nullptr), *hpt_g(nullptr); //leading jets
        std::vector<TH1F*> vectptt; std::vector<TH1F*> vectgen;
        std::vector<TH1F*> vectptt1; std::vector<TH1F*> vectptt2; std::vector<TH1F*> vectptt3; std::vector<TH1F*> vectpttN;
        
        const int n3 = 500;
        vector<double> v3(n3+1);
        for (unsigned int i = 0; i != n3+1; ++i) v3[i] = -5. + i*0.02;
        for ( int k = 0; k<netabins; ++k) //histogramı tanımladıgım yer
        {
         
         /*TH1F *hpt_0   = new TH1F(Form("hpt_0_%d",k),"pT",ny[k],&y[k][0]); vectptt.push_back (hpt_0); hpt_0->Sumw2();
         TH1F *hpt_1   = new TH1F(Form("hpt_1_%d",k),"pT",ny[k],&y[k][0]); vectptt1.push_back (hpt_1); hpt_1->Sumw2();
         TH1F *hpt_2   = new TH1F(Form("hpt_2_%d",k),"pT",ny[k],&y[k][0]); vectptt2.push_back (hpt_2); hpt_2->Sumw2();
         TH1F *hpt_3   = new TH1F(Form("hpt_3_%d",k),"pT",ny[k],&y[k][0]); vectptt3.push_back (hpt_3); hpt_3->Sumw2();
         TH1F *hpt_N   = new TH1F(Form("hpt_N_%d",k),"pT",ny[k],&y[k][0]); vectpttN.push_back (hpt_N); hpt_N->Sumw2();*/
         //TH1F *hpt_g   = new TH1F(Form("hpt_g_%d",k),"pT",nx[k],&x[k][0]);vectgen.push_back (hpt_g); hpt_g->Sumw2();
            
         /*TProfile *pchf= new TProfile(Form("pchf_%d",k),"pchf",ny[k],&y[k][0]); vpchf.push_back (pchf); pchf->Sumw2();
         TProfile *pnhf= new TProfile(Form("pnhf_%d",k),"pnhf",ny[k],&y[k][0]); vpnhf.push_back (pnhf); pnhf->Sumw2();
         TProfile *pnef= new TProfile(Form("pnef_%d",k),"pnef",ny[k],&y[k][0]); vpnef.push_back (pnef); pnef->Sumw2();
         TProfile *pcef= new TProfile(Form("pcef_%d",k),"pcef",ny[k],&y[k][0]); vpcef.push_back (pcef); pcef->Sumw2();
         TProfile *pmuf= new TProfile(Form("pmuf_%d",k),"pmuf",ny[k],&y[k][0]); vpmuf.push_back (pmuf); pmuf->Sumw2();
         TProfile *ppuf= new TProfile(Form("ppuf_%d",k),"ppuf",ny[k],&y[k][0]); vppuf.push_back (ppuf); ppuf->Sumw2();
         TProfile *phhf= new TProfile(Form("phhf_%d",k),"phhf",ny[k],&y[k][0]); vphhf.push_back (phhf); phhf->Sumw2();
         TProfile *phef= new TProfile(Form("phef_%d",k),"phef",ny[k],&y[k][0]); vphef.push_back (phef); phef->Sumw2();
         
         TProfile *pchftp= new TProfile(Form("pchftp_%d",k),"pchftp",ny[k],&y[k][0]); vpchftp.push_back (pchftp); pchftp->Sumw2();
         TProfile *pnhftp= new TProfile(Form("pnhftp_%d",k),"pnhftp",ny[k],&y[k][0]); vpnhftp.push_back (pnhftp); pnhftp->Sumw2();
         TProfile *pneftp= new TProfile(Form("pneftp_%d",k),"pneftp",ny[k],&y[k][0]); vpneftp.push_back (pneftp); pneftp->Sumw2();
         TProfile *pceftp= new TProfile(Form("pceftp_%d",k),"pceftp",ny[k],&y[k][0]); vpceftp.push_back (pceftp); pceftp->Sumw2();
         TProfile *pmuftp= new TProfile(Form("pmuftp_%d",k),"pmuftp",ny[k],&y[k][0]); vpmuftp.push_back (pmuftp); pmuftp->Sumw2();
         TProfile *ppuftp= new TProfile(Form("ppuftp_%d",k),"ppuftp",ny[k],&y[k][0]); vppuftp.push_back (ppuftp); ppuftp->Sumw2();
         TProfile *phhftp= new TProfile(Form("phhftp_%d",k),"phhftp",ny[k],&y[k][0]); vphhftp.push_back (phhftp); phhftp->Sumw2();
         TProfile *pheftp= new TProfile(Form("pheftp_%d",k),"pheftp",ny[k],&y[k][0]); vpheftp.push_back (pheftp); pheftp->Sumw2();
          */
         TH2D *h2mpf= new TH2D (Form("h2mpf_%d",k),"h2mpf",ny[k],&y[k][0],n3,&v3[0]); vmpf.push_back (h2mpf); h2mpf->Sumw2();
         TH2D *h2mpfx= new TH2D (Form("h2mpfx_%d",k),"h2mpfx",ny[k],&y[k][0],n3,&v3[0]); vmpfx.push_back (h2mpfx); h2mpfx->Sumw2();
          
        }
        
        TLorentzVector p4pf;
        TLorentzVector p4gen;
        TLorentzVector w_p4gen;
        //=================Define the variables for Correction========================
        
        float pf_jtpt[100],pf_jt_eta[100],pf_jt_phi[100];    //corr definition (disarda kalabilir)
        vector<float>pfchs_ptcorr;
        
        
        //---------------------------events dongusu icine girdik---------
        
        for (Long64_t jentry=0; jentry<nentries;jentry++)   //MakeClass
        {
            if(jentry==539355 || jentry==2069933) continue;
            Long64_t ientry = LoadTree(jentry);             //MakeClass
            if (ientry < 0) break;                          //MakeClass
            nb = fChain->GetEntry(jentry);   nbytes += nb;  //MakeClass
            cout<<"\r"<<"event number: "<<jentry<<"/ "<<nentries<<flush;
            
            int pf_njt =PFJetsCHS__;  //==========corr definition
            pfchs_ptcorr.resize(pf_njt);    //======corr definition
            
            //weight hesaplama
            double w =EvtHdr__mWeight;
            //double w0 =w;
            w *=1/binmcweight(fChain->GetTreeNumber());
          
           /* //===================== GenPt uzerinden Reweight yapiyoruzzz===========================
           // double w2= 1.;
            double w2[9]={1.,1.,1.,1.,1.,1.,1.,1.,1.};
            
            w_p4gen.SetPxPyPzE(GenJets__fCoordinates_fX[0], GenJets__fCoordinates_fY[0],GenJets__fCoordinates_fZ[0], GenJets__fCoordinates_fT[0]);
            double w_p4gen_pt= (w_p4gen.Pt()>3000. ? 3000. : w_p4gen.Pt());
            double w_pt= (GenJets__!=0 ? max(10.,w_p4gen_pt) : 10.);        //orjinalinde 10,2500 idi
            
            for ( int k = 0; k<netabins; ++k){
                w2[k]= fit_crossfit[k]->Eval(w_pt);
            }
                
             //   w *=w2;
            */
            //===================== GenPt uzerinden Reweight yapiyoruzzz===========================
            
            
            // tag and probe  jets variables
            double jet0_pt=0; double jet1_pt=0; double jet2_pt=0;
            double jet0_phi=0; double jet1_phi=0; double jet2_phi=0;
            double jet0_eta=99; double jet1_eta=99; double jet2_eta=99;
            //double dphi=0; double ptave=0; double alpha=1;
         
            //find leading jets (JEC may change ordering)
            //less then 3 jets we input -1 on the place of the index
            Int_t           jt3leads[3]; // The three leading jets
            Double_t        jtpt[PFJetsCHS__];
            Double_t        jteta[PFJetsCHS__];
            Double_t        jtphi[PFJetsCHS__];
            for (int i = 0; i<3; ++i) jt3leads[i] = -1;
            vector<int> qcdpass;
            
            // Variables for MPF
            double metsumet = PFMet__sumEt_;
            double met    = PFMet__et_;
            double metphi = PFMet__phi_;
            double mex    = met*cos(metphi);
            double mey    = met*sin(metphi);
            double metsumet1 = metsumet;
            
            //====================================RecoPt hesaplamasi basliyor===============================
            for (int j=0; j<PFJetsCHS__; j++) {
                p4pf.SetPxPyPzE(PFJetsCHS__P4__fCoordinates_fX[j],PFJetsCHS__P4__fCoordinates_fY[j],
                                PFJetsCHS__P4__fCoordinates_fZ[j],PFJetsCHS__P4__fCoordinates_fT[j]);
                
                ////////////Calculation of correction ///
                
                pf_jtpt[j]   = p4pf.Pt(); pf_jt_eta[j] = p4pf.Eta(); pf_jt_phi[j] = p4pf.Phi();
                double pf_pt_old=pf_jtpt[j];
                float rho=EvtHdr__mPFRho; double area=PFJetsCHS__area_[j]; double ptraw = pf_jtpt[j]/PFJetsCHS__cor_[j]; //undo islemi
            
                pfchs_jec->setJetEta(pf_jt_eta[j]); pfchs_jec->setJetPt(ptraw);
                pfchs_jec->setRho(rho); pfchs_jec->setJetA(area);
                pfchs_ptcorr[j]= ptraw*pfchs_jec->getCorrection();
                
                double pf_pt=pfchs_ptcorr[j]; double pf_jteta=pf_jt_eta[j]; double pf_phi=pf_jt_phi[j];
                jtpt[j] = pf_pt; jteta[j]= pf_jteta; jtphi[j]= pf_phi; //corr pt for tag and probe
                
                ///////=====================the end of correction////////
                
                //====================mpf type1met=============
                //only use jets with corrected pt>Recopt to equalize DATA and MC thresholds
                if(fabs(pf_jteta)<4.7){
                    if(pf_pt_old > 15.){
                        double dpt= -pf_pt +pf_pt_old;
                        mex += dpt * cos(pf_phi);
                        mey += dpt * sin(pf_phi);
                        metsumet1 += pf_pt -pf_pt_old;
                    }
                }
                
                //==============finding leading jets----tag and probe=================
                 if (jt3leads[0]==-1 or jtpt[jt3leads[0]]<jtpt[j]) {
                   jt3leads[2] = jt3leads[1];
                   jt3leads[1] = jt3leads[0];
                   jt3leads[0] = j;
                 } else if (jt3leads[1]==-1 or jtpt[jt3leads[1]]<jtpt[j]) {
                   jt3leads[2] = jt3leads[1];
                   jt3leads[1] = j;
                 } else if (jt3leads[2]==-1 or jtpt[jt3leads[2]]<jtpt[j]) {
                   jt3leads[2] = j;
                 }
                //==============finding leading jets----tag and probe=================
               
                
                if ( PFJetsCHS__tightID_[j] && (PFMet__et_<0.3*PFMet__sumEt_) && (PFJetsCHS__muf_[j]<0.9) && (PFJetsCHS__cemf_[j]<0.9))
                {
                    
                    for ( int k=0;k<netabins;++k){
                        
                        double etamin=etabins[k];double etamax=etabins[k+1];
                        if (k==8){etamin=etabins[0]; etamax=etabins[9];}
                        if (fabs(pf_jteta)>= etamin && fabs(pf_jteta)< etamax) {
                        //if (fabs(pf_jteta) >= etabins[k] && fabs(pf_jteta) < etabins[k+1] ) {
                            
                            int hot_coldzonebin=hot_coldzonemap->FindBin(pf_jteta,pf_phi);
                            if((hot_coldzonemap->GetBinContent(hot_coldzonebin)==0.0))
                            {
                                double w3= 0.;
                              //  w3=w*w2[k];
                            /*if((pf_pt>114) &&(pf_pt<123) && (k==0) && (0==fChain->GetTreeNumber())){cout<<endl<<jentry<<endl; }*/
                            //  if((pf_pt>346) &&(pf_pt<362) && (k==2) && (0==fChain->GetTreeNumber())){cout<<endl<<jentry<<endl; }
                                
                            /*vectptt[k]->Fill(pf_pt,w3);
                            vpchf[k]->Fill(pf_pt,PFJetsCHS__chf_[j],w);
                            vpnhf[k]->Fill(pf_pt,PFJetsCHS__nhf_[j],w);
                            vpnef[k]->Fill(pf_pt,PFJetsCHS__nemf_[j],w);
                            vpcef[k]->Fill(pf_pt,PFJetsCHS__cemf_[j],w);
                            vpmuf[k]->Fill(pf_pt,PFJetsCHS__muf_[j],w);
                            vppuf[k]->Fill(pf_pt,PFJetsCHS__betaPrime_[j],w);
                            vphhf[k]->Fill(pf_pt,PFJetsCHS__hf_hf_[j],w);
                            vphef[k]->Fill(pf_pt,PFJetsCHS__hf_phf_[j],w);
                            
                            //===================stack icin====================
                            if (j==0) vectptt1[k]->Fill(pf_pt,w);
                            if (j==1) vectptt2[k]->Fill(pf_pt,w);
                            if (j==2) vectptt3[k]->Fill(pf_pt,w);
                            if (j>2)  vectpttN[k]->Fill(pf_pt,w); */
                                
                            if(k==8) qcdpass.push_back(j);
                                
                            }//hotzone
                        }
                    }
                }
              
            } //jetler bitiyor
    
    //======================met1 and metphi=====0 calculation for mpf and mpfx
    
    double met1   = oplus(mex,mey);
    double metphi1= atan2(mey,mex);
    
    bool pass1=false;
    bool pass2=false;
    bool pass3=false;
    for(int jj=0; jj<qcdpass.size(); ++jj){
        if(jj==jt3leads[0]) pass1=true;
        if(jj==jt3leads[1]) pass2=true;
        if(jj==jt3leads[2]) pass3=true;
        
    }
    //-----------------------------------Tag and probe-----------
    //the leading indices
    
    int i0 = jt3leads[0]; int i1 = jt3leads[1]; int i2 = jt3leads[2];
    if (i0 < 0. ) return; //This should not happen
    
    if(pass1 and pass2){
        
    double ptave = (i1>=0 ? 0.5 * (jtpt[i0] + jtpt[i1]) : jtpt[i0]);
    double dphi = (i1>=0 ? DPhi(jtphi[i0], jtphi[i1]) : 0.);
    double dpt = (i1>=0 ? fabs(jtpt[i0]-jtpt[i1])/(2*ptave) : 0.999);
    // If the jetID is bad for the third jet (and the third jet is visible), we set pt3 to ptave (alpha = 1)
    double pt3 = ((i1>=0 and i2>=0 and pass3 and jtpt[i2]>15) ? ((PFJetsCHS__>=3) ? jtpt[i2] : ptave) : 0.);
    double alpha = pt3/ptave;
    

    if ( i0>=0 and jtpt[i0]>15 and fabs(jteta[i0])<1.3 ) { // First leading jet
        if (i1>=0 and jtpt[i1]>15 and fabs(jteta[i1])<1.3 ) { // Second leading jet
           if (alpha < 1.0) { //first loose alpha cut for mpf
              for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
                   int itag = jt3leads[itag_lead];
                   int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
                   double etatag = jteta[itag];
                   double etaprobe = jteta[iprobe];
                   double pttag = jtpt[itag];
                   double ptprobe = jtpt[iprobe];
                   double phiprobe = jtphi[iprobe];
                    
                    // Dijet balance
                    if ((alpha < 0.3) && (dphi > 2.7)) {
                        /*vpchftp[8]->Fill(pttag,PFJetsCHS__chf_[iprobe],w);
                        vpnhftp[8]->Fill(pttag,PFJetsCHS__nhf_[iprobe],w);
                        vpneftp[8]->Fill(pttag,PFJetsCHS__nemf_[iprobe],w);
                        vpceftp[8]->Fill(pttag,PFJetsCHS__cemf_[iprobe],w);
                        vpmuftp[8]->Fill(pttag,PFJetsCHS__muf_[iprobe],w);
                        vppuftp[8]->Fill(pttag,PFJetsCHS__betaPrime_[iprobe],w);
                        vphhftp[8]->Fill(pttag,PFJetsCHS__hf_hf_[iprobe],w);
                        vpheftp[8]->Fill(pttag,PFJetsCHS__hf_phf_[iprobe],w);
                        */
                        }
                  //define mpf and mpfx
                       double mpf    = met1*cos(metphi1 - jtphi[itag])/2;
                       double mpfx   = met1*sin(metphi1 - jtphi[itag])/2;
                  
                       mpf /= ptave;
                       mpfx /= ptave;
                       //fill mpf and mpfx
                       vmpf[8]->Fill(ptave,mpf,1);
                       vmpfx[8]->Fill(ptave,mpfx,1);
                    }
                }
            }
         }
    }
            
            
        } //event bitiyor
        
        myFile.cd();
        myFile.Write();
        myFile.Close();
        
    }
