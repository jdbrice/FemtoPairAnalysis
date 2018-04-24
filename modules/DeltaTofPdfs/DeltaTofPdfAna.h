#ifndef DELTA_TOF_PDFS_H
#define DELTA_TOF_PDFS_H


#include "TreeAnalyzer.h"
#include "XmlRange.h"
#include "RooPlotLib.h"
#include "Reporter.h"
#include "FitConfidence.h"


#include "vendor/loguru.h"

#include "TNamed.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TGaxis.h"


#include "FemtoDstFormat/BranchReader.h"
#include "FemtoPairFormat/FemtoPair.h"

#include <algorithm>




class DeltaTofPdfAna : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

	float nSigmaDeltaY = 3.0;
	float nSigmaDeltaZ = 3.0;

public:

	const int DEBUG = 1;
	DeltaTofPdfAna() {}
	~DeltaTofPdfAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		nSigmaDeltaY = config.getFloat( "p.nSigmaDeltaY", 3.0 );
		nSigmaDeltaZ = config.getFloat( "p.nSigmaDeltaZ", 3.0 );

	}
protected:


	bool signalMass( float m ){

		if ( m > 3.0 && m < 3.2 )
			return true;
		return false;
	}

	bool bgMass( float m ){

		if ( m < 1.5 )
			return true;
		return false;
	}

	bool passDeltaY( float pt, float dy ){
		double sigy = -17.6867 + 18.4528*exp(0.637142/pt);
		double marg = 0.;
		if(pt>3.) marg = 0.5;

		if( fabs( dy ) <= ( nSigmaDeltaY + marg ) * sigy ) 
			return true;
		
		return false;
	}

	bool passDeltaZ( float pt, float dz ){
		double sigz = -32.6793 + 32.6034 * exp( 0.444217 / pt );
		double marg = 0.;
		if(pt>3.) 
			marg = 0.5;
		if( fabs( dz ) <= ( nSigmaDeltaZ + marg ) * sigz ) 
			return true;
		return false;
	}

	virtual void analyzeEvent(){
		FemtoPair * pair = this->_fpr.get();

		TLorentzVector lv1, lv2, lv;
		lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.105 );
		lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.105 );
		lv = lv1 + lv2;

		// TRIGGER FLAG
		if ( pair->d1_mTriggerFlag <= 0 )
			return;
		if ( pair->d2_mTriggerFlag <= 0 )
			return;

		if ( false == passDeltaY( pair->d1_mPt, pair->d1_mDeltaY ) || false == passDeltaY( pair->d2_mPt, pair->d2_mDeltaY ) )
			return;
		if ( false == passDeltaZ( pair->d1_mPt, pair->d1_mDeltaZ ) || false == passDeltaZ( pair->d2_mPt, pair->d2_mDeltaZ ) )
			return;

		if ( pair->d1_mNSigmaPion < -1 || pair->d1_mNSigmaPion > 3 )
			return;
		if ( pair->d2_mNSigmaPion < -1 || pair->d2_mNSigmaPion > 3 )
			return;

		if ( pair->d1_mPt < 1.5 && pair->d2_mPt < 1.5 )
			return;


		////////////////////////////////////////////////////////////////////////
		/// Like-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 != pair->mChargeSum ){
			book->fill( "ls_mass", lv.M() );
			if ( signalMass( lv.M() ) ){
				book->fill( "ls_jpsi_mass", lv.M() );
				
				book->fill( "ls_dtof_vs_pt", pair->d1_mPt, pair->d1_mDeltaTimeOfFlight );
				book->fill( "ls_dtof_vs_pt", pair->d2_mPt, pair->d2_mDeltaTimeOfFlight );
			}

			if ( bgMass( lv.M() ) ){
				book->fill( "ls_bg_mass", lv.M() );
				
				book->fill( "bg_dtof_vs_pt", pair->d1_mPt, pair->d1_mDeltaTimeOfFlight );
				book->fill( "bg_dtof_vs_pt", pair->d2_mPt, pair->d2_mDeltaTimeOfFlight );
			}
		} // 0 != pair->mChargeSum


		////////////////////////////////////////////////////////////////////////
		/// Opposite-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 == pair->mChargeSum ){
			book->fill( "uls_mass", lv.M() );
			if ( signalMass( lv.M() ) ){
				book->fill( "uls_jpsi_mass", lv.M() );
				
				book->fill( "uls_dtof_vs_pt", pair->d1_mPt, pair->d1_mDeltaTimeOfFlight );
				book->fill( "uls_dtof_vs_pt", pair->d2_mPt, pair->d2_mDeltaTimeOfFlight );
			}

			
		} // 0 == pair->mChargeSum


	} //analyzeEvent


	virtual void postMake(){
		TreeAnalyzer::postMake();

		RooPlotLib rpl;
		rpl.link( book );

		gStyle->SetOptStat(0);
		gStyle->SetOptFit( 1111 );
		Reporter rp( "rpdeltaTof.pdf", 1200, 700 );
		
		rp.newPage();
		rp.margins( 0.05, 0.02, 0.13, 0.15 );
		TGaxis::SetMaxDigits(3);


		TH1 * hls = book->get("ls_mass");
		TH1 * hlsjpsi = book->get("ls_jpsi_mass");
		float lsScale = config.getFloat( "p.lsScale", 1.0 );
		hls->Scale( lsScale );
		hlsjpsi->Scale( lsScale );

		rpl.style( "uls_mass" ).set( config, "style.uls_mass" ).set( "yto", 1.0 ).set( "logy", 0 ).draw();
		rpl.style( "ls_mass" ).set( config, "style.ls_mass" ).set( "draw", "same" ).draw();
		rpl.style( "uls_jpsi_mass" ).set( config, "style.uls_mass" ).set( "fillcoloralpha", "#00F", "0.3" ).set("draw", "same").draw();
		rpl.style( "ls_jpsi_mass" ).set( config, "style.ls_mass" ).set( "fillcoloralpha", "#000", "0.3" ).set("draw", "same").draw();

		TLatex lx;
		lx.SetTextSize( 14 / 360.0 );
		lx.DrawLatexNDC( 0.6, 0.8, " -1 < n#sigma_{#pi} < 3" );
		lx.DrawLatexNDC( 0.6, 0.75, " dY, dZ < 3#sigma (+0.5 p_{T} > 3 GeV/c))" );
		lx.DrawLatexNDC( 0.6, 0.7, " p_{T}^{leading} > 1.5 (GeV/c)" );
		lx.DrawLatexNDC( 0.6, 0.65, TString::Format( "bg scale = %0.3f", lsScale ) );
		lx.DrawLatexNDC( 0.6, 0.60, "Signal Region : 3.0 < M < 3.2 (GeV/c^{2})" );

		TLegend * leg = new TLegend( 0.6, 0.85, 0.9, 0.95 );
		leg->AddEntry( book->get("uls_mass"), "unlike sign" );
		leg->AddEntry( hls, "like sign" );
		leg->SetBorderSize( 0 );
		leg->Draw();

		rp.saveImage( "dtof_inv_mass.pdf" );
		rp.next();
		rp.margins( 0.1, 0.15, 0.15, 0.15 );


		TH2 * dtofs = (TH2*)book->get( "uls_dtof_vs_pt" );
		TH2 * dtofbg = (TH2*)book->get( "ls_dtof_vs_pt" );
		dtofbg->Scale( lsScale );

		TH2 * dtofdiff = (TH2*)dtofs->Clone( "dtof_diff" );
		dtofdiff->Add( dtofbg, -1 );

		dtofdiff->GetYaxis()->SetRangeUser( -2, 2 );  
		dtofdiff->SetTitle( "(+-) - [(++)+(--)]" );
		dtofdiff->Draw("colz");

		rp.saveImage( "dtof_diff_2d.pdf" );
		rp.next();
		rp.margins( 0.05, 0.03, 0.13, 0.15 );

		TF1 * fgaus = new TF1( "fgaus", "gaus" );
		fgaus->SetRange( -1, 1 );
		fgaus->SetParameters( 1, 0, 1 );
		dtofdiff->FitSlicesY( fgaus, 0, -1, 2);

		TH1 * dtof_mean = (TH1*)gDirectory->Get( "dtof_diff_1" );
		TH1 * dtof_sigma = (TH1*)gDirectory->Get( "dtof_diff_2" );

		dtof_mean->GetYaxis()->SetTitle( "#mu( #Delta TOF ) (ns)" );
		dtof_sigma->GetYaxis()->SetTitle( "#sigma( #Delta TOF ) (ns)" );

		
		TF1 * fmean = new TF1( "fmean", "[0] + [1]*exp([2]/x)" );
		// dtof_mean->Fit( fmean, "SR", "", 0, 3 );
		float ptswitch = 2.2;
		TFitResultPtr fmr = dtof_mean->Fit( fmean, "SR", "", 1.0, ptswitch );
		dtof_mean->GetYaxis()->SetRangeUser( -0.3, 0.3 );
		rpl.style( dtof_mean ).set( "yto", 1.2 );
		
		lx.DrawLatex( 1, -.2, TString::Format("f = [0] + [1]*exp([2]/x) : p_{T} < %0.2f", ptswitch) );
		lx.DrawLatex( 1, -.25, TString::Format("f = %0.3f : p_{T} >= %0.2f", fmean->Eval( ptswitch ), ptswitch) );

		TF1 * fmean2 = new TF1( "fmean2", "[0]" );
		fmean2->SetParameter( 0, fmean->Eval( ptswitch ) );
		fmean2->SetRange(ptswitch, 10);
		fmean2->Draw("same");
		rp.saveImage( "dtof_mean_fit.pdf" );
		rp.next();

		TF1 * fsigma = new TF1( "fsigma", "[0] + [1]*exp([2]/x)" );

		TFitResultPtr fsr = dtof_sigma->Fit( fsigma, "SR", "", 0, 10 );
		dtof_sigma->GetYaxis()->SetRangeUser( 0, 0.5 );
		rpl.style( dtof_sigma ).set( "yto", 1.5 );

		TGraphErrors * fsge =  FitConfidence::choleskyBands( fsr, fsigma ,100, 100, nullptr, 0, 10 );
		fsge->SetFillColorAlpha( kRed, 0.3 );
		fsigma->Draw("same");
		fsge->Draw("same e3");

		lx.DrawLatex( 3, 0.25, "f = [0] + [1]*exp([2]/x)" );

		rp.saveImage( "dtof_sigma_fit.pdf" );
		rp.next();

		TH1 * hmst = dtofdiff->ProjectionY( "dtof_shape", 5, -1 );
		hmst->Scale( 1.0 / hmst->Integral() );
		hmst->SetLineColor( kRed+1 );
		hmst->DrawClone();
		for ( int i = 5; i < 8; i++ ){
			TH1 * htmp = dtofdiff->ProjectionY( "h0", i, i );
			htmp->Scale( 1.0 / htmp->Integral() );
			htmp->SetLineColor( i - 4  );
			htmp->SetTitle( "Signal Muons" );
			
			htmp->DrawClone("same");
		}

		gPad->SetLogy(1);


		rp.savePage();


		gStyle->SetOptFit(111 );
		hmst->Fit("gaus");
		rp.savePage();


		fmean->Write();
		fsigma->Write();
	}



};


#endif

