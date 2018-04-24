#ifndef DCA_PDFS_H
#define DCA_PDFS_H


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




class DcaPdfAna : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

	float nSigmaDeltaY = 3.0;
	float nSigmaDeltaZ = 3.0;
	float nSigmaDeltaTof = 3.0;

	TH2 * signal_pdf;

public:

	const int DEBUG = 1;
	DcaPdfAna() {}
	~DcaPdfAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		nSigmaDeltaY = config.getFloat( "p.nSigmaDeltaY", 3.0 );
		nSigmaDeltaZ = config.getFloat( "p.nSigmaDeltaZ", 3.0 );
		nSigmaDeltaTof = config.getFloat( "p.nSigmaDeltaTof", 3.0 );

		LOG_F( INFO, "nSigmaDeltaY = %f", nSigmaDeltaY );
		LOG_F( INFO, "nSigmaDeltaZ = %f", nSigmaDeltaZ );
		LOG_F( INFO, "nSigmaDeltaTof = %f", nSigmaDeltaTof );

		signal_pdf = nullptr;

		book->cd();

	}
protected:

	TH1 * subsample( string name, float pt, TH1 * hin, float mean, float sigma, float factor = 1.0, size_t n = 10000 ){
		// TH1 * hinsmooth = (TH1*)hin->Clone( TString::Format( "%s_smooth", hin->GetName() ) );
		TH1 * hinsmooth = new TH1F( name.c_str(), /*TString::Format( "%s_smooth", hin->GetName() ),*/ "", 250, 0, 3.4 );
		float bw = hin->GetXaxis()->GetBinWidth(10);
		hinsmooth->Reset();


		TRandom3 r;
		r.SetSeed(0);
		for ( int i = 0; i < n; i ++ ){
			float v = hin->GetRandom();
			v = v - hin->GetMean();
			v = v / hin->GetStdDev();

			v = v * sigma;
			v = v + mean;
			v = v + r.Gaus( 0, bw * factor );
			if ( v > 3.0  )
				continue;
			hinsmooth->Fill( v );
		}

		return hinsmooth;
	}


	bool signalMass( float m ){

		if ( m > 3.0 && m < 3.2 )
			return true;
		return false;
	}

	bool bgMass( float m ){
		return (!signalMass( m ) && m > 1.0 );
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

	bool passDeltaTof( float pt, float dt ){
		double sigt = 0.1009 + 0.009116 * exp( 3.474 / pt );
		if( fabs( dt ) <= ( nSigmaDeltaTof ) * sigt ) 
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

		if ( false == passDeltaY( pair->d1_mPt, pair->d1_mDeltaY ) || false == passDeltaZ( pair->d1_mPt, pair->d1_mDeltaZ )  )
			return;
		if ( false == passDeltaY( pair->d2_mPt, pair->d2_mDeltaY ) || false == passDeltaZ( pair->d2_mPt, pair->d2_mDeltaZ ) )
			return;

		if ( false == passDeltaTof( pair->d1_mPt, pair->d1_mDeltaTimeOfFlight )  )
			return;
		if ( false == passDeltaTof( pair->d2_mPt, pair->d2_mDeltaTimeOfFlight )  )
			return;

		if ( pair->d1_mNSigmaPion < -1 || pair->d1_mNSigmaPion > 3 )
			return;
		if ( pair->d2_mNSigmaPion < -1 || pair->d2_mNSigmaPion > 3 )
			return;

		// if ( pair->d1_mDCA > 0.75 ) return;
		// if ( pair->d2_mDCA > 0.75 ) return;


		////////////////////////////////////////////////////////////////////////
		/// Like-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 != pair->mChargeSum ){
			
			if ( signalMass( lv.M() ) ){
				book->fill( "ls_jpsi_mass", lv.M() );
			}

			book->fill( "ls_mass", lv.M() );
			book->fill( "ls_dca_vs_pt", pair->d1_mPt, pair->d1_mDCA );
			book->fill( "ls_dca_vs_pt", pair->d2_mPt, pair->d2_mDCA );

			if ( bgMass( lv.M() ) ){
				book->fill( "ls_bg_mass", lv.M() );

				
			}
		} // 0 != pair->mChargeSum


		////////////////////////////////////////////////////////////////////////
		/// Opposite-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 == pair->mChargeSum ){
			book->fill( "uls_mass", lv.M() );
			if ( signalMass( lv.M() ) ){
				book->fill( "uls_jpsi_mass", lv.M() );

				book->fill( "uls_dca_vs_pt", pair->d1_mPt, pair->d1_mDCA );
				book->fill( "uls_dca_vs_pt", pair->d2_mPt, pair->d2_mDCA );
			}

			
		} // 0 == pair->mChargeSum

		return;

	} //analyzeEvent


	virtual void postMake(){
		TreeAnalyzer::postMake();
		book->cd();
		RooPlotLib rpl;
		rpl.link( book );

		gStyle->SetOptStat(0);
		gStyle->SetOptFit( 1111 );
		Reporter rp( "rpDCA.pdf", 1200, 1200 );
		
		rp.newPage();
		rp.margins( 0.05, 0.02, 0.13, 0.15 );
		TGaxis::SetMaxDigits(3);


		TH1 * hls = book->get("ls_mass");
		TH1 * hlsjpsi = book->get("ls_jpsi_mass");
		float lsScale = config.getFloat( "p.lsScale", 1.0 );
		hls->Scale( lsScale );
		hlsjpsi->Scale( lsScale );

		rpl.style( "uls_mass" ).set( config, "style.uls_mass" ).set( "yto", 1.5 ).set( "logy", 0 ).draw();
		rpl.style( "ls_mass" ).set( config, "style.ls_mass" ).set( "draw", "same" ).draw();
		rpl.style( "uls_jpsi_mass" ).set( config, "style.uls_mass" ).set( "fillcoloralpha", "#00F", "0.3" ).set("draw", "same").draw();
		rpl.style( "ls_jpsi_mass" ).set( config, "style.ls_mass" ).set( "fillcoloralpha", "#000", "0.3" ).set("draw", "same").draw();

		TLatex lx;
		lx.SetTextSize( 14 / 360.0 );
		lx.DrawLatex( 3.2, 2300, " -1 < n#sigma_{#pi} < 3" );
		lx.DrawLatex( 3.2, 2100, " dY, dZ < 3#sigma (+0.5 p_{T} > 3 GeV/c)" );
		lx.DrawLatex( 3.2, 1900, " dTof < 3#sigma" );
		// lx.DrawLatex( 3.2, 1900, " p_{T}^{leading} > 1.5 (GeV/c)" );
		lx.DrawLatex( 3.2, 1700, TString::Format( "bg scale = %0.3f", lsScale ) );

		rp.saveImage( "export/dca_jpsi.pdf" );
		rp.next();

		rp.margins( 0.1, 0.15, 0.15, 0.15 );


		TH1 * jpsils = (TH1*)book->get( "ls_jpsi_mass" );
		TH1 * mls = (TH1*)book->get( "ls_mass" );

		TH2 * dcauls = (TH2*)book->get( "uls_dca_vs_pt" );
		TH2 * dcals = (TH2*)book->get( "ls_dca_vs_pt" );
		
		LOG_F( INFO, "bg Scale = %f", lsScale );
		dcals->Scale( jpsils->Integral() / dcals->Integral()  );

		TH2 * dcadiff = (TH2*)dcauls->Clone( "dcadiff" );
		dcadiff->Add( dcals, -1 );


		// Ensure that there are no negative entries
		for ( size_t i = 1; i <= dcadiff->GetXaxis()->GetNbins(); i++ ){
			for ( size_t j = 1; j <= dcadiff->GetYaxis()->GetNbins(); j++ ){
				if ( dcadiff->GetBinContent( i, j ) < 0 )
					dcadiff->SetBinContent( i, j, 0 );
			}
		}

		dcadiff->SetTitle( "" );
		dcadiff->Draw("colz");
		rpl.style(dcadiff).set("logz", 1).set("xr", 0, 5);
		
		rp.saveImage( "export/dca_2d.pdf" );
		rp.next();
		rp.margins( 0.05, 0.02, 0.13, 0.15 );

		dcadiff->SetLineColor(kBlack);
		dcauls->SetLineColor(kBlue);
		dcals->SetLineColor(kRed);

		TH1 * hydcauls = dcauls->ProjectionY("dcaulsy");
		hydcauls->SetTitle(";DCA (cm); dN/dDCA (cm)^{-1}");
		rpl.style( hydcauls ).set( "yto", 1.3 );
		hydcauls->Draw();
		TH1 * hydcals = dcals->ProjectionY("dcalsy");
		hydcals->Draw("same");
		TH1 * hydcadiff = dcadiff->ProjectionY("dcadiffy");
		hydcadiff->Draw("same");

		TLatex lt;
		lt.SetTextSize( 12 / 360.0 );
		lt.SetTextFont(42);
		lx.DrawLatexNDC( 0.5, 0.9, TString::Format( "%0.2f < p_{T} < %0.2f (GeV/c)", 0.0, 15.0) );
		float effs = 1.0;
		float effbg = 1.0;

		effs = hydcadiff->Integral( 1, hydcadiff->GetXaxis()->FindBin(1.0) ) / hydcadiff->Integral();
		effbg = hydcals->Integral( 1, hydcals->GetXaxis()->FindBin(1.0) ) / hydcals->Integral();
		lt.DrawLatexNDC( 0.5, 0.85, TString::Format("DCA < 1.0 : #varepsilon_{sig}=%0.2f, #varepsilon_{bg}=%0.2f", effs, effbg ) );

		effs = hydcadiff->Integral( 1, hydcadiff->GetXaxis()->FindBin(0.5) ) / hydcadiff->Integral();
		effbg = hydcals->Integral( 1, hydcals->GetXaxis()->FindBin(0.5) ) / hydcals->Integral();
		lt.DrawLatexNDC( 0.5, 0.8, TString::Format("DCA < 0.5 : #varepsilon_{sig}=%0.2f, #varepsilon_{bg}=%0.2f", effs, effbg ) );

		TLegend *leg = new TLegend( 0.5, 0.7, 0.9, 0.79 );
		leg->AddEntry( hydcauls, "unlike-sign" );
		leg->AddEntry( hydcals, "like-sign" );
		leg->AddEntry( hydcadiff, "signal" );

		leg->SetNColumns(3);
		leg->SetBorderSize(0);
		leg->Draw();


		gPad->SetLogy(1);
		rp.saveImage( "export/dca_all_pt.pdf" );
		rp.next();

		delete leg;
		

		// dcabg->SetLineColor(kRed);
		for ( size_t i = 5; i <= dcadiff->GetXaxis()->GetNbins(); i++ ){
			
			TH1 * hls = dcals->ProjectionY( "pyls", i, i );
			TH1 * hsig = dcadiff->ProjectionY( "py", i, i );
			TH1 * hrsig = dcauls->ProjectionY( "rpy", i, i );

			if ( hsig->Integral() < 0 )
				continue;
			// rpl.style( hls ).set( "min", 1e-4);
			
			// hsig->SetMinimum( 1e-4 );

			book->get( "param_mean" )->SetBinContent( i, hsig->GetMean() );
			book->get( "param_mean" )->SetBinError( i, hsig->GetMeanError() );
			book->get( "param_sigma" )->SetBinContent( i, hsig->GetStdDev() );
			book->get( "param_sigma" )->SetBinError( i, hsig->GetStdDevError() );
			LOG_F( INFO, "mean=%f, sigma=%f", hsig->GetMean(), hsig->GetStdDev() ) ;

			hrsig->SetTitle("; DCA (cm); dN/dDCA (cm)^{-1}");
			rpl.style( hrsig ).set( "min", 0.1 ).set("yto", 1.3);;
			hrsig->Draw("");
			hsig->Draw( "same" );
			hls->Draw("same");

			float pt = (dcadiff->GetXaxis()->GetBinLowEdge(i)+ dcadiff->GetXaxis()->GetBinUpEdge(i)) / 2.0;
			TH1 * hPdf = subsample( "pdf_" + ts((int)i), pt, hydcadiff, hsig->GetMean(), hsig->GetStdDev(), 1.0, 500000 );
			hPdf->Scale( (hsig->Integral() * ( hsig->GetXaxis()->GetBinWidth(1) / hPdf->GetXaxis()->GetBinWidth(1) ) ) / hPdf->Integral() );
			rpl.style( hPdf ).set( "lc", "#AA0000" );

			lx.DrawLatexNDC( 0.5, 0.9, TString::Format( "%0.2f < p_{T} < %0.2f (GeV/c)", dcadiff->GetXaxis()->GetBinLowEdge(i), dcadiff->GetXaxis()->GetBinUpEdge(i) ) );
			
			effs = hsig->Integral( 1, hsig->GetXaxis()->FindBin(1.0) ) / hsig->Integral();
			effbg = hls->Integral( 1, hls->GetXaxis()->FindBin(1.0) ) / hls->Integral();
			lt.DrawLatexNDC( 0.5, 0.85, TString::Format("DCA < 1.0 : #varepsilon_{sig}=%0.2f, #varepsilon_{bg}=%0.2f", effs, effbg ) );

			effs = hsig->Integral( 1, hsig->GetXaxis()->FindBin(0.5) ) / hsig->Integral();
			effbg = hls->Integral( 1, hls->GetXaxis()->FindBin(0.5) ) / hls->Integral();
			lt.DrawLatexNDC( 0.5, 0.8, TString::Format("DCA < 0.5 : #varepsilon_{sig}=%0.2f, #varepsilon_{bg}=%0.2f", effs, effbg ) );

			TLegend *leg = new TLegend( 0.5, 0.7, 0.9, 0.79 );
			leg->AddEntry( hrsig, "unlike-sign" );
			leg->AddEntry( hls, "like-sign" );
			leg->AddEntry( hsig, "signal" );

			leg->SetNColumns(3);
			leg->SetBorderSize(0);
			leg->Draw();

			gPad->SetLogy(1);
			
			rp.saveImage( "export/dca_pt_" + ts( (int)i ) + ".pdf" );
			rp.next();
			delete leg;


			//  Now draw the sub-sampled PDF shape and the raw signal shape
			hsig->Draw();
			hPdf->Draw("same");


			hsig->SetTitle("; DCA (cm); dN/dDCA (cm)^{-1}");
			rpl.style( hsig ).set("yto", 1.3);
			lx.DrawLatexNDC( 0.5, 0.9, TString::Format( "%0.2f < p_{T} < %0.2f (GeV/c)", dcadiff->GetXaxis()->GetBinLowEdge(i), dcadiff->GetXaxis()->GetBinUpEdge(i) ) );

			effs = hsig->Integral( 1, hsig->GetXaxis()->FindBin(2.0) ) / hsig->Integral();
			effbg = hPdf->Integral( 1, hPdf->GetXaxis()->FindBin(2.0) ) / hPdf->Integral();
			lt.DrawLatexNDC( 0.5, 0.85, TString::Format("DCA < 2.0 : #varepsilon_{sig}=%0.2f, #varepsilon_{PDF}=%0.2f", effs, effbg ) );

			effs = hsig->Integral( 1, hsig->GetXaxis()->FindBin(1.0) ) / hsig->Integral();
			effbg = hPdf->Integral( 1, hPdf->GetXaxis()->FindBin(1.0) ) / hPdf->Integral();
			lt.DrawLatexNDC( 0.5, 0.8, TString::Format("DCA < 1.0 : #varepsilon_{sig}=%0.2f, #varepsilon_{PDF}=%0.2f", effs, effbg ) );

			effs = hsig->Integral( 1, hsig->GetXaxis()->FindBin(0.5) ) / hsig->Integral();
			effbg = hPdf->Integral( 1, hPdf->GetXaxis()->FindBin(0.5) ) / hPdf->Integral();
			lt.DrawLatexNDC( 0.5, 0.75, TString::Format("DCA < 0.5 : #varepsilon_{sig}=%0.2f, #varepsilon_{PDF}=%0.2f", effs, effbg ) );

			leg = new TLegend( 0.5, 0.6, 0.9, 0.74 );
			leg->AddEntry( hsig, "signal" );
			leg->AddEntry( hPdf, "template PDF" );

			leg->SetNColumns(2);
			leg->SetBorderSize(0);
			leg->Draw();

			rp.saveImage( "export/dca_pt_" + ts( (int)i ) + "_compare.pdf" );
			rp.next();
		}

		TF1 * fmean = new TF1( "fmean", "[0] + [1]*exp([2]/x)" );
		rpl.style( "param_mean" ).set("draw", "pe").set( "logy", 0 ).draw();
		book->get("param_mean")->Fit( fmean );
		// book->get( "param_mean" )->Draw();
		
		rp.saveImage( "export/dca_mean.pdf" );
		rp.next();

		TF1 * fsigma = new TF1( "fsigma", "[0] + [1]*exp([2]/x)" );
		rpl.style( "param_sigma" ).set("draw", "pe").draw();
		book->get("param_sigma")->Fit( fsigma );
		rp.saveImage( "export/dca_sigma.pdf" );
		rp.next();


		signal_pdf = (TH2*)book->get( "dca_vs_pt" );
		for ( size_t i = 1; i < signal_pdf->GetXaxis()->GetNbins(); i++ ){
			if ( signal_pdf->GetXaxis()->GetBinLowEdge( i ) < 0.5 ) continue;

			float pt = signal_pdf->GetXaxis()->GetBinCenter( i );
			float m = fmean->Eval( pt );
			float s = fsigma->Eval( pt );
			float bw = hydcadiff->GetXaxis()->GetBinWidth( 10 );

			TRandom3 r;
			r.SetSeed(0);
			for ( int i = 0; i < 100000; i ++ ){
				float v = hydcadiff->GetRandom();
				v = v - hydcadiff->GetMean();
				v = v / hydcadiff->GetStdDev();

				v = v * s;
				v = v + m;
				v = v + r.Gaus( 0, bw );
				if ( v > 3.0  )
					continue;
				signal_pdf->Fill( pt, v );
			}
		}

		signal_pdf->Draw("colz");
		rp.next();

	}





};


#endif

