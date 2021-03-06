#ifndef SIGNAL_MAKER_H
#define SIGNAL_MAKER_H


#include "HistoAnalyzer.h"
#include "vendor/loguru.h"


class SignalMaker : public HistoAnalyzer {
protected:
	HistoBins massBins;
public:

	SignalMaker(){}
	~SignalMaker(){}


	virtual void initialize(){
		LOG_F( INFO, "" );
		book->cd();

		massBins.load( config, "bins.mass" );
	}


	float integral( TH1 * h ){
		float I = 0;
		for ( size_t i = 1; i <= h->GetXaxis()->GetNbins(); i++ ){
			I += h->GetBinContent(i) * h->GetBinWidth(i);
		}
		return I;
	}


	float linterp( float x, float x0, float y0, float x1, float y1 ){
		return y0 + (x - x0) * (y1 - y0) / ( x1 - x0 );
	}

	float corr_factor( float x, TH1 * hbgshape ){

		TAxis * xA = hbgshape->GetXaxis();
		int iBin = xA->FindBin( x );
		float x0 = xA->GetBinCenter( iBin );
		float y0 = hbgshape->GetBinContent( iBin );
		int iBin1 = iBin;
		if ( 1 == iBin && x < x0 ){
			return y0;
		} else if ( xA->GetNbins() == iBin && x > x0 ){
			return y0;
		} else {
			if ( x < x0 )
				iBin1 = iBin - 1;
			else 
				iBin1 = iBin + 1;
		}

		float x1 = xA->GetBinCenter( iBin1 );
		float y1 = hbgshape->GetBinContent( iBin1 );

		return linterp( x, x0, y0, x1, y1 );

	}

	virtual void make(){
		LOG_F( INFO, "" );

		gStyle->SetOptStat(0);

		TH1 * hbg = get<TH1>( "hmw_1_301", "mixed_ls" );
		TH1 * hbgshape = get<TH1>( "tmpfactor", "d" );
		TH1 * hls = get<TH1>( "ls", "d" );
		TH1 * huls = get<TH1>( "uls", "d" );
		TH1 * hmuls = get<TH1>( config.getString("name:mixed_uls", "hmw_1_301"), "mixed_uls" );
		LOG_F( INFO, "like_sign base = %p", hls );
		LOG_F( INFO, "raw signal = %p", huls );
		LOG_F( INFO, "corr factor = %p", hbgshape );

		assert( hbg != 0 );
		assert( hbgshape != 0 );
		assert( hls != 0 );
		assert( huls != 0 );



		
		hbg->Scale( hls->Integral() / hbg->Integral()  );
		// hbg->Scale( hls->GetXaxis()->GetBinWidth( 1 ) / hbg->GetXaxis()->GetBinWidth( 1 )  );

		LOG_F( INFO, "Mass Bins: %s", massBins.toString().c_str() );
		// LOG_F( INFO, "Rebin bgshape" );
		// TH1 * hbgshaperb = hbgshape->Rebin( massBins.nBins()-5, "bgshaperb", massBins.getBins().data() );
		
		// hbgshaperb->Scale( hbgshape->GetXaxis()->GetBinWidth(1), "width" );


		// if ( nullptr == hbgshaperb ){
		// 	LOG_F( INFO, "Abort Rebin bgshape" );
		// 	hbgshaperb = hbgshape;
		// }

		LOG_F( INFO, "Rebin uls" );
		TH1 * hulsrb = huls->Rebin( massBins.nBins(), "ulsrb", massBins.getBins().data() );
		hulsrb->Scale(1.0, "width");
		hbg->Scale(1.0, "width");
		hls->Scale(1.0, "width");
		huls->Scale(1.0, "width");

		TH1 * hcbg = (TH1*)hbg->Clone("cbg");
		hcbg->Reset();

		for ( size_t i = 1; i <= hcbg->GetXaxis()->GetNbins(); i++ ){

			float v = hbg->GetBinContent( i );
			float x = hbg->GetBinCenter( i );

			int iBin = hbgshape->GetXaxis()->FindBin( x );

			float cfactor = corr_factor( x, hbgshape );
			//hbgshape->GetBinContent( hbgshape->GetXaxis()->FindBin( x ) );

			hcbg->SetBinContent( i, v * cfactor );
			hcbg->SetBinError( i, sqrt( v * cfactor ) );

		}

		LOG_F( INFO, "Rebin corbg" );
		TH1 * hcbgrb = hcbg->Rebin( massBins.nBins(), "cbgrb", massBins.getBins().data() );
		hcbgrb->Scale( hcbg->GetXaxis()->GetBinWidth( 1), "width" );


		TH1 * hmulsrb = (TH1*)hmuls->Rebin( massBins.nBins(), "mulsrb", massBins.getBins().data() );
		

		float normx1 = config.get<float>( nodePath + ".normRange:min" );
		float normx2 = config.get<float>( nodePath + ".normRange:max" );
		TAxis * xcbg = hcbg->GetXaxis();
		TAxis * xmuls = hmuls->GetXaxis();
		LOG_F( INFO, "Normaliz in %f < M < %f", normx1, normx2 );
		float factor = hcbg->Integral( xcbg->FindBin( normx1 ), xcbg->FindBin( normx2 )  ) / hmuls->Integral( xmuls->FindBin( normx1 ), xmuls->FindBin( normx2 ) );
		// hmulsrb->Scale(  );
		hmuls->Scale( factor );
		

		LOG_F( INFO, "factor = %f", factor );
	
		hmulsrb->Scale( 1.0, "width" );
		hmulsrb->Scale( (integral(hmuls)) / integral(hmulsrb) );

		Reporter rp( "rpSignal.pdf", 1200, 900 );
		
		rp.newPage();
		rp.margins( 0.05, 0.02, 0.13, 0.15 );

		hbgshape->Draw();
		// hbgshaperb->Draw("same h");

		rp.next();

		huls->Draw();
		hulsrb->Draw("same");
		hls->Draw("same pe");
		hbg->Draw("same");
		hcbg->SetLineColor(kRed);
		hcbg->Draw("same h ");
		hmuls->SetLineColor( kBlack );
		hmuls->Draw("same");
		hmulsrb->SetLineColor( kBlack );
		hmulsrb->Draw("same");

		hcbgrb->SetLineColor(kRed);
		hcbgrb->Draw("same h ");

		gPad->SetLogy(1);
		rp.next();

		TH1 * hsigrb = (TH1*)hulsrb->Clone("hsigrb");
		hsigrb->Add( hcbgrb, -1 );

		hsigrb->Draw("hpe");

		TH1 * hsig = (TH1*)huls->Clone("hsig");
		hsig->Add( hcbg, -1 );

		// hsig->Draw("same hpe");

		rp.next();

	}
};

#endif	