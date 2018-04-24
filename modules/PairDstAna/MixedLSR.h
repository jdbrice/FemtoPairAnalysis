#ifndef MIXED_LSR_H
#define MIXED_LSR_H


#include "HistoAnalyzer.h"
#include "vendor/loguru.h"

#include "RooPlotLib.h"
#include "FitConfidence.h"

#include "TString.h"

class MixedLSR : public HistoAnalyzer {
protected:
public:

	MixedLSR(){}
	~MixedLSR(){}


	virtual void initialize(){
		LOG_F( INFO, "" );

		book->cd();
	}


	virtual void make(){
		LOG_F( INFO, "" );
		book->makeAll( config, nodePath + ".histograms" );
		LOG_F( INFO, "Made Histos" );

		RooPlotLib rpl;


		gStyle->SetOptStat(0);
		gStyle->SetOptFit(111);
		TCanvas * can = new TCanvas( "can", "can", 1200, 900 );


		float m1 = config.get<float>("m1");
		float m2 = config.get<float>("m2");

		TH1 * hlsr = book->get("lsr");

		for ( size_t i = 0; i < 12; i++ ){
			LOG_F( INFO, "Loop %lu", i );
			TH1 * huls_mw = get<TH1>( "uls_weighted_" + ts((int)i) );
			TH1 * hls_mw = get<TH1>( "ls_weighted_" + ts((int)i) );

			TAxis * xuls = huls_mw->GetXaxis();
			TAxis * xls = hls_mw->GetXaxis();

			LOG_F( INFO, "uls_mw=%p, ls_mw=%p", huls_mw, hls_mw );

			double euls = 1;
			double els = 1;
			float vuls = huls_mw->IntegralAndError( xuls->FindBin( m1 ), xuls->FindBin( m2 ), euls );
			float vls = hls_mw->IntegralAndError( xls->FindBin( m1 ), xls->FindBin( m2 ), els );
			euls = sqrt( vuls );
			els = sqrt( vls );

			hlsr->SetBinContent( i+2,  vuls / vls );
			hlsr->SetBinError( i+2,  (vuls / vls) * sqrt( pow(euls / vuls, 2) + pow(els / vls, 2) )  );

			rpl.style( huls_mw ).set("lc", "#000" ).set("draw", "hpe").set("logy", 1).draw();
			rpl.style( hls_mw ).set("lc", "#F00" ).set("draw", "same hpe").draw();

			can->Print( TString::Format( "mlsr_%lu.pdf", i ) );

		}

		TH1 * hlsrd = book->get("lsrd");

		for ( size_t i = 0; i < 12; i++ ){
			LOG_F( INFO, "Loop %lu", i );
			TH1 * huls = get<TH1>( "uls_data_" + ts((int)i) );
			TH1 * hls = get<TH1>( "ls_data_" + ts((int)i) );

			TAxis * xuls = huls->GetXaxis();
			TAxis * xls = hls->GetXaxis();

			LOG_F( INFO, "uls=%p, ls=%p", huls, hls );

			double euls = 1;
			double els = 1;
			float vuls = huls->IntegralAndError( xuls->FindBin( m1 ), xuls->FindBin( m2 ), euls );
			float vls = hls->IntegralAndError( xls->FindBin( m1 ), xls->FindBin( m2 ), els );

			hlsrd->SetBinContent( i+2,  vuls / vls );
			hlsrd->SetBinError( i+2,  (vuls / vls) * sqrt( pow(euls / vuls, 2) + pow(els / vls, 2) )  );

			rpl.style( huls ).set("lc", "#000" ).set("draw", "hpe").set("logy", 1).draw();
			rpl.style( hls ).set("lc", "#F00" ).set("draw", "same hpe").draw();

			can->Print( TString::Format( "dmlsr_%lu.pdf", i ) );

		}

		
		
		TF1 * fitFunc = new TF1( "fitfunc", "pol1" );
		TFitResultPtr fitPointer = hlsrd->Fit( fitFunc, "SR", "", config.get<float>("fit:x1", 0), config.get<float>("fit:x2", 1.3) );
		gPad->SetLogy(0);
		// rpl.style( hlsrdd ).set( "lc", "#000" ).set("draw", "same pe").set("logy", 0).draw();

		fitFunc->SetRange( -0.1, 1.4 );
		TGraphErrors * band = FitConfidence::choleskyBands( fitPointer, fitFunc, 100, 400 );

		float v = hlsrd->GetBinContent( 5 );
		rpl.style(hlsrd).set( "yr", v/2, v*2 );

		band->SetFillColorAlpha( kRed, 0.7 );
		band->Draw( "same e3" );


		can->Print( "lsr.pdf" );


	}
};

#endif	