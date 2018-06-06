#ifndef BACKGROUND_MAKER_H
#define BACKGROUND_MAKER_H


#include "HistoAnalyzer.h"
#include "vendor/loguru.h"


class BackgroundMaker : public HistoAnalyzer {
protected:
	map<string, HistoBins> bins;
public:

	BackgroundMaker(){}
	~BackgroundMaker(){}


	virtual void initialize(){
		LOG_F( INFO, "" );
		book->cd();

		vector<string> paths = config.childrenOf( "bins" );
		for ( auto p : paths ){
			string tn = config.tagName(p);
			LOG_F(  INFO, "Loading Bins: %s", tn.c_str() );
			bins[ tn ].load( config, p );
			LOG_F( INFO, "%s", bins[tn].toString().c_str() );
		}

		book->cd();
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

	size_t how_many_bins( TH1 * h, double x1, double x2 ){

		size_t n = 0;
		for ( size_t i = 1; i <= h->GetXaxis()->GetNbins(); i++ ){

			double x = h->GetBinCenter( i );

			if ( x > x1 && x < x2 )
				n++;

			if ( x > x2 ) break;
		}

		return n;

	}

	void avg_binning( TH1 * hsmall, TH1 * hrb ){

		TAxis * ax = hrb->GetXaxis();
		for ( size_t i = 1; i <= ax->GetNbins(); i++ ){

			double x1 = ax->GetBinLowEdge( i );
			double x2 = ax->GetBinUpEdge( i );
			size_t n_old_bins = how_many_bins( hsmall, x1, x2 );
			// LOG_F( INFO, "Counting # of bins in: (%f->%f) = %lu", x1, x2, n_old_bins );

			// rescale the content by the factor
			float v = hrb->GetBinContent( i );
			if ( 0 == n_old_bins ) continue;


			hrb->SetBinContent( i, v / (double)n_old_bins );
		}
	}

	virtual void make(){
		LOG_F( INFO, "" );
		gStyle->SetOptStat(0);


		TH1 * hraw_uls = get<TH1>( "uls", "data" );
		TH1 * hraw_ls  = get<TH1>( "ls", "data" );

		TH1 * hraw_mix  = get<TH1>( config.get<string>("p.mixedName"), "mixed" );

		TH1 * hraw_corr = get<TH1>( "R_mass", "corr" );

		assert( hraw_uls );
		assert( hraw_ls );
		assert( hraw_mix );
		assert( hraw_corr );

		hraw_mix->Scale( 1.0/5.0 );

		TH1 * hrb_uls  = HistoBins::rebin1D( "rb_uls" , hraw_uls , bins["mass"] );
		TH1 * hrb_ls   = HistoBins::rebin1D( "rb_ls"  , hraw_ls  , bins["mass"] );
		TH1 * hrb_mix  = HistoBins::rebin1D( "rb_mix" , hraw_mix , bins["mass"] );
		TH1 * hrb_corr = HistoBins::rebin1D( "rb_corr", hraw_corr, bins["corrmass"] );

		hrb_uls->Scale( 1.0, "width" );
		hrb_ls ->Scale( 1.0, "width" );
		hrb_mix->Scale( 1.0, "width" );

		avg_binning( hraw_corr, hrb_corr );

		TH1 * hraw_bg = (TH1*)hraw_mix->Clone( "raw_bg" );

		for ( size_t i = 1; i < hraw_bg->GetXaxis()->GetNbins(); i++ ){
			float x  = hraw_bg->GetBinCenter( i );
			float cf = corr_factor( x, hrb_corr );
			float v = hraw_bg->GetBinContent( i );
			hraw_bg->SetBinContent( i, v * cf );
		}

		TH1 * hrb_bg = HistoBins::rebin1D( "rb_bg", hraw_bg, bins["mass"] );
		hrb_bg->Scale( 1.0, "width" );

		Reporter rp( "rpBackground-" + config["sys"] + "-" + config["mod"] + config["debug"] + ".pdf", 1900, 1200 );
		rp.newPage();
		rp.margins( 0.05, 0.02, 0.13, 0.12 );

		RooPlotLib rpl;
		rpl.link( book );

		rpl.style( hrb_ls ).set( config, "style.ls" ).draw();
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");

		rp.next();

		rpl.style( hraw_corr ).set( config, "style.corr" ).draw();

		rp.next();

		rpl.style( hrb_corr ).set( config, "style.corr" ).draw();

		rp.next();

		hraw_mix->Scale(1.0, "width" );
		hraw_bg->Scale(1.0, "width" );

		vector<float> legPos = config.getFloatVector( "TLegend[0]:pos" );
		TLegend *leg0 = new TLegend( legPos[0], legPos[1], legPos[2], legPos[3] );
		leg0->SetBorderSize(0);
		leg0->SetNColumns( config.get<int>( "TLegend[0]:ncol", 2 ) );
		leg0->SetTextSize( config.get<float>( "TLegend[0]:point" ) / 360.0 );
		leg0->AddEntry( hraw_mix, "Mixed LS", "lp" );
		leg0->AddEntry( hraw_bg, "Corr BG", "lp" );




		rpl.style( hraw_bg ).set( config, "style.bg" ).draw();
		rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hrb_ls ).set( config, "style.ls" ).draw("same");
		
		leg0->Draw("same");

		TLatex lx;
		lx.SetTextSize( 14 / 360.0 );
		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );

		rp.next();

		hraw_uls->Scale( 1.0, "width" );

		rpl.style( hraw_uls ).set( config, "style.uls" ).draw();
		rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hrb_ls ).set( config, "style.ls" ).draw("same");
		
		leg0->AddEntry( hraw_uls, "Data #mu^{+}#mu^{-}", "lp" );

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );


		leg0->Draw("same");

		rp.next();


		rpl.style( hraw_uls ).set( config, "style.uls" ).draw();
		rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hrb_ls ).set( config, "style.ls" ).draw("same");
		gPad->SetLogy(1);
		leg0->Draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );

		rp.next();

		TH1 * hraw_sig = (TH1*)hraw_uls->Clone( "hraw_sig" );
		TH1 * hrb_sig  = (TH1*)hrb_uls->Clone( "hrb_sig" );

		hraw_sig->Add( hraw_bg, -1 );
		hrb_sig->Add( hrb_bg, -1 );

		rpl.style( hraw_uls ).set( config, "style.uls" ).draw();
		rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hraw_sig ).set( config, "style.sig" ).draw("same");
		rpl.style( hrb_sig ).set( config, "style.sig" ).draw("same");
		gPad->SetLogy(0);

		leg0->AddEntry( hrb_sig, "Signal", "lp" );
		leg0->Draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );

		rp.next();

		rpl.style( hraw_uls ).set( config, "style.uls" ).draw();
		rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hraw_sig ).set( config, "style.sig" ).draw("same");
		rpl.style( hrb_sig ).set( config, "style.sig" ).draw("same");
		gPad->SetLogy(1);

		leg0->Draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );

		rp.next();



	}
};

#endif	