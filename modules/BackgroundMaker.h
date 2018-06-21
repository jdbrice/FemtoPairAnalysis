#ifndef BACKGROUND_MAKER_H
#define BACKGROUND_MAKER_H


#include "HistoAnalyzer.h"
#include "vendor/loguru.h"


class BackgroundMaker : public HistoAnalyzer {
protected:
	map<string, HistoBins> bins;

	//parts of the cocktail
	map<string,TH2 * > h2d;
	map<string,TH1 * > h1d;
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


	virtual void make_cocktail(){
		// get the 2D histograms
		RooPlotLib rpl;
		float minPt = config.get<float>("p.minPt");

		vector<string> components = { "eta_mumu", "ccbar_mumu", "jpsi_mumu", "psi_mumu", "omega_mumu", "omega_pi0mumu_dalitz", "rho_mumu", "phi_mumu", "eta_gammamumu_dalitz", "etaprime_gammamumu_dalitz" };
		for ( string n : components ){
			h2d[ n ] = get<TH2>( "Scaled_pT_" + n );
			LOG_F( INFO, "Cocktail: %s=%p", n.c_str(), h2d[n]  );
			if ( nullptr == h2d[n] ) continue;

			int ipt1 = h2d[n]->GetYaxis()->FindBin( minPt );
			LOG_F( INFO, "minPt = %f, bin = %d", minPt, ipt1 );
			h1d[ n ] = h2d[ n ]->ProjectionX( (n + "_mass").c_str(), ipt1, -1 );
			if ( n == "eta_mumu" ){
				h1d[ "sum" ] = (TH1*)h1d[ n ]->Clone( "cocktail_sum" );
				h1d[ "sum" ]->Reset();
			}
			if ( n != "ccbar_mumu" ){
				h1d[ n ]->Scale( 3.56e11 * config.get<float>( "p.eff", 1.0 ), "width" );
			} else {
				h1d[ n ]->Scale( 3.56e11 * config.get<float>( "p.eff", 1.0 ) );
			}
			h1d[ "sum" ]->Add( h1d[ n ] );

			// rpl.style( h1d[ n ] ).set( config, "style.cocktail" ).set( config, "style.cocktail_" + n ).draw("same");
		}

		// rpl.style( h1d["sum"] ).set( config, "style.cocktail" ).set( config, "style.cocktail_sum" ).draw( "same" );

	}

	virtual void draw_cocktail(){
		// get the 2D histograms
		RooPlotLib rpl;
		vector<string> components = { "eta_mumu", "ccbar_mumu", "jpsi_mumu", "psi_mumu", "omega_mumu", "omega_pi0mumu_dalitz", "rho_mumu", "phi_mumu", "eta_gammamumu_dalitz", "etaprime_gammamumu_dalitz" };
		for ( string n : components ){
			if ( nullptr == h1d[ n ] ) continue;
			rpl.style( h1d[ n ] ).set( config, "style.cocktail.common" ).set( config, "style.cocktail." + n ).draw("same");
		}
		rpl.style( h1d["sum"] ).set( config, "style.cocktail.common" ).set( config, "style.cocktail.sum" ).draw( "same" );
	}

	virtual void make(){
		LOG_F( INFO, "" );
		gStyle->SetOptStat(0);


		TH1 * hraw_uls = get<TH1>( "uls", "data" );
		TH1 * hraw_ls  = get<TH1>( "ls", "data" );

		TH1 * hraw_mix  = get<TH1>( config.get<string>("p.mixedName"), "mixed" );

		TH1 * hraw_corr = get<TH1>( "R_mass", "corr" );

		hraw_uls->Sumw2();
		hraw_ls->Sumw2();

		assert( hraw_uls );
		assert( hraw_ls );
		assert( hraw_mix );
		assert( hraw_corr );

		hraw_mix->Scale( 1.0/10.0 );

		TH1 * hrb_uls  = HistoBins::rebin1D( "rb_uls" , hraw_uls , bins["mass"] );
		TH1 * hrb_ls   = HistoBins::rebin1D( "rb_ls"  , hraw_ls  , bins["mass"] );
		TH1 * hrb_mix  = HistoBins::rebin1D( "rb_mix" , hraw_mix , bins["mass"] );
		TH1 * hrb_corr = HistoBins::rebin1D( "rb_corr", hraw_corr, bins["corrmass"] );

		hrb_uls->Scale( 1.0, "width" );
		hrb_ls ->Scale( 1.0, "width" );
		hrb_mix->Scale( 1.0, "width" );

		avg_binning( hraw_corr, hrb_corr );

		TH1 * hraw_bg = nullptr;
		if ( config.get<string>( "p.bgsource", "mix" ) == "ls" ){
			hraw_bg = (TH1*)hraw_ls->Clone( "raw_bg" );
		} else {
			hraw_bg = (TH1*)hraw_mix->Clone( "raw_bg" );
		}

		for ( size_t i = 1; i < hraw_bg->GetXaxis()->GetNbins(); i++ ){
			float x  = hraw_bg->GetBinCenter( i );
			float cf = corr_factor( x, hrb_corr );
			float v = hraw_bg->GetBinContent( i );
			hraw_bg->SetBinContent( i, v * cf );
		}

		

		TH1 * hrb_bg = HistoBins::rebin1D( "rb_bg", hraw_bg, bins["mass"] );
		hrb_bg->Scale( 1.0, "width" );
		
		hraw_bg->Scale( config.get<float>("p.bgscaler", 1.0) );
		hrb_bg->Scale( config.get<float>("p.bgscaler", 1.0) );


		Reporter rp( config[nodePath + ".output.Report:url" ] , 1900, 1200 );
		rp.newPage();
		rp.margins( 0.05, 0.02, 0.13, 0.12 );

		RooPlotLib rpl;
		rpl.link( book );
		rpl.link( &config );


		hraw_ls->Scale( 1.0, "width" );
		hraw_mix->Scale( 1.0, "width" );

		/*************************************************************************************************************************
		* Page 1
		* Like-Sign vs. Mixed Like-Sign 
		* ******************************/
		rpl.style( hraw_ls ).set( config, "style.ls" ).draw();
		rpl.style( hrb_ls ).set( config, "style.ls" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");

		rp.next();

		/*************************************************************************************************************************
		* Page 2
		* Raw binning correction factor
		* ******************************/
		rpl.style( hraw_corr ).set( config, "style.corr" ).draw();
		rp.next();

		/*************************************************************************************************************************
		* Page 3
		* Rebinned correction factor
		* ******************************/
		rpl.style( hrb_corr ).set( config, "style.corr" ).draw();
		rp.next();

		hraw_bg->Scale(1.0, "width" );

		vector<float> legPos = config.getFloatVector( "TLegend[0]:pos" );
		TLegend *leg0 = new TLegend( legPos[0], legPos[1], legPos[2], legPos[3] );
		leg0->SetBorderSize(0);
		leg0->SetNColumns( config.get<int>( "TLegend[0]:ncol", 2 ) );
		leg0->SetTextSize( config.get<float>( "TLegend[0]:point" ) / 360.0 );
		leg0->AddEntry( hraw_mix, "Mixed LS", "lp" );
		leg0->AddEntry( hraw_bg, "Corr BG", "lp" );


		/*************************************************************************************************************************
		* Page 4
		* Mixed + correction factor applied BG
		* ******************************/

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

		/*************************************************************************************************************************
		* Page 5
		* Unlike-Sign + Background
		* ******************************/

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

		/*************************************************************************************************************************
		* Page 5b
		* Unlike-Sign + Background LMR
		* ******************************/

		rpl.style( hraw_uls ).set( config, "style.uls" ).set( "x-range", 0.2, 1.2 ).set( "grid-y", 1 ).draw();
		rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hrb_ls ).set( config, "style.ls" ).draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );


		leg0->Draw("same");

		rp.next();

		/*************************************************************************************************************************
		* Page 6
		* Unlike-Sign + Background semilogy
		* ******************************/

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


		/*************************************************************************************************************************
		* Page 7
		* Signal linear y
		* ******************************/
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


		/*************************************************************************************************************************
		* Page 8
		* Signal semi-log y
		* ******************************/
		rpl.style( hraw_uls ).set( config, "style.uls" ).set( "min", 1 ).draw();
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


		/*************************************************************************************************************************
		* Page 9
		* Signal / Background
		* ******************************/

		TH1 * hraw_SoverB = (TH1*)hraw_sig->Clone( "raw_SoverB" );
		TH1 * hrb_SoverB  = (TH1*)hrb_sig->Clone( "rb_SoverB" );

		hraw_SoverB->Divide( hraw_bg );
		hrb_SoverB->Divide( hrb_bg );

		rpl.style( hrb_SoverB ).set( "style.SoverB" ).draw();
		// rpl.style( hraw_uls ).set( config, "style.uls" ).set( "min", 1 ).draw();
		// rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		// rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		// rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		// rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// // rpl.style( hraw_sig ).set( config, "style.sig" ).draw("same");
		// rpl.style( hrb_sig ).set( config, "style.sig" ).draw("same");
		// gPad->SetLogy(1);

		// leg0->Draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );

		rp.next();


		/*************************************************************************************************************************
		* Page 10
		* Signal + Cocktail
		* ******************************/
		rpl.style( hraw_uls ).set( config, "style.uls" ).set( "min", 10 ).draw();
		rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hraw_sig ).set( config, "style.sig" ).draw("same");
		rpl.style( hrb_sig ).set( config, "style.sig" ).draw("same");

		make_cocktail();
		draw_cocktail();

		gPad->SetLogy(1);

		leg0->Draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );
		rp.next();


		/*************************************************************************************************************************
		* Page 11
		* Signal + Cocktail
		* ******************************/
		// rpl.style( hraw_uls ).set( config, "style.uls" ).set( "min", 10 ).draw();
		// rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		// rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		// rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		// rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hraw_sig ).set( config, "style.sig" ).draw("same");
		rpl.style( hrb_sig ).set( config, "style.sig" ).set("min", 1).set("draw", "hpe").draw();

		// make_cocktail();
		draw_cocktail();

		rpl.style( hrb_sig ).set( config, "style.sig" ).set("min", 10).draw("same");

		gPad->SetLogy(1);

		// leg0->Draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );
		rp.next();


		/*************************************************************************************************************************
		* Page 12
		* Signal + Cocktail: LMR
		* ******************************/
		// rpl.style( hraw_uls ).set( config, "style.uls" ).set( "min", 10 ).draw();
		// rpl.style( hrb_uls ).set( config, "style.uls" ).draw("same");
		// rpl.style( hraw_bg ).set( config, "style.bg" ).draw("same");
		// rpl.style( hrb_bg ).set( config, "style.bg" ).draw("same");
		// rpl.style( hraw_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hrb_mix ).set( config, "style.mixed" ).draw("same");
		// rpl.style( hraw_sig ).set( config, "style.sig" ).draw("same");
		rpl.style( hrb_sig ).set( config, "style.sig" ).set("min", 1).set("draw", "hpe").set("x-range", 0, 1.2).set("max", 5e5).draw();

		// make_cocktail();
		draw_cocktail();

		rpl.style( hrb_sig ).set( config, "style.sig" ).set("min", 10).draw("same");

		gPad->SetLogy(1);

		// leg0->Draw("same");

		lx.DrawLatexNDC( 0.6, 0.9, TString::Format( "Run15 %s at #sqrt{s} = 200 GeV", config["sys"].c_str() ) );
		lx.DrawLatexNDC( 0.2, 0.90, "p_{T}^{#mu} > 1.1 (GeV/c), |#eta^{#mu}| < 0.5" );
		lx.DrawLatexNDC( 0.2, 0.85, "p_{T}^{#mu#mu} > 0.0 (GeV/c), |y_{#mu#mu}| < 0.5" );
		rp.next();

		hraw_uls->Write();
		hraw_ls->Write();
		hraw_corr->Write();


		/*************************************************************************************************************************
		* Page 13
		* Signal / Cocktail
		* ******************************/

		TH1 * hrb_ratio = (TH1*)hrb_sig->Clone( "hrb_ratio" );
		hrb_ratio->Reset();
		for ( int i = 1; i < hrb_sig->GetXaxis()->GetNbins(); i++ ){
			double v = hrb_sig->GetBinContent( i );
			double x = hrb_sig->GetBinCenter( i );
			double cv = h1d["sum"]->GetBinContent( h1d["sum"]->GetXaxis()->FindBin( x ) );

			if ( cv == 0 ) continue;
			if (v == 0) continue;
			hrb_ratio->SetBinContent( i, v / cv );
			LOG_F( INFO, "v=%f, x=%f, cv=%f", v, x, cv );
		}
		gPad->SetLogy(0);
		rpl.style( hrb_ratio ).set( config, "style.sig" ).set("min", 0).set("draw", "hpe").set("x-range", 0, 6).set("max", 10).draw();

		rp.next();
	}
};

#endif	