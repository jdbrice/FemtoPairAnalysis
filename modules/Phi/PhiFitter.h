#ifndef PHI_FITTER_H
#define PHI_FITTER_H

#include "HistoAnalyzer.h"

#include "TRandom3.h"
#include "TLatex.h"
#include "TF1.h"

TF1 * phi_fpol = nullptr;
TF1 * phi_fg1 = nullptr;
TF1 * phi_fg2 = nullptr;

Double_t phi_evaluate(Double_t *x, Double_t *p){

	phi_fpol->SetParameters( p[0], p[1], p[2], p[3], p[4], p[5] );
	phi_fg1->SetParameters( p[6], p[7], p[8] );

	return phi_fpol->Eval( x[0] ) + phi_fg1->Eval( x[0] );
}


class PhiFitter : public HistoAnalyzer {
protected:
	TRandom3 r;
	
public:


	

	virtual void initialize(){
		HistoAnalyzer::initialize();

		phi_fpol = new TF1( "phi_fpol", "pol6" );
		phi_fg1 = new TF1( "phi_fg1", "gaus" );
		phi_fg2 = new TF1( "phi_fg2", "gaus" );
		

		r.SetSeed(0);
	}



	virtual void phiFit(){
		RooPlotLib rpl;
		rpl.link( book );
		rpl.link( &config );
		gStyle->SetOptStat(0);

		string resoStr = config.get<string>( "p.reso", "phi" );
		const char * reso = resoStr.c_str();

		const char * name = config.get<string>( "name", "" ).c_str();
		const char * title = config.get<string>( "title", "" ).c_str();

		book->cd();

		TH2 * huls = get<TH2>( "uls" );
		TH2 * hls = get<TH2>( "ls" );
		LOG_F( INFO, "mass vs pT histo : %p", huls );
		float pt1 = config.get<float>( "p.pt1", 0 );
		float pt2 = config.get<float>( "p.pt2", 1000 );

		Reporter rp( TString::Format("%s%s_%0.2f_to_%0.2f.pdf", reso, name, pt1, pt2).Data(), 800, 800 );
		rp.margins( 0.05, 0.05, 0.13, 0.15 );
		rp.newPage();
		

		int ipt1 = huls->GetYaxis()->FindBin( pt1 );
		int ipt2 = huls->GetYaxis()->FindBin( pt2 );

		LOG_F( INFO, "%s pT ( %0.3f -> %0.3f ), bins ( %d -> %d )", reso, pt1, pt2, ipt1, ipt2 );

		TH1 * hmass   = huls->ProjectionX( TString::Format( "mass_%d_to_%d", ipt1, ipt2 ), ipt1, ipt2 );
		TH1 * hmassls = hls ->ProjectionX( TString::Format( "massls_%d_to_%d", ipt1, ipt2 ), ipt1, ipt2 );

		hmass->Sumw2();
		hmassls->Sumw2();


		HistoBins sigMassBins;
		sigMassBins.load( config, "bins.mass" );
		TH1 * hmassrb = hmass->Rebin( 
										sigMassBins.nBins(), 
										TString::Format( "mass_%d_to_%d_rb", ipt1, ipt2 ),
										sigMassBins.bins.data() 
									);

		hmassrb->Scale( 1.0, "width" );

		TH1 * hmassrbls = hmassls->Rebin( 
										sigMassBins.nBins(), 
										TString::Format( "massls_%d_to_%d_rb", ipt1, ipt2 ),
										sigMassBins.bins.data() 
									);

		hmassrbls->Scale( 1.0, "width" );


		if ( config.get<bool>( "p.rmls", false ) == true ){
			hmassrb->Add( hmassrbls, -1 );
		}

		TF1 * f = new TF1( "phif", phi_evaluate, 2, 5, 9 );
		rpl.style( f ).set( config, "style.fit" );
		f->SetLineColor(kBlack);
		
		float mu = 1.013;
		if ( "omega" == resoStr ){
			mu = 0.78;
		} else if ( "k0s" == resoStr ){
			mu = 0.47;
		}

		f->SetParameters( 1.2e5, -1.2e5, -1.2e4, 3.7e4, 2.5e4, -1.9e4, 610, mu, 0.015 );
		
		f->SetNpx( 500 );

		f->SetParLimits( 6, 0, 1e9 );
		f->SetParLimits( 7, mu - 0.1, mu + 0.1);
		f->SetParLimits( 8, 0.001, 0.06 );

		rpl.style( hmassrb ).set( config, "style.mass" ).draw();
		
		float fmin = 0.85;
		float fmax = 1.5;

		if ( "omega" == resoStr ){
			fmin = 0.55;
			fmax = 0.9;
		} else if ( "k0s" == resoStr ){
			fmin = 0.2;
			fmax = 0.9;
		}
		// if ( pt1 < 1.0 ){
		// 	fmax = 4.5;
		// 	fmin = 2.6;
		// }


		// hmassrb->Fit( f, "QR", "", fmin, fmax );
		hmassrb->Fit( f, "RLE", "", fmin, fmax );

		float fmu  = f->GetParameter( 7 );
		float fsig = f->GetParameter( 8 );

		phi_fg1->SetLineColor(kBlue);
		rpl.style( phi_fg1 ).set( config, "style.fit" );
		phi_fg1->SetNpx( 500 );
		phi_fg1->SetRange( fmin, fmax );
		phi_fg1->Draw("same");

		phi_fpol->SetLineColor(kRed);
		rpl.style( phi_fpol ).set( config, "style.fit" );
		phi_fpol->SetRange( fmin, fmax );
		phi_fpol->SetNpx( 500 );
		phi_fpol->DrawClone("same");

		rpl.style( phi_fpol ).set( "style.fitext" );
		phi_fpol->SetRange( 0.2, 2.5 );
		phi_fpol->SetNpx( 500 );
		phi_fpol->Draw("same");

		

		// float bw = hmassrb->GetBinWidth( 5 );
		double Ns = f->Integral( fmu - 3*fsig, fmu + 3*fsig ) ;
		double Nbg = phi_fpol->Integral( fmu - 3*fsig, fmu + 3*fsig ) ;

		double Nse = f->IntegralError( fmu - 3*fsig, fmu + 3*fsig ) ;
		double Nbge = phi_fpol->IntegralError( fmu - 3*fsig, fmu + 3*fsig ) ;

		// Ns *= 1.0/ hmassrb->GetBinWidth(4);
		// Nbg *= 1.0/ hmassrb->GetBinWidth(4);

		LOG_F( INFO, "Ns=%f", Ns );
		LOG_F( INFO, "Nbg=%f", Nbg );


		Ns = Ns - Nbg;
		Nse = Nse + Nbge;

		double SoBe = (Ns/Nbg) * sqrt( pow( sqrt(Ns) / Ns,2) + pow( sqrt(Nbg) / Nbg, 2 ) );

		LOG_F( INFO, "S/B=%0.3e #pm %0.3e", Ns/Nbg, SoBe );

		LOG_F( INFO, "Ns-Nbg=%f", Ns );
		
		double N = Ns + Nbg;
		double Ne = sqrt( N );
		double sig = Ns / Ne;
		double sige = sig * (Nse / Ns);

		TLatex lx;
		lx.SetTextSize( 12.0 / 380.0);
		lx.DrawLatexNDC( .18, 0.9, TString::Format("%0.2f < %s < %0.2f %s", pt1, config.get<string>("p.y-title", "p_{T}").c_str(), pt2, config.get<string>("p.y-units", "(GeV/c)").c_str() ) );
		lx.DrawLatexNDC( .18, 0.85, TString::Format("N_{#phi}^{raw}=%0.3f #pm %0.3f #pm %0.3f", Ns, sqrt(Ns), Nse) );
		lx.DrawLatexNDC( .18, 0.8, TString::Format("S/B=%0.3f", Ns/Nbg) );
		lx.DrawLatexNDC( .18, 0.75, TString::Format("S/#sqrt{S + B}=%0.3f", sig) );

		lx.DrawLatexNDC( .18, 0.6, TString::Format("p0=%0.3f #pm %0.2e", f->GetParameter(6), f->GetParError(6) ) );
		lx.DrawLatexNDC( .18, 0.55, TString::Format("#mu=%0.3f #pm %0.2e", f->GetParameter(7), f->GetParError(7)) );
		lx.DrawLatexNDC( .18, 0.5, TString::Format("#sigma=%0.3f #pm %0.2e", f->GetParameter(8), f->GetParError(8)) );

		lx.DrawLatexNDC( .68, 0.9, TString::Format( "%s", title ) );

		rp.savePage();

		float mpt = (pt1 + pt2)/2.0;
		if ( true == config.get<bool>( "p.stepMinOnly", false ) )
			mpt = pt1 + 1e-5;
		if ( book->get( "yield" ) ){
			int ibin = book->get("yield")->GetXaxis()->FindBin( mpt );
			book->setBin( "yield", ibin, Ns, sqrt(Ns) );
		}
		if ( book->get( "mass" ) ){
			int ibin = book->get("mass")->GetXaxis()->FindBin( mpt );
			book->get("mass")->SetBinContent( ibin, f->GetParameter( 7 ) );
			book->get("mass")->SetBinError( ibin, f->GetParError( 7 ) );
			
		}
		if ( book->get( "width" ) ){
			int ibin = book->get("width")->GetXaxis()->FindBin( mpt );
			book->get("width")->SetBinContent( ibin, f->GetParameter( 8 ) );
			book->get("width")->SetBinError( ibin, f->GetParError( 8 ) );
			
		}
		if ( book->get( "SoverB" ) ){
			int ibin = book->get("SoverB")->GetXaxis()->FindBin( mpt );
			book->get("SoverB")->SetBinContent( ibin, Ns/Nbg );
			book->get("SoverB")->SetBinError( ibin, SoBe );
			
		}
		if ( book->get( "significance" ) ){
			int ibin = book->get("significance")->GetXaxis()->FindBin( mpt );
			book->get("significance")->SetBinContent( ibin, sig );
			book->get("significance")->SetBinError( ibin, sige );
			
		}


		f->Write();
	}

	virtual void make(){
		LOG_F( INFO, "" );

		book->cd();
		book->makeAll( nodePath + ".histograms" );

		for ( string n : {"mass", "width", "yield", "SoverB", "significance"} ){
			if ( book->get(n) )
				book->get(n)->Sumw2();
		}

		HistoBins fit_bins( config, "bins.fit" );
		LOG_F( INFO, "Found %d slices to fit", fit_bins.nBins() );
		for ( size_t i = 0; i < fit_bins.nBins(); i++ ){
			LOG_F( INFO, "fit %lu", i );
			float pt1 = fit_bins.getBins()[i];
			float pt2 = fit_bins.getBins()[i+1];
			LOG_F( INFO, "Setting fit y=(%f, %f)", pt1, pt2 );

			config.set( "p.pt1", ts(pt1) );
			if ( false == config.get<bool>( "p.stepMinOnly", false ) )
				config.set( "p.pt2", ts(pt2) );

			phiFit();
		}
	}

	virtual void postMake(){
		
		Reporter rp( "yield_phi_summary.pdf", 1600, 800 );
		rp.margins( 0.03, 0.05, 0.15, 0.15 );
		rp.newPage();

		RooPlotLib rpl;
		rpl.link( book );
		rpl.link( &config );

		if ( book->get( "yield" ) ){
			if ( true == config.get<bool>( "p.scalebw", true ) )
				book->get( "yield" )->Scale( 1.0, "width" );
			rpl.style( "yield" ).set( "style.yield" ).draw();
			gPad->SetLogy(1);
			rp.next();
		}
		gPad->SetLogy(0);
		rpl.style( "mass" ).set( "style.massfit" ).draw();
		rp.next();
		rpl.style( "width" ).set( "style.width" ).draw();
		rp.next();
		rpl.style( "SoverB" ).set( "style.SoverB" ).draw();
		rp.next();
		rpl.style( "significance" ).set( "style.significance" ).draw();
		rp.next();
	}

};



#endif