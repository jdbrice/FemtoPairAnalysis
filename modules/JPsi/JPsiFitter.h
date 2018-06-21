#ifndef JPSI_FITTER_H
#define JPSI_FITTER_H

#include "HistoAnalyzer.h"

#include "TRandom3.h"
#include "TLatex.h"
#include "TF1.h"

TF1 * fpol = nullptr;
TF1 * fg1 = nullptr;
TF1 * fg2 = nullptr;
bool rejectPsi2S = true;


Double_t evaluate(Double_t *x, Double_t *p){
	
	if ( rejectPsi2S &&  x[0] > 3.6 && x[0] < 3.8 ){
		TF1::RejectPoint();
		return 0;
	}


	fpol->SetParameters( p[0], p[1], p[2], p[3], p[4], p[5] );
	fg1->SetParameters( p[6], p[7], p[8] );

	return fpol->Eval( x[0] ) + fg1->Eval( x[0] );
}


class JPsiFitter : public HistoAnalyzer {
protected:
	TRandom3 r;
	
public:


	

	virtual void initialize(){
		HistoAnalyzer::initialize();

		fpol = new TF1( "fpol", "pol6" );
		fg1 = new TF1( "fg1", "gaus" );
		fg2 = new TF1( "fg2", "gaus" );
		

		r.SetSeed(0);
	}



	virtual void jpsiFit(){
		RooPlotLib rpl;
		rpl.link( book );
		gStyle->SetOptStat(0);

		

		book->cd();

		TH2 * huls = get<TH2>( "uls" );
		LOG_F( INFO, "mass vs pT histo : %p", huls );
		float pt1 = config.get<float>( "p.pt1", 0 );
		float pt2 = config.get<float>( "p.pt2", 1000 );

		Reporter rp( TString::Format("yield_jpsi_%0.2f_to_%0.2f.pdf", pt1, pt2).Data(), 800, 800 );
		rp.margins( 0.01, 0.05, 0.15, 0.15 );
		rp.newPage();
		

		int ipt1 = huls->GetYaxis()->FindBin( pt1 );
		int ipt2 = huls->GetYaxis()->FindBin( pt2 );

		LOG_F( INFO, "JPsi pT ( %0.3f -> %0.3f ), bins ( %d -> %d )", pt1, pt2, ipt1, ipt2 );

		TH1 * hmass = huls->ProjectionX( TString::Format( "mass_%d_to_%d", ipt1, ipt2 ), ipt1, ipt2 );

		HistoBins sigMassBins;
		sigMassBins.load( config, "bins.mass" );
		TH1 * hmassrb = hmass->Rebin( 
										sigMassBins.nBins(), 
										TString::Format( "mass_%d_to_%d_rb", ipt1, ipt2 ),
										sigMassBins.bins.data() 
									);

		hmassrb->Scale( 1.0, "width" );

		TF1 * f = new TF1( "jpsif", evaluate, 2, 5, 9 );
		rpl.style( f ).set( config, "style.fit" );
		f->SetLineColor(kBlack);
		f->SetParameters( 100, 3.1, 0.3, 4, 5, 6, 610, 3.1, 0.05 );
		f->SetParameter( 9, 100 );
		f->SetParameter( 10, 3.6 );
		f->SetParameter( 11, 0.5 );
		f->SetNpx( 500 );

		f->SetParLimits( 6, 0, 1e9 );
		f->SetParLimits( 7, 2.9, 3.3);
		f->SetParLimits( 8, 0.01, 0.1 );

		// f->SetParLimits( 9, 0, 10000 );
		// f->SetParLimits( 10, 3.5, 3.8);
		// f->SetParLimits( 11, 0, 5.0 );
		// hmassrb->GetXaxis()->SetRangeUser( 2.5, 3.4 );
		// hmassrb->SetMinimum( 0 );
		rpl.style( hmassrb ).set( config, "style.mass" ).draw();
		
		float fmax = 4.5;
		float fmin = 2.6;
		if ( pt1 < 1.0 ){
			fmax = 4.5;
			fmin = 2.6;
		}


		hmassrb->Fit( f, "RL", "", fmin, fmax );
		hmassrb->Fit( f, "RL", "", fmin, fmax );
		// hmassrb->Fit( f, "RL", "", fmin, fmax );

		LOG_F( INFO, "Mass = %0.3f +/- %0.3f", f->GetParameter( 7 ), f->GetParError( 7 ) );
		LOG_F( INFO, "Width = %0.3f +/- %0.3f", f->GetParameter( 8 ), f->GetParError( 8 ) );

		fg1->SetLineColor(kBlue);
		rpl.style( fg1 ).set( config, "style.fit" );
		fg1->SetNpx( 500 );
		fg1->SetRange( 2.5, 4.5 );
		fg1->Draw("same");

		fpol->SetLineColor(kRed);
		rpl.style( fpol ).set( config, "style.fit" );
		fpol->SetNpx( 500 );
		fpol->SetRange( 2.5, 4.5 );
		fpol->Draw("same");

		// float bw = hmassrb->GetBinWidth( 5 );
		double Ns = f->Integral( 2.9, 3.3 ) ;
		double Nbg = fpol->Integral( 2.9, 3.3 ) ;

		double Nse = f->IntegralError( 2.9, 3.3 ) ;
		double Nbge = fpol->IntegralError( 2.9, 3.3 ) ;


		// Ns *= 1.0/ hmassrb->GetBinWidth(4);
		// Nbg *= 1.0/ hmassrb->GetBinWidth(4);

		LOG_F( INFO, "Ns=%f", Ns );
		LOG_F( INFO, "Nbg=%f", Nbg );

		Ns = Ns - Nbg;
		Nse = Nse + Nbge;
		
		double sig = Ns / sqrt( Ns + Nbg );

		TLatex lx;
		lx.SetTextSize( 12.0 / 380.0);
		lx.DrawLatexNDC( .18, 0.95, TString::Format("%0.2f < p_{T} < %0.2f (GeV/c)", pt1, pt2) );
		lx.DrawLatexNDC( .18, 0.9, TString::Format("N_{J/#Psi}^{raw}=%0.3f #pm %0.3f #pm %0.3f", Ns, sqrt(Ns), Nse) );
		lx.DrawLatexNDC( .18, 0.85, TString::Format("S/B=%0.3f", Ns/Nbg) );
		lx.DrawLatexNDC( .18, 0.8, TString::Format("S/#sqrt{S + B}=%0.3f", sig) );

		lx.DrawLatexNDC( .18, 0.6, TString::Format("p0=%0.3f", fg1->GetParameter(0) ) );
		lx.DrawLatexNDC( .18, 0.55, TString::Format("#mu=%0.3f", fg1->GetParameter(1)) );
		lx.DrawLatexNDC( .18, 0.5, TString::Format("#sigma=%0.3f", fg1->GetParameter(2)) );

		rp.savePage();

		float mpt = (pt1 + pt2)/2.0;
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


		f->Write();
	}

	virtual void make(){
		LOG_F( INFO, "" );

		book->cd();
		book->makeAll( nodePath + ".histograms" );

		for ( string n : {"mass", "width", "yield"} ){
			if ( book->get(n) )
				book->get(n)->Sumw2();
		}

		HistoBins fit_bins( config, "bins.fit" );
		LOG_F( INFO, "Found %d slices to fit", fit_bins.nBins() );
		for ( size_t i = 0; i < fit_bins.nBins(); i++ ){
			float pt1 = fit_bins.getBins()[i];
			float pt2 = fit_bins.getBins()[i+1];
			LOG_F( INFO, "Setting pt=(%f, %f)", pt1, pt2 );

			config.set( "p.pt1", ts(pt1) );
			config.set( "p.pt2", ts(pt2) );

			jpsiFit();
		}
	}

	virtual void postMake(){
		
		Reporter rp( "yield_jpsi_summary.pdf", 1600, 800 );
		rp.margins( 0.03, 0.05, 0.15, 0.15 );
		rp.newPage();

		RooPlotLib rpl;
		rpl.link( book );
		rpl.link( &config );

		if ( book->get( "yield" ) ){
			book->get( "yield" )->Scale( 1.0, "width" );
			rpl.style( "yield" ).set( "style.yield" ).draw();
			gPad->SetLogy(1);
			rp.next();
		}
		gPad->SetLogy(0);
		rpl.style( "mass" ).set( "style.massfit" ).draw();
		rp.next();
		rpl.style( "width" ).set( "style.width" ).draw();
		rp.savePage();
	}



	virtual void psi2sFit(){
		RooPlotLib rpl;
		rpl.link( book );
		gStyle->SetOptStat(0);

		rejectPsi2S = false;

		book->cd();

		TH2 * huls = get<TH2>( "uls" );
		LOG_F( INFO, "mass vs pT histo : %p", huls );
		float pt1 = config.get<float>( "p.pt1", 0 );
		float pt2 = config.get<float>( "p.pt2", 1000 );

		Reporter rp( TString::Format("psi2s_%0.2f_to_%0.2f.pdf", pt1, pt2).Data(), 800, 800 );
		rp.margins( 0.01, 0.05, 0.15, 0.15 );
		rp.newPage();
		

		int ipt1 = huls->GetYaxis()->FindBin( pt1 );
		int ipt2 = huls->GetYaxis()->FindBin( pt2 );

		LOG_F( INFO, "Psi2s pT ( %0.3f -> %0.3f ), bins ( %d -> %d )", pt1, pt2, ipt1, ipt2 );

		TH1 * hmass = huls->ProjectionX( TString::Format( "mass_%d_to_%d", ipt1, ipt2 ), ipt1, ipt2 );

		hmass->Sumw2();

		HistoBins sigMassBins;
		sigMassBins.load( config, "bins.mass" );
		TH1 * hmassrb = hmass->Rebin( 
										sigMassBins.nBins(), 
										TString::Format( "mass_%d_to_%d_rb", ipt1, ipt2 ),
										sigMassBins.bins.data() 
									);

		TF1 * f = new TF1( "psi2sf", evaluate, 2, 5, 9 );
		rpl.style( f ).set( config, "style.fit" );
		f->SetLineColor(kBlack);
		f->SetParameters( 100, 3.1, 0.3, 4, 5, 6, 40, 3.6, 0.05 );
		f->SetNpx( 500 );

		f->SetParLimits( 6, 0, 10000 );
		f->SetParLimits( 7, 3.5, 3.9);
		f->SetParLimits( 8, 0.01, 1.1 );

		rpl.style( hmassrb ).set( config, "style.mass" ).draw();
		
		float fmin = 3.2;
		float fmax = 4.5;

		hmassrb->Fit( f, "RL", "", fmin, fmax );
		hmassrb->Fit( f, "RL", "", fmin, fmax );
		hmassrb->Fit( f, "RL", "", fmin, fmax );

		fg1->SetLineColor(kBlue);
		rpl.style( fg1 ).set( config, "style.fit" );
		fg1->SetNpx( 500 );
		fg1->SetRange( 3.3, 4.5 );
		fg1->Draw("same");

		fpol->SetLineColor(kRed);
		rpl.style( fpol ).set( config, "style.fit" );
		fpol->SetNpx( 500 );
		fpol->SetRange( 3.3, 4.5 );
		fpol->Draw("same");

		double Ns = f->Integral( 3.3, 4.0);
		double Nbg = fpol->Integral(3.3, 4.0 );

		Ns *= 1.0/ hmassrb->GetBinWidth(4);
		Nbg *= 1.0/ hmassrb->GetBinWidth(4);

		LOG_F( INFO, "Ns=%f", Ns );
		LOG_F( INFO, "Nbg=%f", Nbg );

		Ns = Ns - Nbg;
		
		double sig = Ns / sqrt( Ns + Nbg );

		TLatex lx;
		lx.SetTextSize( 12.0 / 380.0);
		lx.DrawLatexNDC( .18, 0.95, TString::Format("%0.2f < p_{T} < %0.2f (GeV/c)", pt1, pt2) );
		lx.DrawLatexNDC( .18, 0.9, TString::Format("N_{#Psi(2S)}^{raw}=%0.3f", Ns) );
		lx.DrawLatexNDC( .18, 0.85, TString::Format("S/B=%0.3f", Ns/Nbg) );
		lx.DrawLatexNDC( .18, 0.8, TString::Format("S/#sqrt{S + B}=%0.3f", sig) );

		lx.DrawLatexNDC( .18, 0.6, TString::Format("p0=%0.3f", fg1->GetParameter(0) ) );
		lx.DrawLatexNDC( .18, 0.55, TString::Format("#mu=%0.3f", fg1->GetParameter(1)) );
		lx.DrawLatexNDC( .18, 0.5, TString::Format("#sigma=%0.3f", fg1->GetParameter(2)) );

		rp.savePage();


		f->Write();
	}

};



#endif