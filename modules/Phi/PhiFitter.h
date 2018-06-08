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
		
		float mu = 1.05;
		if ( "omega" == resoStr ){
			mu = 0.78;
		} else if ( "k0s" == resoStr ){
			mu = 0.47;
		}

		f->SetParameters( 100, 1.04, 0.3, 4, 5, 6, 610, mu, 0.05 );
		
		f->SetNpx( 500 );

		f->SetParLimits( 6, 0, 1e9 );
		f->SetParLimits( 7, mu - 0.1, mu + 0.1);
		f->SetParLimits( 8, 0.01, 0.1 );

		rpl.style( hmassrb ).set( config, "style.mass" ).draw();
		
		float fmin = 0.85;
		float fmax = 1.3;

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


		hmassrb->Fit( f, "RL", "", fmin, fmax );
		hmassrb->Fit( f, "RL", "", fmin, fmax );
		hmassrb->Fit( f, "RL", "", fmin, fmax );

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

		// Ns *= 1.0/ hmassrb->GetBinWidth(4);
		// Nbg *= 1.0/ hmassrb->GetBinWidth(4);

		LOG_F( INFO, "Ns=%f", Ns );
		LOG_F( INFO, "Nbg=%f", Nbg );


		Ns = Ns - Nbg;

		LOG_F( INFO, "Ns-Nbg=%f", Ns );
		
		double sig = Ns / sqrt( Ns + Nbg );

		TLatex lx;
		lx.SetTextSize( 12.0 / 380.0);
		lx.DrawLatexNDC( .18, 0.9, TString::Format("%0.2f < p_{T} < %0.2f (GeV/c)", pt1, pt2) );
		lx.DrawLatexNDC( .18, 0.85, TString::Format("N_{#phi}^{raw}=%0.3f", Ns) );
		lx.DrawLatexNDC( .18, 0.8, TString::Format("S/B=%0.3f", Ns/Nbg) );
		lx.DrawLatexNDC( .18, 0.75, TString::Format("S/#sqrt{S + B}=%0.3f", sig) );

		lx.DrawLatexNDC( .18, 0.6, TString::Format("p0=%0.3f", phi_fg1->GetParameter(0) ) );
		lx.DrawLatexNDC( .18, 0.55, TString::Format("#mu=%0.3f", phi_fg1->GetParameter(1)) );
		lx.DrawLatexNDC( .18, 0.5, TString::Format("#sigma=%0.3f", phi_fg1->GetParameter(2)) );

		lx.DrawLatexNDC( .68, 0.9, TString::Format( "%s", title ) );

		rp.savePage();


		f->Write();
	}

	virtual void make(){
		LOG_F( INFO, "" );

		phiFit();
	}

};



#endif