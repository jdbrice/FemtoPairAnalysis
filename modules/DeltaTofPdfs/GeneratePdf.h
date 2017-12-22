#ifndef GENERATE_PDF_H
#define GENERATE_PDF_H

#include "HistoAnalyzer.h"

#include "TRandom3.h"


class GeneratePdf : public HistoAnalyzer {
protected:
	TRandom3 r;
	TH2 * hdeltaTOF;
	TH1 * hshape;
	TF1 * fmean;
	TF1 * fsigma;
	float binwidth = 1;
public:

	virtual void initialize(){
		HistoAnalyzer::initialize();

		hdeltaTOF = get<TH2>( "dtof_diff" );
		hshape = get<TH1>( "dtof_shape" );

		// non-negative the shape bins so that it can be treated as a PDF
		for ( int i = 1; i <= hshape->GetXaxis()->GetNbins(); i++ ){
			if ( hshape->GetBinContent( i ) < 0 )
				hshape->SetBinContent( i, 0 );
		}

		fmean = get<TF1>( "fmean" );
		fsigma = get<TF1>( "fsigma" );

		binwidth = hshape->GetXaxis()->GetBinWidth( 10 );

		r.SetSeed(0);
	}

	float linterp( float x, float x0, float y0, float x1, float y1 ){
		return y0 + (x - x0) * (y1 - y0) / ( x1 - x0 );
	}

	float meanAtPt( float pt ){
		if ( pt < 1.0 )
			return 0.15;

		float m = fmean->Eval( pt );
		int ptbin = hdeltaTOF->GetXaxis()->FindBin( pt );
		float ptbincent = hdeltaTOF->GetXaxis()->GetBinCenter( ptbin );
		
		if ( pt > 2.2 )
			m = fmean->Eval( 2.2 );
		else if ( pt < 1.0 )
		{
			TH1 * h1 = hdeltaTOF->ProjectionY( "tmp", ptbin, ptbin );
			float m0 = h1->GetMean();
			float pt0 = ptbincent;

			delete h1;

			int ptbin2 = ptbin+1;
			if ( pt < ptbincent && ptbin > 5  ){
				ptbin2 = ptbin - 1;
			}
			h1 = hdeltaTOF->ProjectionY( "tmp", ptbin2, ptbin2 );
			float m1 = h1->GetMean();
			float pt1 = hdeltaTOF->GetXaxis()->GetBinCenter( ptbin2 );

			m = linterp( pt, pt0, m0, pt1, m1 );
			delete h1;
		}

		return m;
	}


	TH1 * pdfAtPt( 	float pt, int nSamples = 100000 ){

		float mean = hshape->GetMean();
		float sigma = hshape->GetRMS() / 1.0; // tuned to match data

		float m = meanAtPt( pt );
		float s = fsigma->Eval( pt );

		cout << "pT=" << pt << " : m=" << m << ", s=" << s << endl;

		TH1F * hsample = new TH1F( "hsample", "", 1250, -5, 5 );
		for ( int i = 0; i < nSamples; i ++ ){
			float v = (hshape->GetRandom() - mean) / (sigma );
			
			v *= s;
			v += m;
			v += r.Gaus( 0, binwidth );

			hsample->Fill( v  );
		}
		hsample->SetLineColor( kRed+1 );
		return hsample;
	}

	virtual void make_comparison(float pt){
		gStyle->SetOptStat( 0 );

		int ptbin = hdeltaTOF->GetXaxis()->FindBin( pt );
		float pt1 = hdeltaTOF->GetXaxis()->GetBinLowEdge( ptbin );
		float pt2 = hdeltaTOF->GetXaxis()->GetBinUpEdge( ptbin );
		
		cout << "ptbin = " << ptbin << endl;
		TH1 * h1 = hdeltaTOF->ProjectionY( "hdata", ptbin, ptbin );
		TH1 * hsample = pdfAtPt( pt, config.getInt( "p.nSamples", 1000 ) );

		h1->Scale( 1.0 / h1->Integral() );
		hsample->Scale( 1.0 / hsample->Integral() );
		hsample->Scale( 1.0, "width" );
		h1->Scale( 1.0, "width" );
		
		hsample->SetTitle( TString::Format( "%0.3f < p_{T} < %0.3f; #DeltaTOF (ns); dN/d #DeltaTOF (ns)^{-1}", pt1, pt2 ) );
		hsample->Draw("");
		h1->Draw("same");
		gPad->SetLogy(1);

		TLegend * leg = new TLegend( 0.6, 0.8, 0.9, 0.9 );
		leg->AddEntry( h1, "data (J/#Psi #mus)" );
		leg->AddEntry( hsample, "Sampled" );
		leg->Draw("same");



		gPad->Print( TString::Format("compare_sample_ptbin_%d.pdf", ptbin) );
	}


	virtual void make(){
		LOG_F( INFO, "" );

		book->cd();

		for ( float pt : {1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 5.0, 10.0} ){
			make_comparison( pt );
		}


		TH1 * hshape = get<TH1>( "dtof_shape" );
		float mean = hshape->GetMean();
		float sigma = hshape->GetRMS();

		int nSamples = config.getInt( "p.nSamples", 1000 );
		TH2F * hsignalPdf = new TH2F( "hsignalPdf", "", 1500, 0, 15, 500, -5, 5 );
		for ( int ix = 1; ix < hsignalPdf->GetXaxis()->GetNbins()+1; ix++ ){
			float pt = hsignalPdf->GetXaxis()->GetBinCenter( ix );
			if ( pt < 0.8 ) continue;
			TH1 * hsample = pdfAtPt( pt, nSamples );

			for ( int is = 0; is < nSamples; is++ ){
				float dtof = hsample->GetRandom();
				hsignalPdf->Fill( pt, dtof );
			}

			delete hsample;
		}
	}

};



#endif