



void sample_test( double pt = 2.0 ){

	gStyle->SetOptStat( 0 );

	TFile *f = new TFile( "deltaTofPdfs.root" );

	TH1 * hshape = (TH1*)f->Get( "dtof_shape" );
	float mean = hshape->GetMean();
	float sigma = hshape->GetRMS();//hshape->GetStdDev();
	// TH1 * hsample = (TH1*)hshape->Clone( "hsample" );
	// hsample->Reset();
	TH1F * hsample = new TH1F( "hsample", "", 500, -2, 2 );
	hsample->SetLineColor( kRed+1 );

	TH2 * h2 = (TH2*)f->Get( "dtof_diff" );
	int ptbin = h2->GetXaxis()->FindBin( pt );
	float pt1 = h2->GetXaxis()->GetBinLowEdge( ptbin );
	float pt2 = h2->GetXaxis()->GetBinUpEdge( ptbin );
	float binwidth = hshape->GetXaxis()->GetBinWidth( 10 );
	cout << "ptbin = " << ptbin << endl;
	TH1 * h1 = h2->ProjectionY( "hdata", ptbin, ptbin );
	// hshape = h1;


	TF1 * fmean = (TF1*)f->Get( "fmean" );
	TF1 * fsigma = (TF1*)f->Get( "fsigma" );

	float s = fsigma->Eval( pt );
	float m = fmean->Eval( pt );
	if ( pt > 2.0 )
		m = fmean->Eval( 2.0 );
	else if ( pt < 1.75 ){
		m = h1->GetMean();
	}

	cout << "pT=" << pt << " : m=" << m << ", s=" << s << endl;

	TRandom3 r;
	r.SetSeed(0);

	// non-negative the shape bins
	for ( int i = 1; i <= hshape->GetXaxis()->GetNbins(); i++ ){
		if ( hshape->GetBinContent( i ) < 0 )
			hshape->SetBinContent( i, 0 );
	}

	for ( int i = 0; i < 100000; i ++ ){
		float v = (hshape->GetRandom() - mean) / (sigma / 1.5);
		
		v *= s;
		v += m;
		v += r.Gaus( 0, binwidth );

		hsample->Fill( v  );
	}
	

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