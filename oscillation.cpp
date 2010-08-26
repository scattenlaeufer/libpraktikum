#include "oscillation.h"

void Oscillation::init() {
		oscillationPad = NULL;
		frequencyPad = NULL;

		oscillationFunction = NULL;
		expFunction = NULL;
		oscillationGraph = NULL;
		expStatistics = NULL;
		frequencyGraph = NULL;

		oscillationVisible = false;
		frequencyVisible = false;
		doFit = false;

		xFreq = NULL;
		yFreq = NULL;
		iterations = 10000;
		x_low = 800;
		x_high = 1600;
}

void Oscillation::draw() {
	oscillationVisible = true;
	frequencyVisible = true;
	drawOscillation();
	drawFrequency();
}

void Oscillation::draw(double x_1, double x_2){

	x_low = x_1;
	x_high = x_2;
	draw();
}

void Oscillation::drawOs(){

	oscillationVisible = true;
	frequencyVisible = false;
	drawOscillation();
}

void Oscillation::drawOscillation() {
	canvas->cd();
	if (frequencyVisible)
		oscillationPad = new TPad("oscillationPad", "Oscillation", 0, 0.5, 1, 1);
	else
		oscillationPad = new TPad("oscillationPad", "Oscillation", 0, 0, 1, 1);
	oscillationPad->SetGrid();
	oscillationPad->Draw();
	oscillationPad->cd();

	// Create a graph of the oscillation
	if (hasErrors)
		oscillationGraph = new TGraphErrors(length, x, y, xErrors, yErrors);
	else
		oscillationGraph = new TGraph(length, x, y);
	oscillationGraph->SetName("oscillationGraph");
	oscillationGraph->SetTitle("Oscillation");
	oscillationGraph->Draw("APC");
	if (doFit) {
		//TODO do fit to a damped harmonic oscillator here
		fit();
	}
}

void Oscillation::drawFrequency() {
	canvas->cd();
	if (oscillationVisible)
		frequencyPad = new TPad("frequencyPad", "Frequency", 0, 0, 1, 0.5);
	else
		frequencyPad = new TPad("frequencyPad", "Frequency", 0, 0, 1, 1);
	frequencyPad->SetGrid();
	frequencyPad->Draw();
	frequencyPad->cd();
	runFFT();

	frequencyGraph = new TGraph(iterations, xFreq, yFreq);
	frequencyGraph->SetName("frequencyGraph");
	frequencyGraph->SetTitle("Frequency");
	frequencyGraph->Draw("APL");
	
	printFreqData();	
}

void Oscillation::runFFT() {
	xFreq = new double[iterations];
	yFreq = new double[iterations];
	fourier(length, x, y, iterations, x_low, x_high, xFreq, yFreq, false);
	/*
	TVirtualFFT *fft = TVirtualFFT::FFT(1, (Int_t*)(&length), "R2C M K");
	if (!fft) return;
	fft->SetPoints(y);
	fft->Transform();
	TH1 *hr = NULL;
	hr = TH1::TransformHisto(fft, hr, "MAG");
	hr->Draw();*/
	/*
	// Code copied from praktlib.h
	TH1::AddDirectory(kFALSE);
	Double_t xStart = x[0];
	Double_t xStop = x[length-1];
	TH1D *hfft = new TH1D("hfft", "hfft", length, xStart, xStop);
	for (unsigned int i=0; i<length; i++){
		hfft->SetBinContent(i+1, y[i]);
	}
	TH1 *hm =0;
	TVirtualFFT::SetTransform(0);
	hm = hfft->FFT(hm, "MAG");
	hm->SetStats(kFALSE);
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	Double_t *im_full = new Double_t[length];
	fft->GetPointsComplex(yFreq, im_full);
	int maxPoint=0;
	for (unsigned int i=0;i<length/2;i++)
	{
		if (yFreq[i]>yFreq[maxPoint])
		{
			maxPoint=i;
		}
	}
	// Scale x
	for (unsigned int i=0;i<length;i++)
	{
		xFreq[i] = i/TMath::Abs((x[length-1]-x[0]));
	}
	*/
}



int Oscillation::fourier(int n_datasets, const double* data_t, const double* data_a, int n_f, double f_min, double f_max, double* f_out, double* amp_out, bool progress){

	//calculate distance between two values of frequency
	double df=((double)(f_max-f_min))/((double)(n_f-1));

	double sumsin=0;
	double sumcos=0;
	double omega=0;

	//calculate amplitudes for the requested values of frequency
	for(int i=0; i < n_f; i++){
		sumsin=0;
		sumcos=0;

		//calculate the currently used frequency
		f_out[i]=f_min+(i*df);

		//calculate omega once
		omega=f_out[i]*2*TMath::Pi();

		//Do the magic ...
		for(int j=0; j< n_datasets; j++)
		{
			sumcos +=data_a[j]*cos(omega*data_t[j]);
		}
		for(int j=0; j< n_datasets; j++)
		{
			sumsin +=data_a[j]*sin(omega*data_t[j]);
		}
		amp_out[i]= sqrt(pow(sumcos,2)+pow(sumsin,2));
		//... magic done

		//Debugoutput, calculations may take a while, so you can see the progress
		if((progress==true) && (i % 32 ==0)) 
			cout <<f_out[i] << " " << amp_out[i]<< endl;
	}

	return 0;
}



double Oscillation::peakfinderSchwerpunkt(double * x, double * y, int n_datasets, double &sigm_peak, int start, double val){
	int i0=start; //index or left border of relevant section of peak
	int i1=n_datasets; //index of right border of relevant section of peak
	
	cout << i0 << " " << i1 <<endl;
	
	//if left border is greater than the right one, use the whole array of datapoints
	if(i0>=n_datasets) i0=0;
	
	double x1; //x-values of border of relevant section of peak
	
	double ymax=y[i0]; //maximum of peak
	
	
	//Find maximum of peak
	for(int i=i0; i<n_datasets; i++){
		if(y[i]>ymax){
			ymax=y[i];
		}
		
	}
	
	//find left border of relevant section of peak
	//all values have to be greate than val*ymax
	for(int i=i0; i<n_datasets; i++){
		if(y[i]>ymax*val){ 
			i0=i;
			break;
		}
	}
	
	//find right border of relevant section of peak
	//all values have to be greate than val*ymax
	for(int i=i0+1; i<n_datasets; i++){
		if(y[i]<ymax*val){
			x1=x[i];
			i1=i;
			break;
		}
		x1=x[n_datasets];
	}
	
	
	double sum_y=0;
	double sum_xy=0;
	double sum_xy2=0;
	
	//calculate some values for the relevant section of the peak
	for(int i=i0; i<i1; i++){
		sum_y+=y[i];
		sum_xy+=x[i]*y[i];
		sum_xy2+=x[i]*x[i]*y[i];
	}
	
	
	//do some magic; Not fully understood
	double x_peak=sum_xy/sum_y;
	sigm_peak=sqrt(fabs(sum_xy2/sum_y - pow(x_peak,2)));
	
	return x_peak;
}

void Oscillation::printFreqData(){

	double peak_error;
	double peak = peakfinderSchwerpunkt(xFreq,yFreq,iterations,peak_error,0,1./sqrt(2.));
	printf("peak: %f, peakerror: %f\n",peak,peak_error);
	freqStat = new TPaveText(0.65, 0.7, 0.85, 0.8, "NDC" );
	char buf[64];
	snprintf(buf, sizeof(buf), "f = (%s)Hz",
			utils::printNumber(peak, peak_error).c_str());
	freqStat -> AddText(buf);
	freqStat -> Draw();
}

void Oscillation::fit(){

	printf("start fitting\n");

	// number of maximas
	int n_max = 0;				

	// Arrays needed by TSpectrum
	Float_t *dest = new Float_t[length];
	Float_t *source = new Float_t[length];

	// "converting" double to Float_t, because TSpectrum cannot deal with double
	for (int i = 0; i < length; ++i){
		source[i] = y[i];
	}

	TSpectrum *s = new TSpectrum();

	n_max = s -> SearchHighRes(source, dest, length, 5,5, false, 2, false, 3);

	// get the position in the array of the maxima
	Float_t *xpos = s -> GetPositionX();

	// get the actual x-/y-values of the maxima
	double x_max[n_max];
	double y_max[n_max];
	double x_max_error[n_max];
	double y_max_error[n_max];

	if (hasErrors){
		for (int i = 0; i < n_max; ++i){
			x_max[i] = x[(int) xpos[i]];
			y_max[i] = y[(int) xpos[i]];
			x_max_error[i] = xErrors[(int) xpos[i]];
			y_max_error[i] = yErrors[(int) xpos[i]];
		}

	}
	else{
		for (int i = 0; i < n_max; ++i){
			x_max[i] = x[(int) xpos[i]];
			y_max[i] = y[(int) xpos[i]];
		}
	}

	// get the minima, in principle the same as getting the maxima only with inverted y-array
	
	for (int i = 0; i < length; ++i){
		source[i] = -1 * y[i];
	}

	int n_min = s -> SearchHighRes(source, dest, length, 5,5, false, 2, false, 3);
	
//	delete [] xpos;
	xpos = s -> GetPositionX();
	if (n_min);

	double x_min[n_min];
	double y_min[n_min];
	double x_min_error[n_min];
	double y_min_error[n_min];

	if (hasErrors){
		for (int i = 0; i < n_min; ++i){
			x_min[i] = x[(int) xpos[i]];
			y_min[i] = y[(int) xpos[i]];
			x_min_error[i] = xErrors[(int) xpos[i]];
			y_min_error[i] = yErrors[(int) xpos[i]];
		}

	}
	else{
		for (int i = 0; i < n_min; ++i){
			x_min[i] = x[(int) xpos[i]];
			y_min[i] = y[(int) xpos[i]];
		}
	}

	// combine maxima and minima
	
	int n = n_max + n_min;
	double x_ext[n];
	double y_ext[n];
	double x_ext_error[n];
	double y_ext_error[n];

	if (hasErrors){
		for (int i = 0; i < n_max; ++i){
			x_ext[i] = x_max[i];
			y_ext[i] = y_max[i];
			x_ext_error[i] = y_max_error[i];
			y_ext_error[i] = y_max_error[i];
		}
		for (int i = 0; i < n_min; ++i){
			x_ext[n_max+i] = x_min[i];
			y_ext[n_max+i] = y_min[i];
			x_ext_error[n_max+i] = y_min_error[i];
			y_ext_error[n_max+i] = y_min_error[i];
		}
	}
	else{
		for (int i = 0; i < n_max; ++i){
			x_ext[i] = x_max[i];
			y_ext[i] = y_max[i];
		}
		for (int i = 0; i < n_min; ++i){
			x_ext[n_max+i] = x_min[i];
			y_ext[n_max+i] = -1 * y_min[i];
		}
	}


	// printing values, just for debugging
	for (int i = 0; i < n; ++i){
		printf("%d:   %f  ->   %f - %f\n",i,xpos[i], x_ext[i],y_ext[i]);
	}

	printf("number of minima: %d\n",n_max);

	TGraph *g_max;
	TGraph *g_min;
	if (hasErrors){
		g_max = new TGraphErrors(n_max,x_max,y_max,x_max_error,y_max_error);
		g_min = new TGraphErrors(n_min,x_min,y_min,x_min_error,y_min_error);
	}
	else{
		g_max = new TGraph(n_max,x_max,y_max);
		g_min = new TGraph(n_min,x_min,y_min);
	}

	// fitting exponetial function to maxima
	expFunction = new TF1("envalope","[0]*exp(-[1]*x) + [2]",x[0],x[length-1]);
	expFunction -> SetLineColor(4);

	TF1 *expFunction_min = new TF1("envalope_min","[0]*exp(-[1]*x) + [2]",x[0],x[length-1]);
	expFunction_min -> SetLineColor(4);

	g_max -> Fit("envalope");
	g_min -> Fit("envalope_min");

	oscillationPad -> cd();
//	g_ein -> Draw("*SAME");
	expFunction -> Draw("SAME");
	expFunction_min -> Draw("SAME");
	printf("%d: %f\n",length,x[length-1]);

	// print fitting parameters
	double gamma_max = expFunction -> GetParameter(1);
	double gamma_min = expFunction_min -> GetParameter(1);
	double gamma_max_error = expFunction -> GetParError(1);
	double gamma_min_error = expFunction_min -> GetParError(1);

	if (hasErrors){
		expStatistics = new TPaveText(0.65, 0.6, 0.85, 0.8, "NDC" );
	}
	else{
		expStatistics = new TPaveText(0.65, 0.7, 0.85, 0.8, "NDC" );
	}
	char buf[64];
	snprintf(buf, sizeof(buf), "#gamma_{1} = (%s)#frac{1}{s}",
			utils::printNumber(gamma_max, gamma_max_error).c_str());
	expStatistics -> AddText(buf);
	if (hasErrors){
		snprintf(buf, sizeof(buf), "#frac{#chi^{2}_{1}}{ndf_{1}} = %s",
				utils::toString(expFunction -> GetChisquare() / expFunction ->GetNDF(), 2).c_str());
		expStatistics -> AddText(buf);
	}
	snprintf(buf, sizeof(buf), "#gamma_{2} = (%s)#frac{1}{s}",
			utils::printNumber(gamma_min,gamma_min_error).c_str());
	expStatistics -> AddText(buf);
	if (hasErrors){
		snprintf(buf, sizeof(buf), "#frac{#chi^{2}_{2}}{ndf_{2}} = %s",
				utils::toString(expFunction_min -> GetChisquare() / expFunction_min ->GetNDF(), 2).c_str());
		expStatistics -> AddText(buf);
	}

	double gamma[2] = {gamma_min, gamma_max};
	double gamma_error[2] = {gamma_min_error, gamma_max_error};
	double gamma_mean_error;
	double rms;
	double gamma_mean = utils::weightedMean(gamma,gamma_error,2,gamma_mean_error,rms);

	snprintf(buf, sizeof(buf), "#bar{#gamma} = (%s)#frac{1}{s}",
			utils::printNumber(gamma_mean,gamma_mean_error).c_str());
	expStatistics -> AddText(buf);
	expStatistics -> Draw();
	printf("fitting done\n");

	delete g_max;
	delete s;
	delete [] dest;
	delete [] source;
}


