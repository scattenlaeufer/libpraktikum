#include "dataanalysis.h"
#include "options.h"
#include "utils.h"

#include <iostream>
#include <stdio.h>
#include <TStyle.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TList.h>
#include <TVirtualFFT.h>
#include <TSpectrum.h>
#include <TMath.h>


/** \class Oscillation
 * \brief Analysis of harmonic oscillations
 *
 */

class Oscillation : DataAnalysis {
	public:

		/** Creates an oscillation object with errors.
		 * \param[in] *_x Pointer to the array of x-values
		 * \param[in] *_y Pointer to the array of y-values
		 * \param[in] *_xErrors Pointer to the array of x-errors
		 * \param[in] *_yErrors Pointer to the array of y-errors
		 * \param[in] _length Length of the arrays
		 */
		inline Oscillation(const double *_x, const double *_y, const double *_xErrors, const double *_yErrors, const unsigned int _length)
			: DataAnalysis(_x, _y, _xErrors, _yErrors, _length)
		{
			length = _length;
			init();
		}

		/** Creates an oscillation object without errors.
		 * \param[in] *_x Pointer to the array of x-values
		 * \param[in] *_y Pointer to the array of y-values
		 * \param[in] _length Length of the arrays
		 */
		inline Oscillation(const double *_x, const double *_y, const unsigned int _length)
			: DataAnalysis(_x, _y, _length)
		{
			length = _length;
			init();
		}

		/** Shuold in some time perform a fast fourier tranformation, but at the moment peformce 
		 * a fourier transformation with 10000 itterations as default.
		 */
		void runFFT();

		/** Starts the analysis of the given data. 
		 * Sets oscillationVisible and frequencyVisible to true and peformce drawOscillation() 
		 * and drawFrequency().
		 */
		void draw();

		/** Performce draw() with predifined borders of the fft
		 * \param[in] x_1 lower border for fft
		 * \param[in] x_2 upper border for fft
		 */
		void draw(double x_1,double x_2);

		/** Draws only the oscillationPad.
		 */
		void drawOs();

		/** Does alot of magic. Especialy creating and filling of the osillationPad.
		 */
		void drawOscillation();

		/** Does alot of magic. Especialy creating and filling of the frequencyPad.
		 */
		void drawFrequency();

		inline TPad *getOscillationPad() const {
			return oscillationPad;
		}

		inline TPad *getFrequencyPad() const {
			return frequencyPad;
		}


		inline TF1 *getExpFunction() const {
			return expFunction;
		}

		inline TGraph *getOscillationGraph() const {
			return oscillationGraph;
		}

		inline TGraph *getFrequencyGraph() const {
			return frequencyGraph;
		}

		inline TPaveText *getExpStats() const {
			return expStatistics;
		}
		/** Enables the fitting of an envalope to a dampd harmonic oscillation
		 */
		inline void enableFitting() {
			doFit = true;
		}


	protected:
		TPad *oscillationPad;
		TPad *frequencyPad;

		TF1 *oscillationFunction;
		TF1 *expFunction;
		TGraph *oscillationGraph;
		TPaveText *expStatistics;
		TPaveText *freqStat;
		TGraph *frequencyGraph;

		bool oscillationVisible;
		bool frequencyVisible;
		bool doFit;

		// the x values of the frequency spectrum
		double *xFreq;
		// the y values of the frequency spectrum, i.e. the amplitude
		double *yFreq;
		unsigned int iterations;
		// borders for fourier
		double x_low;
		double x_high;

		int length;

	private:
		// This is called by the constructors
		void init();

		/**
		 * \brief simple discrete fourier transformation
		 * 
		 * Fourier transformation as in the maple library of Anfängerpraktikum at RWTH \n
		 * \a f_out and \a amp_out have to be allocated with size>=\a n_f 
		 * 
		 * \param[in] n_datasets numbers of values in \a data_t and \a data_a
		 * \param[in] data_t time of each datapoint
		 * \param[in] data_a input value of each datapoint
		 * \param[in] n_f number of discrete frequencies to calculate amplitude for
		 * \param[in] f_min minimum frequency
		 * \param[in] f_max maximum frequency
		 * \param[out] f_out pointer to double, frequency value of each generated datapoint (array of \a n_f double values has to be allocated)
		 * \param[out] amp_out pointer to double, amplitude of each generated datapoint (array of \a n_f double values has to be allocated)
		 * \param[in] progress (optional) weather to show generated values or not (can be useful to monitor progress which can take a while) (default: true)
		 * \return allways 0
		 */
		int fourier(int n_datasets, const double* data_t, const double* data_a, int n_f, double f_min, double f_max, double* f_out, double* amp_out, bool progress);

		/**
		 * \brief Peakfinder with integrated calculation of RMS of peak
		 * 
		 * Algorithm inspired by MAPLE Library for Anfängerpraktikum \n
		 * Only finds one single Peak
		 * 
		 * \param[in] x Array with frequency belonging to the amplitude data (generated as f_out by fourier, eg.)
		 * \param[in] y Array with the amplitude data (generated as amp_out by fourier, eg.)
		 * \param[in] n_datasets Number of datapoints in amp and f (if setting lower than real size of the array, the rest of the points are ignored)
		 * \param[out] sigm_peak RMS of the found peak
		 * \param[in] start (optional) datapoint to start with (ignorred if greater than n_freqs) (default: 0)
		 * \param[in] val (optional) only datapoints greater than ymax*val in the peak are used (default: 1/sqrt(2))
		 * 
		 * \returns frequency of the found peak (as double, NOT as index of a data point)
		*/
		double peakfinderSchwerpunkt(double * x, double * y, int n_datasets, double &sigm_peak, int start, double val);

		void printFreqData();

		void fit();
};


