/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2016)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef _PROG_RES_DIR
#define _PROG_RES_DIR

#include <iostream>
#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
//#include <data/symmetries.h>
#include <data/sampling.h>
#include <data/metadata.h>
//#include <data/matrix2d.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_fftw.h>
#include <math.h>
#include <limits>
#include <complex>
#include "fourier_filter.h"
#include <data/filters.h>
#include <string>
#include "symmetrize.h"

/**@defgroup Monogenic Resolution
   @ingroup ReconsLibrary */
//@{
/** SSNR parameters. */

class ProgResDir : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnOut, fnVol, fnVol2, fnMask, fnchim, fnSym, fnMeanVol, fnMaskOut, fnMaxVol,
	fnMinVol, fnMd, fnVar, fnDoA, fnSph, fnDirections;

	/** sampling rate, minimum resolution, and maximum resolution */
	double sampling, minRes, maxRes, R, ang_sampling, N_points, pepep, N_directions;

	/** Is the volume previously masked?*/
	int NVoxelsOriginalMask, Nvoxels;

	/** Step in digital frequency */
	double N_freq, trimBound, significance;

	/** The search for resolutions is linear or inverse**/
	bool exactres, noiseOnlyInHalves, nonmanual_mask;

public:

    void defineParams();
    void readParams();
    void produceSideInfo();

    /* Mogonogenid amplitud of a volume, given an input volume,
     * the monogenic amplitud is calculated and low pass filtered at frequency w1*/
    void amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
    		double w1, double w1l, double wH, MultidimArray<double> &amplitude,
    		int count, int dir, FileName fnDebug,
    		double rot, double tilt);

    void amplitudeMonogenicSignal3D2(MultidimArray< std::complex<double> > &myfftV,
    		double w1, double w1l, double wH, MultidimArray<double> &amplitude,
    		int count, int dir, FileName fnDebug,
    		double rot, double tilt);

    void inertiaMatrix(MultidimArray<double> &resolutionVol,
			   MultidimArray<double> &Inertia_11,
			   MultidimArray<double> &Inertia_12,
			   MultidimArray<double> &Inertia_13,
			   MultidimArray<double> &Inertia_22,
			   MultidimArray<double> &Inertia_23,
			   MultidimArray<double> &Inertia_33,
			   MultidimArray<double> &SumRes,
			   double rot, double tilt);

    void diagSymMatrix3x3(Matrix2D<double> A,
			Matrix1D<double> &eigenvalues, Matrix2D<double> &P);

    void degreeOfAnisotropy(Matrix1D<double> eigenvalues, Matrix2D<double> eigenvectors,
			double &doa, double &direction_x, double &direction_y, double &direction_z,
			int &counter);

    void resolution2eval(int &count_res, double step,
			double &resolution, double &last_resolution,
			double &freq, double &freqL, double &freqH,
			int &last_fourier_idx,
			bool &continueIter, bool &breakIter, bool &doNextIteration);

    void generateGridProjectionMatching(FileName fnVol_, double smprt,
    		Matrix2D<double> &angles);

    void defineDirection(Matrix1D<double> &r0, Matrix1D<double> &rF,
			Matrix2D<double> &direction, double &eigenvalue, double &eigenvalue_max, int eigdir,
			int k , int i, int j);

    void defineSegment(Matrix1D<double> &r0, Matrix1D<double> &rF,
			MultidimArray<int> &arrows, double &elongation);


    void run();

public:
    Image<int> mask;
    MultidimArray<double> iu, VRiesz; // Inverse of the frequency
	MultidimArray< std::complex<double> > fftV, *fftN; // Fourier transform of the input volume
	FourierTransformer transformer_inv, transformer_direct;
	MultidimArray< std::complex<double> > fftVRiesz, fftVRiesz_aux;
	FourierFilter lowPassFilter, FilterBand;
	bool halfMapsGiven;
	int N_smoothing;
	Sampling mysampling;
	Matrix2D<double> angles;
	Matrix1D<double> freq_fourier;

};
//@}
#endif
