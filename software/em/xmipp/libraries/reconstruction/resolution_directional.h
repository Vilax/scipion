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
	FileName fnOut, fnVol, fnMask, fnDoA, fnDirections, fnradial, fnazimuthal, fnMonoRes, fnAniRes,
	fnMDradial, fnMDazimuthal, fnMeanResolution, fnMDThr, fnLowestResolution, fnHighestResolution;

	/** sampling rate, minimum resolution, and maximum resolution */
	double sampling, minRes, maxRes, R, ang_sampling, N_directions, Rparticle, res_step;

	/** Is the volume previously masked?*/
	int NVoxelsOriginalMask, Nvoxels, Nthr;

	/** Step in digital frequency */
	double N_freq, significance;

	/** Analyze radial and azimuthal resolutoin */
	bool checkellipsoids;

public:

    void defineParams();
    void readParams();
    void produceSideInfo();

    /* Mogonogenid amplitud of a volume, given an input volume,
     * the monogenic amplitud is calculated and low pass filtered at frequency w1*/
    void amplitudeMonogenicSignal3D_fast(const MultidimArray< std::complex<double> > &myfftV,
    		double w1, double w1l, double wH, MultidimArray<double> &amplitude,
    		int count, int dir, FileName fnDebug,
    		double rot, double tilt);

    void defineCone(MultidimArray< std::complex<double> > &myfftV,
    		MultidimArray< std::complex<double> > &conefilter, double rot, double tilt);

    void diagSymMatrix3x3(Matrix2D<double> A,
			Matrix1D<double> &eigenvalues, Matrix2D<double> &P);

    void resolution2eval_(int &fourier_idx, double min_step,
			double &resolution, double &last_resolution,
			int &last_fourier_idx,
			double &freq, double &freqL, double &freqH,
			bool &continueIter, bool &breakIter, bool &doNextIteration);

    double firstMonoResEstimation(MultidimArray< std::complex<double> > &myfftV,
    		double w1, double w1l, MultidimArray<double> &amplitude);

    void generateGridProjectionMatching(FileName fnVol_, double smprt,
    		Matrix2D<double> &angles);

    void removeOutliers(Matrix2D<double> &anglesMat, Matrix2D<double> &resolutionMat);

    void ellipsoidFitting(Matrix2D<double> &anglesMat,
			Matrix2D<double> &resolutionMat,
			Matrix2D<double> &axis);

    void radialAzimuthalResolution(Matrix2D<double> &resolutionMat,
    		MultidimArray<int> &pmask,
    		MultidimArray<double> &radial,
    		MultidimArray<double> &azimuthal,
    		MultidimArray<double> &meanResolution,
    		MultidimArray<double> &lowestResolution,
    		MultidimArray<double> &highestResolution,
    		double &radial_Thr, double &azimuthal_Thr,
			MetaData &mdprefDirs);

    void radialAverageInMask(MultidimArray<int> &mask,
    						MultidimArray<double> &inputVol, MetaData &md);

    void run();

public:
    Image<int> mask;
    MultidimArray<double> iu, VRiesz; // Inverse of the frequency
	MultidimArray< std::complex<double> > fftV, conefilter; // Fourier transform of the input volume
	FourierTransformer transformer_inv;
	MultidimArray< std::complex<double> > fftVRiesz, fftVRiesz_aux;
	int N_smoothing;
	Matrix2D<double> angles, resolutionMatrix, maskMatrix, trigProducts;
	Matrix1D<double> freq_fourier;

};
//@}
#endif
