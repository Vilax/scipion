/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 				jlvilas@cnb.csic.es (2018)
 * 			   Javier Vargas            		jvargas@cnb.csic.es (2018)
 * 			   Carlos Oscar S. Sorzano          coss@cnb.csic.es (2018)
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

#ifndef _PROG_LOC_SHARP
#define _PROG_LOC_SHARP

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

/**@defgroup Local Sharpening
   @ingroup ReconsLibrary */
//@{
/** SSNR parameters. */

class ProgLocSharp : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnOut, fnVol, fnVol2, fnMask, fnchim, fnSym, fnMeanVol, fnMaskOut, fnMaxVol,
	fnMinVol, fnMd, fnVar, fnSrp, fnSph, fnDirections;

	/** sampling rate, minimum resolution, and maximum resolution */
	double sampling, minRes, maxRes, R, ang_sampling, N_points, pepep, N_directions;

	/** Is the volume previously masked?*/
	int NVoxelsOriginalMask, Nvoxels;

	/** Step in digital frequency */
	double N_freq, trimBound, significance;

	/** The search for resolutions is linear or inverse**/
	bool nonmanual_mask;

public:

    void defineParams();
    void readParams();
    void produceSideInfo();

    /* Mogonogenid amplitud of a volume, given an input volume,
     * the monogenic amplitud is calculated and low pass filtered at frequency w1*/
    void amplitudeMonogenicSignal3D_fast(MultidimArray< std::complex<double> > &myfftV,
    		double w1, double w1l, double wH, MultidimArray<double> &amplitude,
    		int count, FileName fnDebug);

    void resolution2eval(int &fourier_idx, double min_step,
			double &resolution, double &last_resolution,
			int &last_fourier_idx,
			double &freq, double &freqL, double &freqH,
			bool &continueIter, bool &breakIter, bool &doNextIteration);

    void run();

public:
    Image<int> mask;
    MultidimArray<double> iu, VRiesz; // Inverse of the frequency
	MultidimArray< std::complex<double> > fftV, *fftN; // Fourier transform of the input volume
	FourierTransformer transformer_inv, transformer_direct;
	MultidimArray< std::complex<double> > fftVRiesz, fftVRiesz_aux;
	FourierFilter lowPassFilter, FilterBand;
	int N_smoothing;
	Sampling mysampling;
	Matrix2D<double> angles;
	Matrix1D<double> freq_fourier;
	Image<double> Vfiltered, VresolutionFiltered;

};
//@}
#endif
