/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2018)
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

#ifndef _PROG_LOC_SHARPENING
#define _PROG_LOC_SHARPENING

#include <iostream>
#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
#include <data/metadata.h>
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

class ProgLocSharpening : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnOut, fnVol, fnRes;

	/** sampling rate, minimum resolution, and maximum resolution */
	double sampling, maxRes, minRes, lambda;
	int Niter;

public:

    void defineParams();
    void readParams();
    void produceSideInfo();

    /* Mogonogenid amplitud of a volume, given an input volume,
     * the monogenic amplitud is calculated and low pass filtered at frequency w1*/
    void lowPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
    		double w, double wL, MultidimArray<double> &filteredVol, int count);

    void bandPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
    		double w, double wL, MultidimArray<double> &filteredVol, int count);

    void resolution2eval(int &count_res, double step,
    								double &resolution, double &last_resolution,
    								double &freq, double &freqL,
    								int &last_fourier_idx,
    								bool &continueIter,	bool &breakIter,
    								bool &doNextIteration);

    void maxMinResolution(MultidimArray<double> &resVol,
			double &maxRes, double &minRes);

    void richardsonLucy(MultidimArray<double> &vol, MultidimArray<double> &filteredVol,
			MultidimArray<double> &outVol);

    void resVol2IdxVol(const MultidimArray<double> &resVol, double &maxres, MultidimArray<int> &idxVol,
    		std::vector<int> &idxList);

    void run();

public:
    std::vector<int> idxList;
    MultidimArray<int> idxVol;
    MultidimArray<double> Vorig;
    MultidimArray<double> resVol;
    MultidimArray<double> iu, sharpenedMap; // Inverse of the frequency
	MultidimArray< std::complex<double> > fftV, fftVfilter; // Fourier transform of the input volume
	FourierTransformer transformer_inv, transformer;
	Image<double> Vfiltered, VresolutionFiltered;
	Matrix1D<double> freq_fourier;
	Matrix2D<double> resolutionMatrix, maskMatrix;
};
//@}
#endif
