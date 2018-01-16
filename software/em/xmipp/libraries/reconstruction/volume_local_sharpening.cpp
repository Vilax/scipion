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

#include "volume_local_sharpening.h"
//#define DEBUG
//#define DEBUG_MASK
//#define DEBUG_DIR
//define DEBUG_FILTER
//#define MONO_AMPLITUDE
//define DEBUG_SYMMETRY

void ProgLocSharp::readParams()
{
	fnVol = getParam("--vol");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	sampling = getDoubleParam("--sampling_rate");
	R = getDoubleParam("--volumeRadius");
	fnSym = getParam("--sym");
	significance = getDoubleParam("--significance");
	fnSrp = getParam("--sharpened_vol");
}


void ProgLocSharp::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">        : Input volume");
	addParamsLine("  [--mask <vol_file=\"\">]     : Mask defining the macromolecule");
	addParamsLine("                               :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                               :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--sym <symmetry>]			  : Symmetry (c1, c2, c3,..d1, d2, d3,...)");
	addParamsLine("  [--sampling_rate <s=1>]      : Sampling rate (A/px)");
	addParamsLine("  [--volumeRadius <s=100>]     : This parameter determines the radius of a sphere where the volume is");
	addParamsLine("  [--significance <s=0.95>]    : The level of confidence for the hypothesis test.");
	addParamsLine("  [--md_resdir <file=\".\">]   : Metadata with mean resolution by direction.");
}

void ProgLocSharp::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;

	Image<double> V;
	V.read(fnVol);

	V().setXmippOrigin();

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);
	N_freq = ZSIZE(inputVol);
	maxRes = ZSIZE(inputVol);
	minRes = 2*sampling;

	transformer.FourierTransform(inputVol, fftV);
	iu.initZeros(fftV);

	// Frequency volume
	double uz, uy, ux, uz2, u2, uz2y2;
	long n=0;

	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
		uz2=uz*uz;

		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
			uz2y2=uz2+uy*uy;

			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
				u2=uz2y2+ux*ux;
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
				++n;
			}
		}
	}


	// Prepare low pass filter
	lowPassFilter.FilterShape = RAISED_COSINE;
	lowPassFilter.raised_w = 0.01;
	lowPassFilter.do_generate_3dmask = false;
	lowPassFilter.FilterBand = LOWPASS;

	// Prepare mask
	MultidimArray<int> &pMask=mask();

	if (fnMask != "")
	{
		mask.read(fnMask);
		mask().setXmippOrigin();
	}
	else
	{
		std::cout << "Error: a mask ought to be provided" << std::endl;
		exit(0);
	}

	//use the mask for preparing resolution volumes
	Image<double> AvgResoltion, VarianzeResolution;
	AvgResoltion().resizeNoCopy(inputVol);
	VarianzeResolution().resizeNoCopy(inputVol);
	AvgResoltion().initZeros();
	VarianzeResolution().initZeros();

	AvgResoltion.write(fnOut);
	VarianzeResolution.write(fnVar);

	AvgResoltion.clear();
	VarianzeResolution.clear();

	N_smoothing = 15;
	NVoxelsOriginalMask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
		if (A3D_ELEM(pMask, k, i, j) == 1)
			NVoxelsOriginalMask++;
		if (i*i+j*j+k*k > (R-N_smoothing)*(R-N_smoothing))
			A3D_ELEM(pMask, k, i, j) = -1;
	}

	#ifdef DEBUG_MASK
	mask.write("mask.vol");
	#endif

	fftN=&fftV;

	freq_fourier.initZeros(ZSIZE(inputVol));
	int size = ZSIZE(inputVol);
	maxRes = size;
	minRes = 1;
	V.clear();


	double u;
	int size_fourier(ZSIZE(fftV));

	for(size_t k=0; k<size_fourier; ++k)
	{
		FFT_IDX2DIGFREQ(k,size, u);
		VEC_ELEM(freq_fourier,k) = u;
		std::cout << "freq_fourier = " << sampling/u  << std::endl;
	}
}


void ProgLocSharp::amplitudeMonogenicSignal3D_fast(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, FileName fnDebug)
{
	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);
	amplitude.resizeNoCopy(VRiesz);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(freq-freqH);

	double uz, uy, ux;
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);

				double un=1.0/iun;
				if (freqH<=un && un<=freq)
				{
					//TODO: Check if fftVRiesz can be made equal to myfftV
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
				} else if (un>freq)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
				}
				++n;
			}
		}
	}


	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

//	#ifdef DEBUG_DIR
//	if (count == 0)
//	{
//		Image<double> filteredvolume;
//		filteredvolume = VRiesz;
//		filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
//	}
//	#endif


	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				ux = VEC_ELEM(freq_fourier,j);
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second and third component of Riesz vector
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		uz = VEC_ELEM(freq_fourier,k);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier,i);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uz*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+= DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	transformer_inv.inverseFourierTransform(fftVRiesz_aux, VRiesz);

	amplitude.setXmippOrigin();
	int z_size = ZSIZE(amplitude);
	int x_size = XSIZE(amplitude);
	int y_size = YSIZE(amplitude);

	double limit_radius = (z_size*0.5-N_smoothing);
	n=0;
	for(int k=0; k<z_size; ++k)
	{
		uz = (k - z_size*0.5);
		for(int i=0; i<y_size; ++i)
		{
			uy = (i - y_size*0.5);
			for(int j=0; j<x_size; ++j)
			{
				ux = (j - x_size*0.5);
				DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
				DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));

				double radius = sqrt(ux*ux + uy*uy + uz*uz);
				if ((radius>=limit_radius) && (radius<=(z_size*0.5)))
					DIRECT_MULTIDIM_ELEM(amplitude, n) *= 0.5*(1+cos(PI*(limit_radius-radius)/(N_smoothing)));
				else if (radius>(0.5*z_size))
					DIRECT_MULTIDIM_ELEM(amplitude, n) = 0;
				++n;
			}
		}
	}
	//TODO: change (k - z_size*0.5)
	//TODO: check the number of un=1.0/DIRECT_MULTIDIM_ELEM(iu,n); and optimize
	//TODO: Use the square of the monogenic amplitude


//		#ifdef MONO_AMPLITUDE
//		Image<double> saveImg2;
//		saveImg2 = amplitude;
//		if (fnDebug.c_str() != "")
//		{
//			FileName iternumber = formatString("smoothed_volume_%i.vol", count);
//			saveImg2.write(fnDebug+iternumber);
//		}
//		saveImg2.clear();
//		#endif


	//amplitude.setXmippOrigin();

	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);

    double raised_w = PI/(freqL-freq);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
	{
		double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
		if (freqL>=un && un>=freq)
			DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-freq)));
		else
			if (un>freqL)
				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
	}
	transformer_inv.inverseFourierTransform();

//	#ifdef MONO_AMPLITUDE
//
//	if (fnDebug.c_str() != "")
//	{
//		saveImg2 = amplitude;
//		FileName iternumber = formatString("_Filtered_Amplitude_%i.vol", count);
//		saveImg2.write(fnDebug+iternumber);
//	}
//	saveImg2.clear();
////	#endif // DEBUG

}


void ProgLocSharp::resolution2eval(int &fourier_idx, double min_step,
								double &resolution, double &last_resolution,
								int &last_fourier_idx,
								double &freq, double &freqL, double &freqH,
								bool &continueIter, bool &breakIter, bool &doNextIteration)
{
	int volsize = ZSIZE(VRiesz);

	FFT_IDX2DIGFREQ(fourier_idx, volsize, freq);

	resolution = sampling/freq;
	std::cout << "res = " << resolution << std::endl;
	std::cout << "min_step = " << min_step << std::endl;

	//TODO: I am sure that the abs can be removed
	if ( fabs(resolution - last_resolution)<min_step )
	{
		std::cout << "entro last_resolution = "  << last_resolution << "res = " << resolution  << std::endl;
		freq = sampling/(last_resolution-min_step);
		DIGFREQ2FFT_IDX(freq, volsize, fourier_idx);
		FFT_IDX2DIGFREQ(fourier_idx, volsize, freq);

		if (fourier_idx == last_fourier_idx)
		{
			continueIter = true;
			++fourier_idx;
			return;
		}
	}

	resolution = sampling/freq;
	last_resolution = resolution;

	double step = 0.05*resolution;

	double resolution_L, resolution_H;

	if ( step < min_step)
	{
		resolution_L = resolution - min_step;
		resolution_H = resolution + min_step;
	}
	else
	{
		resolution_L = 0.95*resolution;
		resolution_H = 1.05*resolution;
	}

	freqH = sampling/(resolution_H);
	freqL = sampling/(resolution_L);

	std::cout << "freq_H = " << freqH << std::endl;
	std::cout << "freq_L = " << freqL << std::endl;

	if (freqH>0.5 || freqH<0)
		freqH = 0.5;

	if (freqL>0.5 || freqL<0)
		freqL = 0.5;
	int fourier_idx_H, fourier_idx_L;

	DIGFREQ2FFT_IDX(freqH, volsize, fourier_idx_H);
	DIGFREQ2FFT_IDX(freqL, volsize, fourier_idx_L);

	if (fourier_idx_H == fourier_idx)
		fourier_idx_H = fourier_idx - 1;

	if (fourier_idx_L == fourier_idx)
		fourier_idx_L = fourier_idx + 1;

	FFT_IDX2DIGFREQ(fourier_idx_H, volsize, freqH);
	FFT_IDX2DIGFREQ(fourier_idx_L, volsize, freqL);

	std::cout << "freq_H = " << freqH << std::endl;
	std::cout << "freq_L = " << freqL << std::endl;

	if (freq>0.49 || freq<0)
	{
		std::cout << "Nyquist limit reached" << std::endl;
		breakIter = true;
		doNextIteration = false;
		return;
	}
	else
	{
		breakIter = false;
		doNextIteration = true;
	}
	std::cout << "resolution = " << resolution << "  resolutionL = " <<
				sampling/(freqL) << "  resolutionH = " << sampling/freqH
				<< "  las_res = " << last_resolution << std::endl;
	last_fourier_idx = fourier_idx;
	++fourier_idx;
}


void ProgLocSharp::run()
{
	produceSideInfo();

	bool continueIter = false, breakIter = false;
	double criticalZ=icdf_gauss(significance);

	double range = maxRes-minRes;
	double step = range/N_freq;

	if (step<0.3)
		step=0.3;

	step = 0.3;

	std::cout << "maxRes = " << maxRes << std::endl;
	std::cout << "minRes = " << minRes << std::endl;
	std::cout << "step = " << step << std::endl;

	Image<double> outputResolution;
	outputResolution().initZeros(VRiesz);
	MultidimArray<double> &pOutputResolution = outputResolution();
	MultidimArray<double> amplitudeMS, amplitudeMN;
	MultidimArray<int> mask_aux = mask();
	MultidimArray<int> &pMask = mask_aux;
	std::vector<double> list;
	double resolution;  //A huge value for achieving last_resolution < resolution
	double freq, freqL, freqH, counter, resolution_2;
	double max_meanS = -1e38;
	double cut_value = 0.025;

	bool doNextIteration=true;
	bool lefttrimming = false;

	int fourier_idx = 2, last_fourier_idx = -1, iter = 0, fourier_idx_2;
	int count_res = 0;
	double criticalW=-1;
	std::cout << "--------------NEW DIRECTION--------------" << std::endl;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		if (DIRECT_MULTIDIM_ELEM(pMask, n) == 1)
			DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = maxRes;

	std::vector<double> noiseValues;
	FileName fnDebug;
	double last_resolution = 0;

	do
	{
		continueIter = false;
		breakIter = false;
		//std::cout << "--------------Frequency--------------" << std::endl;

		resolution2eval(fourier_idx, step,
						resolution, last_resolution, last_fourier_idx,
						freq, freqL, freqH,
						continueIter, breakIter, doNextIteration);

		if (breakIter)
			break;

		if (continueIter)
			continue;

		std::cout << "resolution = " << resolution << "  resolutionL = " << sampling/freqL << "  resolutionH = " << sampling/freqH << std::endl;
		std::cout << "resolution = " << freq 	   << "  resolutionL = " << freqL 		   << "  resolutionH = " << freqH << std::endl;


		list.push_back(resolution);

		if (iter<2)
			resolution_2 = list[0];
		else
			resolution_2 = list[iter - 2];

		fnDebug = "Signal";

		amplitudeMonogenicSignal3D_fast(fftV, freq, freqH, freqL, amplitudeMS, iter, fnDebug);

		double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
		noiseValues.clear();

		double amplitudeValue;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
			if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
			{
				sumS  += amplitudeValue;
				sumS2 += amplitudeValue*amplitudeValue;
				++NS;
			}
			else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
			{
				sumN  += amplitudeValue;
				sumN2 += amplitudeValue*amplitudeValue;
				++NN;
			}
		}


		if ( (NS/NVoxelsOriginalMask)<cut_value ) //when the 2.5% is reached then the iterative process stops
		{
			std::cout << "Search of resolutions stopped due to mask has been completed" << std::endl;
			doNextIteration =false;
			Nvoxels = 0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
			{
			  if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n) > 0)
				DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
			}
		#ifdef DEBUG_MASK
		mask.write("partial_mask.vol");
		#endif
		lefttrimming = true;
		}
		else
		{
		if (NS == 0)
		{
			std::cout << "There are no points to compute inside the mask" << std::endl;
			std::cout << "If the number of computed frequencies is low, perhaps the provided"
					"mask is not enough tight to the volume, in that case please try another mask" << std::endl;
			break;
		}

		double meanS=sumS/NS;
//			double sigma2S=sumS2/NS-meanS*meanS;
		double meanN=sumN/NN;
		double sigma2N=sumN2/NN-meanN*meanN;

		if (meanS>max_meanS)
			max_meanS = meanS;

		if (meanS<0.00001*max_meanS)
		{
			//std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
			//std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
			std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
			std::cout << "\n"<< std::endl;
			doNextIteration = false;
		}
		else
		{
			// Check local resolution
			double thresholdNoise;
			thresholdNoise = meanN+criticalZ*sqrt(sigma2N);

			#ifdef DEBUG
			  std::cout << "Iteration = " << iter << ",   Resolution= " << resolution << ",   Signal = " << meanS << ",   Noise = " << meanN << ",  Threshold = " << thresholdNoise <<std::endl;
			#endif

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
			{
				if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
					if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
					{
						DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
						DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution;//sampling/freq;
					}
					else
					{
						DIRECT_MULTIDIM_ELEM(pMask, n) += 1;
						if (DIRECT_MULTIDIM_ELEM(pMask, n) >2)
						{
							DIRECT_MULTIDIM_ELEM(pMask, n) = -1;
							DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution_2; //resolution + counter*step;
						}
					}
			}

			#ifdef DEBUG_MASK
			FileName fnmask_debug;
			fnmask_debug = formatString("maske_%i.vol", iter);
			mask.write(fnmask_debug);
			#endif

			//#ifdef DEBUG
				std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
//					std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
				std::cout << "  meanS= " << meanS << " NS= " << NS << std::endl;
				std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
			//#endif

			if (doNextIteration)
				if (resolution <= (minRes-0.001))
					doNextIteration = false;
			}
		}

		++iter;
		last_resolution = resolution;
	}while(doNextIteration);


	if (lefttrimming == false)
	{
	  Nvoxels = 0;
	  FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
	  {
		if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n) == 0)
		{
		  DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
		}
		else
		{
		  Nvoxels++;
		  DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
		}
	  }
	#ifdef DEBUG_MASK
	  //mask.write(fnMaskOut);
	#endif
	}
	amplitudeMN.clear();
	amplitudeMS.clear();
	fftVRiesz.clear();


	double last_resolution_2 = resolution;



}
