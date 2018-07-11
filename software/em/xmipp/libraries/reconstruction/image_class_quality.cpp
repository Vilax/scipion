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

#include "image_class_quality.h"
//#define DEBUG
//#define DEBUG_MASK

void ProgClassQuality::readParams()
{
	fnImg = getParam("-i");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	sampling = getDoubleParam("--sampling");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	freq_step = getDoubleParam("--step");
	significance = getDoubleParam("--significance");
	fnMd = getParam("--md_outputdata");
	nthrs = getIntParam("--threads");
}


void ProgClassQuality::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  -i <img_file=\"\">   : Input volume");
	addParamsLine("  --mask <img_file=\"\">  : Mask defining the macromolecule");
	addParamsLine("                          :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                          :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--sampling <s=1>]   : Sampling rate (A/px)");
	addParamsLine("                            : Use -1 to disable this option");
	addParamsLine("  [--step <s=0.25>]       : The resolution is computed at a number of frequencies between mininum and");
	addParamsLine("                            : maximum resolution px/A. This parameter determines that number");
	addParamsLine("  [--minRes <s=30>]         : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]          : Maximum resolution (A)");
	addParamsLine("  [--significance <s=0.95>]    : The level of confidence for the hypothesis test.");
	addParamsLine("  [--md_outputdata <file=\".\">]  : It is a control file. The provided mask can contain voxels of noise.");
	addParamsLine("                                  : Moreover, voxels inside the mask cannot be measured due to an unsignificant");
	addParamsLine("                                  : SNR. Thus, a new mask is created. This metadata file, shows, the number of");
	addParamsLine("                                  : voxels of the original mask, and the created mask");
	addParamsLine("  [--threads <s=4>]               : Number of threads");
}


void ProgClassQuality::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;
	Image<double> img;
    img.read(fnImg);

	MultidimArray<double> &inputImg_aux = img();

	int sizeVol;
	sizeVol = XSIZE(img());

	MultidimArray<double> inputImg;
	inputImg.setXmippOrigin();
//	if(sizeVol<256)
//	{
//		inputImg_aux.setXmippOrigin();
//		inputImg_aux.window(inputImg, -128, -128, 128, 128, 0.0);
//
////		int N_smoothing;
////		N_smoothing = floor( 0.5*(256 - YSIZE(inputImg)) );
////
////		int y_size = YSIZE(inputImg);
////		int siz = y_size*0.5;
////
////		double limit_radius = (siz-N_smoothing);
////		long n=0;
////
////		int ux, uy;
////		for(int i=0; i<y_size; ++i)
////		{
////			uy = (i - siz);
////			uy *= uy;
////			for(int j=0; j<y_size; ++j)
////			{
////				ux = (j - siz);
////				ux *= ux;
////				double radius = sqrt(ux + uy);
////				if ((radius>=limit_radius))
////					DIRECT_MULTIDIM_ELEM(inputImg, n) *= 0.5*(1+cos(PI*(limit_radius-radius)/(N_smoothing)));
////				++n;
////			}
////		}
//
//	}
//	else
//	{
		inputImg = inputImg_aux;

//		int N_smoothing =60;
//		inputImg.setXmippOrigin();
//		FOR_ALL_ELEMENTS_IN_ARRAY2D(inputImg)
//		{
//
//			if (i*i+j*j < N_smoothing*N_smoothing)
//				A2D_ELEM(inputImg, i, j) = cos(0.5*sqrt(i*i + j*j));
//			else
//				A2D_ELEM(inputImg, i, j) = 0;
//		}

//	}

//		MultidimArray<double> inputImg = img();
//		inputImg.setXmippOrigin();

//	Image<double> saveimg;
//	saveimg() = inputImg;
//	saveimg.write("circul.xmp");



//	transformer_inv.setThreadsNumber(nthrs);

	FourierTransformer transformer;

	VRiesz.resizeNoCopy(inputImg);

	transformer.FourierTransform(inputImg, fftV);
	iu.initZeros(fftV);


	// Calculate u and first component of Riesz vector
	double uz, uy, ux, uy2, u2;
	long n=0;
	for(size_t i=0; i<YSIZE(fftV); ++i)
	{
		FFT_IDX2DIGFREQ(i,YSIZE(inputImg),uy);
		uy2=uy*uy;
		for(size_t j=0; j<XSIZE(fftV); ++j)
		{
			FFT_IDX2DIGFREQ(j,XSIZE(inputImg),ux);
			u2=uy2+ux*ux;
			if ((i != 0) || (j != 0))
				DIRECT_A2D_ELEM(iu,i,j) = 1.0/sqrt(u2);
			else
				DIRECT_A2D_ELEM(iu,i,j) = 1e38;
			++n;
		}
	}



//	#ifdef DEBUG
	Image<double> saveiu;
	saveiu = 1/iu;
	saveiu.write("iu.vol");
//	#endif

	//Prepare mask
	MultidimArray<int> &pMask=mask();

	if (fnMask != "")
	{
		mask.read(fnMask);
		mask().setXmippOrigin();

		R = floor(XSIZE(inputImg)*0.4);


		std::cout << "R= " << R << std::endl;
		NVoxelsOriginalMask = 0;

		FOR_ALL_ELEMENTS_IN_ARRAY2D(pMask)
		{
			if (A2D_ELEM(pMask, i, j) == 1)
				++NVoxelsOriginalMask;
			if (i*i+j*j > R*R)
				A2D_ELEM(pMask, i, j) = -1;
		}
	}
	else
	{
		std::cout << "Error: a mask ought to be provided" << std::endl;
		exit(0);
	}


//	saveiu() = inputImg;
//	saveiu.write("circulo.xmp");

//	#ifdef DEBUG_MASK
	mask.write("mask.vol");
//	#endif

	fftN=&fftV;



	double u;
	int size_fourier = YSIZE(fftV);
	freq_fourier.initZeros(size_fourier);

	int size = XSIZE(inputImg);
	img.clear();
	std::cout << "XSIZE= " << size << "  YSIZE= " << size_fourier << std::endl;

	VEC_ELEM(freq_fourier,0) = 1e-38;
	for(size_t k=1; k<size_fourier; ++k)
	{
		FFT_IDX2DIGFREQ(k,size, u);
		VEC_ELEM(freq_fourier,k) = u;
		std::cout << k << " " << VEC_ELEM(freq_fourier,k) << std::endl;
	}
}



void ProgClassQuality::amplitudeMonogenicSignal(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, FileName fnDebug)
{
	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(freq-freqH);


	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double iun=DIRECT_MULTIDIM_ELEM(iu,n);
		double un=1.0/iun;
		if (freqH<=un && un<=freq)
		{
			//double H=0.5*(1+cos((un-w1)*ideltal));
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
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);


//	#ifdef DEBUG
	Image<double> filteredvolume;
	filteredvolume = VRiesz;
	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
//	#endif

	amplitude.resizeNoCopy(VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	double uy, ux;
	n=0;
	for(size_t i=0; i<YSIZE(myfftV); ++i)
	{
//		uy = VEC_ELEM(freq_fourier,i);
		for(size_t j=0; j<XSIZE(myfftV); ++j)
		{
			ux = VEC_ELEM(freq_fourier,j);
			A2D_ELEM(fftVRiesz, i, j) = ux*A2D_ELEM(fftVRiesz_aux, i, j);
//			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
//			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
			++n;
		}
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	filteredvolume = VRiesz;
	filteredvolume.write(formatString("RieszX_%i.vol", count));

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	std::cout << "YSIZE(myfftV) = "<<  YSIZE(myfftV) << "  XSIZE(myfftV) = " << XSIZE(myfftV) << std::endl;

	n=0;
	fftVRiesz.initZeros(myfftV);
	for(size_t i=0; i<YSIZE(myfftV); ++i)
	{
		uy = VEC_ELEM(freq_fourier,i);
		for(size_t j=0; j<XSIZE(myfftV); ++j)
		{
			A2D_ELEM(fftVRiesz, i, j) = uy*A2D_ELEM(fftVRiesz_aux, i, j);
//			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
			++n;
		}
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	filteredvolume = VRiesz;
	filteredvolume.write(formatString("RieszY_%i.vol", count));

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
	}


//	#ifdef DEBUG

	FileName iternumber;
	if (fnDebug.c_str() != "")
	{
	Image<double> saveImg;
	saveImg = amplitude;

	iternumber = formatString("_Amplitude_%i.vol", count);
	saveImg.write(fnDebug+iternumber);
	saveImg.clear();
	}
//	#endif // DEBUG
//


	// Low pass filter the monogenic amplitude
	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);
	double raised_w = PI/(freqL-freq);

	n=0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
	{
		double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//		std::cout << "un = " << un << "  freqL = " << freqL << " freq = " << freq << std::endl;
		if ((freqL)>=un && un>=freq)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-freq)));
		}
		else
		{
			if (un>freqL)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
			}
		}
	}
	transformer_inv.inverseFourierTransform();


//	#ifdef DEBUG
	Image<double> saveImg2;
//	FileName iternumber;
	saveImg2 = amplitude;
	FileName fnSaveImg2;
	if (fnDebug.c_str() != "")
	{
		iternumber = formatString("_Filtered_Amplitude_%i.vol", count);
		saveImg2.write(fnDebug+iternumber);
	}
	saveImg2.clear();
//	#endif // DEBUG
}


void ProgClassQuality::resolution2eval(int &count_res, double step,
								double &resolution, double &last_resolution,
								double &freq, double &freqL,
								int &last_fourier_idx,
								bool &continueIter,	bool &breakIter,
								bool &doNextIteration)
{
	resolution = maxRes - count_res*step;
	freq = sampling/resolution;
	std::cout << "Res = " << resolution << " " << freq << std::endl;
	++count_res;

	double Nyquist = 2*sampling;
	double aux_frequency;
	int fourier_idx;

	DIGFREQ2FFT_IDX(freq, YSIZE(VRiesz), fourier_idx);

//	std::cout << "Resolution = " << resolution << "   iter = " << count_res-1 << std::endl;
//	std::cout << "freq = " << freq << "   Fourier index = " << fourier_idx << std::endl;

	FFT_IDX2DIGFREQ(fourier_idx, YSIZE(VRiesz), aux_frequency);

	freq = aux_frequency;

	if (fourier_idx == last_fourier_idx)
	{
//		std::cout << "entro en el if"  << std::endl;
		continueIter = true;
		return;
	}

	last_fourier_idx = fourier_idx;
	resolution = sampling/aux_frequency;


	if (count_res == 0)
		last_resolution = resolution;

	if ( ( resolution<Nyquist ))// || (resolution > last_resolution) )
	{
		breakIter = true;
		return;
	}


	freqL = sampling/(resolution + step);

	int fourier_idx_2;

	DIGFREQ2FFT_IDX(freqL, XSIZE(VRiesz), fourier_idx_2);

	if (fourier_idx_2 == fourier_idx)
	{
		if (fourier_idx > 0){
			//std::cout << " index low =  " << (fourier_idx - 1) << std::endl;
			FFT_IDX2DIGFREQ(fourier_idx - 1, XSIZE(VRiesz), freqL);
		}
		else{
			freqL = sampling/(resolution + step);
		}
	}


}


void ProgClassQuality::run()
{
	produceSideInfo();


	Image<double> outputResolution;
	outputResolution().initZeros(VRiesz);

	MultidimArray<int> &pMask = mask();
	MultidimArray<double> amplitudeMS, amplitudeMN;

	std::cout << "Looking for maximum frequency ..." << std::endl;
	double criticalZ=icdf_gauss(significance);
	double criticalW=-1;
	double resolution, resolution_2, last_resolution = 10000;  //A huge value for achieving
												//last_resolution < resolution
	double freq, freqH, freqL, resVal, counter;
	double max_meanS = -1e38;
	double cut_value = 0.025;

	double R_ = freq_step;

	if (R_<0.25)
		R_=0.25;

	double Nyquist = 2*sampling;
	if (minRes<2*sampling)
		minRes = Nyquist;

	bool doNextIteration=true;
	bool lefttrimming = false;
	int fourier_idx, last_fourier_idx = -1, fourier_idx_2;

	//A first MonoRes estimation to get an accurate mask

	double mean_Signal, mean_noise;


	int count_res = 0;
	FileName fnDebug;

	int iter=0;
	std::vector<double> list, passTest;

	std::cout << "Analyzing frequencies" << std::endl;
	std::vector<double> noiseValues;

	MetaData mdPoints;
	size_t objId;


	do
	{
		bool continueIter = false;
		bool breakIter = false;

		resolution2eval(count_res, R_,
						resolution, last_resolution,
						freq, freqH,
						last_fourier_idx, continueIter, breakIter, doNextIteration);

		if (continueIter)
			continue;

		if (breakIter)
			break;

		std::cout << "resolution = " << resolution << std::endl;

		list.push_back(resolution);

		long idxPass;
		if (iter <2)
		{
			idxPass = 0;
			resolution_2 = list[0];
		}
		else
		{
			resolution_2 = list[iter - 2];
			idxPass = iter - 2;
		}

		fnDebug = "Signal";

		freqL = freq + 0.01;

		amplitudeMonogenicSignal(fftV, freq, freqH, freqL, amplitudeMS, iter, fnDebug);


		double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
		noiseValues.clear();

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
			if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
			{
				sumS  += amplitudeValue;
				sumS2 += amplitudeValue*amplitudeValue;
				++NS;
			}
			else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
			{
				noiseValues.push_back(amplitudeValue);
				sumN  += amplitudeValue;
				sumN2 += amplitudeValue*amplitudeValue;
				++NN;
			}
		}


//		#ifdef DEBUG
		std::cout << "NS" << NS << std::endl;
		std::cout << "NVoxelsOriginalMask" << NVoxelsOriginalMask << std::endl;
		std::cout << "NS/NVoxelsOriginalMask = " << NS/NVoxelsOriginalMask << std::endl;
//		#endif

//		std::cout << "NS = " << NS << "  NVoxelsOriginalMask" << NVoxelsOriginalMask << std::endl;
		if ( (NS/NVoxelsOriginalMask)<cut_value ) //when the 2.5% is reached then the iterative process stops
		{
			std::cout << "Search of resolutions stopped due to mask has been completed" << std::endl;
			doNextIteration =false;
			Nvoxels = 0;

			#ifdef DEBUG_MASK
			mask.write("partial_mask.vol");
			#endif
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
			double sigma2S=sumS2/NS-meanS*meanS;
			double meanN=sumN/NN;
			double sigma2N=sumN2/NN-meanN*meanN;

			// Check local resolution
			double thresholdNoise;
			std::sort(noiseValues.begin(),noiseValues.end());
			//thresholdNoise = noiseValues[size_t(noiseValues.size()*significance)];
			thresholdNoise = meanN+criticalZ*sqrt(sigma2N);
//			#ifdef DEBUG
			  std::cout << "Iteration = " << iter << ",   Resolution= " << resolution <<
					  ",   Signal = " << meanS << ",   Noise = " << meanN << ",  Threshold = "
					  << thresholdNoise << std::endl;
//			#endif

			double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);

			if (meanS>max_meanS)
				max_meanS = meanS;

			if (meanS<0.001*max_meanS)
			{
				std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
				break;
			}

			long goodPoints = 0;
			long badPoints = 0;

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
			{
				if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
					if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
					{
						DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
						++goodPoints;
					}
					else{
						DIRECT_MULTIDIM_ELEM(pMask, n) += 1;
						if (DIRECT_MULTIDIM_ELEM(pMask, n) >2)
						{
							++badPoints;
							DIRECT_MULTIDIM_ELEM(pMask, n) = -1;
						}
					}
			}

			objId = mdPoints.addObject();
			mdPoints.setValue(MDL_IMAGE, fnOut, objId);
			mdPoints.setValue(MDL_COUNT, (size_t) goodPoints, objId);
			mdPoints.setValue(MDL_RESOLUTION_FREQ, resolution, objId);

			passTest.push_back(goodPoints);

			// Is the mean inside the signal significantly different from the noise?
			z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
			#ifdef DEBUG
				std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
				std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
				std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
				std::cout << "  z=" << z << " (" << criticalZ << ")" << std::endl;
			#endif
			if (z<criticalZ)
			{
				criticalW = freq;
				std::cout << "Search stopped due to z>Z (hypothesis test)" << std::endl;
				doNextIteration=false;
			}
			if (doNextIteration)
			{
				if (resolution <= (minRes-0.001))
					doNextIteration = false;
			}

		}
		iter++;
		last_resolution = resolution;
	} while (doNextIteration);

	amplitudeMN.clear();
	amplitudeMS.clear();

	mdPoints.write(fnOut);

//	MetaData md;
//
//	objId = md.addObject();
//	md.setValue(MDL_IMAGE, fnOut, objId);
//	md.setValue(MDL_COUNT, (size_t) NVoxelsOriginalMask, objId);
//	md.setValue(MDL_COUNT2, (size_t) Nvoxels, objId);
//
//	md.write(fnMd);
}
