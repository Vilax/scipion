
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

#include "resolution_directional.h"
//#define DEBUG
//#define DEBUG_MASK
//#define DEBUG_DIR
//#define DEBUG_FILTER
#define MONO_AMPLITUDE
//define DEBUG_SYMMETRY

void ProgResDir::readParams()
{
	fnVol = getParam("--vol");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	sampling = getDoubleParam("--sampling_rate");
	ang_sampling = getDoubleParam("--angular_sampling");
	R = getDoubleParam("--volumeRadius");
	significance = getDoubleParam("--significance");
	fnDoA = getParam("--doa_vol");
	fnDirections = getParam("--directions");
	fnradial =getParam("--radialRes");
	fnazimuthal =getParam("--azimuthalRes");
	fnMDradial =getParam("--radialAvg");
	fnMDazimuthal =getParam("--azimuthalAvg");
}


void ProgResDir::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">        : Input volume");
	addParamsLine("  [--mask <vol_file=\"\">]     : Mask defining the macromolecule");
	addParamsLine("                               :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                               :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]      : Sampling rate (A/px)");
	addParamsLine("  [--angular_sampling <s=15>]  : Angular Sampling rate (degrees)");
	addParamsLine("  [--volumeRadius <s=100>]     : This parameter determines the radius of a sphere where the volume is");
	addParamsLine("  [--significance <s=0.95>]    : The level of confidence for the hypothesis test.");
	addParamsLine("  [--doa_vol <vol_file=\"\">]  : Output filename with DoA volume");
	addParamsLine("  [--directions <vol_file=\"\">]  : Output preffered directions");
	addParamsLine("  [--radialRes <vol_file=\"\">]  : Output radial resolution map");
	addParamsLine("  [--azimuthalRes <vol_file=\"\">]  : Output azimuthal resolution map");
	addParamsLine("  [--radialAvg <vol_file=\"\">]  : Radial Average of the radial resolution map");
	addParamsLine("  [--azimuthalAvg <vol_file=\"\">]  : Radial Average of the azimuthal resolution map");
}

void ProgResDir::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;

	Image<double> V;
	V.read(fnVol);

	V().setXmippOrigin();

	//Sweeping the projection sphere
	std::cout << "Obtaining angular projections..." << std::endl;
	generateGridProjectionMatching(fnVol, ang_sampling, angles);

	FourierTransformer transformer;

	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);
	N_freq = ZSIZE(inputVol);
	maxRes = ZSIZE(inputVol);
	minRes = 2*sampling;

	transformer_inv.setThreadsNumber(4);

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
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol), ux);
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
	Image<double> AvgResoltion;
	AvgResoltion().resizeNoCopy(inputVol);
	AvgResoltion().initZeros();
	AvgResoltion.write(fnOut);
	AvgResoltion.clear();

	N_smoothing = 7;
	NVoxelsOriginalMask = 0;
	double radius = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
		if (A3D_ELEM(pMask, k, i, j) == 1)
		{
			if ((k*k + i*i + j*j)>radius)
				radius = k*k + i*i + j*j;
		}
//		std::cout << "i j k " << i << " " << j << " " << k << std::endl;

		if (A3D_ELEM(pMask, k, i, j) == 1)
			++NVoxelsOriginalMask;
		if (i*i+j*j+k*k > (R-N_smoothing)*(R-N_smoothing))
			A3D_ELEM(pMask, k, i, j) = -1;
	}
	Rparticle = round(sqrt(radius));
	std::cout << "particle radiues = " << Rparticle << std::endl;
	size_t xrows = angles.mdimx;
//	Matrix2D<double> aaa;
//
//	aaa.initConstant(3, 4, 5);
//	MAT_ELEM(aaa, 1, 2) = 0;
//	std::cout << aaa << std::endl;
//	std::cout << "NVoxelsOriginalMask = " << NVoxelsOriginalMask << std::endl;

	resolutionMatrix.initConstant(xrows, NVoxelsOriginalMask, maxRes);


	#ifdef DEBUG_MASK
	std::cout << "-------------DEBUG-----------" <<std::endl;
	std::cout << "Next number ought to be the same than number of directions"
			<< std::endl;
	std::cout << "number of angles" << xrows << std::endl;
	std::cout << "---------END-DEBUG--------" <<std::endl;
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

	VEC_ELEM(freq_fourier,0) = 1e-38;
	for(size_t k=1; k<size_fourier; ++k)
	{
		FFT_IDX2DIGFREQ(k,size, u);
		VEC_ELEM(freq_fourier,k) = u;
//		std::cout << "freq_fourier = " << sampling/u  << std::endl;
	}
}


void ProgResDir::generateGridProjectionMatching(FileName fnVol_, double smprt,
												Matrix2D<double> &angles)
{
  /*
	FileName fnGallery;
	FileName fnanglesmd = "angles.doc";
	// Generate projections
	fnGallery=formatString("angles.stk");
	String args=formatString("-i %s -o %s --sampling_rate %f",
			fnVol_.c_str(), fnGallery.c_str(), smprt); //fnGallery.c_str(),smprt);
	String cmd=(String)"xmipp_angular_project_library " + args;
	std::cout << cmd << std::endl;
	system(cmd.c_str());
	MetaData md;
	md.read(fnanglesmd);
	size_t md_size = md.size();
	//				if ((fabs(uz) <= 0.1) || (fabs(uy) <= 0.1) )//|| (fabs(uz) <= 0.1))
	//					DIRECT_MULTIDIM_ELEM(iu,n) = ux;
	Matrix2D<double> aux_angles(2,md_size);
	size_t count = 0;
	double rot, tilt;
	FOR_ALL_OBJECTS_IN_METADATA(md)
	{
		md.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
		md.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
		if ( (rot>=0) && (tilt>=0) )
		{
			MAT_ELEM(aux_angles,0, count) = rot;
			if (tilt>90)
				MAT_ELEM(aux_angles,1, count) = tilt-180;
			else
				MAT_ELEM(aux_angles,1, count) = tilt;
			++count;
		}
	}
	N_directions = count;
	angles.initZeros(4,N_directions);
	//The loop beging in 1 instead of 0 avoiding the repeated direction rot=0 tilt=0.
	for (size_t k = 1; k<(N_directions); k++)
	{
		MAT_ELEM(angles, 0, k) = MAT_ELEM(aux_angles,0, k);
		MAT_ELEM(angles, 1, k) = MAT_ELEM(aux_angles,1, k);
		std::cout << "k=" << k << "  rot = "  << MAT_ELEM(angles, 0, k) <<
								  "  tilt = " << MAT_ELEM(angles, 1, k) << std::endl;
	}
	//TODO: check if the angles output are correct
	*/
	angles.initZeros(2,82);
	MAT_ELEM(angles, 0, 0) = 0.000000;	 MAT_ELEM(angles, 1, 0) = 0.000000;
	MAT_ELEM(angles, 0, 1) = 36.000000;	 MAT_ELEM(angles, 1, 1) = 15.858741;
	MAT_ELEM(angles, 0, 2) = 36.000000;	 MAT_ELEM(angles, 1, 2) = 31.717482;
	MAT_ELEM(angles, 0, 3) = 36.000000;	 MAT_ELEM(angles, 1, 3) = 47.576224;
	MAT_ELEM(angles, 0, 4) = 36.000000;	 MAT_ELEM(angles, 1, 4) = 63.434965;
	MAT_ELEM(angles, 0, 5) = 62.494295;	 MAT_ELEM(angles, 1, 5) = -76.558393;
	MAT_ELEM(angles, 0, 6) = 54.000000;	 MAT_ELEM(angles, 1, 6) = 90.000000;
	MAT_ELEM(angles, 0, 7) = 45.505705;	 MAT_ELEM(angles, 1, 7) = 76.558393;
	MAT_ELEM(angles, 0, 8) = 108.000000;	 MAT_ELEM(angles, 1, 8) = 15.858741;
	MAT_ELEM(angles, 0, 9) = 108.000000;	 MAT_ELEM(angles, 1, 9) = 31.717482;
	MAT_ELEM(angles, 0, 10) = 108.000000;	 MAT_ELEM(angles, 1, 10) = 47.576224;
	MAT_ELEM(angles, 0, 11) = 108.000000;	 MAT_ELEM(angles, 1, 11) = 63.434965;
	MAT_ELEM(angles, 0, 12) = 134.494295;	 MAT_ELEM(angles, 1, 12) = -76.558393;
	MAT_ELEM(angles, 0, 13) = 126.000000;	 MAT_ELEM(angles, 1, 13) = 90.000000;
	MAT_ELEM(angles, 0, 14) = 117.505705;	 MAT_ELEM(angles, 1, 14) = 76.558393;
	MAT_ELEM(angles, 0, 15) = 144.000000;	 MAT_ELEM(angles, 1, 15) = -15.858741;
	MAT_ELEM(angles, 0, 16) = 144.000000;	 MAT_ELEM(angles, 1, 16) = -31.717482;
	MAT_ELEM(angles, 0, 17) = 144.000000;	 MAT_ELEM(angles, 1, 17) = -47.576224;
	MAT_ELEM(angles, 0, 18) = 144.000000;	 MAT_ELEM(angles, 1, 18) = -63.434965;
	MAT_ELEM(angles, 0, 19) = 170.494295;	 MAT_ELEM(angles, 1, 19) = 76.558393;
	MAT_ELEM(angles, 0, 20) = 162.000000;	 MAT_ELEM(angles, 1, 20) = 90.000000;
	MAT_ELEM(angles, 0, 21) = 153.505705;	 MAT_ELEM(angles, 1, 21) = -76.558393;
	MAT_ELEM(angles, 0, 22) = 72.000000;	 MAT_ELEM(angles, 1, 22) = -15.858741;
	MAT_ELEM(angles, 0, 23) = 72.000000;	 MAT_ELEM(angles, 1, 23) = -31.717482;
	MAT_ELEM(angles, 0, 24) = 72.000000;	 MAT_ELEM(angles, 1, 24) = -47.576224;
	MAT_ELEM(angles, 0, 25) = 72.000000;	 MAT_ELEM(angles, 1, 25) = -63.434965;
	MAT_ELEM(angles, 0, 26) = 98.494295;	 MAT_ELEM(angles, 1, 26) = 76.558393;
	MAT_ELEM(angles, 0, 27) = 90.000000;	 MAT_ELEM(angles, 1, 27) = 90.000000;
	MAT_ELEM(angles, 0, 28) = 81.505705;	 MAT_ELEM(angles, 1, 28) = -76.558393;
	MAT_ELEM(angles, 0, 29) = 0.000000;	 MAT_ELEM(angles, 1, 29) = -15.858741;
	MAT_ELEM(angles, 0, 30) = 0.000000;	 MAT_ELEM(angles, 1, 30) = -31.717482;
	MAT_ELEM(angles, 0, 31) = 0.000000;	 MAT_ELEM(angles, 1, 31) = -47.576224;
	MAT_ELEM(angles, 0, 32) = 0.000000;	 MAT_ELEM(angles, 1, 32) = -63.434965;
	MAT_ELEM(angles, 0, 33) = 26.494295;	 MAT_ELEM(angles, 1, 33) = 76.558393;
	MAT_ELEM(angles, 0, 34) = 18.000000;	 MAT_ELEM(angles, 1, 34) = 90.000000;
	MAT_ELEM(angles, 0, 35) = 9.505705;	 MAT_ELEM(angles, 1, 35) = -76.558393;
	MAT_ELEM(angles, 0, 36) = 12.811021;	 MAT_ELEM(angles, 1, 36) = 42.234673;
	MAT_ELEM(angles, 0, 37) = 18.466996;	 MAT_ELEM(angles, 1, 37) = 59.620797;
	MAT_ELEM(angles, 0, 38) = 0.000000;	 MAT_ELEM(angles, 1, 38) = 90.000000;
	MAT_ELEM(angles, 0, 39) = 8.867209;	 MAT_ELEM(angles, 1, 39) = 75.219088;
	MAT_ELEM(angles, 0, 40) = 72.000000;	 MAT_ELEM(angles, 1, 40) = 26.565058;
	MAT_ELEM(angles, 0, 41) = 59.188979;	 MAT_ELEM(angles, 1, 41) = 42.234673;
	MAT_ELEM(angles, 0, 42) = 84.811021;	 MAT_ELEM(angles, 1, 42) = 42.234673;
	MAT_ELEM(angles, 0, 43) = 53.533003;	 MAT_ELEM(angles, 1, 43) = 59.620797;
	MAT_ELEM(angles, 0, 44) = 72.000000;	 MAT_ELEM(angles, 1, 44) = 58.282544;
	MAT_ELEM(angles, 0, 45) = 90.466996;	 MAT_ELEM(angles, 1, 45) = 59.620797;
	MAT_ELEM(angles, 0, 46) = 72.000000;	 MAT_ELEM(angles, 1, 46) = 90.000000;
	MAT_ELEM(angles, 0, 47) = 63.132791;	 MAT_ELEM(angles, 1, 47) = 75.219088;
	MAT_ELEM(angles, 0, 48) = 80.867209;	 MAT_ELEM(angles, 1, 48) = 75.219088;
	MAT_ELEM(angles, 0, 49) = 144.000000;	 MAT_ELEM(angles, 1, 49) = 26.565058;
	MAT_ELEM(angles, 0, 50) = 131.188979;	 MAT_ELEM(angles, 1, 50) = 42.234673;
	MAT_ELEM(angles, 0, 51) = 156.811021;	 MAT_ELEM(angles, 1, 51) = 42.234673;
	MAT_ELEM(angles, 0, 52) = 125.533003;	 MAT_ELEM(angles, 1, 52) = 59.620797;
	MAT_ELEM(angles, 0, 53) = 144.000000;	 MAT_ELEM(angles, 1, 53) = 58.282544;
	MAT_ELEM(angles, 0, 54) = 162.466996;	 MAT_ELEM(angles, 1, 54) = 59.620797;
	MAT_ELEM(angles, 0, 55) = 144.000000;	 MAT_ELEM(angles, 1, 55) = 90.000000;
	MAT_ELEM(angles, 0, 56) = 135.132791;	 MAT_ELEM(angles, 1, 56) = 75.219088;
	MAT_ELEM(angles, 0, 57) = 152.867209;	 MAT_ELEM(angles, 1, 57) = 75.219088;
	MAT_ELEM(angles, 0, 58) = 180.000000;	 MAT_ELEM(angles, 1, 58) = -26.565058;
	MAT_ELEM(angles, 0, 59) = 167.188979;	 MAT_ELEM(angles, 1, 59) = -42.234673;
	MAT_ELEM(angles, 0, 60) = 180.000000;	 MAT_ELEM(angles, 1, 60) = -58.282544;
	MAT_ELEM(angles, 0, 61) = 161.533003;	 MAT_ELEM(angles, 1, 61) = -59.620797;
	MAT_ELEM(angles, 0, 62) = 180.000000;	 MAT_ELEM(angles, 1, 62) = 90.000000;
	MAT_ELEM(angles, 0, 63) = 171.132791;	 MAT_ELEM(angles, 1, 63) = -75.219088;
	MAT_ELEM(angles, 0, 64) = 108.000000;	 MAT_ELEM(angles, 1, 64) = -26.565058;
	MAT_ELEM(angles, 0, 65) = 120.811021;	 MAT_ELEM(angles, 1, 65) = -42.234673;
	MAT_ELEM(angles, 0, 66) = 95.188979;	 MAT_ELEM(angles, 1, 66) = -42.234673;
	MAT_ELEM(angles, 0, 67) = 126.466996;	 MAT_ELEM(angles, 1, 67) = -59.620797;
	MAT_ELEM(angles, 0, 68) = 108.000000;	 MAT_ELEM(angles, 1, 68) = -58.282544;
	MAT_ELEM(angles, 0, 69) = 89.533003;	 MAT_ELEM(angles, 1, 69) = -59.620797;
	MAT_ELEM(angles, 0, 70) = 108.000000;	 MAT_ELEM(angles, 1, 70) = 90.000000;
	MAT_ELEM(angles, 0, 71) = 116.867209;	 MAT_ELEM(angles, 1, 71) = -75.219088;
	MAT_ELEM(angles, 0, 72) = 99.132791;	 MAT_ELEM(angles, 1, 72) = -75.219088;
	MAT_ELEM(angles, 0, 73) = 36.000000;	 MAT_ELEM(angles, 1, 73) = -26.565058;
	MAT_ELEM(angles, 0, 74) = 48.811021;	 MAT_ELEM(angles, 1, 74) = -42.234673;
	MAT_ELEM(angles, 0, 75) = 23.188979;	 MAT_ELEM(angles, 1, 75) = -42.234673;
	MAT_ELEM(angles, 0, 76) = 54.466996;	 MAT_ELEM(angles, 1, 76) = -59.620797;
	MAT_ELEM(angles, 0, 77) = 36.000000;	 MAT_ELEM(angles, 1, 77) = -58.282544;
	MAT_ELEM(angles, 0, 78) = 17.533003;	 MAT_ELEM(angles, 1, 78) = -59.620797;
	MAT_ELEM(angles, 0, 79) = 36.000000;	 MAT_ELEM(angles, 1, 79) = 90.000000;
	MAT_ELEM(angles, 0, 80) = 44.867209;	 MAT_ELEM(angles, 1, 80) = -75.219088;
	MAT_ELEM(angles, 0, 81) = 27.132791;	 MAT_ELEM(angles, 1, 81) = -75.219088;
}

//	amplitude.setXmippOrigin();
//
//
//    double raised_w = PI/(freqL-freq);
//
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
//	{
//			double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//			if (freqL>=un && un>=freq)
//				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-freq)));
//			else
//				if (un>freqL)
//					DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
//	}
//	transformer_inv.inverseFourierTransform();
//
//
//	if (fnDebug.c_str() != "")
//	{
//
//		saveImg2 = amplitude;
//		FileName iternumber = formatString("_Filtered_Amplitude_%i_%i.vol", dir, count);
//		saveImg2.write(fnDebug+iternumber);
//	}
//	saveImg2.clear();


void ProgResDir::amplitudeMonogenicSignal3D_fast(const MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, int dir, FileName fnDebug,
		double rot, double tilt)
{
	fftVRiesz.initZeros(myfftV);
//	MultidimArray<double> coneVol;
//	coneVol.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(freq-freqH);

	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
//				double iun = *ptriun;
				double un=1.0/iun;
				if (freqH<=un && un<=freq)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= DIRECT_MULTIDIM_ELEM(conefilter, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
//					DIRECT_MULTIDIM_ELEM(coneVol, n) = DIRECT_MULTIDIM_ELEM(conefilter, n);
				} else if (un>freq)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= DIRECT_MULTIDIM_ELEM(conefilter, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
//					DIRECT_MULTIDIM_ELEM(coneVol, n) = DIRECT_MULTIDIM_ELEM(conefilter, n);
				}
				++n;
			}
		}
	}

//	#ifdef DEBUG_DIR
////	if ( (count == 0) )
////	{
//		Image<double> direction;
//		direction = coneVol;
//		direction.write(formatString("cone_%i_%i.vol", dir+1, count));
////	}
//	#endif

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	#ifdef DEBUG_DIR
		Image<double> filteredvolume;
		filteredvolume = VRiesz;
		filteredvolume.write(formatString("Volumen_filtrado_%i_%i.vol", dir+1,count));
	#endif




//	amplitude.initZeros(VRiesz);
	amplitude.resizeNoCopy(VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	}

	// Calculate first component of Riesz vector
	double ux;
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
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
	}

	// Calculate second and third component of Riesz vector
	n=0;
	double uy, uz;
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
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+= DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	}

	transformer_inv.inverseFourierTransform(fftVRiesz_aux, VRiesz);


//	amplitude.setXmippOrigin();
	int z_size = ZSIZE(amplitude);
	int siz = z_size*0.5;

	double limit_radius = (siz-N_smoothing);
	n=0;
	for(int k=0; k<z_size; ++k)
	{
		uz = (k - siz);
		uz *= uz;
		for(int i=0; i<z_size; ++i)
		{
			uy = (i - siz);
			uy *= uy;
			for(int j=0; j<z_size; ++j)
			{
				ux = (j - siz);
				ux *= ux;
				DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
				DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
				double radius = sqrt(ux + uy + uz);
				if ((radius>=limit_radius) && (radius<=(z_size*0.5)))
					DIRECT_MULTIDIM_ELEM(amplitude, n) *= 0.5*(1+cos(PI*(limit_radius-radius)/(N_smoothing)));
				else if (radius>(0.5*z_size))
					DIRECT_MULTIDIM_ELEM(amplitude, n) = 0;
				++n;
			}
		}
	}
	//TODO: change (k - z_size*0.5)

		#ifdef MONO_AMPLITUDE
		Image<double> saveImg2;
		saveImg2 = amplitude;
		if (fnDebug.c_str() != "")
		{
			FileName iternumber = formatString("smoothed_volume_%i_%i.vol", dir+1, count);
			saveImg2.write(fnDebug+iternumber);
		}
		saveImg2.clear();
		#endif

	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);

	double raised_w = PI/(freqL-freq);

	n=0;
	std::cout << "freqL = " << freqL << std::endl;
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
			if (un>(freqL))
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
			}
		}
	}
	transformer_inv.inverseFourierTransform();

	#ifdef MONO_AMPLITUDE

//	if (fnDebug.c_str() != "")
//	{
		saveImg2 = amplitude;
		FileName iternumber = formatString("_Filtered_Amplitude_%i_%i.vol", dir+1, count);
		saveImg2.write(fnDebug+iternumber);
//	}
	#endif // DEBUG

		exit(0);
}


//void ProgResDir::amplitudeMonogenicSignal3D_fast(const MultidimArray< std::complex<double> > &myfftV,
//		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, int dir, FileName fnDebug,
//		double rot, double tilt)
//{
//	fftVRiesz.initZeros(myfftV);
////	MultidimArray<double> coneVol;
////	coneVol.initZeros(myfftV);
//	fftVRiesz_aux.initZeros(myfftV);
//	std::complex<double> J(0,1);
//
//	// Filter the input volume and add it to amplitude
//	long n=0;
//	double ideltal=PI/(freq-freqH);
//
//
//	std::complex<double> *ptrfftVRiesz=&DIRECT_MULTIDIM_ELEM(fftVRiesz, 0);
//	std::complex<double> *ptrmyfftV=&DIRECT_MULTIDIM_ELEM(myfftV, 0);
//	std::complex<double> *ptrfftVRiesz_aux=&DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, 0);
//	double *ptriun=&DIRECT_MULTIDIM_ELEM(iu, 0);
//	std::complex<double> aux1;
//
//	for(size_t k=0; k<ZSIZE(myfftV); ++k)
//	{
//		for(size_t i=0; i<YSIZE(myfftV); ++i)
//		{
//			for(size_t j=0; j<XSIZE(myfftV); ++j)
//			{
////				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
//				double iun = *ptriun;
//				double un=1.0/iun;
////				if (freqH<=un && un<=freq)
////				{
////					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//////					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= DIRECT_MULTIDIM_ELEM(conefilter, n);
////					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
////					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
////					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
////					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
//////					DIRECT_MULTIDIM_ELEM(coneVol, n) = DIRECT_MULTIDIM_ELEM(conefilter, n);
////				} else if (un>freq)
////				{
////					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//////					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= DIRECT_MULTIDIM_ELEM(conefilter, n);
////					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
////					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
////					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
//////					DIRECT_MULTIDIM_ELEM(coneVol, n) = DIRECT_MULTIDIM_ELEM(conefilter, n);
////				}
//
//				if (freqH<=un && un<=freq)
//				{
//					aux1 = *ptrmyfftV;
//					aux1 *= 0.5*(1+cos((un-freq)*ideltal));
//					*ptrfftVRiesz *= aux1;
//					aux1 *= -J;
//					aux1 *= iun;
//					*ptrfftVRiesz_aux = aux1;
//				} else if (un>freq)
//				{
//					*ptrfftVRiesz = *ptrmyfftV;
//					aux1 = -J;
//					aux1 *= iun;
//					aux1 *= *ptrmyfftV;
//					*ptrfftVRiesz_aux = aux1;
//				}
//				ptrfftVRiesz++;
//				ptrmyfftV++;
//				ptrfftVRiesz_aux++;
//				ptriun++;
////				++n;
//			}
//		}
//	}
//
////	#ifdef DEBUG_DIR
//////	if ( (count == 0) )
//////	{
////		Image<double> direction;
////		direction = coneVol;
////		direction.write(formatString("cone_%i_%i.vol", dir+1, count));
//////	}
////	#endif
//
//	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
//
//	#ifdef DEBUG_DIR
//		Image<double> filteredvolume;
//		filteredvolume = VRiesz;
//		filteredvolume.write(formatString("Volumen_filtrado_%i_%i.vol", dir+1,count));
//	#endif
//
//
//
//
////	amplitude.initZeros(VRiesz);
//	amplitude.resizeNoCopy(VRiesz);
//	double *ptrVRiesz=&DIRECT_MULTIDIM_ELEM(VRiesz, 0);
//	double *ptramplitude=&DIRECT_MULTIDIM_ELEM(amplitude, 0);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
//	{
////		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
//		*ptramplitude++ = (*ptrVRiesz)*(*ptrVRiesz);
//		ptrVRiesz++;
////		ptramplitude++;
//	}
//
//	// Calculate first component of Riesz vector
//	double ux;
////	n=0;
//	ptrfftVRiesz=&DIRECT_MULTIDIM_ELEM(fftVRiesz, 0);
//	ptrfftVRiesz_aux=&DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, 0);
//	for(size_t k=0; k<ZSIZE(myfftV); ++k)
//	{
//		for(size_t i=0; i<YSIZE(myfftV); ++i)
//		{
//			for(size_t j=0; j<XSIZE(myfftV); ++j)
//			{
////				ux = VEC_ELEM(freq_fourier,j);
////				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
////				++n;
//				*ptrfftVRiesz = *ptrfftVRiesz_aux;
//				*ptrfftVRiesz++ *= VEC_ELEM(freq_fourier,j);
////				ptrfftVRiesz++;
//				ptrfftVRiesz_aux++;
//			}
//		}
//	}
//
//	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
//
//	ptrVRiesz=&DIRECT_MULTIDIM_ELEM(VRiesz, 0);
//	ptramplitude=&DIRECT_MULTIDIM_ELEM(amplitude, 0);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
//	{
////		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
//		*ptramplitude++ += (*ptrVRiesz)*(*ptrVRiesz);
//		ptrVRiesz++;
////		ptramplitude++;
//	}
///*
//	// Calculate second and third component of Riesz vector
////	n=0;
//	ptrfftVRiesz=&DIRECT_MULTIDIM_ELEM(fftVRiesz, 0);
//	ptrfftVRiesz_aux=&DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, 0);
//	double uy, uz;
//	for(size_t k=0; k<ZSIZE(myfftV); ++k)
//	{
//		uz = VEC_ELEM(freq_fourier,k);
//		for(size_t i=0; i<YSIZE(myfftV); ++i)
//		{
//			uy = VEC_ELEM(freq_fourier,i);
//			for(size_t j=0; j<XSIZE(myfftV); ++j)
//			{
//				*ptrfftVRiesz = uz;
//				*ptrfftVRiesz++ *= (*ptrfftVRiesz_aux);
//				*ptrfftVRiesz_aux++ *= uy;
////				ptrfftVRiesz++;
////				ptrfftVRiesz_aux++;
////				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uz*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
////				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
////				++n;
//			}
//		}
//	}
//	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
//
//	ptrVRiesz=&DIRECT_MULTIDIM_ELEM(VRiesz, 0);
//	ptramplitude=&DIRECT_MULTIDIM_ELEM(amplitude, 0);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
//	{
////		DIRECT_MULTIDIM_ELEM(amplitude,n)+= DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
//		*ptramplitude++ += (*ptrVRiesz)*(*ptrVRiesz);
//		ptrVRiesz++;
////		ptramplitude++;
//	}
//
//	transformer_inv.inverseFourierTransform(fftVRiesz_aux, VRiesz);
//*/
//	ptrfftVRiesz=&DIRECT_MULTIDIM_ELEM(fftVRiesz, 0);
//	ptrfftVRiesz_aux=&DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, 0);
//	double uy, uz;
//	for(size_t k=0; k<ZSIZE(myfftV); ++k)
//	{
//		for(size_t i=0; i<YSIZE(myfftV); ++i)
//		{
//			uy = VEC_ELEM(freq_fourier,i);
//			for(size_t j=0; j<XSIZE(myfftV); ++j)
//			{
//				*ptrfftVRiesz = *ptrfftVRiesz_aux;
//				*ptrfftVRiesz++ *= uy;
////				ptrfftVRiesz++;
//			}
//		}
//	}
//	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
//
//	ptrVRiesz=&DIRECT_MULTIDIM_ELEM(VRiesz, 0);
//	ptramplitude=&DIRECT_MULTIDIM_ELEM(amplitude, 0);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
//	{
////		DIRECT_MULTIDIM_ELEM(amplitude,n)+= DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
//		*ptramplitude++ += (*ptrVRiesz)*(*ptrVRiesz);
//		ptrVRiesz++;
////		ptramplitude++;
//	}
//
//	ptrfftVRiesz=&DIRECT_MULTIDIM_ELEM(fftVRiesz, 0);
//	ptrfftVRiesz_aux=&DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, 0);
//
//	for(size_t k=0; k<ZSIZE(myfftV); ++k)
//	{
//		uz = VEC_ELEM(freq_fourier,k);
//		for(size_t i=0; i<YSIZE(myfftV); ++i)
//		{
//
//			for(size_t j=0; j<XSIZE(myfftV); ++j)
//			{
//				*ptrfftVRiesz = *ptrfftVRiesz_aux;
//				*ptrfftVRiesz++ *= uy;
////				ptrfftVRiesz++;
//			}
//		}
//	}
//	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
//
//	amplitude.setXmippOrigin();
//	int z_size = ZSIZE(amplitude);
////	int x_size = XSIZE(amplitude);
////	int y_size = YSIZE(amplitude);
//	int siz = z_size*0.5;
//
//	double limit_radius = (siz-N_smoothing);
//	n=0;
//	ptrVRiesz=&DIRECT_MULTIDIM_ELEM(VRiesz, 0);
//	ptramplitude=&DIRECT_MULTIDIM_ELEM(amplitude, 0);
//	for(int k=0; k<z_size; ++k)
//	{
//		uz = (k - siz);
//		uz *= uz;
//		for(int i=0; i<z_size; ++i)
//		{
//			uy = (i - siz);
//			uy *= uy;
//			for(int j=0; j<z_size; ++j)
//			{
//				ux = (j - siz);
//				ux *= ux;
////				DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
////				DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
//				*ptramplitude += (*ptrVRiesz)*(*ptrVRiesz);
//				*ptramplitude = sqrt(*ptramplitude);
//				double radius = sqrt(ux + uy + uz);
////				if ((radius>=limit_radius) && (radius<=(z_size*0.5)))
////					DIRECT_MULTIDIM_ELEM(amplitude, n) *= 0.5*(1+cos(PI*(limit_radius-radius)/(N_smoothing)));
////				else if (radius>(0.5*z_size))
////					DIRECT_MULTIDIM_ELEM(amplitude, n) = 0;
//				if ((radius>=limit_radius) && (radius<=(siz)))
//					*ptramplitude *= 0.5*(1+cos(PI*(limit_radius-radius)/(N_smoothing)));
//				else if (radius>(siz))
//					*ptramplitude = 0;
////				++n;
//				ptrVRiesz++;
//				ptramplitude++;
//			}
//		}
//	}
//	//TODO: change (k - z_size*0.5)
//
//		#ifdef MONO_AMPLITUDE
//		Image<double> saveImg2;
//		saveImg2 = amplitude;
//		if (fnDebug.c_str() != "")
//		{
//			FileName iternumber = formatString("smoothed_volume_%i_%i.vol", dir+1, count);
//			saveImg2.write(fnDebug+iternumber);
//		}
//		saveImg2.clear();
//		#endif
//
//	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);
//
//	double raised_w = PI/(freqL-freq);
//
//	ptrfftVRiesz=&DIRECT_MULTIDIM_ELEM(fftVRiesz, 0);
//	ptriun=&DIRECT_MULTIDIM_ELEM(iu, 0);
//	n=0;
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
//	{
//		double un = 1.0/(*ptriun);
////		double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//		if ((freqL)>=un && un>=freq)
//		{
////			DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-freq)));
//			*ptrfftVRiesz *= 0.5*(1 + cos(raised_w*(un-freq)));
//		}
//		else
//		{
//			if (un>(freqL))
//			{
////				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
//				*ptrfftVRiesz = 0;
//			}
//		}
//		ptriun++;
//		ptrfftVRiesz++;
////		n++;
//	}
//	transformer_inv.inverseFourierTransform();
//
//	#ifdef MONO_AMPLITUDE
//
////	if (fnDebug.c_str() != "")
////	{
//		saveImg2 = amplitude;
//		FileName iternumber = formatString("_Filtered_Amplitude_%i_%i.vol", dir+1, count);
//		saveImg2.write(fnDebug+iternumber);
////	}
//	#endif // DEBUG
//
//}

void ProgResDir::defineCone(MultidimArray< std::complex<double> > &myfftV,
		MultidimArray< std::complex<double> > &conefilter, double rot, double tilt)
{
//	conefilter.initZeros(myfftV);
	conefilter = myfftV;
	// Filter the input volume and add it to amplitude

//	MultidimArray<double> conetest;
//	conetest.initZeros(myfftV);
//	#ifdef DEBUG_DIR
//	MultidimArray<double> coneVol;
//	coneVol.initZeros(iu);
//	#endif

	double x_dir, y_dir, z_dir;

	x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir = cos(tilt*PI/180);

	double uz, uy, ux;
	long n = 0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		uz = VEC_ELEM(freq_fourier,k);
		uz *= z_dir;
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier,i);
			uy *= y_dir;
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				ux = VEC_ELEM(freq_fourier,j);
				ux *= x_dir;

				//BE CAREFULL with the order
				//double dotproduct = (uy*x_dir + ux*y_dir + uz*z_dir)*iun;
				iun *= (ux + uy + uz);
				double acosine = acos(fabs(iun));
				//TODO: remove fabs

				//4822.53 mean a smoothed cone angle of 20 degrees
				double arg_exp = acosine*acosine*acosine*acosine*4822.53;
				DIRECT_MULTIDIM_ELEM(conefilter, n) *= exp(-arg_exp);
//				DIRECT_MULTIDIM_ELEM(conetest, n) = exp(-arg_exp);
				++n;
			}
		}
	}

//	Image<double> saveImg2;
//	saveImg2 = conetest;
//	saveImg2.write("cono.vol");

}

void ProgResDir::diagSymMatrix3x3(Matrix2D<double> A,
					Matrix1D<double> &eigenvalues, Matrix2D<double> &eigenvectors)
{
	Matrix2D<double> B;
	B.initZeros(3,3);

	MAT_ELEM(B,0,0) = 1;
	MAT_ELEM(B,1,1) = 1;
	MAT_ELEM(B,2,2) = 1;

	generalizedEigs(A, B, eigenvalues, eigenvectors);
}

void ProgResDir::resolution2eval_(int &fourier_idx, double min_step,
								double &resolution, double &last_resolution,
								int &last_fourier_idx,
								double &freq, double &freqL, double &freqH,
								bool &continueIter, bool &breakIter, bool &doNextIteration)
{
	int volsize = ZSIZE(VRiesz);

	FFT_IDX2DIGFREQ(fourier_idx, volsize, freq);

	resolution = sampling/freq;
//	std::cout << "res = " << resolution << std::endl;
//	std::cout << "min_step = " << min_step << std::endl;

	//TODO: I am sure that the abs can be removed
	if ( fabs(resolution - last_resolution)<min_step )
	{
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

//	std::cout << "freq_H = " << freqH << std::endl;
//	std::cout << "freq_L = " << freqL << std::endl;

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
//	std::cout << "resolution = " << resolution << "  resolutionL = " <<
//				sampling/(freqL) << "  resolutionH = " << sampling/freqH
//				<< "  las_res = " << last_resolution << std::endl;
	last_fourier_idx = fourier_idx;
	++fourier_idx;
}

void ProgResDir::removeOutliers(Matrix2D<double> &anglesMat, Matrix2D<double> &resolutionMat)
{
	double x1, y1, z1, x2, y2, z2, distance, resolution, sigma,
				rot, tilt, threshold, sigma2, lastMinDistance;
	double meandistance = 0, distance_2 = 0;
	int xrows = angles.mdimx, N=0;

//	std::cout << "xrows = " << xrows << std::endl;

	double criticalZ = icdf_gauss(significance);

	for (int k = 0; k<NVoxelsOriginalMask; ++k)
	{
		meandistance = 0;
		distance_2 = 0;

		for (int i = 0; i<xrows; ++i)
		{
			resolution = MAT_ELEM(resolutionMat, i, k);
//			rot = MAT_ELEM(anglesMat,0, i)*PI/180;
//			tilt = MAT_ELEM(anglesMat,1, i)*PI/180;
			x1 = resolution*MAT_ELEM(trigProducts, 0, i);
			y1 = resolution*MAT_ELEM(trigProducts, 1, i);
			z1 = resolution*MAT_ELEM(trigProducts, 2, i);
			lastMinDistance = 1e38;
			for (int j = 0; j<xrows; ++j)
			{
//				rot = MAT_ELEM(anglesMat,0, j)*PI/180;
//				tilt = MAT_ELEM(anglesMat,1, j)*PI/180;
				resolution = MAT_ELEM(resolutionMat, j, k);
				x2 = resolution*MAT_ELEM(trigProducts, 0, j);
				y2 = resolution*MAT_ELEM(trigProducts, 1, j);
				z2 = resolution*MAT_ELEM(trigProducts, 2, j);

				if (i != j)
				{
					distance = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
//					std::cout << "distance = " << distance << std::endl;
					if (distance < lastMinDistance)
						lastMinDistance = distance;
				}
			}
			meandistance += lastMinDistance;
			distance_2 += lastMinDistance*lastMinDistance;
		}

		meandistance = (meandistance)/((double) xrows);
		sigma2 = distance_2/((double) xrows) - meandistance*meandistance;

		threshold = meandistance + criticalZ*sqrt(sigma2);

		for (int i = 0; i<xrows; ++i)
		{
			resolution = MAT_ELEM(resolutionMat, i, k);
//			rot = MAT_ELEM(anglesMat,0, i)*PI/180;
//			tilt = MAT_ELEM(anglesMat,1, i)*PI/180;
			x1 = resolution*MAT_ELEM(trigProducts, 0, i);
			y1 = resolution*MAT_ELEM(trigProducts, 1, i);
			z1 = resolution*MAT_ELEM(trigProducts, 2, i);
			lastMinDistance = 1e38;
			for (int j = 0; j<xrows; ++j)
			{
//				rot = MAT_ELEM(anglesMat,0, j)*PI/180;
//				tilt = MAT_ELEM(anglesMat,1, j)*PI/180;
				resolution = MAT_ELEM(resolutionMat, j, k);
				x2 = resolution*MAT_ELEM(trigProducts, 0, j);
				y2 = resolution*MAT_ELEM(trigProducts, 1, j);
				z2 = resolution*MAT_ELEM(trigProducts, 2, j);

				if (i != j)
				{
					distance = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
					if (distance < lastMinDistance)
						lastMinDistance = distance;
				}
			}
			if (k<70)
			{
				std::cout << k << " " << MAT_ELEM(resolutionMat, i, k) << " " << MAT_ELEM(trigProducts, 0, i) << " " << MAT_ELEM(trigProducts, 1, i) << " " << MAT_ELEM(trigProducts, 2, i) << ";"<< std::endl;
			}
			if (lastMinDistance>=threshold)
			{
				MAT_ELEM(resolutionMat, i, k) = -1;
			}
		}
//		std::vector<double> xx(xrows), yy(xrows), zz(xrows);
//		for (int i = 0; i<xrows; ++i)
//		{
//			resolution = MAT_ELEM(resolutionMat, i, k);
//
//			if (resolution>0)
//			{
//				if (k<50)
//				{
//					std::cout << k << " " << resolution << " " << MAT_ELEM(trigProducts, 0, i) << " " << MAT_ELEM(trigProducts, 1, i) << " " << MAT_ELEM(trigProducts, 2, i) << ";"<< std::endl;
//				}
//				xx.push_back(fabs(resolution*MAT_ELEM(trigProducts, 0, i)));
//				yy.push_back(fabs(resolution*MAT_ELEM(trigProducts, 1, i)));
//				zz.push_back(fabs(resolution*MAT_ELEM(trigProducts, 2, i)));
//
//			}
//		}
//		std::sort(xx.begin(),xx.end());
//		std::sort(yy.begin(),yy.end());
//		std::sort(zz.begin(),zz.end());
//
//		double xh = xx[(int) floor(0.95*((double) xrows))];
//		double yh = yy[(int) floor(0.95*((double) xrows))];
//		double zh = zz[(int) floor(0.95*((double) xrows))];
//
//		for (int i = 0; i<xrows; ++i)
//		{
//			resolution = MAT_ELEM(resolutionMat, i, k);
//
//			if (resolution>0)
//			{
//				x1 = fabs(resolution*MAT_ELEM(trigProducts, 0, i));
//				y1 = fabs(resolution*MAT_ELEM(trigProducts, 1, i));
//				z1 = fabs(resolution*MAT_ELEM(trigProducts, 2, i));
//				if ((x1<xh) || (y1<yh) || (z1<zh))
//					MAT_ELEM(resolutionMat, i, k) = -1;
//			}
//		}
	}
}

void ProgResDir::ellipsoidFitting(Matrix2D<double> &anglesMat,
									Matrix2D<double> &resolutionMat,
									Matrix2D<double> &axis)
{
	double x, y, z, a, b, c, resolution, rot, tilt;
	int xrows = angles.mdimx;
	std::vector<double> list_distances;

	//std::cout << "xrows = " << xrows << std::endl;

	//MAT_ELEM(resolutionMat, direccion, resolucion)

	Matrix2D<double> ellipMat;
	int dimMatrix = 0;
	size_t mycounter;
	Matrix2D<double> pseudoinv, quadricMatrix;
	Matrix1D<double> onesVector, leastSquares;
	Matrix2D<double> eigenvectors;
	Matrix1D<double> eigenvalues;
	//rows 1 2 3 (length axis a b c- where a is the smallest one)
	//rows 4 5 6 x y z coordinates of the first eigenvector
	//rows 7 8 9 x y z coordinates of the second eigenvector
	//rows 10 11 12 x y z coordinates of the third eigenvector
	axis.initZeros(12, NVoxelsOriginalMask);

	quadricMatrix.initZeros(3,3);
	for (int k = 0; k<NVoxelsOriginalMask; ++k)
	{
		dimMatrix = 0;
		for (int i = 0; i<xrows; ++i)
		{
			if (MAT_ELEM(resolutionMat, i, k) > 0)
				++dimMatrix;
		}

		ellipMat.initZeros(dimMatrix, 6);
		mycounter = 0; //It is required to store the matrix ellipMat
		for (int i = 0; i<xrows; ++i)
		{
			resolution = MAT_ELEM(resolutionMat, i, k);

			if (resolution>0)
			{
				x = resolution*MAT_ELEM(trigProducts, 0, i);
				y = resolution*MAT_ELEM(trigProducts, 1, i);
				z = resolution*MAT_ELEM(trigProducts, 2, i);

				MAT_ELEM(ellipMat, mycounter, 0) = x*x;
				MAT_ELEM(ellipMat, mycounter, 1) = y*y;
				MAT_ELEM(ellipMat, mycounter, 2) = z*z;
				MAT_ELEM(ellipMat, mycounter, 3) = 2*x*y;
				MAT_ELEM(ellipMat, mycounter, 4) = 2*x*z;
				MAT_ELEM(ellipMat, mycounter, 5) = 2*y*z;
				++mycounter;
			}
		}

		ellipMat.inv(pseudoinv);

		onesVector.initConstant(mycounter, 1.0);
		leastSquares = pseudoinv*onesVector;

		MAT_ELEM(quadricMatrix, 0, 0) = VEC_ELEM(leastSquares, 0);
		MAT_ELEM(quadricMatrix, 0, 1) = VEC_ELEM(leastSquares, 3);
		MAT_ELEM(quadricMatrix, 0, 2) = VEC_ELEM(leastSquares, 4);
		MAT_ELEM(quadricMatrix, 1, 0) = VEC_ELEM(leastSquares, 3);
		MAT_ELEM(quadricMatrix, 1, 1) = VEC_ELEM(leastSquares, 1);
		MAT_ELEM(quadricMatrix, 1, 2) = VEC_ELEM(leastSquares, 5);
		MAT_ELEM(quadricMatrix, 2, 0) = VEC_ELEM(leastSquares, 4);
		MAT_ELEM(quadricMatrix, 2, 1) = VEC_ELEM(leastSquares, 5);
		MAT_ELEM(quadricMatrix, 2, 2) = VEC_ELEM(leastSquares, 2);

		diagSymMatrix3x3(quadricMatrix, eigenvalues, eigenvectors);

		if (VEC_ELEM(eigenvalues, 0)<0){//This is de the vectorial product
			VEC_ELEM(eigenvalues, 0) = VEC_ELEM(eigenvalues, 1);
			MAT_ELEM(axis,3, k) = MAT_ELEM(eigenvectors,0,1)*MAT_ELEM(eigenvectors,2,2)-
									MAT_ELEM(eigenvectors,2,1)*MAT_ELEM(eigenvectors,1,2);
			MAT_ELEM(axis,4, k) = MAT_ELEM(eigenvectors,2,1)*MAT_ELEM(eigenvectors,0,2)-
									MAT_ELEM(eigenvectors,0,1)*MAT_ELEM(eigenvectors,2,2);
			MAT_ELEM(axis,5, k) = MAT_ELEM(eigenvectors,0,1)*MAT_ELEM(eigenvectors,1,2)-
									MAT_ELEM(eigenvectors,1,1)*MAT_ELEM(eigenvectors,0,2);
		}

		a = 1/sqrt(VEC_ELEM(eigenvalues, 0));
		b = 1/sqrt(VEC_ELEM(eigenvalues, 1));
		c = 1/sqrt(VEC_ELEM(eigenvalues, 2));

//		std::cout << "a = " << a << std::endl;
//		std::cout << "b = " << b << std::endl;
//		std::cout << "c = " << c << std::endl;

		MAT_ELEM(axis,0, k) = a;
		MAT_ELEM(axis,1, k) = b;
		MAT_ELEM(axis,2, k) = c;

		MAT_ELEM(axis,3, k) = MAT_ELEM(eigenvectors,0,0);
		MAT_ELEM(axis,4, k) = MAT_ELEM(eigenvectors,1,0);
		MAT_ELEM(axis,5, k) = MAT_ELEM(eigenvectors,2,0);

		MAT_ELEM(axis,6, k) = MAT_ELEM(eigenvectors,0,1);
		MAT_ELEM(axis,7, k) = MAT_ELEM(eigenvectors,1,1);
		MAT_ELEM(axis,8, k) = MAT_ELEM(eigenvectors,2,1);

		MAT_ELEM(axis,9, k) = MAT_ELEM(eigenvectors,0,2);
		MAT_ELEM(axis,10, k) = MAT_ELEM(eigenvectors,1,2);
		MAT_ELEM(axis,11, k) = MAT_ELEM(eigenvectors,2,2);
	}
}

void ProgResDir::radialAverageInMask(MultidimArray<int> &mask, MultidimArray<double> &inputVol, MetaData &md)
{
		double u_inf, u_sup, u;

		MultidimArray<int> &pMask = mask;
		int step = 1;

		double N;
		MultidimArray<double> radialAvg(XSIZE(inputVol)*0.5);
		MultidimArray<double> test_ring;
		test_ring.initZeros(inputVol);
		//DIRECT_MULTIDIM_ELEM(radialAvg,0) = sqrt(real(conj(A3D_ELEM(fftV, 0,0,0))*A3D_ELEM(fftV, 0,0,0)));
		std::cout << "XSIZE = " << XSIZE(inputVol) << std::endl;

		int uk, uj, ui;

		int siz = XSIZE(inputVol);
		size_t objId;

		inputVol.setXmippOrigin();
		pMask.setXmippOrigin();

		for(size_t kk=1; kk<siz*0.5; ++kk)
		{
			double cum_mean = 0;
			N = 0;
			u_sup = kk + step;
			u_inf = kk - step;

			FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
			{
				 if (A3D_ELEM(pMask, k, i, j)>0)
				{
				  //std::cout << "entro " << std::endl;
					u = sqrt(k*k + i*i + j*j);
					if ((u<u_sup) && (u>=u_inf))
					{
						cum_mean += A3D_ELEM(inputVol, k, i, j);
						N = N + 1;
					}

				 }

			}

			objId = md.addObject();
			if (cum_mean==0)
			{
				md.setValue(MDL_IDX, kk, objId);
				md.setValue(MDL_AVG, cum_mean, objId);
			}
			else
			{
				md.setValue(MDL_IDX, kk, objId);
				md.setValue(MDL_AVG, (cum_mean/N), objId);
			}
		}
}

void ProgResDir::radialAzimuthalResolution(Matrix2D<double> &resolutionMat,
		MultidimArray<int> &pmask,
		MultidimArray<double> &radial,
		MultidimArray<double> &azimuthal)
{

	radial.initZeros(pmask);
	azimuthal.initZeros(pmask);
	double radial_angle = 15*PI/180;
	double azimuthal_resolution = 0;
	double radial_resolution = 0;
	double azimuthal_angle = 75*PI/180;
	double resolution, dotproduct, x, y, z, iu, arcos;
	int xrows = angles.mdimx;
	int idx;
	idx = 0;
	double count_radial, count_azimuthal;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(pmask)
	{
		//i defines the direction and k the voxel
		if (A3D_ELEM(pmask,k,i,j) > 0 )
		{
			iu = 1/sqrt(i*i + j*j + k*k);
			count_radial = 0;
			count_azimuthal = 0;
			for (int ii = 0; ii<xrows; ++ii)
			{
				resolution = MAT_ELEM(resolutionMat, ii, idx);

				if (resolution>0)
				{

					x = MAT_ELEM(trigProducts, 0, ii);
					y = MAT_ELEM(trigProducts, 1, ii);
					z = MAT_ELEM(trigProducts, 2, ii);

					dotproduct = (x*i + y*j + z*k)*iu;
					arcos = acos(fabs(dotproduct));
					if (arcos<=azimuthal_angle)
					{
						count_azimuthal = count_azimuthal + 1;
						azimuthal_resolution += resolution;
					}
					if (arcos>=radial_angle)
					{
						count_radial = count_radial + 1;
						radial_resolution += resolution;
					}

				}
			}
			std::cout << "count_radial = " << count_radial << std::endl;
			std::cout << "count_azimuthal = " << count_azimuthal << std::endl;
			std::cout << "  " << std::endl;
			++idx;
		}
		A3D_ELEM(radial,k,i,j) = radial_resolution/count_radial;
		A3D_ELEM(azimuthal,k,i,j) = azimuthal_resolution/count_azimuthal;
		azimuthal_resolution = 0;
		radial_resolution = 0;
	}


}

double ProgResDir::firstMonoResEstimation(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, MultidimArray<double> &amplitude)
{
	fftVRiesz.initZeros(myfftV);
	amplitude.resizeNoCopy(VRiesz);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double iw=1.0/freq;
	double iwl=1.0/freqH;
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
		} else if (un>freq)
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	double uz, uy, ux;
	n=0;

	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				ux = VEC_ELEM(freq_fourier,j);
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				double un=1.0/iun;
				if (freqH<=un && un<=freq)
				{
					//double H=0.5*(1+cos((un-w1)*ideltal));
					//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-J*ux*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
					//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *=0.5*(1+cos((un-w1)*ideltal));//H;
					//Next lines are an optimization of the commented ones
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= -ux*iun*0.5*(1+cos((un-freq)*ideltal));//H;
				} else if (un>freq)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-ux*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
				}
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;

	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier,i);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				double un=1.0/iun;
				if (freqH<=un && un<=freq)
				{
					//double H=0.5*(1+cos((un-w1)*ideltal));
					//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-J*uy*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
					//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *=0.5*(1+cos((un-w1)*ideltal));//H;
					//Next lines are an optimization of the commented ones
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= -uy*iun*0.5*(1+cos((un-freq)*ideltal));//H;
				} else if (un>freq)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-uy*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
				}
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate third component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		uz = VEC_ELEM(freq_fourier,k);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				double un=1.0/iun;
				if (freqH<=un && un<=freq)
				{
					//double H=0.5*(1+cos((un-w1)*ideltal));
					//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-J*uz*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
					//DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-w1)*ideltal));//H;
					//Next lines are an optimization of the commented ones
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= -uz*iun*0.5*(1+cos((un-freq)*ideltal));//H;
				} else if (un>freq)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (-uz*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= J;
				}
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
	}

	// Low pass filter the monogenic amplitude
	lowPassFilter.w1 = freq;
	amplitude.setXmippOrigin();
	lowPassFilter.applyMaskSpace(amplitude);

	double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
	MultidimArray<int> &pMask = mask();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitude, n);
		if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
		{
			sumN  += amplitudeValue;
			sumN2 += amplitudeValue*amplitudeValue;
			++NN;
		}
	}

	double mean_noise = sumN/NN;
	return mean_noise;

}

void ProgResDir::run()
{
	produceSideInfo();

	bool continueIter = false, breakIter = false;
	double criticalZ=icdf_gauss(significance);

	double range = maxRes-minRes;
	double step = range/N_freq;

	if (step<0.3)
		step=0.3;

	step = 0.3;

	std::cout << "Analyzing directions " << std::endl;
	std::cout << "maxRes = " << maxRes << std::endl;
	std::cout << "minRes = " << minRes << std::endl;
	std::cout << "N_freq = " << N_freq << std::endl;
	std::cout << "step = " << step << std::endl;
	std::cout << "criticalZ = " << criticalZ << std::endl;

	Image<double> outputResolution;
	MultidimArray<double> amplitudeMS, amplitudeMN;

	double w, wH;
	int volsize = ZSIZE(VRiesz);
//	FFT_IDX2DIGFREQ(10, volsize, w);
//	FFT_IDX2DIGFREQ(11, volsize, wH); //Frequency chosen for a first estimation

	//Checking with MonoRes at 50A;
	int aux_idx;
	double aux_freq;
	aux_freq = sampling/50;
	if (maxRes>30)
	{
		DIGFREQ2FFT_IDX(sampling/30, volsize, aux_idx);
		FFT_IDX2DIGFREQ(aux_idx, volsize, w);
		FFT_IDX2DIGFREQ(aux_idx+1, volsize, wH); //Frequency chosen for a first estimation
	}
	else
	{
		FFT_IDX2DIGFREQ(3, volsize, w);
		FFT_IDX2DIGFREQ(4, volsize, w);
		aux_idx = 3;
	}
	std::cout << "Calling MonoRes core as a first estimation at " << sampling/w << "A." << std::endl;


	double AvgNoise;
	AvgNoise = firstMonoResEstimation(fftV, w, wH, amplitudeMS)/9.0;

	N_directions=angles.mdimx;

	std::cout << "N_directions = " << N_directions << std::endl;

	double cone_angle = 20.0; //(degrees)
	cone_angle = PI*cone_angle/180;

	trigProducts.initZeros(3, N_directions);

	for (size_t dir=0; dir<N_directions; dir++)
	{
		outputResolution().initZeros(VRiesz);
//		MultidimArray<double> &pOutputResolution = outputResolution();
		double freq, freqL, freqH, counter, resolution_2;
		MultidimArray<int> mask_aux = mask();
		MultidimArray<int> &pMask = mask_aux;
		std::vector<double> list;
		double resolution;  //A huge value for achieving last_resolution < resolution

		double max_meanS = -1e38;
		double cut_value = 0.025;

		bool doNextIteration=true;

		int fourier_idx, last_fourier_idx = -1, iter = 0, fourier_idx_2;
		fourier_idx = aux_idx;
		int count_res = 0;
		double rot = MAT_ELEM(angles, 0, dir);
		double tilt = MAT_ELEM(angles, 1, dir);
		MAT_ELEM(trigProducts, 0, dir) = sin(tilt*PI/180)*cos(rot*PI/180);
		MAT_ELEM(trigProducts, 1, dir) = sin(tilt*PI/180)*sin(rot*PI/180);
		MAT_ELEM(trigProducts, 2, dir) = cos(tilt*PI/180);
		std::cout << "--------------NEW DIRECTION--------------" << std::endl;
		std::cout << "direction = " << dir+1 << "   rot = " << rot << "   tilt = " << tilt << std::endl;


		std::vector<double> noiseValues;
		FileName fnDebug;
		double last_resolution = 0;

		defineCone(fftV, conefilter, rot, tilt);
		maskMatrix.initConstant(1, NVoxelsOriginalMask, 1);
		do
		{
			continueIter = false;
			breakIter = false;
			//std::cout << "--------------Frequency--------------" << std::endl;

			resolution2eval_(fourier_idx, step,
							resolution, last_resolution, last_fourier_idx,
							freq, freqL, freqH,
							continueIter, breakIter, doNextIteration);

			if (breakIter)
				break;

			if (continueIter)
				continue;

			std::cout << "resolution = " << resolution << "  resolutionL = " << sampling/freqL << "  resolutionH = " << sampling/freqH << " iter = " << iter << std::endl;
//			std::cout << "resolution = " << freq 	   << "  resolutionL = " << freqL 		   << "  resolutionH = " << freqH << std::endl;


			list.push_back(resolution);

			if (iter<2)
				resolution_2 = list[0];
			else
				resolution_2 = list[iter - 2];

			fnDebug = "Signal";

			amplitudeMonogenicSignal3D_fast(conefilter, freq, freqH, freqL, amplitudeMS, iter, dir, fnDebug, rot, tilt);

			double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
			noiseValues.clear();

			double amplitudeValue;

			double x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
			double y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
			double z_dir = cos(tilt*PI/180);

			double uz, uy, ux;

			//TODO: check if can be taken out side the loop
//			MultidimArray<double> coneVol;
//			coneVol.initZeros(amplitudeMS);
			int n=0;
			int z_size = ZSIZE(amplitudeMS);
			int x_size = XSIZE(amplitudeMS);
			int y_size = YSIZE(amplitudeMS);

			for(int k=0; k<z_size; ++k)
			{
//				std::cout << " k = " << k  <<std::endl;
				for(int i=0; i<y_size; ++i)
				{
					for(int j=0; j<x_size; ++j)
					{
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						{
							amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
							sumS  += amplitudeValue;
							++NS;
						}
						else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
						{
							uz = (k - z_size*0.5);
//							std::cout << " uz = " << uz  <<std::endl;
							ux = (j - x_size*0.5);
							uy = (i - y_size*0.5);

							double rad = sqrt(ux*ux + uy*uy + uz*uz);
							double iun = 1/rad;

							//BE CAREFULL with the order
							double dotproduct = (uy*y_dir + ux*x_dir + uz*z_dir)*iun;

							double acosine = acos(dotproduct);

							//TODO: change efficienty the if condition by means of fabs(cos(angle))
							if (((acosine<(cone_angle)) || (acosine>(PI-cone_angle)) )
									&& (rad>Rparticle))
							{

//								DIRECT_MULTIDIM_ELEM(coneVol, n) = 1;
								amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
								sumN  += amplitudeValue;
								sumN2 += amplitudeValue*amplitudeValue;
								++NN;
							}
						}
						++n;
					}
				}
			}

//			#ifdef DEBUG_DIR
//				if (iter == 0)
//				{
//				Image<double> img;
//
//				FileName iternumber;
//				iternumber = formatString("cone_noise_%i_%i.vol", dir, iter);
//				img = coneVol;
//				img.write(iternumber);
//				}
//			#endif

//				std::cout << "NS = " << NS << std::endl;
			if ( (NS/(double) NVoxelsOriginalMask)<cut_value ) //when the 2.5% is reached then the iterative process stops
			{
				std::cout << "Search of resolutions stopped due to mask has been completed" << std::endl;
				doNextIteration =false;
				Nvoxels = 0;
//				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
//				{
//				  if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n) > 0)
//					DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
//				}
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
	//			double sigma2S=sumS2/NS-meanS*meanS;
				double meanN=sumN/NN;
				double sigma2N=sumN2/NN-meanN*meanN;

				if (meanS>max_meanS)
					max_meanS = meanS;

				if (meanS<0.001*AvgNoise)//0001*max_meanS)
				{
					//std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS	= " << NS << std::endl;
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

					size_t maskPos = 0;
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
					{
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						{
							if (MAT_ELEM(maskMatrix, 0, maskPos) >=1)
							{
								if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
								{
	//								DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution;//sampling/freq;
									MAT_ELEM(resolutionMatrix, dir, maskPos) = resolution;
									MAT_ELEM(maskMatrix, 0, maskPos) = 1;
								}
								else
								{
									MAT_ELEM(maskMatrix, 0, maskPos) += 1;
									if (MAT_ELEM(maskMatrix, 0, maskPos) >2)
									{
										MAT_ELEM(maskMatrix, 0, maskPos) = 0;
										MAT_ELEM(resolutionMatrix, dir, maskPos) = resolution_2;
	//									DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution_2; //resolution + counter*step;
									}
								}
							}
							++maskPos;
						}
					}

//					#ifdef DEBUG_MASK
//					FileName fnmask_debug;
//					fnmask_debug = formatString("maske_%i.vol", iter);
//					mask.write(fnmask_debug);
//					#endif

					//#ifdef DEBUG
//						std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
//						std::cout << "  meanS= " << meanS << " NS= " << NS << std::endl;
//						std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
					//#endif

					if (doNextIteration)
						if (resolution <= (minRes-0.001))
							doNextIteration = false;
					}
			}
			++iter;
			last_resolution = resolution;
		}while(doNextIteration);



//		amplitudeMN.clear();
//		amplitudeMS.clear();
//		fftVRiesz.clear();
		
		int NVoxelsOriginalMask_bis = 0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mask())
		{
			if (DIRECT_MULTIDIM_ELEM(mask(), n) == 1)
				++NVoxelsOriginalMask_bis;
		}
		//////////////////
		//INERTIA MOMENT//
		//////////////////
//		inertiaMatrixNew(resolutionMatrix, inertiaMatrixVariable, rot, tilt, dir);
//		#ifdef DEBUG_DIR

		size_t maskPos=0;
		Image<double> ResolutionVol;
		MultidimArray<double> &pResolutionVol = ResolutionVol();

		pResolutionVol.initZeros(amplitudeMS);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pResolutionVol)
		{
			if (DIRECT_MULTIDIM_ELEM(mask(), n) == 1)
			{
				double myres = MAT_ELEM(resolutionMatrix, dir, maskPos);
				DIRECT_MULTIDIM_ELEM(pResolutionVol, n) = myres;
				++maskPos;
			}
		}
		//#endif
//		#ifdef DEBUG_DIR
		Image<double> saveImg;
		saveImg = pResolutionVol;
		FileName fnres = formatString("resolution_dir_%i.vol", dir+1);
		saveImg.write(fnres);
		saveImg.clear();
//		#endif
		pResolutionVol.clear();
		list.clear();

		std::cout << "----------------direction-finished----------------" << std::endl;
	}

	//Remove outliers
	removeOutliers(trigProducts, resolutionMatrix);
//	removeOutliers(angles, resolutionMatrix);

	//Ellipsoid fitting
	Matrix2D<double> axis;
	ellipsoidFitting(trigProducts, resolutionMatrix, axis);
//	ellipsoidFitting(angles, resolutionMatrix, axis);

	Image<double> doaVol;
	MultidimArray<double> &pdoaVol = doaVol();

	pdoaVol.initZeros(NSIZE(mask()),ZSIZE(mask()), YSIZE(mask()), XSIZE(mask()));



	int idx = 0;
//	std::cout << "antes del for = " << MAT_ELEM(axis, 0, 0) << std::endl;
	std::cout << "NVoxelsOriginalMask = " << NVoxelsOriginalMask << std::endl;


	double niquist;

	niquist = 2*sampling;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdoaVol)
	{
		if (DIRECT_MULTIDIM_ELEM(mask(), n) >0 ) //before ==1
		{

			double a = MAT_ELEM(axis, 0, idx);
			double c = MAT_ELEM(axis, 2, idx);
//			if (idx<100)
//				std::cout << c << " " << a << ";" << std::endl;
//			std::cout << "a = " << a << std::endl;
//			std::cout << "c = " << c << std::endl;
			DIRECT_MULTIDIM_ELEM(pdoaVol, n) = (c)/(a);
			++idx;
		}
	}



	Image<double> imgdoa;
	imgdoa = pdoaVol;
	imgdoa.write(fnDoA);

	MultidimArray<double> radial, azimuthal;
	radialAzimuthalResolution(resolutionMatrix, mask(), radial, azimuthal);

	imgdoa = radial;
	imgdoa.write(fnradial);
	imgdoa = azimuthal;
	imgdoa.write(fnazimuthal);

	std::cout << "Calculating the radial and azimuthal resolution " << std::endl;


	MetaData mdRadial, mdAzimuthal;

//	Image<double> V;
//	V.read(fnVol);
//	MultidimArray<double> &inputVol = V();
//
//	radialAverageInMask(mask(), inputVol, mdAzimuthal);
//	mdAzimuthal.write(fnMDazimuthal);

	radialAverageInMask(mask(), azimuthal, mdAzimuthal);
	radialAverageInMask(mask(), radial, mdRadial);

	mdAzimuthal.write(fnMDazimuthal);
	mdRadial.write(fnMDradial);


///////////////////////

	double lambda_1, lambda_2, lambda_3, doa;
	double direction_x, direction_y, direction_z;
	int counter = 0;
	Matrix2D<double> eigenvectors;
	Matrix1D<double> eigenvalues, r0_1(3), rF_1(3), r0_2(3), rF_2(3), r0_3(3), rF_3(3), r(3);
	MultidimArray<int> arrows;
	arrows.initZeros(mask());
	const int gridStep=10;
	size_t n=0;
	int maskPos=0;
///////////////////////

	idx = 0;
	int siz;
	siz = XSIZE(arrows);
	double xcoor, ycoor, zcoor, rad, rot, tilt;
	MetaData md;
	size_t objId;
	FileName fn_md;

	FOR_ALL_ELEMENTS_IN_ARRAY3D(arrows)
	{
			if (A3D_ELEM(mask(),k,i,j) > 0 ) //before ==1
			{

				//lambda_3 is assumed as the least eigenvalue
				if ( (i%gridStep==0) && (j%gridStep==0) && (k%gridStep==0) )
				{
					double lambda_1 = MAT_ELEM(axis, 0, idx);
					double lambda_3 = MAT_ELEM(axis, 2, idx);

					xcoor = MAT_ELEM(axis, 3, idx);
					ycoor = MAT_ELEM(axis, 4, idx);
					zcoor = MAT_ELEM(axis, 5, idx);

					rot = atan2(ycoor, xcoor)*180/PI;
					tilt = acos(zcoor)*180/PI;


//					rotation3DMatrix(double ang, const Matrix1D<double> &axis,
//					                      Matrix2D<double> &result, bool homogeneous)

					double sc;
					sc = lambda_1/8.0;
//					std::cout << "a = " << lambda_3 << "  c= " << lambda_1 << std::endl;
//					std::cout << "sc = " << sc << "  c/sc= " << lambda_1/sc << std::endl;

					//write md with values!
					objId = md.addObject();
					md.setValue(MDL_ANGLE_ROT, rot, objId);
					md.setValue(MDL_ANGLE_TILT, tilt, objId);
					md.setValue(MDL_XCOOR, (int) j, objId);
					md.setValue(MDL_YCOOR, (int) i, objId);
					md.setValue(MDL_ZCOOR, (int) k, objId);
					md.setValue(MDL_MAX, 7.0, objId);
					md.setValue(MDL_MIN, lambda_3/sc, objId);
					md.setValue(MDL_INTSCALE, lambda_3/lambda_1, objId);
				}
				++idx;
			}
			++n;
	}

	md.write(fnDirections);
}

