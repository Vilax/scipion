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
//define DEBUG_FILTER
#define MONO_AMPLITUDE
//define DEBUG_SYMMETRY

void ProgResDir::readParams()
{
	fnVol = getParam("--vol");
	fnVol2 = getParam("--vol2");
	fnMeanVol = getParam("--meanVol");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	fnMaskOut = getParam("--mask_out");
	fnchim = getParam("--chimera_volume");
	sampling = getDoubleParam("--sampling_rate");
	ang_sampling = getDoubleParam("--angular_sampling");
	R = getDoubleParam("--volumeRadius");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	fnVar = getParam("--varVol");
	fnSym = getParam("--sym");
	N_freq = getDoubleParam("--number_frequencies");
	noiseOnlyInHalves = checkParam("--noiseonlyinhalves");
	significance = getDoubleParam("--significance");
	fnMd = getParam("--md_resdir");
	fnDoA = getParam("--doa_vol");
}


void ProgResDir::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">        : Input volume");
	addParamsLine("  [--mask <vol_file=\"\">]     : Mask defining the macromolecule");
	addParamsLine("                               :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                               :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [--mask_out <vol_file=\"\">] : sometimes the provided mask is not perfect, and contains voxels out of the particle");
	addParamsLine("                               :+ Thus the algorithm calculated a tight mask to the volume");
	addParamsLine("  [--vol2 <vol_file=\"\">]     : Half volume 2");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--meanVol <vol_file=\"\">]  : Mean volume of half1 and half2 (only it is neccesary the two haves are used)");
	addParamsLine("  --sym <symmetry>: Symmetry (c1, c2, c3,..d1, d2, d3,...)");
	addParamsLine("  [--chimera_volume <output=\"Chimera_resolution_volume.vol\">]: Local resolution volume for chimera viewer (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]      : Sampling rate (A/px)");
	addParamsLine("  [--angular_sampling <s=15>]  : Angular Sampling rate (degrees)");
	addParamsLine("  [--volumeRadius <s=100>]     : This parameter determines the radius of a sphere where the volume is");
	addParamsLine("  [--number_frequencies <w=50>]       : The resolution is computed at a number of frequencies between mininum and");
	addParamsLine("                               : maximum resolution px/A. This parameter determines that number");
	addParamsLine("  [--minRes <s=30>]            : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]             : Maximum resolution (A)");
	addParamsLine("  [--maxVol <vol_file=\"\">]   : Output filename with maximum resolution volume");
	addParamsLine("  [--minVol <vol_file=\"\">]   : Output filename with minimum resolution volume");
	addParamsLine("  [--varVol <vol_file=\"\">]   : Output filename with varianze resolution volume");
	addParamsLine("  [--noiseonlyinhalves]        : The noise estimation is only performed inside the mask");
	addParamsLine("  [--significance <s=0.95>]    : The level of confidence for the hypothesis test.");
	addParamsLine("  [--md_resdir <file=\".\">]   : Metadata with mean resolution by direction.");
	addParamsLine("  [--doa_vol <vol_file=\"\">]  : Output filename with DoA volume");
	addParamsLine("  [--sphericity <vol_file=\"\">]  : Output filename with DoA volume");
}

void ProgResDir::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;

	Image<double> V;
	if ((fnVol !="") && (fnVol2 !=""))
	{
		Image<double> V1, V2;
		V1.read(fnVol);
		V2.read(fnVol2);
		V()=0.5*(V1()+V2());
		V.write(fnMeanVol);
		halfMapsGiven = true;
	}
	else{
	    V.read(fnVol);
	    halfMapsGiven = false;
	}
	V().setXmippOrigin();

	//Sweeping the projection sphere
	std::cout << "Obtaining angular projections..." << std::endl;
	generateGridProjectionMatching(fnVol, ang_sampling, angles);

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);

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



//	#ifdef DEBUG_DIR
//	Image<double> freqimg;
//	freqimg = iu;
//	freqimg.write(formatString("frequencys.vol"));
//	#endif



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

	if (halfMapsGiven)
	{
		Image<double> V1, V2;
		V1.read(fnVol);
		V2.read(fnVol2);

		V1()-=V2();
		V1()/=2;
		fftN=new MultidimArray< std::complex<double> >;
		FourierTransformer transformer2;
		#ifdef DEBUG
		  V1.write("diff_volume.vol");
		#endif
		transformer2.FourierTransform(V1(), *fftN);
	}
	else
	{
		fftN=&fftV;
	}


	freq_fourier.initZeros(ZSIZE(inputVol));
	int size = ZSIZE(inputVol);
	V.clear();


	double u;
	int size_fourier(ZSIZE(fftV));

	for(size_t k=0; k<size_fourier; ++k)
	{
		FFT_IDX2DIGFREQ(k,size, u);
		VEC_ELEM(freq_fourier,k) = u;
//		std::cout << "freq_fourier = " << u  << std::endl;
	}
}


void ProgResDir::generateGridProjectionMatching(FileName fnVol_, double smprt,
												Matrix2D<double> &angles)
{
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
//			std::cout << "k = "<< count << "  rot = "  << rot <<
//					  "  tilt = " << tilt << std::endl;
			count++;
		}
	}

	N_directions = count;
	angles.initZeros(4,N_directions);

	for (size_t k = 0; k<N_directions; k++)
	{
		MAT_ELEM(angles, 0, k) = MAT_ELEM(aux_angles,0, k);
		MAT_ELEM(angles, 1, k) = MAT_ELEM(aux_angles,1, k);
		std::cout << "k=" << k << "  rot = "  << MAT_ELEM(angles, 0, k) <<
								  "  tilt = " << MAT_ELEM(angles, 1, k) << std::endl;
	}
}


void ProgResDir::amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
		double w1, double w1l, double wH, MultidimArray<double> &amplitude, int count, int dir, FileName fnDebug,
		double angle_cone, double rot, double tilt)
{
	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);
	amplitude.resizeNoCopy(VRiesz);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(w1-w1l);

	#ifdef DEBUG_DIR
	MultidimArray<double> coneVol;
	coneVol.initZeros(iu);
	#endif


	double x_dir, y_dir, z_dir;

	x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir = cos(tilt*PI/180);

	double uz, uy, ux, uxxuyy;
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		uz = VEC_ELEM(freq_fourier,k);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier,i);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);

				ux = VEC_ELEM(freq_fourier,j);
				double dotproduct;

				//BE CAREFULL with the order
				dotproduct = (uy*x_dir + ux*y_dir + uz*z_dir)*iun;

				double acosine = acos(fabs(dotproduct));

				double arg_exp = acosine*acosine*acosine*acosine/(0.12*0.12*0.12*0.12);
				double un=1.0/iun;

//				DIRECT_MULTIDIM_ELEM(coneVol, n) = exp(-arg_exp);

				if (w1l<=un && un<=w1)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= exp(-arg_exp);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-w1)*ideltal));//H;
				} else if (un>w1)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= exp(-arg_exp);
				}

				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J*iun*DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
				++n;
			}
		}
	}

	#ifdef DEBUG_DIR
	if ( (count == 0) )
	{
		Image<double> direction;
		direction = coneVol;
		direction.write(formatString("cone_%i.vol", dir));
	}
	#endif

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	//#ifdef DEBUG_DIR

//		Image<double> filteredvolume;
//		filteredvolume = VRiesz;
//		filteredvolume.write(formatString("Volumen_filtrado_%i_%i.vol", dir,count));

	//#endif

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
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
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
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
			++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+= DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	//transformer_inv.inverseFourierTransform(fftVRiesz_aux, VRiesz);



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
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= uz*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
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

	Image<double> saveImg2;
//	saveImg2 = amplitude;
	FileName fnSaveImg2;
	FileName iternumber;

//	#ifdef MONO_AMPLITUDE

//	if (fnDebug.c_str() != "")
//	{
//		saveImg2 = amplitude;
//		iternumber = formatString("No_Filtered_Amplitude_%i_%i.vol", dir, count);
//		saveImg2.write(fnDebug+iternumber);
//	}
//	#endif


	size_t Ydim, Xdim, Zdim, Ndim;
	amplitude.getDimensions(Xdim, Ydim, Zdim, Ndim);


	n=0;
//	//Creating mask for smoothing
//		for (size_t i = 0; i< Xdim; i++)
//		{
//			for (size_t j = 0; j< Ydim; j++)
//			{
//				for (size_t k = 0; k< Zdim; k++)
//				{
//					if (i<= N_smoothing)
//						DIRECT_MULTIDIM_ELEM(amplitude, n) = DIRECT_MULTIDIM_ELEM(amplitude, n)*0.5*(1+cos(PI*(N_smoothing-i)/(N_smoothing)));
//					if (j<= N_smoothing)
//						DIRECT_MULTIDIM_ELEM(amplitude, n) = DIRECT_MULTIDIM_ELEM(amplitude, n)*0.5*(1+cos(PI*(N_smoothing-j)/(N_smoothing)));
//					if (k<= N_smoothing)
//						DIRECT_MULTIDIM_ELEM(amplitude, n) = DIRECT_MULTIDIM_ELEM(amplitude, n)*0.5*(1+cos(PI*(N_smoothing-k)/(N_smoothing)));
//					if (i>= (Xdim - N_smoothing))
//						DIRECT_MULTIDIM_ELEM(amplitude, n) = DIRECT_MULTIDIM_ELEM(amplitude, n)*0.5*(1+cos(PI*(N_smoothing+Xdim-i)/(N_smoothing)));
//					if (j>= (Ydim - N_smoothing))
//						DIRECT_MULTIDIM_ELEM(amplitude, n) = DIRECT_MULTIDIM_ELEM(amplitude, n)*0.5*(1+cos(PI*(N_smoothing+Ydim-j)/(N_smoothing)));
//					if (k>= (Zdim - N_smoothing))
//						DIRECT_MULTIDIM_ELEM(amplitude, n) = DIRECT_MULTIDIM_ELEM(amplitude, n)*0.5*(1+cos(PI*(N_smoothing+Zdim-k)/(N_smoothing)));
//					++n;
//				}
//			}
//		}

		amplitude.setXmippOrigin();
		int z_size = ZSIZE(amplitude);
		int x_size = XSIZE(amplitude);
		int y_size = YSIZE(amplitude);

		std::cout << "z_size = " << z_size << std::endl;
		double limit_radius = (z_size/2-N_smoothing);
		for(int k=0; k<z_size; ++k)
		{
			uz = (k - z_size*0.5);
			for(int i=0; i<y_size; ++i)
			{
				uy = (i - y_size*0.5);
				for(int j=0; j<x_size; ++j)
				{
					ux = (j - x_size*0.5);
					double radius = sqrt(ux*ux + uy*uy + uz*uz);
					if ((radius>=limit_radius) && (radius<=(z_size*0.5)))
						DIRECT_MULTIDIM_ELEM(amplitude, n) *= 0.5*(1+cos(PI*(limit_radius-radius)/(N_smoothing)));
					else if (radius>(0.5*z_size))
						DIRECT_MULTIDIM_ELEM(amplitude, n) = 0;
					++n;
				}
			}
		}


//		#ifdef MONO_AMPLITUDE
//		saveImg2 = amplitude;
//		if (fnDebug.c_str() != "")
//		{
//			iternumber = formatString("smoothed_volume_%i_%i.vol", dir, count);
//			saveImg2.write(fnDebug+iternumber);
//		}
////		saveImg2.clear();
//		#endif


	amplitude.setXmippOrigin();

	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);

    double raised_w = PI/(wH-w1);
//    w1 = 0.4;

//	for (size_t k=0; k<ZSIZE(fftVRiesz); k++)
//	{
//		for (size_t i=0; i<YSIZE(fftVRiesz); i++)
//		{
//			for (size_t j=0; j<XSIZE(fftVRiesz); j++)
//			{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
	{
			double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//			double un=1.0/DIRECT_A3D_ELEM(iu,k,i,j);

			if ((wH)>=un && un>=w1)
				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-w1)));
			else
				if (un>(wH))
					DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
//			n++;
//			}
//		}
//	}
	}
	transformer_inv.inverseFourierTransform();


	// Low pass filter the monogenic amplitude
//	amplitude.setXmippOrigin();
//	lowPassFilter.w1 = w1;
//
//	lowPassFilter.applyMaskSpace(amplitude);

//	#ifdef MONO_AMPLITUDE

//	saveImg2 = amplitude;
//
//	if (fnDebug.c_str() != "")
//	{
//		iternumber = formatString("_Filtered_Amplitude_%i_%i.vol", dir, count);
//		saveImg2.write(fnDebug+iternumber);
//	}
//	saveImg2.clear();
////	#endif // DEBUG

}


void ProgResDir::postProcessingLocalResolutions(MultidimArray<double> &resolutionVol,
				std::vector<double> &list, MultidimArray<double> &resolutionChimera,
				double &cut_value, MultidimArray<int> &pMask)
{
	MultidimArray<double> resolutionVol_aux = resolutionVol;
	double last_resolution_2 = list[(list.size()-1)];

	// Count number of voxels with resolution
	size_t N=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>(last_resolution_2-0.001)) //the value 0.001 is a tolerance
			++N;

	// Get all resolution values
	MultidimArray<double> resolutions(N);
	size_t N_iter=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n)>(last_resolution_2-0.001))
			DIRECT_MULTIDIM_ELEM(resolutions,N_iter++)=DIRECT_MULTIDIM_ELEM(resolutionVol, n);

	// Sort value and get threshold
	std::sort(&A1D_ELEM(resolutions,0),&A1D_ELEM(resolutions,N));
	double filling_value = A1D_ELEM(resolutions, (int)(0.5*N)); //median value
	double trimming_value = A1D_ELEM(resolutions, (int)((1-cut_value)*N));
	double freq, res, init_res, last_res;



	init_res = list[0];
	last_res = list[(list.size()-1)];
	
	resolutionChimera = resolutionVol;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
	{
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n) < last_res)
		{
			if (DIRECT_MULTIDIM_ELEM(pMask, n) >=1)
			{
				DIRECT_MULTIDIM_ELEM(resolutionVol, n) = filling_value;
			}
			else
			{
				DIRECT_MULTIDIM_ELEM(resolutionVol, n) = 0;
				DIRECT_MULTIDIM_ELEM(pMask,n) = 0;
			}

		}
		if (DIRECT_MULTIDIM_ELEM(resolutionVol, n) > trimming_value)
		{
		  DIRECT_MULTIDIM_ELEM(pMask,n) = 0;
		  DIRECT_MULTIDIM_ELEM(resolutionVol, n) = filling_value;
		}
	}
	//#ifdef DEBUG_MASK
	Image<int> imgMask;
	imgMask = pMask;
	imgMask.write(fnMaskOut);
	//#endif
}


void ProgResDir::inertiaMatrix(MultidimArray<double> &resolutionVol,
							   MultidimArray<double> &Inertia_11,
							   MultidimArray<double> &Inertia_12,
							   MultidimArray<double> &Inertia_13,
							   MultidimArray<double> &Inertia_22,
							   MultidimArray<double> &Inertia_23,
							   MultidimArray<double> &Inertia_33,
							   MultidimArray<double> &SumRes,
							   double rot, double tilt)
{
	double x_dir, y_dir, z_dir, resVal, x_dir_sym, y_dir_sym, z_dir_sym, r_xyz2;
	x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir = cos(tilt*PI/180);
	x_dir_sym = -sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir_sym = -sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir_sym = -cos(tilt*PI/180);

	//std::cout << "x_dir = " << x_dir << "  y_dir = " << y_dir<< "  z_dir = " << z_dir << std::endl;

	int idx = 0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resolutionVol)
	{
		if (DIRECT_MULTIDIM_ELEM(mask(), n) == 1)
		{
			resVal = DIRECT_MULTIDIM_ELEM(resolutionVol,n);//*DIRECT_MULTIDIM_ELEM(resolutionVol,n);
			r_xyz2 = resVal*resVal;
//			r_xyz2 = 10;

			DIRECT_MULTIDIM_ELEM(Inertia_11,n) += r_xyz2*(1-x_dir*x_dir);
			DIRECT_MULTIDIM_ELEM(Inertia_12,n) -= r_xyz2*x_dir*y_dir;
			DIRECT_MULTIDIM_ELEM(Inertia_13,n) -= r_xyz2*x_dir*z_dir;
			DIRECT_MULTIDIM_ELEM(Inertia_22,n) += r_xyz2*(1-y_dir*y_dir);
			DIRECT_MULTIDIM_ELEM(Inertia_23,n) -= r_xyz2*y_dir*z_dir;
			DIRECT_MULTIDIM_ELEM(Inertia_33,n) += r_xyz2*(1-z_dir*z_dir);

			DIRECT_MULTIDIM_ELEM(Inertia_11,n) += r_xyz2*(1-x_dir_sym*x_dir_sym);
			DIRECT_MULTIDIM_ELEM(Inertia_12,n) -= r_xyz2*x_dir_sym*y_dir_sym;
			DIRECT_MULTIDIM_ELEM(Inertia_13,n) -= r_xyz2*x_dir_sym*z_dir_sym;
			DIRECT_MULTIDIM_ELEM(Inertia_22,n) += r_xyz2*(1-y_dir_sym*y_dir_sym);
			DIRECT_MULTIDIM_ELEM(Inertia_23,n) -= r_xyz2*y_dir_sym*z_dir_sym;
			DIRECT_MULTIDIM_ELEM(Inertia_33,n) += r_xyz2*(1-z_dir_sym*z_dir_sym);

			DIRECT_MULTIDIM_ELEM(SumRes,n) += 2;
			idx++;

			if (idx == 34)
			{
				MetaData md;
				md.read("res.xmd");
				size_t objId;
				objId = md.addObject();
				md.setValue(MDL_ANGLE_ROT, rot, objId);
				md.setValue(MDL_ANGLE_TILT, tilt, objId);
				md.setValue(MDL_RESOLUTION_FREQREAL, resVal, objId);


				md.write("res.xmd");

				std::cout << "resVal = " << resVal << std::endl;
			}

			//resVal = 1;

		}
	}
}

void ProgResDir::diagSymMatrix3x3(Matrix2D<double> A, int Ndirections,
					double &lambda_1, double &lambda_2, double &lambda_3)
{
	double b, c, d, p, q, Delta;

	b = MAT_ELEM(A, 0, 0) + MAT_ELEM(A, 1, 1) + MAT_ELEM(A, 2, 2);
	c = MAT_ELEM(A, 0, 0)*MAT_ELEM(A, 1, 1) +
		MAT_ELEM(A, 0, 0)*MAT_ELEM(A, 2, 2) +
		MAT_ELEM(A, 1, 1)*MAT_ELEM(A, 2, 2) -
		MAT_ELEM(A, 0, 1)*MAT_ELEM(A, 0, 1) -
		MAT_ELEM(A, 0, 2)*MAT_ELEM(A, 0, 2) -
		MAT_ELEM(A, 1, 2)*MAT_ELEM(A, 1, 2);

	d = MAT_ELEM(A, 0, 0)*MAT_ELEM(A, 1, 2)*MAT_ELEM(A, 1, 2) +
		MAT_ELEM(A, 1, 1)*MAT_ELEM(A, 0, 2)*MAT_ELEM(A, 0, 2) +
		MAT_ELEM(A, 2, 2)*MAT_ELEM(A, 0, 1)*MAT_ELEM(A, 0, 1) -
		MAT_ELEM(A, 0, 0)*MAT_ELEM(A, 1, 1)*MAT_ELEM(A, 2, 2) -
		2*MAT_ELEM(A, 0, 1)*MAT_ELEM(A, 0, 2)*MAT_ELEM(A, 1, 2);

	p = b*b - 3*c;
	q = 2*b*b*b - 9*b*c - 27*d;

	Delta = acos(q/(2*sqrt(p*p*p)));

	lambda_1 =(1.0/3)*(b+2.0*sqrt(p)*cos(Delta/3));
	lambda_2 =(1.0/3)*(b+2.0*sqrt(p)*cos((Delta+2.0*PI)/3));
	lambda_3 =(1.0/3)*(b+2.0*sqrt(p)*cos((Delta-2.0*PI)/3));
}

void ProgResDir::degreeOfAnisotropy(double &lambda_1, double &lambda_2, double &lambda_3,
							double &doa)
{
	//This equation comes from Thomsen's formula, with an error lesser than 1.061%
	/*double p = 1.6075;
	double A_ellip, V_ellip;
	double aux = (pow((lambda_1*lambda_2),p) + pow((lambda_1*lambda_3),p) +
				pow((lambda_2*lambda_3),p));
	A_ellip = 4.0*PI*pow(aux/3.0,1.0/p);
	V_ellip = (4.0/3.0)*PI*lambda_1*lambda_2*lambda_3;
//	std::cout << "lambda_1  = " << lambda_1 << "  lambda_2  = " << lambda_2 << "  lambda_3  = " << lambda_3 << std::endl;
//	std::cout << "A_ellip  = " << A_ellip << "  V_ellip  = " << V_ellip << std::endl;
	sph = pow(PI,(1.0/3.0))*pow(6.0*V_ellip,2.0/3.0)/A_ellip;
	*/
	MultidimArray<double> eigs(3);
	A1D_ELEM(eigs,0) = lambda_1;
	A1D_ELEM(eigs,1) = lambda_2;
	A1D_ELEM(eigs,2) = lambda_3;

	// Sort value and get threshold
	std::sort(&A1D_ELEM(eigs,0),&A1D_ELEM(eigs,2));

	double max = 1/sqrt(A1D_ELEM(eigs,0));
	double min = 1/sqrt(A1D_ELEM(eigs,2));


	//Defining DoA
	doa = (max-min)/(max+min);


	/*

	Matrix2D<double> invA;
	Matrix1D<double> inertia_vector(3), ellipsoid_axes(3);
	invA.initConstant(3,3,0.5);

	MAT_ELEM(invA, 0, 0) = -0.5;
	MAT_ELEM(invA, 1, 1) = -0.5;
	MAT_ELEM(invA, 2, 2) = -0.5;

	VEC_ELEM(inertia_vector, 0) = 3.0*lambda_1;
	VEC_ELEM(inertia_vector, 1) = 3.0*lambda_2;
	VEC_ELEM(inertia_vector, 2) = 3.0*lambda_3;

	ellipsoid_axes = invA*inertia_vector;


	lambda_1 = fabs(sqrt(round(VEC_ELEM(ellipsoid_axes, 0)*100.0)/100.0));
	lambda_2 = fabs(sqrt(round(VEC_ELEM(ellipsoid_axes, 1)*100.0)/100.0));
	lambda_3 = fabs(sqrt(round(VEC_ELEM(ellipsoid_axes, 2)*100.0)/100.0));

	if ((lambda_3>lambda_2) && (lambda_3>lambda_1))
		sph = (lambda_1*lambda_2)/(lambda_3*lambda_3);
	else if ((lambda_2>lambda_3) && (lambda_2>lambda_1))
		sph = (lambda_1*lambda_3)/(lambda_2*lambda_2);
		else if ((lambda_1>lambda_2) && (lambda_1>lambda_3))
			sph = (lambda_2*lambda_3)/(lambda_1*lambda_1);

*/
}


void ProgResDir::resolution2eval(int &count_res, double step,
								double &resolution, double &last_resolution,
								double &freq, double &freqL, double &freqH,
								int &last_fourier_idx,
								bool &continueIter,	bool &breakIter, bool &doNextIteration)
{
	resolution = maxRes - count_res*step;
	freq = sampling/resolution;


	++count_res;

	double Nyquist = 2*sampling;
	double aux_frequency, aux_frequency_next;
	int fourier_idx;

	DIGFREQ2FFT_IDX(freq, ZSIZE(VRiesz), fourier_idx);


//	std::cout << "Resolution = " << resolution << "   iter = " << count_res-1 << std::endl;
//	std::cout << "freq = " << freq << "   Fourier index = " << fourier_idx << std::endl;

	FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(VRiesz), aux_frequency);

	/////////////////////
	FFT_IDX2DIGFREQ((fourier_idx+1), ZSIZE(VRiesz), aux_frequency_next);
	/////////////////////

	freq = aux_frequency;
	freqH = aux_frequency_next;


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
		//std::cout << "Nyquist limit reached" << std::endl;
		breakIter = true;
		return;
	}
	if ( ( freq>0.49))// || (resolution > last_resolution) )
	{
		std::cout << "Nyquist limit reached" << std::endl;
		doNextIteration = false;
		return;
	}

	freqL = sampling/(resolution + step);

	int fourier_idx_2;

	DIGFREQ2FFT_IDX(freqL, ZSIZE(VRiesz), fourier_idx_2);
	double caca, caca2;

	if (fourier_idx_2 == fourier_idx)
	{
		if (fourier_idx > 0){
			//std::cout << " index low =  " << (fourier_idx - 1) << std::endl;
			FFT_IDX2DIGFREQ(fourier_idx - 1, ZSIZE(VRiesz), freqL);
		}
		else{
			freqL = sampling/(resolution + step);
		}
	}
}


void ProgResDir::createVectorField(Image<double> volume, MultidimArray<int> &mask)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	volume.getDimensions(Xdim, Ydim, Zdim, Ndim);

	int height_cylinder = 8;
	int height_cone = 4;
	int radius_cylinder = 2;
	int radius_cone = 4;

	std::vector<int> x_coor, y_coor, z_coor;

	int step_x = Xdim/(height_cone + height_cylinder);
	int step_y = Ydim/(height_cone + height_cylinder);
	int step_z = Zdim/(height_cone + height_cylinder);

	std::ofstream descrfile ("arrows.descr");
	descrfile << Xdim << " " << Ydim << " " << Zdim << " " << Ndim << "\n" << std::endl;

	//cyl = 1 0 0 0 5 5 20 0 0 0

	for (size_t k = 0; k< Zdim; k+step_z)
	{
		if (k != 0 )
			k = k - 0.5*Zdim;
		else
			k = 0.5*Zdim;
		for (size_t i = 0; i< Ydim; i+step_y)
		{
			if (i != 0 )
				i = i - 0.5*Ydim;
			else
				i = 0.5*Ydim;
			for (size_t j = 0; j< Zdim; j+step_x)
			{
				if (DIRECT_A3D_ELEM(mask, k,i,j) == 1)
				{
					if (j != 0 )
						j = j - 0.5*Xdim;
					else
						j = 0.5*Xdim;

					x_coor.push_back(i);
					y_coor.push_back(j);
					z_coor.push_back(k);
					descrfile << "my text here!" << std::endl;
				}
			}
		}
	}

	descrfile.close();

}


void ProgResDir::run()
{
	produceSideInfo();

	MetaData md;
	bool continueIter = false, breakIter = false;
	double criticalZ=icdf_gauss(significance);
	Image<double> Inertia_00, Inertia_01, Inertia_02, Inertia_11, Inertia_12, Inertia_22,
					SumRes;

	MultidimArray<double> &pInertia_00 = Inertia_00(), &pInertia_01 = Inertia_01(),
						  &pInertia_02 = Inertia_02(), &pInertia_11 = Inertia_11(),
						  &pInertia_12 = Inertia_12(), &pInertia_22 = Inertia_22(),
						  &pSumRes = SumRes();

	FileName fnInertia_00 = "Inertia_00.vol", fnInertia_01 = "Inertia_01.vol",
			 fnInertia_02 = "Inertia_02.vol", fnInertia_11 = "Inertia_11.vol",
			 fnInertia_12 = "Inertia_12.vol", fnInertia_22 = "Inertia_22.vol";

	double Nyquist = 2*sampling;
	double range = maxRes-minRes;
	double step = range/N_freq;

	if (step<0.1)
		step=0.1;

		std::cout << "Analyzing directions" << std::endl;

//	N_directions=1;

	for (size_t dir=0; dir<N_directions; dir++)
	{
		Image<double> outputResolution;
		outputResolution().resizeNoCopy(VRiesz);
		outputResolution().initZeros();
		MultidimArray<double> &pOutputResolution = outputResolution();
		MultidimArray<double> amplitudeMS, amplitudeMN;
		MultidimArray<int> mask_aux = mask();
		MultidimArray<int> &pMask = mask_aux;
		std::vector<double> list;
		double resolution, last_resolution = maxRes;  //A huge value for achieving last_resolution < resolution
		double freq, freqL, freqH, resVal, counter, resolution_2;
		double max_meanS = -1e38;
		double cut_value = 0.025;

		bool doNextIteration=true;
		bool lefttrimming = false;


		int fourier_idx, last_fourier_idx = -1, iter = 0, fourier_idx_2;
		int count_res = 0;
		double criticalW=-1;
		double angle_cone = ang_sampling;
		double rot = MAT_ELEM(angles, 0, dir);
		double tilt = MAT_ELEM(angles, 1, dir);
		std::cout << "--------------NEW DIRECTION--------------" << std::endl;
		std::cout << "direction = " << dir << "   rot = " << rot << "   tilt = " << tilt << std::endl;

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
			if (DIRECT_MULTIDIM_ELEM(pMask, n) == 1)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = maxRes;

		std::vector<double> noiseValues;
		FileName fnDebug;
		do
		{
			continueIter = false;
			breakIter = false;
			//std::cout << "--------------Frequency--------------" << std::endl;

			resolution2eval(count_res, step,
							resolution, last_resolution,
							freq, freqL, freqH,
							last_fourier_idx, continueIter, breakIter, doNextIteration);

			if (continueIter)
				continue;

			if (breakIter)
				break;

			std::cout << "resolution = " << resolution << "  resolutionL = " <<
					sampling/(freqL) << "  resolutionH = " << sampling/freqH << "  " << freqH << std::endl;

			std::cout << "\n"<< std::endl;


			list.push_back(resolution);

			if (iter <2)
				resolution_2 = list[0];
			else
				resolution_2 = list[iter - 2];

			fnDebug = "Signal";

			amplitudeMonogenicSignal3D(fftV, freq, freqL, freqH, amplitudeMS, iter, dir, fnDebug, angle_cone, rot, tilt);
			if (halfMapsGiven)
			{
				fnDebug = "Noise";
				amplitudeMonogenicSignal3D(*fftN, freq, freqL, freqH, amplitudeMN, iter, dir, fnDebug, angle_cone, rot, tilt);
			}


			double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
			noiseValues.clear();

			double amplitudeValue;
			if (halfMapsGiven)
			{
				if (noiseOnlyInHalves)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
					{
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						{
							amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
							double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
							sumS  += amplitudeValue;
							sumS2 += amplitudeValue*amplitudeValue;
							++NS;
							sumN  += amplitudeValueN;
							sumN2 += amplitudeValueN*amplitudeValueN;
							++NN;
						}
					}
				}
				else
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
					{
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						{
							amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
							sumS  += amplitudeValue;
							sumS2 += amplitudeValue*amplitudeValue;
							++NS;
						}
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=0)
						{
							double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
							sumN  += amplitudeValueN;
							sumN2 += amplitudeValueN*amplitudeValueN;
							++NN;
						}
					}
				}
			}
			else
			{
//				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
//				{
//
//					if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
//					{
//						amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
//						sumS  += amplitudeValue;
//						sumS2 += amplitudeValue*amplitudeValue;
//						++NS;
//					}
//					else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
//					{
//						amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
//						sumN  += amplitudeValue;
//						sumN2 += amplitudeValue*amplitudeValue;
//						++NN;
//					}
//				}
				double x_dir, y_dir, z_dir;

				x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
				y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
				z_dir = cos(tilt*PI/180);

				double uz, uy, ux;

				amplitudeMS.setXmippOrigin();


				MultidimArray<double> coneVol;
				coneVol.initZeros(amplitudeMS);
				int n=0;
				int z_size = ZSIZE(amplitudeMS);
				int x_size = XSIZE(amplitudeMS);
				int y_size = YSIZE(amplitudeMS);

				for(int k=0; k<z_size; ++k)
				{
					for(int i=0; i<y_size; ++i)
					{
						for(int j=0; j<x_size; ++j)
						{
							if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
							{
								amplitudeValue=DIRECT_A3D_ELEM(amplitudeMS, k,i,j);
								sumS  += amplitudeValue;
								sumS2 += amplitudeValue*amplitudeValue;
								++NS;
							}
							else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
							{
								uz = (k - z_size*0.5);
								ux = (j - x_size*0.5);
								uy = (i - y_size*0.5);

								double iun = 1/sqrt(ux*ux + uy*uy + uz*uz);

								//BE CAREFULL with the order
								double dotproduct;
								dotproduct = (ux*y_dir + uy*x_dir + uz*z_dir)*iun;

								double acosine = acos(fabs(dotproduct));

								//DIRECT_A3D_ELEM(coneVol, k,i,j) = 1;

								if (acosine<(PI*20/180))
								{
									DIRECT_A3D_ELEM(coneVol, k,i,j) = 1;
									amplitudeValue=DIRECT_A3D_ELEM(amplitudeMS, k,i,j);
									sumN  += amplitudeValue;
									sumN2 += amplitudeValue*amplitudeValue;
									++NN;
								}
							}
							n++;
						}
					}
				}

//				if (iter == 0)
//				{
//				Image<double> img;
//
//				FileName iternumber;
//				iternumber = formatString("cone_%i.vol", dir);
//				img = coneVol;
//				img.write(iternumber);
//				}
			}



			if ( (NS/NVoxelsOriginalMask)<cut_value ) //when the 2.5% is reached then the iterative process stops
			{
				std::cout << "Search of resolutions stopped due to mask has been completed" << std::endl;
				doNextIteration =false;
				Nvoxels = 0;
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
				{
//				  if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n) == 0)
//					DIRECT_MULTIDIM_ELEM(pMask, n) = 0;
//				  else
//				  {
//					Nvoxels++;
//					DIRECT_MULTIDIM_ELEM(pMask, n) = 1;
//				  }
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
			double sigma2S=sumS2/NS-meanS*meanS;
			double meanN=sumN/NN;
			double sigma2N=sumN2/NN-meanN*meanN;

			if (meanS>max_meanS)
				max_meanS = meanS;

			if (meanS<0.00001*max_meanS)
			{
				std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
				std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
				std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
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

				// Is the mean inside the signal significantly different from the noise?
				double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
				//#ifdef DEBUG
					std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
					std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
					std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
					std::cout << "  z= " << z << " (" << criticalZ << ")" << std::endl;
				//#endif
				if (z<criticalZ)
				{
					criticalW = freq;
					std::cout << "Search stopped due to z>Z (hypothesis test)" << std::endl;
					doNextIteration=false;
				}
				if (doNextIteration)
					if (resolution <= (minRes-0.001))
						doNextIteration = false;
				}
			}

			iter++;
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
		if (fnSym!="c1")
		{
			SymList SL;
			SL.readSymmetryFile(fnSym);
			MultidimArray<double> VSimetrized;
			symmetrizeVolume(SL, pOutputResolution, VSimetrized, LINEAR, DONT_WRAP);
			outputResolution() = VSimetrized;
			VSimetrized.clear();
		}

		#ifdef DEBUG_SYMMETRY
			outputResolution.write("resolution_simple_simmetrized.vol");

		Image<double> saveImg;
		saveImg = pOutputResolution;
		FileName fnres;
		fnres = formatString("resolution_simmetrized_%i.vol", dir);
		saveImg.write(fnres);
		saveImg.clear();
		#endif

		//MultidimArray<double> resolutionFiltered, resolutionChimera;
		//postProcessingLocalResolutions(pOutputResolution, list, resolutionChimera, cut_value, pMask);


		//////////////////
		//INERTIA MOMENT//
		//////////////////
		if (dir == 0)
		{
			pInertia_00.initZeros(pOutputResolution);
			pInertia_01.initZeros(pOutputResolution);
			pInertia_02.initZeros(pOutputResolution);
			pInertia_11.initZeros(pOutputResolution);
			pInertia_12.initZeros(pOutputResolution);
			pInertia_22.initZeros(pOutputResolution);
			pSumRes.initZeros(pOutputResolution);

		}
		else
		{
			Inertia_00.read(fnInertia_00);
			Inertia_01.read(fnInertia_01);
			Inertia_02.read(fnInertia_02);
			Inertia_11.read(fnInertia_11);
			Inertia_12.read(fnInertia_12);
			Inertia_22.read(fnInertia_22);

			SumRes.read("sumRes.vol");
		}

		inertiaMatrix(pOutputResolution, pInertia_00, pInertia_01,
					pInertia_02, pInertia_11, pInertia_12, pInertia_22,
					pSumRes,	rot, tilt);

		SumRes.write("sumRes.vol");
		Inertia_00.write(fnInertia_00);
		Inertia_01.write(fnInertia_01);
		Inertia_02.write(fnInertia_02);
		Inertia_11.write(fnInertia_11);
		Inertia_12.write(fnInertia_12);
		Inertia_22.write(fnInertia_22);

		Inertia_00.clear();
		Inertia_01.clear();
		Inertia_02.clear();
		Inertia_11.clear();
		Inertia_12.clear();
		Inertia_22.clear();
		SumRes.clear();
		//////////////////////////////

//		std::cout << "after postprocessing" << std::endl;
		Image<double> VarianzeResolution, AvgResolution;

		//They were initialized in producesideinfo()
		AvgResolution.read(fnOut);
		VarianzeResolution.read(fnVar);

		MultidimArray<double> &pAvgResolution = AvgResolution();
		MultidimArray<double> &pVarianzeResolution = VarianzeResolution();

		// Calculating min and max resolution volumes
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		{
			if (DIRECT_MULTIDIM_ELEM(mask(), n) == 1)
			{
			double Resvoxel = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
			DIRECT_MULTIDIM_ELEM(pAvgResolution, n) += Resvoxel;
			DIRECT_MULTIDIM_ELEM(pVarianzeResolution, n) += Resvoxel*Resvoxel;
			}
			else
			{
				DIRECT_MULTIDIM_ELEM(pAvgResolution, n) = 0;
				DIRECT_MULTIDIM_ELEM(pVarianzeResolution, n) = 0;
			}
		}

		AvgResolution.write(fnOut);
		VarianzeResolution.write(fnVar);

		double sumS=0, sumS2=0, N_elem=0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		{
			double value=DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
			if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>0)
			{
				sumS  += value;
				sumS2 += value*value;
				++N_elem;
			}
		}

		double avg = sumS/N_elem;
		MAT_ELEM(angles, 2, dir) = avg;
		MAT_ELEM(angles, 3, dir) = sumS2/N_elem-avg*avg;

		/////////////////////////
		size_t objId;
		objId = md.addObject();
		md.setValue(MDL_ANGLE_ROT, rot*180/PI, objId);
		md.setValue(MDL_ANGLE_TILT, tilt*180/PI, objId);
		md.setValue(MDL_RESOLUTION_FREQREAL, avg, objId);
		md.setValue(MDL_STDDEV, MAT_ELEM(angles, 3, dir), objId);

		md.write(fnMd);

		/////////////////////////
//		#ifdef DEBUG_DIR
//		Image<double> saveImg;
//		saveImg = pOutputResolution;
//		FileName fnres = formatString("resolution_dir_%i.vol", dir);
//		saveImg.write(fnres);
//		saveImg.clear();
//		#endif

		pOutputResolution.clear();
		pVarianzeResolution.clear();
		list.clear();

//		exit(0);
		std::cout << "----------------direction-finished----------------" << std::endl;
	}


	Inertia_00.read(fnInertia_00);
	Inertia_01.read(fnInertia_01);
	Inertia_02.read(fnInertia_02);
	Inertia_11.read(fnInertia_11);
	Inertia_12.read(fnInertia_12);
	Inertia_22.read(fnInertia_22);

	SumRes.read("sumRes.vol");

	Matrix2D<double> InertiaMatrix;
	InertiaMatrix.initZeros(3,3);

	Image<double> AvgResolution;
	AvgResolution.read(fnOut);
	//MultidimArray<double> &pAvgResolution = AvgResolution();

	double lambda_1, lambda_2, lambda_3, sph;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pInertia_11)
	{
		if (DIRECT_MULTIDIM_ELEM(mask(),n) == 1 )
		{
			double val = 1/DIRECT_MULTIDIM_ELEM(pSumRes,n);

			MAT_ELEM(InertiaMatrix, 0, 0) = DIRECT_MULTIDIM_ELEM(pInertia_00,n)*val;
			MAT_ELEM(InertiaMatrix, 0, 1) = DIRECT_MULTIDIM_ELEM(pInertia_01,n)*val;
			MAT_ELEM(InertiaMatrix, 0, 2) = DIRECT_MULTIDIM_ELEM(pInertia_02,n)*val;
			MAT_ELEM(InertiaMatrix, 1, 0) = DIRECT_MULTIDIM_ELEM(pInertia_01,n)*val;//MAT_ELEM(InertiaMatrix, 0, 1);
			MAT_ELEM(InertiaMatrix, 1, 1) = DIRECT_MULTIDIM_ELEM(pInertia_11,n)*val;
			MAT_ELEM(InertiaMatrix, 1, 2) = DIRECT_MULTIDIM_ELEM(pInertia_12,n)*val;
			MAT_ELEM(InertiaMatrix, 2, 0) = DIRECT_MULTIDIM_ELEM(pInertia_02,n)*val;//MAT_ELEM(InertiaMatrix, 0, 2)
			MAT_ELEM(InertiaMatrix, 2, 1) = DIRECT_MULTIDIM_ELEM(pInertia_12,n)*val;//MAT_ELEM(InertiaMatrix, 1, 2)
			MAT_ELEM(InertiaMatrix, 2, 2) = DIRECT_MULTIDIM_ELEM(pInertia_22,n)*val;


			//std::cout << "InertiaMatrix = " << InertiaMatrix << std::endl;

			diagSymMatrix3x3(InertiaMatrix, N_directions, lambda_1, lambda_2, lambda_3);


//			std::cout << "lambda_1 = " << lambda_1 <<
//					   "  lambda_2 = " << lambda_2 <<
//					   "  lambda_3 = " << lambda_3 << std::endl;


			degreeOfAnisotropy(lambda_1, lambda_2, lambda_3, sph);

//			std::cout << "LAMBDA_1 = " << lambda_1 <<
//					   "  LAMBDA_2 = " << lambda_2 <<
//					   "  LAMBDA_3 = " << lambda_3 << std::endl;

			DIRECT_MULTIDIM_ELEM(pInertia_00,n) = lambda_1;
			DIRECT_MULTIDIM_ELEM(pInertia_01,n) = lambda_2;
			DIRECT_MULTIDIM_ELEM(pInertia_02,n) = lambda_3;

			DIRECT_MULTIDIM_ELEM(pInertia_11,n) = sph;
//			DIRECT_MULTIDIM_ELEM(pInertia_12,n) = (lambda_1-lambda_3)/(lambda_1+lambda_3);

		}
		else
		{
			DIRECT_MULTIDIM_ELEM(pInertia_11,n) = 0;
			DIRECT_MULTIDIM_ELEM(pInertia_12,n) =0;
		}
	}
	Inertia_00.write("lambda_1.vol");
	Inertia_01.write("lambda_2.vol");
	Inertia_02.write("lambda_3.vol");
	Inertia_11.write(fnDoA);


	Image<double> VarianzeResolution;
	AvgResolution.read(fnOut);
	MultidimArray<double> &pAvgResolution = AvgResolution();

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pAvgResolution)
	{
		DIRECT_MULTIDIM_ELEM(pAvgResolution, n) /= N_directions;
		DIRECT_MULTIDIM_ELEM(pAvgResolution, n) *= DIRECT_MULTIDIM_ELEM(mask(), n);
	}
	AvgResolution.write(fnOut);
	#ifdef DEBUG
	Image<int> saveImg_int;
	FileName fnm;
	saveImg_int = mask();
	saveImg_int.write("maskfinalk.vol");
	saveImg_int.clear();
	#endif
	VarianzeResolution.read(fnVar);

	MultidimArray<double> &pVarianzeResolution = VarianzeResolution();

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pAvgResolution)
	{
		if (DIRECT_MULTIDIM_ELEM(mask(), n) == 1)
		{
			double resmean = DIRECT_MULTIDIM_ELEM(pAvgResolution, n);
			DIRECT_MULTIDIM_ELEM(pVarianzeResolution, n) = (DIRECT_MULTIDIM_ELEM(pVarianzeResolution, n)/N_directions) - (resmean*resmean/(N_directions*N_directions));
		}
	}
//	AvgResolution.write(fnDoA);
	VarianzeResolution.write(fnVar);
}
