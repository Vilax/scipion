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
//#define MONO_AMPLITUDE
//define DEBUG_SYMMETRY

void ProgResDir::readParams()
{
	fnVol = getParam("--vol");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	sampling = getDoubleParam("--sampling_rate");
	ang_sampling = getDoubleParam("--angular_sampling");
	R = getDoubleParam("--volumeRadius");
	Rparticle = getDoubleParam("--particleRadius");
	fnSym = getParam("--sym");
	significance = getDoubleParam("--significance");
	fnDoA = getParam("--doa_vol");
	fnDirections = getParam("--directions");
}


void ProgResDir::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file=\"\">        : Input volume");
	addParamsLine("  [--mask <vol_file=\"\">]     : Mask defining the macromolecule");
	addParamsLine("                               :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                               :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--sym <symmetry>]			  : Symmetry (c1, c2, c3,..d1, d2, d3,...)");
	addParamsLine("  [--sampling_rate <s=1>]      : Sampling rate (A/px)");
	addParamsLine("  [--angular_sampling <s=15>]  : Angular Sampling rate (degrees)");
	addParamsLine("  [--volumeRadius <s=100>]     : This parameter determines the radius of a sphere where the volume is");
	addParamsLine("  [--particleRadius <s=100>]     : This parameter determines the radius of the particle");
	addParamsLine("  [--varVol <vol_file=\"\">]   : Output filename with varianze resolution volume");
	addParamsLine("  [--significance <s=0.95>]    : The level of confidence for the hypothesis test.");
	addParamsLine("  [--doa_vol <vol_file=\"\">]  : Output filename with DoA volume");
	addParamsLine("  [--directions <vol_file=\"\">]  : Output preffered directions");
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
	Image<double> AvgResoltion;
	AvgResoltion().resizeNoCopy(inputVol);
	AvgResoltion().initZeros();
	AvgResoltion.write(fnOut);
	AvgResoltion.clear();

	N_smoothing = 7;
	NVoxelsOriginalMask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
		if (A3D_ELEM(pMask, k, i, j) == 1)
			++NVoxelsOriginalMask;
		if (i*i+j*j+k*k > (R-N_smoothing)*(R-N_smoothing))
			A3D_ELEM(pMask, k, i, j) = -1;
	}

	size_t xrows = angles.mdimx;
//	Matrix2D<double> aaa;
//
//	aaa.initConstant(3, 4, 5);
//	MAT_ELEM(aaa, 1, 2) = 0;
//	std::cout << aaa << std::endl;
//	std::cout << "NVoxelsOriginalMask = " << NVoxelsOriginalMask << std::endl;

	resolutionMatrix.initConstant(xrows, NVoxelsOriginalMask, maxRes);
	inertiaMatrixVariable.initZeros(7, NVoxelsOriginalMask);
	maskMatrix.initConstant(1, NVoxelsOriginalMask, 1);

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
		MAT_ELEM(angles, 0, k-1) = MAT_ELEM(aux_angles,0, k);
		MAT_ELEM(angles, 1, k-1) = MAT_ELEM(aux_angles,1, k);
		std::cout << "k=" << k << "  rot = "  << MAT_ELEM(angles, 0, k-1) <<
								  "  tilt = " << MAT_ELEM(angles, 1, k-1) << std::endl;
	}
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


void ProgResDir::amplitudeMonogenicSignal3D_fast(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, int dir, FileName fnDebug,
		double rot, double tilt)
{
	fftVRiesz.initZeros(myfftV);
	MultidimArray<double> coneVol;
	coneVol.initZeros(myfftV);
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
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= DIRECT_MULTIDIM_ELEM(conefilter, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
//					DIRECT_MULTIDIM_ELEM(coneVol, n) = DIRECT_MULTIDIM_ELEM(conefilter, n);
				} else if (un>freq)
				{
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= DIRECT_MULTIDIM_ELEM(conefilter, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
					DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
//					DIRECT_MULTIDIM_ELEM(coneVol, n) = DIRECT_MULTIDIM_ELEM(conefilter, n);
				}
				++n;
			}
		}
	}

	#ifdef DEBUG_DIR
//	if ( (count == 0) )
//	{
//		Image<double> direction;
//		direction = coneVol;
//		direction.write(formatString("cone_%i_%i.vol", dir+1, count));
//	}
	#endif

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	#ifdef DEBUG_DIR
//	if (count == 0)
//	{
//		Image<double> filteredvolume;
//		filteredvolume = VRiesz;
//		filteredvolume.write(formatString("Volumen_filtrado_%i_%i.vol", dir+1,count));
//	}
	#endif


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


		#ifdef MONO_AMPLITUDE
//		Image<double> saveImg2;
//		saveImg2 = amplitude;
//		if (fnDebug.c_str() != "")
//		{
//			FileName iternumber = formatString("smoothed_volume_%i_%i.vol", dir+1, count);
//			saveImg2.write(fnDebug+iternumber);
//		}
//		saveImg2.clear();
		#endif


	//amplitude.setXmippOrigin();

	//transformer_inv.FourierTransform(amplitude, fftVRiesz, false);


	// Low pass filter the monogenic amplitude
	lowPassFilter.w1 = freq;
	amplitude.setXmippOrigin();
	lowPassFilter.applyMaskSpace(amplitude);

//    double raised_w = PI/(freqL-freq);
//
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
//	{
//		double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//
//	}
//	transformer_inv.inverseFourierTransform();

	#ifdef MONO_AMPLITUDE
//
//	if (fnDebug.c_str() != "")
//	{
//		saveImg2 = amplitude;
//		iternumber = formatString("_Filtered_Amplitude_%i_%i.vol", dir+1, count);
//		saveImg2.write(fnDebug+iternumber);
//	}
//	saveImg2.clear();
	#endif // DEBUG

}

void ProgResDir::defineCone(MultidimArray< std::complex<double> > &myfftV,
		MultidimArray<double> &conefilter, double rot, double tilt)
{
	conefilter.initZeros(myfftV);

	// Filter the input volume and add it to amplitude
	long n=0;

	#ifdef DEBUG_DIR
	MultidimArray<double> coneVol;
	coneVol.initZeros(iu);
	#endif

	double x_dir, y_dir, z_dir;

	x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir = cos(tilt*PI/180);

	double uz, uy, ux;
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

				//BE CAREFULL with the order
				//double dotproduct = (uy*x_dir + ux*y_dir + uz*z_dir)*iun;
				double dotproduct = (ux*x_dir + uy*y_dir + uz*z_dir)*iun;
				double acosine = acos(fabs(dotproduct));

				//4822.53 mean a smoothed cone angle of 20 degrees
				double arg_exp = acosine*acosine*acosine*acosine*4822.53;
				DIRECT_MULTIDIM_ELEM(conefilter, n) = exp(-arg_exp);
				++n;
			}
		}
	}
}


void ProgResDir::inertiaMatrix(MultidimArray<double> &resolutionVol,
							   MultidimArray<double> &Inertia_11,
							   MultidimArray<double> &Inertia_12,
							   MultidimArray<double> &Inertia_13,
							   MultidimArray<double> &Inertia_22,
							   MultidimArray<double> &Inertia_23,
							   MultidimArray<double> &Inertia_33,
							   MultidimArray<double> &SumRes,
							   double rot, double tilt, size_t dir)
{
	double x_dir, y_dir, z_dir, resVal, x_dir_sym, y_dir_sym, z_dir_sym, r_xyz2;
	x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir = cos(tilt*PI/180);
	x_dir_sym = -sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir_sym = -sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir_sym = -cos(tilt*PI/180);

	//std::cout << "x_dir = " << x_dir << "  y_dir = " << y_dir<< "  z_dir = " << z_dir << std::endl;

	size_t idx = 0;



	int nn=0;
	for(int k=0; k<ZSIZE(resolutionVol); ++k)
	{
		for(int i=0; i<YSIZE(resolutionVol); ++i)
		{
			for(int j=0; j<XSIZE(resolutionVol); ++j)
			{
				if (DIRECT_MULTIDIM_ELEM(mask(), nn) == 1)
				{
					resVal = DIRECT_MULTIDIM_ELEM(resolutionVol,nn);//*DIRECT_MULTIDIM_ELEM(resolutionVol,n);
					double aux_val = DIRECT_A3D_ELEM(resolutionVol,k,i,j);
					r_xyz2 = resVal*resVal;
		//			r_xyz2 = 10;

					if ( ( (k == 30)  || (k == 35) || (k == 40) || (k == 45) || (k == 50) ||
							(k == 55) || (k == 60) || (k == 65) || (k == 70) || (k == 75)
							|| (k == 80) || (k == 85)) && (i==52) && (j==62) )
					{
						MetaData md;
						size_t objId;
						FileName fn_md;
						fn_md = formatString("res_%i.xmd", k);
						md.read(fn_md);

						objId = md.addObject();
						md.setValue(MDL_IDX, dir, objId);
						md.setValue(MDL_XCOOR, i, objId);
						md.setValue(MDL_YCOOR, j, objId);
						md.setValue(MDL_ZCOOR, k, objId);
						md.setValue(MDL_ANGLE_ROT, rot, objId);
						md.setValue(MDL_ANGLE_TILT, tilt, objId);
						md.setValue(MDL_RESOLUTION_FREQREAL, resVal, objId);
						md.setValue(MDL_RESOLUTION_FREQ, aux_val, objId);
						md.write(fn_md);
					}
					DIRECT_MULTIDIM_ELEM(Inertia_11,nn) += r_xyz2*(1-x_dir*x_dir);
					DIRECT_MULTIDIM_ELEM(Inertia_12,nn) -= r_xyz2*x_dir*y_dir;
					DIRECT_MULTIDIM_ELEM(Inertia_13,nn) -= r_xyz2*x_dir*z_dir;
					DIRECT_MULTIDIM_ELEM(Inertia_22,nn) += r_xyz2*(1-y_dir*y_dir);
					DIRECT_MULTIDIM_ELEM(Inertia_23,nn) -= r_xyz2*y_dir*z_dir;
					DIRECT_MULTIDIM_ELEM(Inertia_33,nn) += r_xyz2*(1-z_dir*z_dir);

					DIRECT_MULTIDIM_ELEM(Inertia_11,nn) += r_xyz2*(1-x_dir_sym*x_dir_sym);
					DIRECT_MULTIDIM_ELEM(Inertia_12,nn) -= r_xyz2*x_dir_sym*y_dir_sym;
					DIRECT_MULTIDIM_ELEM(Inertia_13,nn) -= r_xyz2*x_dir_sym*z_dir_sym;
					DIRECT_MULTIDIM_ELEM(Inertia_22,nn) += r_xyz2*(1-y_dir_sym*y_dir_sym);
					DIRECT_MULTIDIM_ELEM(Inertia_23,nn) -= r_xyz2*y_dir_sym*z_dir_sym;
					DIRECT_MULTIDIM_ELEM(Inertia_33,nn) += r_xyz2*(1-z_dir_sym*z_dir_sym);

					DIRECT_MULTIDIM_ELEM(SumRes,nn) += 2;
					++idx;
				}
				++nn;
			}
		}
	}

}

void ProgResDir::inertiaMatrixNew(Matrix2D<double> &resolutionMatrix,
							   Matrix2D<double> &inertiaMatrix,
							   double rot, double tilt, size_t dir)
{
	double x_dir, y_dir, z_dir, resVal, x_dir_sym, y_dir_sym, z_dir_sym, r_xyz2;
	x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir = cos(tilt*PI/180);
	x_dir_sym = -sin(tilt*PI/180)*cos(rot*PI/180);
	y_dir_sym = -sin(tilt*PI/180)*sin(rot*PI/180);
	z_dir_sym = -cos(tilt*PI/180);

	int nn=0;
	int maskpos = 0;
	std::cout << "antes del for" << std::endl;
	for(int k=0; k<ZSIZE(mask()); ++k)
	{
		for(int i=0; i<YSIZE(mask()); ++i)
		{
			for(int j=0; j<XSIZE(mask()); ++j)
			{
				if (DIRECT_MULTIDIM_ELEM(mask(), nn) == 1)
				{
					resVal = MAT_ELEM(resolutionMatrix, dir, maskpos);
					r_xyz2 = resVal*resVal;

					if ( ( (k == 30)  || (k == 35) || (k == 40) || (k == 45) || (k == 50) ||
							(k == 55) || (k == 60) || (k == 65) || (k == 70) || (k == 75)
							|| (k == 80) || (k == 85)) && (i==52) && (j==62) )
					{
						MetaData md;
						size_t objId;
						FileName fn_md;
						fn_md = formatString("res_%i.xmd", k);
						md.read(fn_md);

						objId = md.addObject();
						md.setValue(MDL_IDX, dir, objId);
						md.setValue(MDL_XCOOR, i, objId);
						md.setValue(MDL_YCOOR, j, objId);
						md.setValue(MDL_ZCOOR, k, objId);
						md.setValue(MDL_ANGLE_ROT, rot, objId);
						md.setValue(MDL_ANGLE_TILT, tilt, objId);
						md.setValue(MDL_RESOLUTION_FREQREAL, resVal, objId);
						md.write(fn_md);
					}

					MAT_ELEM(inertiaMatrix, 0, maskpos) += (r_xyz2*(1-x_dir*x_dir) + r_xyz2*(1-x_dir_sym*x_dir_sym));
					MAT_ELEM(inertiaMatrix, 1, maskpos) -= (r_xyz2*x_dir*y_dir + r_xyz2*x_dir_sym*y_dir_sym);
					MAT_ELEM(inertiaMatrix, 2, maskpos) -= (r_xyz2*x_dir*z_dir + r_xyz2*x_dir_sym*z_dir_sym);
					MAT_ELEM(inertiaMatrix, 3, maskpos) += (r_xyz2*(1-y_dir*y_dir) + r_xyz2*(1-y_dir_sym*y_dir_sym));
					MAT_ELEM(inertiaMatrix, 4, maskpos) -= (r_xyz2*y_dir*z_dir + r_xyz2*y_dir_sym*z_dir_sym);
					MAT_ELEM(inertiaMatrix, 5, maskpos) += (r_xyz2*(1-z_dir*z_dir) + r_xyz2*(1-z_dir_sym*z_dir_sym));
					MAT_ELEM(inertiaMatrix, 6, maskpos) += 2;
					++maskpos;
				}
				++nn;
			}
		}
	}
	std::cout << "despues del for" << std::endl;

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

void ProgResDir::degreeOfAnisotropy(Matrix1D<double> eigenvalues,
									Matrix2D<double> eigenvectors,
									double &doa,
							double &direction_x, double &direction_y, double &direction_z,
							int &counter)
{
	++counter;

	direction_x = MAT_ELEM(eigenvectors, 0, 1);
	direction_y = MAT_ELEM(eigenvectors, 1, 1);
	direction_z = MAT_ELEM(eigenvectors, 2, 1);

	Matrix2D<double> invA;
	Matrix1D<double> inertia_vector(3), ellipsoid_axes(3);
	invA.initConstant(3,3,0.5);

	MAT_ELEM(invA, 0, 0) = -0.5;
	MAT_ELEM(invA, 1, 1) = -0.5;
	MAT_ELEM(invA, 2, 2) = -0.5;

	VEC_ELEM(inertia_vector, 0) = 3.0*VEC_ELEM(eigenvalues,0);
	VEC_ELEM(inertia_vector, 1) = 3.0*VEC_ELEM(eigenvalues,1);
	VEC_ELEM(inertia_vector, 2) = 3.0*VEC_ELEM(eigenvalues,2);

	ellipsoid_axes = invA*inertia_vector;

	//lambda_1, lambda_2, lambda_3 are the square of lambda_1, lambda_2, lambda_3
	double lambda_1 = sqrt(VEC_ELEM(ellipsoid_axes, 0));
	double lambda_2 = sqrt(VEC_ELEM(ellipsoid_axes, 1));
	double lambda_3 = sqrt(VEC_ELEM(ellipsoid_axes, 2));

	double lambda_values[3] = {lambda_1, lambda_2, lambda_3};

	std::sort(lambda_values, lambda_values+3);

	double aux = (lambda_values[0]/lambda_values[2]);
	doa = 1 - aux*aux;

	if (counter == 34)
	{
		std::cout << "lambda_1 = " << lambda_1 << std::endl;
		std::cout << "lambda_2 = " << lambda_2 << std::endl;
		std::cout << "lambda_3 = " << lambda_3 << std::endl;

		std::cout << "doa = " << doa << std::endl;
		std::cout << "lambda_list[2] = " << lambda_values[2] << std::endl;
		std::cout << "lambda_list[1] = " << lambda_values[1] << std::endl;
		std::cout << "lambda_list[0] = " << lambda_values[0] << std::endl;
	}
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

//	std::cout << "freq_H = " << freqH << std::endl;
//	std::cout << "freq_L = " << freqL << std::endl;

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

void ProgResDir::defineDirection(Matrix1D<double> &r0, Matrix1D<double> &rF,
							Matrix2D<double> &direction, double &eigenvalue, double &eigenvalue_max,  int eigdir,
							int k, int i, int j)
{
//	std::cout << "direction = " << MAT_ELEM(direction,0,eigdir) <<
//			" " << MAT_ELEM(direction,1,eigdir) << " " << MAT_ELEM(direction,2,eigdir) << std::endl;
//	std::cout << "eigenvalue = " << eigenvalue << std::endl;

	double aux = eigenvalue/eigenvalue_max;

	VECTOR_R3(r0,j-5.0*aux*MAT_ELEM(direction,0,eigdir),
				 i-5.0*aux*MAT_ELEM(direction,1,eigdir),
				 k-5.0*aux*MAT_ELEM(direction,2,eigdir));

	VECTOR_R3(rF,j+5.0*aux*MAT_ELEM(direction,0,eigdir),
				 i+5.0*aux*MAT_ELEM(direction,1,eigdir),
				 k+5.0*aux*MAT_ELEM(direction,2,eigdir));
	std::cout << " " << std::endl;
}

void ProgResDir::defineSegment(Matrix1D<double> &r0, Matrix1D<double> &rF,
							MultidimArray<int> &arrows, double &elongation)
{
	Matrix1D<double> r(3);
	XX(r)=(1-elongation)*XX(r0)+elongation*XX(rF);
	YY(r)=(1-elongation)*YY(r0)+elongation*YY(rF);
	ZZ(r)=(1-elongation)*ZZ(r0)+elongation*ZZ(rF);

	//A3D_ELEM(arrows,(int)round(XX(r)),(int)round(YY(r)),(int)round(ZZ(r)))=1;

	DIRECT_A3D_ELEM(arrows,(int)round(ZZ(r)),(int)round(YY(r)),(int)round(XX(r)))=1;

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
//	std::cout << "maxRes = " << maxRes << std::endl;
//	std::cout << "minRes = " << minRes << std::endl;
//	std::cout << "N_freq = " << N_freq << std::endl;
//	std::cout << "step = " << step << std::endl;

	Image<double> outputResolution;
	MultidimArray<double> amplitudeMS, amplitudeMN;

	double w, wH;
	int volsize = ZSIZE(VRiesz);
	FFT_IDX2DIGFREQ(10, volsize, w)
	FFT_IDX2DIGFREQ(11, volsize, wH)

	double AvgNoise;
	AvgNoise = firstMonoResEstimation(fftV, w, wH, amplitudeMS)/9.0;


//	N_directions=1;

	std::cout << "N_directions = " << N_directions << std::endl;

	double cone_angle = 20.0; //(degrees)


	for (size_t dir=0; dir<N_directions; dir++)
	{
		outputResolution().initZeros(VRiesz);
		MultidimArray<double> &pOutputResolution = outputResolution();
		double freq, freqL, freqH, counter, resolution_2;
		MultidimArray<int> mask_aux = mask();
		MultidimArray<int> &pMask = mask_aux;
		std::vector<double> list;
		double resolution;  //A huge value for achieving last_resolution < resolution

		double max_meanS = -1e38;
		double cut_value = 0.025;

		bool doNextIteration=true;

		int fourier_idx = 5, last_fourier_idx = -1, iter = 0, fourier_idx_2;
		int count_res = 0;
		double criticalW=-1;
		double rot = MAT_ELEM(angles, 0, dir);
		double tilt = MAT_ELEM(angles, 1, dir);
		std::cout << "--------------NEW DIRECTION--------------" << std::endl;
		std::cout << "direction = " << dir+1 << "   rot = " << rot << "   tilt = " << tilt << std::endl;

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
			if (DIRECT_MULTIDIM_ELEM(pMask, n) == 1)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = maxRes;

		std::vector<double> noiseValues;
		FileName fnDebug;
		double last_resolution = 0;

		defineCone(fftV, conefilter, rot, tilt);

		int aa = 0;
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

			amplitudeMonogenicSignal3D_fast(fftV, freq, freqH, freqL, amplitudeMS, iter, dir, fnDebug, rot, tilt);

			double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
			noiseValues.clear();

			double amplitudeValue;

			double x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
			double y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
			double z_dir = cos(tilt*PI/180);

			double uz, uy, ux;

			amplitudeMS.setXmippOrigin();

			//TODO: check if can be taken out side the loop
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
							amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
							sumS  += amplitudeValue;
//							sumS2 += amplitudeValue*amplitudeValue;
							++NS;
						}
						else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
						{
							uz = (k - z_size*0.5);
							ux = (j - x_size*0.5);
							uy = (i - y_size*0.5);

							double rad = sqrt(ux*ux + uy*uy + uz*uz);
//							std::cout << "rad = " << rad << std::endl;
							double iun = 1/rad;

							//BE CAREFULL with the order
							double dotproduct = (uy*y_dir + ux*x_dir + uz*z_dir)*iun;

							double acosine = acos(dotproduct);

							if (((acosine<(PI*cone_angle/180)) || (acosine>(PI*(180-cone_angle)/180)) )
									&& (rad>Rparticle))
							{

								DIRECT_MULTIDIM_ELEM(coneVol, n) = 1;
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

			#ifdef DEBUG_DIR
				if (iter == 0)
				{
				Image<double> img;

				FileName iternumber;
				iternumber = formatString("cone_noise_%i_%i.vol", dir, iter);
				img = coneVol;
				img.write(iternumber);
				}
			#endif

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

					size_t maskPos = 0;
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
					{
						if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
						{
							if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
							{
								DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution;//sampling/freq;
								MAT_ELEM(resolutionMatrix, dir, maskPos) = resolution;
								MAT_ELEM(maskMatrix, dir, maskPos) = 1;
							}
							else
							{
								MAT_ELEM(maskMatrix, dir, maskPos) += 1;
								if (MAT_ELEM(maskMatrix, dir, maskPos) >2)
								{
									MAT_ELEM(maskMatrix, dir, maskPos) = 0;
									MAT_ELEM(resolutionMatrix, dir, maskPos) = resolution_2;
									DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = resolution_2; //resolution + counter*step;
								}
							}
							++maskPos;
						}
					}
					#ifdef DEBUG_MASK
					FileName fnmask_debug;
					fnmask_debug = formatString("maske_%i.vol", iter);
					mask.write(fnmask_debug);
					#endif

					//#ifdef DEBUG
						std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
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
		std::cout << "Finalizo" << std::endl;

//		amplitudeMN.clear();
//		amplitudeMS.clear();
//		fftVRiesz.clear();
		
		int NVoxelsOriginalMask_bis = 0;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(mask())
		{
			if (A3D_ELEM(mask(), k, i, j) == 1)
				++NVoxelsOriginalMask_bis;
		}
		std::cout << "NVoxelsOriginalMask_bis = " << NVoxelsOriginalMask_bis << std::endl;
		//////////////////
		//INERTIA MOMENT//
		//////////////////

		std::cout << "antes del inertia matrix New" << std::endl;
		inertiaMatrixNew(resolutionMatrix, inertiaMatrixVariable, rot, tilt, dir);
		std::cout << "despues del inertia matrix New" << std::endl;
		#ifdef DEBUG_DIR
		size_t maskPos=0;
		Image<double> ResolutionVol;
		MultidimArray<double> &pResolutionVol = ResolutionVol();

		pResolutionVol.initZeros(pOutputResolution);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		{
			if (DIRECT_MULTIDIM_ELEM(mask(), n) == 1)
			{
				double myres = MAT_ELEM(resolutionMatrix, dir, maskPos);
				DIRECT_MULTIDIM_ELEM(pResolutionVol, n) = myres;
				++maskPos;
			}
		}

		Image<double> saveImg;
		saveImg = pOutputResolution;
		FileName fnres = formatString("resolution_dir_%i.vol", dir+1);
		saveImg.write(fnres);
		saveImg.clear();
		#endif
		pOutputResolution.clear();
		list.clear();

//		exit(0);

		std::cout << "----------------direction-finished----------------" << std::endl;
	}

	Matrix2D<double> InertiaMatrix;
	InertiaMatrix.initZeros(3,3);

	//MultidimArray<double> &pAvgResolution = AvgResolution();

	double lambda_1, lambda_2, lambda_3, doa;
	double direction_x, direction_y, direction_z;
	int counter = 0;
	Matrix2D<double> eigenvectors;
	Matrix1D<double> eigenvalues, r0_1(3), rF_1(3), r0_2(3), rF_2(3), r0_3(3), rF_3(3), r(3);
	MultidimArray<int> arrows;
	MultidimArray<double> doaVol;
	doaVol.initZeros(arrows);
	arrows.initZeros(mask());
	const int gridStep=10;
	size_t n=0;
	int maskPos=0;

	std::cout << "Antes del FOR ALL DIRECT ELEMENTS" << std::endl;
//	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(arrows)
//	{
	for(int k=0; k<ZSIZE(arrows); ++k)
	{
		for(int i=0; i<YSIZE(arrows); ++i)
		{
			for(int j=0; j<XSIZE(arrows); ++j)
			{
				std::cout << "entro en la mascara" << std::endl;
		if (DIRECT_MULTIDIM_ELEM(mask(),n) == 1 )
		{
			std::cout << "entro en la mascara IF" << std::endl;
			double val = 1/MAT_ELEM(inertiaMatrixVariable, 6, maskPos);

			MAT_ELEM(InertiaMatrix, 0, 0) = MAT_ELEM(inertiaMatrixVariable, 0, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 0, 1) = MAT_ELEM(inertiaMatrixVariable, 1, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 0, 2) = MAT_ELEM(inertiaMatrixVariable, 2, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 1, 1) = MAT_ELEM(inertiaMatrixVariable, 3, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 1, 0) = MAT_ELEM(inertiaMatrixVariable, 1, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 1, 2) = MAT_ELEM(inertiaMatrixVariable, 4, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 2, 0) = MAT_ELEM(inertiaMatrixVariable, 2, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 2, 1) = MAT_ELEM(inertiaMatrixVariable, 4, maskPos)*val;
			MAT_ELEM(InertiaMatrix, 2, 2) = MAT_ELEM(inertiaMatrixVariable, 5, maskPos)*val;

			diagSymMatrix3x3(InertiaMatrix, eigenvalues, eigenvectors);

			degreeOfAnisotropy(eigenvalues, eigenvectors, doa,
					direction_x, direction_y, direction_z, counter);

			DIRECT_MULTIDIM_ELEM(doaVol,n) = doa;

			//lambda_1 is assumed as the least eigenvalue
			if ( (i%gridStep==0) && (j%gridStep==0) && (k%gridStep==0) )
			{
//				std::cout << "----------------------------------" << std::endl;
//				std::cout << "eigenvalues = " << eigenvalues << std::endl;
//				std::cout << "eigenvectors = " << eigenvectors << std::endl;

				double lambda_1 = VEC_ELEM(eigenvalues, 0);
				double lambda_2 = VEC_ELEM(eigenvalues, 1);
				double lambda_3 = VEC_ELEM(eigenvalues, 2);

//				std::cout << "LAMBDA_1 = " << lambda_1 <<
//						   "  LAMBDA_2 = " << lambda_2 <<
//						   "  LAMBDA_3 = " << lambda_3 << std::endl;

				defineDirection(r0_1, rF_1, eigenvectors, lambda_1, lambda_3, 0, k, i, j);
				defineDirection(r0_2, rF_2, eigenvectors, lambda_2, lambda_3, 1, k, i, j);
				defineDirection(r0_3, rF_3, eigenvectors, lambda_3, lambda_3, 2, k, i, j);

				for (double t=0; t<1; t+=0.02)
				{
					defineSegment(r0_1, rF_1, arrows, t);
					defineSegment(r0_2, rF_2, arrows, t);
					defineSegment(r0_3, rF_3, arrows, t);
				}
			}
			++maskPos;
		}
		++n;
			}
		}
	}
//	}
	Image<int> preferreddir;
	preferreddir()=arrows;
	preferreddir.write(fnDirections);
	Image<double> imgdoa;
	imgdoa() = doaVol;
	imgdoa.write(fnDoA);
}
