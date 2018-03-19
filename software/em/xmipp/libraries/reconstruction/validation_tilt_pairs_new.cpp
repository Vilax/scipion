/***************************************************************************
 *
 * Authors:    Jose Luis Vilas          (jlvilas@cnb.csic.es)
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
#include "validation_tilt_pairs_new.h"
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <complex>
//#include <cmath>

/*
 * xmipp_validation_tilt_pairs --tilt /home/vilas/ScipionUserData/projects/rct/Runs/001623_XmippProtValidateTilt/extra/tilted/angles_iter001_00.xmd --untilt /home/vilas/ScipionUserData/projects/rct/Runs/001623_XmippProtValidateTilt/extra/untilted/angles_iter001_00.xmd -o caca
 * */

//Define Program parameters
//void ProgValidationTiltPairs::defineParams()
//{
//    //Usage
//    addUsageLine("Takes two coordinates sets and defines the coordinate transformation between them");
//	addUsageLine("First set defines the untilted coordinates, second set defines the tilted coordinates");
//	addParamsLine(" --tilt <metadata> : Metadata with angular assignment for the tilted images");
//	addParamsLine(" --untilt <metadata> : Metadata with angular assignment for the untilted images");
//	addParamsLine(" -o <metadata> : Metadata with matrix transformation");
//}
//
void ProgValidationTiltPairs::defineParams()
{
    //Usage
    addUsageLine("Takes two coordinates sets and defines the coordinate transformation between them");
    addUsageLine("First set defines the untilted coordinates, second set defines the tilted coordinates");
    addParamsLine(" --tilt <md_file=\"\"> : Metadata with angular assignment for the tilted images");
    addParamsLine(" --untilt <md_file=\"\"> : Metadata with angular assignment for the untilted images");
    addParamsLine(" -o <md_file=\"\"> : Metadata with matrix transformation");
    addParamsLine("  [--vol <img_file=\"\">]    : Input reference volume");
    addParamsLine("  [--angular_sampling <s=5>]   : Angular sampling rate in degrees. This sampling "
		    "represents the accuracy for assigning directions");
    addParamsLine("  [--maxshift <s=10>]   : Maximum shift for aligning images (in pixels)");
}


//Read params

void ProgValidationTiltPairs::readParams()
{
    fntilt = getParam("--tilt");  //Set of tilted coordinates
    fnuntilt = getParam("--untilt");
    fnOut = getParam("-o");  //Output file
    fnVol = getParam("--vol");
    smprt = getDoubleParam("--angular_sampling");
    maxshift = getIntParam("--maxshift");
}

void ProgValidationTiltPairs::generateProjections(const FileName &fnVol, double smprt)
{
    //TODO	: Check the behaviour of this function
    FileName fnGallery, fnGalleryMetaData;

    // Generate projections
    fnGallery=formatString("%s/gallery.stk",fnOut.c_str());
    String args=formatString("-i %s -o %s --sampling_rate %f",
		    fnVol.c_str(),fnGallery.c_str(),smprt);
		    //We are considering the psi sampling = angular sampling rate

    std::cout << args << std::endl;
    String cmd=(String)"xmipp_angular_project_library " + args;
    system(cmd.c_str());
}

void ProgValidationTiltPairs::assignParticle(Image<double> ImgUn_exp,
											const Matrix2D<double> &angles_rot_tilt,
											const MultidimArray<double> allGalleryProjection,
											Matrix2D<double> &outAssignment, int Nbest)
{
	MultidimArray<double> imgGallery;
	double corr, bestcorr1;
	Matrix2D<double> transformation_matrix, bestAngularAsignMatrix;
	std::vector<double> rotVector, tiltVector, psiVector, XshiftVector, YshiftVector, corrVector;

	size_t len_p;
	len_p = NSIZE(allGalleryProjection);

	std::cout << "NSIZE(allGalleryProjection) = " << len_p << std::endl;

	int count = 0;


	for (int j = 0; j<len_p; ++j)
	{
		//Untilt assignment
		imgGallery.aliasImageInStack(allGalleryProjection,j);
		MultidimArray<double> imgGallery_orig = imgGallery;
		imgGallery_orig.setXmippOrigin();

		//////////////////////////////////////
		//CORRELATION UNTILT AND PROJECTIONS//
		corr = alignImages(ImgUn_exp(), imgGallery_orig, transformation_matrix, true);

		if ((fabs(MAT_ELEM(transformation_matrix, 0, 2)) > maxshift) || (fabs(MAT_ELEM(transformation_matrix, 1, 2)) > maxshift))
			continue;

		if ((corr <0.7*bestcorr1) || (corr<0))
			continue;

		/////////////////////
		//UNTILT ASSIGNMENT//
		double psi = acos( 0.5*( MAT_ELEM(transformation_matrix,0,0) + MAT_ELEM(transformation_matrix,1,1) ) )*180/PI;

		/////////////////////////////////
		//OUTPUT ANGLES AND ASSIGNMENTS//
		if (corr>bestcorr1)
		{
			rotVector.push_back( MAT_ELEM(angles_rot_tilt, j, 0) );
			tiltVector.push_back( MAT_ELEM(angles_rot_tilt, j, 1) );
			psiVector.push_back(psi);
			XshiftVector.push_back(-MAT_ELEM(transformation_matrix,0,2));
			YshiftVector.push_back(-MAT_ELEM(transformation_matrix,1,2));
			corrVector.push_back(corr);
			bestcorr1 = corr;
			++count;
		}
	}

	outAssignment.initZeros(Nbest, 6);

	size_t len;
	len = rotVector.size();

	for (size_t k=0; k<Nbest; ++k)
	{
		MAT_ELEM(outAssignment, k, 0) = rotVector[len-k];
		MAT_ELEM(outAssignment, k, 1) = tiltVector[len-k];
		MAT_ELEM(outAssignment, k, 2) = psiVector[len-k];
		MAT_ELEM(outAssignment, k, 3) = XshiftVector[len-k];
		MAT_ELEM(outAssignment, k, 4) = YshiftVector[len-k];
		MAT_ELEM(outAssignment, k, 5) = corrVector[len-k];
	}


}

void ProgValidationTiltPairs::run()
{
	std::cout << "Starting..." << std::endl;

	// Generating projection from the volume, fnVol, with an angular sampling rate of smprt degrees
	std::cout << "Generating projections" << std::endl;
	generateProjections(fnVol, smprt);

	//Reading particles from the untilted stack and projections
	MetaData mduntilt_exp, mdtilt_exp, mdproj;
	FileName fnprojection = fnOut+"/gallery.doc";

	mduntilt_exp.read(fnuntilt);
	mdtilt_exp.read(fntilt);
	mdproj.read(fnprojection);


	MultidimArray<double> allGalleryProjection;
	FileName fnprojection_stk = fnOut+"/gallery.stk";
	Image<double> imgstack;
	imgstack.read(fnprojection_stk);
	allGalleryProjection = imgstack();
	std::cout << "Data and projections were loaded"<< std::endl;

	//Reading untilted, tilted particles and the projections of the volume,
	//then they are stored them into string vectors.
	//Thus we are avoiding accessing once and again to metadata
	std::vector<std::string> Untilted_filenames, Tilted_filenames, proj_filenames;
	FileName fnuntilt_exp, fntilt_exp, fnproj;
	double rot, tilt;

	Matrix2D<double> angles_rot_tilt;
	angles_rot_tilt.initZeros(mdproj.size(),2); //first column ->rot, second column->tilt

	FOR_ALL_OBJECTS_IN_METADATA(mduntilt_exp)
	{
		mduntilt_exp.getValue(MDL_IMAGE, fnuntilt_exp, __iter.objId);
		//mduntilt_exp.getValue(MDL_ITEM_ID, fnuntilt_exp, __iter.objId);
		Untilted_filenames.push_back(fnuntilt_exp);
	}
	std::cout << "Untilt MetaData read" << std::endl;

	FOR_ALL_OBJECTS_IN_METADATA(mdtilt_exp)
	{
		mdtilt_exp.getValue(MDL_IMAGE, fntilt_exp, __iter.objId);
		//mdtilt_exp.getValue(MDL_ITEM_ID, fntilt_exp, __iter.objId);
		Tilted_filenames.push_back(fntilt_exp);
	}
	std::cout << "Tilt MetaData read" << std::endl;



	FOR_ALL_OBJECTS_IN_METADATA(mdproj)
	{
		mdproj.getValue(MDL_IMAGE, fnproj, __iter.objId);
		mdproj.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
		mdproj.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
		proj_filenames.push_back(fnproj);

		MAT_ELEM(angles_rot_tilt, (__iter.objId)-1, 0) = rot;
		MAT_ELEM(angles_rot_tilt, (__iter.objId)-1, 1) = tilt;
	}
	std::cout << "Projection MetaData gallery read" << std::endl;
	size_t len_u = Untilted_filenames.size();
	size_t len_t = Tilted_filenames.size();
	size_t len_p = proj_filenames.size();

	//Cleaning memory
	mduntilt_exp.clear();
	mdtilt_exp.clear();
	mdproj.clear();
	imgstack.clear();
		//State: Finished

	Matrix2D<double> ZYZ_u, ZYZ_t, ZYZ_angles, ZYZ_theo, ZYZ_theo2;
	Matrix1D<double> Eu_dir;
	double alpha, beta, gamma;
	ZYZ_theo.initZeros(4,4);
	ZYZ_theo2.initZeros(4,4);
	ZYZ_u.initZeros(4,4);
	ZYZ_t.initZeros(4,4);
	ZYZ_angles.initZeros(4,4);

	MultidimArray<double> imgGallery;
	Image<double> ImgUn_exp, ImgTilted_exp;
	Matrix2D<double> outAssignment_untilted, outAssignment_tilted;
	int Nbest = 4;

	MetaData mddir;
	size_t objId_dir;
	for (int i = 0; i<len_u; ++i)
	{
		std::cout << "particle  " << i+1 << "/" << len_u-1 << std::endl;

		ImgUn_exp.read(Untilted_filenames[i]);	//Reading image
		ImgUn_exp().setXmippOrigin();

		assignParticle(ImgUn_exp, angles_rot_tilt, allGalleryProjection,
				outAssignment_untilted, Nbest);

		ImgUn_exp.read(Tilted_filenames[i]);	//Reading image
		ImgTilted_exp().setXmippOrigin();

		assignParticle(ImgTilted_exp, angles_rot_tilt, allGalleryProjection,
				outAssignment_tilted, Nbest);

		double corr1, corr2, corrTotal;
		double bestTotal = -1;
		size_t idxun, idxt;

		for (size_t j=0; j<Nbest; ++j)
		{
			corr1 = MAT_ELEM(outAssignment_untilted, j, 5);
			for (size_t k=0; k<Nbest; ++k)
			{
				corr2 = MAT_ELEM(outAssignment_tilted, k, 5);
				corrTotal = corr1 + corr2;
				if (corrTotal>bestTotal)
				{
					bestTotal = corrTotal;
					idxun = j;
					idxt = k;
				}
			}
		}

		Euler_angles2matrix(MAT_ELEM(outAssignment_untilted, 0, idxun),
							MAT_ELEM(outAssignment_untilted, 1, idxun),
							MAT_ELEM(outAssignment_untilted, 2, idxun), ZYZ_u);
		Euler_angles2matrix(MAT_ELEM(outAssignment_tilted, 0, idxt),
							MAT_ELEM(outAssignment_tilted, 1, idxt),
							MAT_ELEM(outAssignment_tilted, 2, idxt), ZYZ_t);

//		Euler_angles2matrix(180, 4.879613, 85.217391, ZYZ_u);
//		Euler_angles2matrix(-100.486262, 40.6494, 7.4767, ZYZ_t);
		Euler_angles2matrix(0, 40, 0, ZYZ_theo);

		ZYZ_angles = ZYZ_t* (ZYZ_u.inv());
		std::cout << "rot_u  " << MAT_ELEM(outAssignment_untilted, 0, idxun) << std::endl;
		std::cout << "tilt_u " << MAT_ELEM(outAssignment_untilted, 1, idxun) << std::endl;
		std::cout << "psi_u  " << MAT_ELEM(outAssignment_untilted, 2, idxun) << std::endl;
		std::cout << "corr_u  " << MAT_ELEM(outAssignment_untilted, 5, idxun) << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "rot_t  " << MAT_ELEM(outAssignment_tilted, 0, idxt) << std::endl;
		std::cout << "tilt_t " << MAT_ELEM(outAssignment_tilted, 1, idxt) << std::endl;
		std::cout << "psi_t  " << MAT_ELEM(outAssignment_tilted, 2, idxt) << std::endl;
		std::cout << "corr_t  " << MAT_ELEM(outAssignment_tilted, 5, idxt) << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "ZYZ_matrix = " << ZYZ_angles << std::endl;
		std::cout << "ZYZ_theoretical = " << ZYZ_theo << std::endl;

		Euler_matrix2angles(ZYZ_angles, alpha, beta, gamma);
		alpha *= -1;
		std::cout << "alpha = " << alpha << std::endl;
		std::cout << "betaa = " << beta << std::endl;
		std::cout << "gamma = " << gamma << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "        " << std::endl;


		Euler_direction(alpha, beta, gamma,  Eu_dir);
		objId_dir = mddir.addObject();
		mddir.setValue(MDL_ANGLE_ROT, alpha, objId_dir);
		mddir.setValue(MDL_ANGLE_TILT, beta, objId_dir);
		mddir.setValue(MDL_ANGLE_PSI, gamma, objId_dir);
		mddir.setValue(MDL_X, VEC_ELEM(Eu_dir,0), objId_dir);
		mddir.setValue(MDL_Y, VEC_ELEM(Eu_dir,1), objId_dir);
		mddir.setValue(MDL_Z, VEC_ELEM(Eu_dir,2), objId_dir);
	}

	mddir.write((String)"output_validation.xmd" );
}
