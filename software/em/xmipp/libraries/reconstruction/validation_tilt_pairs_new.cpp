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

void ProgValidationTiltPairs::defineParams()
{
    //Usage
    addUsageLine("Takes two coordinates sets and defines the coordinate transformation between them");
    addUsageLine("First set defines the untilted coordinates, second set defines the tilted coordinates");
    addParamsLine(" --tilt <md_file=\"\">   : Metadata with angular assignment for the tilted images");
    addParamsLine(" --untilt <md_file=\"\"> : Metadata with angular assignment for the untilted images");
    addParamsLine(" -o <md_file=\"\"> : Metadata with matrix transformation");
    addParamsLine("  [--vol <img_file=\"\">]    : Input reference volume");
    addParamsLine("  [--angular_sampling <s=5>]   : Angular sampling rate in degrees. This sampling "
		    "represents the accuracy for assigning directions");
    addParamsLine("  [--maxshift <s=10>]   : Maximum shift for aligning images (in pixels)");
    addParamsLine("  [--nbest <s=4>]   : Number of best correlations to be considered in the alignment");
}

void ProgValidationTiltPairs::readParams()
{
    fntilt = getParam("--tilt");  //Set of tilted coordinates
    fnuntilt = getParam("--untilt");
    fnOut = getParam("-o");  //Output file
    fnVol = getParam("--vol");
    smprt = getDoubleParam("--angular_sampling");
    maxshift = getIntParam("--maxshift");
    Nbest = getIntParam("--nbest");
}

void ProgValidationTiltPairs::generateProjections(const FileName &fnVol, double smprt)
{
    //TODO	: Check the behaviour of this function
    FileName fnGallery, fnGalleryMetaData;

    // Generate projections
    fnGallery=formatString("gallery.stk",fnOut.c_str());
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
		corr = alignImages(ImgUn_exp(), imgGallery_orig, transformation_matrix, false);

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
	std::cout << " len " << len <<  std::endl;
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

	// Generating projection from the volume,
	//fnVol, with an angular sampling rate of smprt degrees
	std::cout << "Generating projections" << std::endl;
	generateProjections(fnVol, smprt);

	//Reading particles from the untilted stack and projections
	MetaData mduntilt_exp, mdtilt_exp, mdproj;
	FileName fnprojection = "gallery.doc", fnprojection_stk = "gallery.stk";
	MultidimArray<double> allGalleryProjection;

	mduntilt_exp.read(fnuntilt);
	mdtilt_exp.read(fntilt);
	mdproj.read(fnprojection);

	Image<double> imgstack;
	imgstack.read(fnprojection_stk);
	allGalleryProjection = imgstack();
	std::cout << "Data and projections were loaded"<< std::endl;

	//Reading untilted, tilted particles and the projections of the volume,
	//then they are stored them into string vectors.
	//Thus we are avoiding accessing once and again to metadata
	std::vector<std::string> Untilted_filenames, Tilted_filenames, proj_filenames;
	FileName fnuntilt_exp, fntilt_exp, fnproj;
	double rot, tilt, psi, shiftX, shiftY;

	Matrix2D<double> angles_rot_tilt;
	angles_rot_tilt.initZeros(mdproj.size(),2); //first column ->rot, second column->tilt


	Matrix2D<double> transformation_matrix;
	transformation_matrix.initZeros(3,3);
	MAT_ELEM(transformation_matrix,2,2) = 1;
	Image<double> ImgUntilted, ImgTilted;
	MultidimArray<double> corrected, avgUntilted, avgTilted;
	mduntilt_exp.getValue(MDL_IMAGE, fnuntilt_exp, 1);
	ImgUntilted.read(fnuntilt_exp);
	avgUntilted.initZeros(ImgUntilted());
	ImgTilted.read(fnuntilt_exp);
	avgTilted.initZeros(ImgTilted());
//	std::cout << "matrix = " << transformation_matrix << std::endl;

//	Image<double> sav;

	int count = 0;
	FOR_ALL_OBJECTS_IN_METADATA(mduntilt_exp)
	{
		mduntilt_exp.getValue(MDL_IMAGE, fnuntilt_exp, __iter.objId);
		mduntilt_exp.getValue(MDL_SHIFT_X, shiftX, __iter.objId);
		mduntilt_exp.getValue(MDL_SHIFT_Y, shiftY, __iter.objId);
		mduntilt_exp.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
		MAT_ELEM(transformation_matrix,0,0) = cos(psi*PI/180);
		MAT_ELEM(transformation_matrix,1,1) = MAT_ELEM(transformation_matrix,0,0);
		MAT_ELEM(transformation_matrix,0,1) = sin(psi*PI/180);
		MAT_ELEM(transformation_matrix,1,0) = -sin(psi*PI/180);
		MAT_ELEM(transformation_matrix,0,2) = shiftX;
		MAT_ELEM(transformation_matrix,1,2) = shiftY;
		ImgUntilted.read(fnuntilt_exp);
//		ImgTilted.read(fnuntilt_exp);
		applyGeometry(LINEAR, corrected, ImgUntilted(), transformation_matrix, IS_NOT_INV, true);
		avgUntilted += corrected;
		MAT_ELEM(transformation_matrix,0,2) = 0;
		MAT_ELEM(transformation_matrix,1,2) = 0;
		mdtilt_exp.getValue(MDL_IMAGE, fntilt_exp, __iter.objId);
		ImgTilted.read(fntilt_exp);
		applyGeometry(LINEAR, corrected, ImgTilted(), transformation_matrix, IS_NOT_INV, true);
		avgTilted += corrected;
		Untilted_filenames.push_back(fnuntilt_exp);
		++count;
	}
	avgUntilted = avgUntilted/((double) count);
	ImgUntilted() = avgUntilted;;
	ImgUntilted.write("avgUntilted.xmp");
	ImgUntilted() = avgTilted;;
	ImgUntilted.write("avgTilted.xmp");

	std::cout << "Class average computed" << std::endl;


	FOR_ALL_OBJECTS_IN_METADATA(mdtilt_exp)
	{
		mdtilt_exp.getValue(MDL_IMAGE, fntilt_exp, __iter.objId);
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

	MultidimArray<double> imgGallery, aux_avg;
	Image<double> ImgUn_exp, ImgTilted_exp;
	Matrix2D<double> outAssignment_untilted, outAssignment_tilted;
	int Nbest = 4;

	MetaData mddir, md_U_assignment, md_T_assignment;
	size_t objId_dir, objId_U_assignment, objId_T_assignment;
	////////new////////
	FileName fnAvgUntilted = "avgUntilted.xmp", fnAvgTilted = "avgTilted.xmp";
	ImgUn_exp.read(fnAvgUntilted);	//Reading image
	ImgUn_exp().setXmippOrigin();
	std::cout << "particle untilted" << std::endl;
	assignParticle(ImgUn_exp, angles_rot_tilt, allGalleryProjection,
					outAssignment_untilted, Nbest);

/*
	ImgTilted_exp.read(fnAvgTilted);	//Reading image
	ImgTilted_exp().setXmippOrigin();
	std::cout << "particle untilted" << std::endl;
	assignParticle(ImgTilted_exp, angles_rot_tilt, allGalleryProjection,
					outAssignment_tilted, Nbest);
*/

	MetaData md_outAssignment_untilted, md_outAssignment_tilted;
	size_t objId;
	for (size_t k = 0; k<Nbest; ++k)
	{
		objId = md_outAssignment_untilted.addObject();
		md_outAssignment_untilted.setValue(MDL_IMAGE, fnAvgUntilted, objId);
		md_outAssignment_untilted.setValue(MDL_ANGLE_ROT, MAT_ELEM(outAssignment_untilted, k, 0), objId);
		md_outAssignment_untilted.setValue(MDL_ANGLE_TILT, MAT_ELEM(outAssignment_untilted, k, 1), objId);
		md_outAssignment_untilted.setValue(MDL_ANGLE_PSI, MAT_ELEM(outAssignment_untilted, k, 2), objId);
		md_outAssignment_untilted.setValue(MDL_SHIFT_X, MAT_ELEM(outAssignment_untilted, k, 3), objId);
		md_outAssignment_untilted.setValue(MDL_SHIFT_Y, MAT_ELEM(outAssignment_untilted, k, 4), objId);
/*
		md_outAssignment_tilted.setValue(MDL_IMAGE, fnAvgTilted, objId);
		md_outAssignment_tilted.setValue(MDL_ANGLE_ROT, MAT_ELEM(outAssignment_tilted, k, 0), objId);
		md_outAssignment_tilted.setValue(MDL_ANGLE_TILT, MAT_ELEM(outAssignment_tilted, k, 1), objId);
		md_outAssignment_tilted.setValue(MDL_ANGLE_PSI, MAT_ELEM(outAssignment_tilted, k, 2), objId);
		md_outAssignment_tilted.setValue(MDL_SHIFT_X, MAT_ELEM(outAssignment_tilted, k, 3), objId);
		md_outAssignment_tilted.setValue(MDL_SHIFT_Y, MAT_ELEM(outAssignment_tilted, k, 4), objId);
*/
	}

	md_outAssignment_untilted.write((String)"myuntiled_direction.xmd" );
	md_outAssignment_tilted.write((String)"mytiled_direction.xmd" );




	for (int i = 0; i<len_t; ++i)
	{
		std::cout << "particle  " << i+1 << "/" << len_t-1 << std::endl;

		ImgTilted_exp.read(Tilted_filenames[i]);	//Reading image
		ImgTilted_exp().setXmippOrigin();
		std::cout << "particle tilted" << std::endl;

		std::cout << " Tilted_filenames[i] " << Tilted_filenames[i] << std::endl;
		Image<double> save;
		save = ImgTilted_exp;
		save.write("PPP.xmp");
		assignParticle(ImgTilted_exp, angles_rot_tilt, allGalleryProjection,
						outAssignment_tilted, Nbest);

		std::cout << "after particle" << std::endl;
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


		objId_U_assignment = md_U_assignment.addObject();
		mduntilt_exp.getValue(MDL_ANGLE_PSI, psi, objId_U_assignment);


		md_U_assignment.setValue(MDL_ANGLE_ROT, MAT_ELEM(outAssignment_untilted, 0, idxun), objId_U_assignment);
		md_U_assignment.setValue(MDL_ANGLE_TILT, MAT_ELEM(outAssignment_untilted, 1, idxun), objId_U_assignment);
		md_U_assignment.setValue(MDL_ANGLE_PSI, psi, objId_U_assignment);
		//		md_U_assignment.setValue(MDL_ANGLE_PSI, MAT_ELEM(outAssignment_untilted, 2, idxun), objId_U_assignment);

		objId_T_assignment = md_T_assignment.addObject();
		md_T_assignment.setValue(MDL_ANGLE_ROT, MAT_ELEM(outAssignment_tilted, 0, idxt), objId_T_assignment);
		md_T_assignment.setValue(MDL_ANGLE_TILT, MAT_ELEM(outAssignment_tilted, 1, idxt), objId_T_assignment);
		md_T_assignment.setValue(MDL_ANGLE_PSI, MAT_ELEM(outAssignment_tilted, 2, idxt), objId_T_assignment);


		Euler_angles2matrix(MAT_ELEM(outAssignment_untilted, 0, idxun),
							MAT_ELEM(outAssignment_untilted, 1, idxun),
							psi, ZYZ_u);
							//							MAT_ELEM(outAssignment_untilted, 2, idxun), ZYZ_u);
		Euler_angles2matrix(MAT_ELEM(outAssignment_tilted, 0, idxt),
							MAT_ELEM(outAssignment_tilted, 1, idxt),
							MAT_ELEM(outAssignment_tilted, 2, idxt), ZYZ_t);

		ZYZ_angles = ZYZ_t* (ZYZ_u.inv());

		Euler_matrix2angles(ZYZ_angles, alpha, beta, gamma);
		alpha *= -1;
		std::cout << "alpha = " << alpha << std::endl;
		std::cout << "betaa = " << beta << std::endl;
		std::cout << "gamma = " << gamma << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "        " << std::endl;

		Euler_direction(alpha, beta, gamma,  Eu_dir);
		std::cout << "x = " << VEC_ELEM(Eu_dir,0) << std::endl;
		std::cout << "y = " << VEC_ELEM(Eu_dir,1) << std::endl;
		std::cout << "z = " << VEC_ELEM(Eu_dir,2) << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "---------------------" << std::endl;
		objId_dir = mddir.addObject();
		mddir.setValue(MDL_ANGLE_ROT, alpha, objId_dir);
		mddir.setValue(MDL_ANGLE_TILT, beta, objId_dir);
		mddir.setValue(MDL_ANGLE_PSI, gamma, objId_dir);
		mddir.setValue(MDL_X, VEC_ELEM(Eu_dir,0), objId_dir);
		mddir.setValue(MDL_Y, VEC_ELEM(Eu_dir,1), objId_dir);
		mddir.setValue(MDL_Z, VEC_ELEM(Eu_dir,2), objId_dir);
	}

/*
	for (int i = 0; i<len_u; ++i)
	{
		std::cout << "particle  " << i+1 << "/" << len_u-1 << std::endl;

		ImgUn_exp.read(Untilted_filenames[i]);	//Reading image
		ImgUn_exp().setXmippOrigin();
		std::cout << "particle untilted" << std::endl;
		assignParticle(ImgUn_exp, angles_rot_tilt, allGalleryProjection,
						outAssignment_untilted, Nbest);

		ImgTilted_exp.read(Tilted_filenames[i]);	//Reading image
		ImgTilted_exp().setXmippOrigin();
		std::cout << "particle tilted" << std::endl;

		std::cout << " Tilted_filenames[i] " << Tilted_filenames[i] << std::endl;
		Image<double> save;
		save = ImgTilted_exp;
		save.write("PPP.xmp");
		assignParticle(ImgTilted_exp, angles_rot_tilt, allGalleryProjection,
						outAssignment_tilted, Nbest);

		std::cout << "after particle" << std::endl;
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

		objId_U_assignment = md_U_assignment.addObject();
		md_U_assignment.setValue(MDL_ANGLE_ROT, MAT_ELEM(outAssignment_untilted, 0, idxun), objId_U_assignment);
		md_U_assignment.setValue(MDL_ANGLE_TILT, MAT_ELEM(outAssignment_untilted, 1, idxun), objId_U_assignment);
		md_U_assignment.setValue(MDL_ANGLE_PSI, MAT_ELEM(outAssignment_untilted, 2, idxun), objId_U_assignment);

		objId_T_assignment = md_T_assignment.addObject();
		md_T_assignment.setValue(MDL_ANGLE_ROT, MAT_ELEM(outAssignment_tilted, 0, idxun), objId_T_assignment);
		md_T_assignment.setValue(MDL_ANGLE_TILT, MAT_ELEM(outAssignment_tilted, 1, idxun), objId_T_assignment);
		md_T_assignment.setValue(MDL_ANGLE_PSI, MAT_ELEM(outAssignment_tilted, 2, idxun), objId_T_assignment);


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
//		std::cout << "rot_u  " << MAT_ELEM(outAssignment_untilted, 0, idxun) << std::endl;
//		std::cout << "tilt_u " << MAT_ELEM(outAssignment_untilted, 1, idxun) << std::endl;
//		std::cout << "psi_u  " << MAT_ELEM(outAssignment_untilted, 2, idxun) << std::endl;
//		std::cout << "corr_u  " << MAT_ELEM(outAssignment_untilted, 5, idxun) << std::endl;
//		std::cout << "--------------------" << std::endl;
//		std::cout << "rot_t  " << MAT_ELEM(outAssignment_tilted, 0, idxt) << std::endl;
//		std::cout << "tilt_t " << MAT_ELEM(outAssignment_tilted, 1, idxt) << std::endl;
//		std::cout << "psi_t  " << MAT_ELEM(outAssignment_tilted, 2, idxt) << std::endl;
//		std::cout << "corr_t  " << MAT_ELEM(outAssignment_tilted, 5, idxt) << std::endl;
//		std::cout << "--------------------" << std::endl;
//		std::cout << "ZYZ_matrix = " << ZYZ_angles << std::endl;
//		std::cout << "ZYZ_theoretical = " << ZYZ_theo << std::endl;

		Euler_matrix2angles(ZYZ_angles, alpha, beta, gamma);
		alpha *= -1;
		std::cout << "alpha = " << alpha << std::endl;
		std::cout << "betaa = " << beta << std::endl;
		std::cout << "gamma = " << gamma << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "        " << std::endl;

		Euler_direction(alpha, beta, gamma,  Eu_dir);
		std::cout << "x = " << VEC_ELEM(Eu_dir,0) << std::endl;
		std::cout << "y = " << VEC_ELEM(Eu_dir,1) << std::endl;
		std::cout << "z = " << VEC_ELEM(Eu_dir,2) << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "---------------------" << std::endl;
		objId_dir = mddir.addObject();
		mddir.setValue(MDL_ANGLE_ROT, alpha, objId_dir);
		mddir.setValue(MDL_ANGLE_TILT, beta, objId_dir);
		mddir.setValue(MDL_ANGLE_PSI, gamma, objId_dir);
		mddir.setValue(MDL_X, VEC_ELEM(Eu_dir,0), objId_dir);
		mddir.setValue(MDL_Y, VEC_ELEM(Eu_dir,1), objId_dir);
		mddir.setValue(MDL_Z, VEC_ELEM(Eu_dir,2), objId_dir);
	}
*/
	mddir.write((String)"output_validation.xmd" );
}
