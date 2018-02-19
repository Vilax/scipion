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
#include "validation_tilt_pairs.h"
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


void ProgValidationTiltPairs::angularAssignment(size_t len_u, size_t len_p, const Matrix2D<double> &angles_rot_tilt,
		Matrix2D<double> &position_u_gallery_and_psi,
		const std::vector<std::string> Untilted_filenames,
		const MultidimArray<double> allGalleryProjection, const FileName particle_type)
{
	Matrix2D<double> transformation_matrix;
	position_u_gallery_and_psi.initZeros(6,len_u);

	MultidimArray<double> imgGallery;
	Image<double> ImgUn_exp;


	for (size_t k=0; k<(len_u); k++)
	{
		std::cout << "Iteration:  " << k << "/" << len_u-1 << std::endl;

		ImgUn_exp.read(Untilted_filenames[k]);	//Reading image
		ImgUn_exp().setXmippOrigin();

		double corr1 = 0, bestcorr1 = 0;

		for (size_t j =0; j<(len_p); j++)
		{
			//Untilt assignment
			imgGallery.aliasImageInStack(allGalleryProjection,j);
			MultidimArray<double> imgGallery_orig = imgGallery;
			imgGallery_orig.setXmippOrigin();

			//////////////////////////////////////////
			//CORRELATION UNTILT AND PROJECTIONS
			corr1 = alignImages(ImgUn_exp(), imgGallery_orig, transformation_matrix, true);

			if ((fabs(MAT_ELEM(transformation_matrix, 0, 2)) > maxshift) || (fabs(MAT_ELEM(transformation_matrix, 1, 2)) > maxshift))
				continue;

			if ((corr1 <0.7*bestcorr1) || (corr1<0))
				continue;

			//////////////////////////////////////////
			//UNTILT ASSIGNMENT
			double psi = acos( 0.5*( MAT_ELEM(transformation_matrix,0,0) + MAT_ELEM(transformation_matrix,1,1) ) )*180/PI;

			//////////////////////////////////////////
			//OUTPUT ANGLES AND ASSIGNMENTS
			if (corr1>bestcorr1)
			{
				bestcorr1 = corr1;
				MAT_ELEM(position_u_gallery_and_psi, 0, k) = MAT_ELEM(angles_rot_tilt, j, 0);
				MAT_ELEM(position_u_gallery_and_psi, 1, k) = MAT_ELEM(angles_rot_tilt, j, 1);
//				if (MAT_ELEM(position_u_gallery_and_psi, 1, k) < 15)
//				{
//					MAT_ELEM(position_u_gallery_and_psi, 2, k) = 0;
//				}
//				else
//				{
//					MAT_ELEM(position_u_gallery_and_psi, 2, k) = psi;
//				}
				MAT_ELEM(position_u_gallery_and_psi, 2, k) = psi;
				MAT_ELEM(position_u_gallery_and_psi, 3, k) = corr1;
				MAT_ELEM(position_u_gallery_and_psi, 4, k) = -MAT_ELEM(transformation_matrix,0,2);  // X Shift  in untilted image
				MAT_ELEM(position_u_gallery_and_psi, 5, k) = -MAT_ELEM(transformation_matrix,1,2);  // Y Shift
				#ifdef DEBUG
					std::cout << "File " << Untilted_filenames[k] << std::endl;
					std::cout << "Gallery " << j << std::endl;
					std::cout << "rot  = " << MAT_ELEM(position_u_gallery_and_psi, 0, k) << std::endl;
					std::cout << "tilt = " << MAT_ELEM(position_u_gallery_and_psi, 1, k) << std::endl;
					std::cout << "psi  = " << MAT_ELEM(position_u_gallery_and_psi, 2, k) << std::endl;
					std::cout << "corr  = " << corr1 << std::endl;
					std::cout << "---------------------" << std::endl;

					Image<double> save;
					save()=ImgUn_exp();
					save.write((String) particle_type);
					save()=imgGallery;
					save.write("Galeria.xmp");
					save()=imgGallery_orig;
					save.write("Galeria_orig.xmp");
					std::cout << "Press any key" << std::endl;
//					int c;
//					c= getchar();

				#endif
			}
		}
	}
	std::cout << "   " << std::endl;
}




void ProgValidationTiltPairs::run()
{
	std::cout << "Starting..." << std::endl;

	// Generating projection from the volume, fnVol, with an angular sampling rate of smprt degrees
	std::cout << "Generating projections" << std::endl;
	generateProjections(fnVol, smprt);
		//State: Finished

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
		//State: Finished
	std::cout << " ------------------- "<< std::endl;

	//Reading untilted, tilted particles and the projections of the volume,then they are stored them into string vectors.
	//Thus we are avoiding accessing once and again to metadata
	std::vector<std::string> Untilted_filenames, Tilted_filenames, proj_filenames;
	FileName fnuntilt_exp, fntilt_exp, fnproj;
	double rot, tilt;

	Matrix2D<double> angles_rot_tilt;
	angles_rot_tilt.initZeros(mdproj.size(),2); //first column is the rot angle and second column is the tilt angle

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

	//For each experimental untilted image, an angular assignment is performed
	std::cout << "Searching correlations" << std::endl;

	Matrix2D<double> position_u_gallery_and_psi, position_t_gallery_and_psi;

	position_u_gallery_and_psi.initZeros(6,len_u);
	position_t_gallery_and_psi.initZeros(6,len_t);

	Image<double> ImgUn_exp;

	FileName tilted_particle="Tilted_par.xmp";
	FileName untilted_particle="Untilted_par.xmp";

	//Untilt Assignment
	std::cout << "Untilt assignment" << std::endl;
	angularAssignment(len_u, len_p, angles_rot_tilt, position_u_gallery_and_psi, Untilted_filenames, allGalleryProjection, untilted_particle);

	//Tilt Assignment
	std::cout << "Tilt assignment" << std::endl;
	angularAssignment(len_u, len_p, angles_rot_tilt, position_t_gallery_and_psi, Tilted_filenames, allGalleryProjection, tilted_particle);

	std::cout << "Correlations ended " << std::endl;
		//State: Finished


	if (len_u!=len_t)
	{
		std::cerr << "ERROR: Mismatch dimensions: Number of untilted and tilted particles is not the same" << std::endl;
		exit(0);
	}

//	//Calculating micrograph tilt angles
//	std::cout << "----------------" << std::endl;
//	std::cout << "UNTILTED" << std::endl;
//	std::cout << position_u_gallery_and_psi << std::endl;
//	std::cout << "----------------" << std::endl;
//	std::cout << "TILTED" << std::endl;
//	std::cout << position_t_gallery_and_psi << std::endl;
//	std::cout << "----------------" << std::endl;

	Matrix2D<double> ZYZ_u, ZYZ_t, ZYZ_angles, ZYZ_theo, ZYZ_theo2;
	Matrix1D<double> Eu_dir;
	ZYZ_theo.initZeros(4,4);
	ZYZ_theo2.initZeros(4,4);
	ZYZ_u.initZeros(4,4);
	ZYZ_t.initZeros(4,4);
	ZYZ_angles.initZeros(4,4);
	double alpha, beta, gamma;
	MetaData mddir;
	size_t objId_dir;
	for (int j = 0; j<len_u; j++)
	{
		Euler_angles2matrix(MAT_ELEM(position_u_gallery_and_psi, 0, j),
				            MAT_ELEM(position_u_gallery_and_psi, 1, j),
				            MAT_ELEM(position_u_gallery_and_psi, 2, j), ZYZ_u);
		Euler_angles2matrix(MAT_ELEM(position_t_gallery_and_psi, 0, j),
							MAT_ELEM(position_t_gallery_and_psi, 1, j),
							MAT_ELEM(position_t_gallery_and_psi, 2, j), ZYZ_t);
//		Euler_angles2matrix(180, 4.879613, 85.217391, ZYZ_u);
//		Euler_angles2matrix(-100.486262, 40.6494, 7.4767, ZYZ_t);
		Euler_angles2matrix(0, 40, 0, ZYZ_theo);

		ZYZ_angles = ZYZ_t* (ZYZ_u.inv());


		std::cout << "rot_u  " << MAT_ELEM(position_u_gallery_and_psi, 0, j) << std::endl;
		std::cout << "tilt_u " << MAT_ELEM(position_u_gallery_and_psi, 1, j) << std::endl;
		std::cout << "psi_u  " << MAT_ELEM(position_u_gallery_and_psi, 2, j) << std::endl;
		std::cout << "corr_u  " << MAT_ELEM(position_u_gallery_and_psi, 3, j) << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "rot_t  " << MAT_ELEM(position_t_gallery_and_psi, 0, j) << std::endl;
		std::cout << "tilt_t " << MAT_ELEM(position_t_gallery_and_psi, 1, j) << std::endl;
		std::cout << "psi_t  " << MAT_ELEM(position_t_gallery_and_psi, 2, j) << std::endl;
		std::cout << "corr_t  " << MAT_ELEM(position_t_gallery_and_psi, 3, j) << std::endl;
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


	//Storing assignments into output metadata
	MetaData mduntilt_output, mdtilt_output;
	size_t objId_un, objId_t;
	for (size_t k=0; k<(len_u); k++)
	{
		objId_un = mduntilt_output.addObject();
		mduntilt_output.setValue(MDL_IMAGE, Untilted_filenames[k], objId_un);
		mduntilt_output.setValue(MDL_PARTICLE_ID, k, objId_un);
		mduntilt_output.setValue(MDL_ANGLE_ROT, MAT_ELEM(position_u_gallery_and_psi, 0, k), objId_un);
		mduntilt_output.setValue(MDL_ANGLE_TILT, MAT_ELEM(position_u_gallery_and_psi, 1, k), objId_un);
		mduntilt_output.setValue(MDL_ANGLE_PSI, MAT_ELEM(position_u_gallery_and_psi, 2, k), objId_un);
		mduntilt_output.setValue(MDL_MAXCC, MAT_ELEM(position_u_gallery_and_psi, 3, k), objId_un);
		mduntilt_output.setValue(MDL_SHIFT_X, MAT_ELEM(position_u_gallery_and_psi, 4, k), objId_un);
		mduntilt_output.setValue(MDL_SHIFT_Y, MAT_ELEM(position_u_gallery_and_psi, 5, k), objId_un);
	}
	for (size_t k=0; k<(len_t); k++)
	{
		objId_t = mdtilt_output.addObject();
		mdtilt_output.setValue(MDL_IMAGE, Tilted_filenames[k], objId_t);
		mdtilt_output.setValue(MDL_PARTICLE_ID, k, objId_t);
		mdtilt_output.setValue(MDL_ANGLE_ROT, MAT_ELEM(position_t_gallery_and_psi, 0, k), objId_t);
		mdtilt_output.setValue(MDL_ANGLE_TILT, MAT_ELEM(position_t_gallery_and_psi, 1, k), objId_t);
		mdtilt_output.setValue(MDL_ANGLE_PSI, MAT_ELEM(position_t_gallery_and_psi, 2, k), objId_t);
		mdtilt_output.setValue(MDL_MAXCC, MAT_ELEM(position_t_gallery_and_psi, 3, k), objId_t);
		mdtilt_output.setValue(MDL_SHIFT_X, MAT_ELEM(position_t_gallery_and_psi, 4, k), objId_t);
		mdtilt_output.setValue(MDL_SHIFT_Y, MAT_ELEM(position_t_gallery_and_psi, 5, k), objId_t);
	}

	mduntilt_output.write((String)"particles@"+fnOut+'/'+fnuntilt.getBaseName() +"_angular_assignment"+ ".xmd" );
	mdtilt_output.write((String)"particles@"+fnOut+'/'+fntilt.getBaseName() +"_angular_assignment"+ ".xmd" );
		//State: Finished
}


//void ProgValidationTiltPairs::readParams()
//{
//    fntiltimage_In = getParam("--tilt");  //Set of tilted coordinates
//    fnuntiltimage_In = getParam("--untilt");
//	fnOut = getParam("-o");  //Output file
//}
//
//
//void ProgValidationTiltPairs::quaternion2Paulibasis(double rot, double tilt, double psi, std::complex<double> (&L)[4])
//{
//	double cr, ct, cp, sr, st, sp;
//
//	cr = cos(rot/2);
//	ct = cos(tilt/2);
//	cp = cos(psi/2);
//	sr = sin(rot/2);
//	st = sin(tilt/2);
//	sp = sin(psi/2);
//
//	L[0] = cr*ct*cp - sr*ct*sp;
//	L[1] = 1i*(sr*st*cp - cr*st*sp);
//	L[2] = 1i*(cr*st*cp + sr*st*sp);
//	L[3] = 1i*(sr*ct*cp+cr*ct*sp);
//}
//
//
//void ProgValidationTiltPairs::matrix2Paulibasis(std::complex<double> M[4],
//		std::complex<double> (&P)[4])
//{
//	//M[0] = m11; M[1]=m12; M[2]=m21; M[3]=m22
//
//	std::complex<double> I=1i;
//	std::complex<double> aux=0.5;
//	P[0]=(M[0]+M[3])*aux;
//	P[1]=(M[1]+M[2])*aux;
//	P[2]=I*(M[1]-M[2])*aux;
//	P[3]=(M[0]-M[3])*aux;
//}
//
//void ProgValidationTiltPairs::InversefromPaulibasis(std::complex<double> Original[4],
//		std::complex<double> (&Inver)[4])
//{
//	//It takes a Pauli expression and returns its inverse expressed in Pauli basis
//
//	//TODO Raise an exception is the Original matrix does not belong to SU(2) group
//
//	std::complex<double> Inver_matrix[4], mat[4];
//	std::complex<double> NOriginal;
//	double aux=0.5;
//
//	Paulibasis2matrix(Original,mat);
//
//	Inver_matrix[0] = mat[3];
//	Inver_matrix[1] = -mat[1];
//	Inver_matrix[2] = -mat[2];
//	Inver_matrix[3] = mat[0];
//
//	matrix2Paulibasis(Inver_matrix,Inver);
//}
//
//void ProgValidationTiltPairs::inverse_matrixSU2(std::complex<double> Original[4],
//		std::complex<double> (&Inver)[4])
//{
//	//It takes a matrix and returns its inverse expressed in Pauli basis
//
//	//TODO Raise an exception is the Original matrix does not belong to SU(2) group
//
//	std::complex<double> Inver_matrix[4];
//
//	Inver_matrix[0] = Original[3];
//	Inver_matrix[1] = -Original[1];
//	Inver_matrix[2] = -Original[2];
//	Inver_matrix[3] = Original[0];
//
//	matrix2Paulibasis(Inver_matrix,Inver);
//}
//
//void ProgValidationTiltPairs::Paulibasis2matrix(std::complex<double> P[4], std::complex<double> (&M)[4])
//{
//	std::complex<double> I=1i;
//	M[0] = (P[0]+P[3]);
//	M[1] = (P[1]-I*P[2]);
//	M[2] = (P[1]+I*P[2]);
//	M[3] = (P[0]-P[3]);
//}
//
//void ProgValidationTiltPairs::Pauliproduct(std::complex<double> A[4], std::complex<double> B[4],
//		std::complex<double> (&P)[4])
//{
//	std::complex<double> A_matrix[4], B_matrix[4], aux[4];
//
//	Paulibasis2matrix(A,A_matrix);
//	Paulibasis2matrix(B,B_matrix);
//
//	aux[0] = A_matrix[0]*B_matrix[0] + A_matrix[1]*B_matrix[2];
//	aux[2] = A_matrix[2]*B_matrix[0] + A_matrix[3]*B_matrix[2];
//	aux[1] = A_matrix[0]*B_matrix[1] + A_matrix[1]*B_matrix[3];
//	aux[3] = A_matrix[2]*B_matrix[1] + A_matrix[3]*B_matrix[3];
//
//	matrix2Paulibasis(aux, P);
//}
//
//void ProgValidationTiltPairs::extrarotationangles(std::complex<double> R[4], double &alpha_x, double &alpha_y)
//{
//	std::complex<double> I=1i;
//	std::complex<double> aux1 = I*R[1]/R[0],  aux2 = I*R[2]/R[0], alpha_aux_x, alpha_aux_y;
//
//	if ((aux1.imag() == 0) && (aux2.imag() == 0))
//	{
//		alpha_aux_x = 2*atan(aux1.real());
//		alpha_aux_y = 2*atan(aux2.real());
//		alpha_x = alpha_aux_x.real()*180/PI;
//		alpha_y = alpha_aux_y.real()*180/PI;
//	}
//	else
//	{
//		std::cout << "Error, argument of atan is complex" << std::endl;
//	}
//}
//
//
//void ProgValidationTiltPairs::angles2tranformation(double untilt_angles[3],
//		double tilt_angles[3], double alpha_x, double alpha_y)
//{
//	double rotu = untilt_angles[0], tiltu = untilt_angles[1], psiu = untilt_angles[2];
//	double rott = tilt_angles[0], tiltt = tilt_angles[1], psit = tilt_angles[2];
//	std::complex<double> qu[4], M[4], Inv_qu[4], qt[4], R[4];
//
//	quaternion2Paulibasis(rotu, tiltu, psiu, qu);
//    Paulibasis2matrix(qu,M);
//    inverse_matrixSU2(M, Inv_qu);
//    quaternion2Paulibasis(rott, tiltt, psit, qt);
//    Pauliproduct(qt, Inv_qu, R);
//
//    extrarotationangles(R, alpha_x, alpha_y);
//    std::cout << "alpha_x = " << alpha_x << std::endl;
//    std::cout << "alpha_y = " << alpha_y << std::endl;
//}
//
//void ProgValidationTiltPairs::run()
//{
//	MetaData MD_tilted, MD_untilted, DF1sorted, DF2sorted, DFweights;
//
//	MD_tilted.read(fntiltimage_In);
//	MD_untilted.read(fnuntiltimage_In);
//
//	DF1sorted.sort(MD_tilted,MDL_ITEM_ID,true);
//	DF2sorted.sort(MD_untilted,MDL_ITEM_ID,true);
//
//	MDIterator iter1(DF1sorted), iter2(DF2sorted);
//	std::vector< Matrix1D<double> > ang1, ang2;
//	Matrix1D<double> rotTiltPsi(3), z(3);
//	size_t currentId;
//	bool anotherImageIn2=iter2.hasNext();
//	size_t id1, id2;
//	bool mirror;
//	Matrix2D<double> Eu, Et, R;
//	double alpha, beta;
//	while (anotherImageIn2)
//	{
//		ang1.clear();
//		ang2.clear();
//
//		// Take current id
//		DF2sorted.getValue(MDL_ITEM_ID,currentId,iter2.objId);
//
//		// Grab all the angles in DF2 associated to this id
//		bool anotherIteration=false;
//		do
//		{
//			DF2sorted.getValue(MDL_ITEM_ID,id2,iter2.objId);
//			anotherIteration=false;
//			if (id2==currentId)
//			{
//				DF2sorted.getValue(MDL_ANGLE_ROT,XX(rotTiltPsi),iter2.objId);
//				DF2sorted.getValue(MDL_ANGLE_TILT,YY(rotTiltPsi),iter2.objId);
//				DF2sorted.getValue(MDL_ANGLE_PSI,ZZ(rotTiltPsi),iter2.objId);
//				DF2sorted.getValue(MDL_FLIP,mirror,iter2.objId);
//				std::cout << "From DF2:" << XX(rotTiltPsi) << " " << YY(rotTiltPsi) << " " << ZZ(rotTiltPsi) << " " << mirror << std::endl;
//				//LINEA ANTERIOR ORIGINAL
//				if (mirror)
//				{
//					double rotp, tiltp, psip;
//					Euler_mirrorX(XX(rotTiltPsi),YY(rotTiltPsi),ZZ(rotTiltPsi), rotp, tiltp, psip);
//					XX(rotTiltPsi)=rotp;
//					YY(rotTiltPsi)=tiltp;
//					ZZ(rotTiltPsi)=psip;
//				}
//				ang2.push_back(rotTiltPsi);
//				iter2.moveNext();
//				if (iter2.hasNext())
//					anotherIteration=true;
//			}
//		} while (anotherIteration);
//
//		// Advance Iter 1 to catch Iter 2
//		double N=0, cumulatedDistance=0;
//		size_t newObjId=0;
//		if (iter1.objId>0)
//		{
//			DF1sorted.getValue(MDL_ITEM_ID,id1,iter1.objId);
//			while (id1<currentId && iter1.hasNext())
//			{
//				iter1.moveNext();
//				DF1sorted.getValue(MDL_ITEM_ID,id1,iter1.objId);
//			}
//
//			// If we are at the end of DF1, then we did not find id1 such that id1==currentId
//			if (!iter1.hasNext())
//				break;
//
//			// Grab all the angles in DF1 associated to this id
//			anotherIteration=false;
//			do
//			{
//				DF1sorted.getValue(MDL_ITEM_ID,id1,iter1.objId);
//				anotherIteration=false;
//				if (id1==currentId)
//				{
//					DF1sorted.getValue(MDL_ANGLE_ROT,XX(rotTiltPsi),iter1.objId);
//					DF1sorted.getValue(MDL_ANGLE_TILT,YY(rotTiltPsi),iter1.objId);
//					DF1sorted.getValue(MDL_ANGLE_PSI,ZZ(rotTiltPsi),iter1.objId);
//					DF1sorted.getValue(MDL_FLIP,mirror,iter1.objId);
//					std::cout << "From DF1:" << XX(rotTiltPsi) << " " << YY(rotTiltPsi) << " " << ZZ(rotTiltPsi) << " " << mirror << std::endl;
//					//LINEA ANTERIOR ORIGINAL
//					if (mirror)
//					{
//						double rotp, tiltp, psip;
//						Euler_mirrorX(XX(rotTiltPsi),YY(rotTiltPsi),ZZ(rotTiltPsi), rotp, tiltp, psip);
//						XX(rotTiltPsi)=rotp;
//						YY(rotTiltPsi)=tiltp;
//						ZZ(rotTiltPsi)=psip;
//					}
//					ang1.push_back(rotTiltPsi);
//					iter1.moveNext();
//					if (iter1.hasNext())
//						anotherIteration=true;
//				}
//			} while (anotherIteration);
//
//			// Process both sets of angles
//			for (size_t i=0; i<ang2.size(); ++i)
//			{
//				const Matrix1D<double> &anglesi=ang2[i];
//				double rotu=XX(anglesi);
//				double tiltu=YY(anglesi);
//				double psiu=ZZ(anglesi);
//				Euler_angles2matrix(rotu,tiltu,psiu,Eu,false);
//				/*std::cout << "------UNTILTED MATRIX------" << std::endl;
//				std::cout << Eu << std::endl;
//				std::cout << "vector" << std::endl;
//				std::cout << Eu(2,0) << "  " << Eu(2,1) << "  " << Eu(2,2) << std::endl;*/
//
//
//				for (size_t j=0; j<ang1.size(); ++j)
//				{
//					const Matrix1D<double> &anglesj=ang1[j];
//					double rott=XX(anglesj);
//					double tiltt=YY(anglesj);
//					double psit=ZZ(anglesj);
//					double alpha_x, alpha_y;
//					Euler_angles2matrix(rott,tiltt,psit,Et,false);
//					//////////////////////////////////////////////////////////////////
//					double untilt_angles[3]={rotu, tiltu, psiu}, tilt_angles[3]={rott, tiltt, psit};
//					angles2tranformation(untilt_angles, tilt_angles, alpha_x, alpha_y);
//					//std::cout << "alpha = " << (alpha_x*alpha_x+alpha_y*alpha_y) << std::endl;
//					//////////////////////////////////////////////////////////////////
//					/*std::cout << "------TILTED MATRIX------" << std::endl;
//					std::cout << Et << std::endl;
//					std::cout << "vector" << std::endl;
//					std::cout << Et(2,0) << "  " << Et(2,1) << "  " << Et(2,2) << std::endl;
//					std::cout << "---------------------------" << std::endl;
//					std::cout << "---------------------------" << std::endl;*/
//					R=Eu*Et.transpose();
//					double rotTransf, tiltTransf, psiTransf;
//					Euler_matrix2angles(R, rotTransf, tiltTransf, psiTransf);
//					std::cout << "Rot_and_Tilt " << rotTransf << " " << tiltTransf << std::endl;
//					//LINEA ANTERIOR ORIGINAL
//
//				XX(z) = Eu(2,0) - Et(2,0);
//				YY(z) = Eu(2,1) - Et(2,1);
//				ZZ(z) = Eu(2,2) - Et(2,2);
//
//				alpha = atan2(YY(z), XX(z));        //Expressed in rad
//				beta = atan2(XX(z)/cos(alpha), ZZ(z));   //Expressed in rad
//				std::cout << "alpha = " << alpha*180/PI << std::endl;
//				std::cout << "beta = " << beta*180/PI << std::endl;
//				}
//			}
//		}
//		else
//			N=0;
//
//		if (N>0)
//		{
//			double meanDistance=cumulatedDistance/ang2.size();
//			DFweights.setValue(MDL_ANGLE_DIFF,meanDistance,newObjId);
//		}
//		else
//			if (newObjId>0)
//				DFweights.setValue(MDL_ANGLE_DIFF,-1.0,newObjId);
//		anotherImageIn2=iter2.hasNext();
//	}
//
//	std::complex<double> qu[4], qt[4], M[4], Inv_qu[4], test[4], P1[4], P2[4], Inv_quu[4];
//	double rotu=34*PI/180, tiltu=10*PI/180, psiu=5*PI/180;
//	double rott=25*PI/180, tiltt=15*PI/180, psit=40*PI/180;
//
//
//    quaternion2Paulibasis(rotu, tiltu, psiu, qu);
//    /*std::cout << "quaternion2Pauli" << std::endl;
//    std::cout << "Untilted " << qu[0] << " " << qu[1] << " " << qu[2] << " " << qu[3] << std::endl;
//    std::cout << "      " << std::endl;*/
//
//    Paulibasis2matrix(qu,M);
//    /*std::cout << "Pauli2matrix" << std::endl;
//    std::cout << "Matriz   " << M[0] << " " << M[1] << " " << M[2] << " " << M[3] << std::endl;
//    std::cout << "      " << std::endl;*/
//
//    inverse_matrixSU2(M, Inv_qu);
//    /*std::cout << "inverse_matrixSU(2)" << std::endl;
//    std::cout << "Inversa  " << Inv_qu[0] << " " << Inv_qu[1] << " " << Inv_qu[2] << " " << Inv_qu[3] << std::endl;
//    std::cout << "      " << std::endl;*/
//
//    quaternion2Paulibasis(rott, tiltt, psit, qt);
//    /*std::cout << "quaternion2Pauli" << std::endl;
//    std::cout << "Tilted " << qt[0] << " " << qt[1] << " " << qt[2] << " " << qt[3] << std::endl;
//    std::cout << "      " << std::endl;*/
//
//    InversefromPaulibasis(qu,Inv_quu);
//
//    Pauliproduct(qt, Inv_qu, P1);
//    /*std::cout << "Pauliproduct" << std::endl;
//    std::cout << "quaternion qt  " << P1[0] << " " << P1[1] << " " << P1[2] << " " << P1[3] << std::endl;
//    std::cout << "      " << std::endl;
//    std::cout << "-----------------------------------" << std::endl;*/
//
//    //double alpha_x, alpha_y;
//    //extrarotationangles(P1, alpha_x, alpha_y);
//    //std::cout << "alpha_x = " << alpha_x << " " << "alpha_y = " << alpha_y << std::endl;
//}
//


