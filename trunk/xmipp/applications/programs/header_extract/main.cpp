/***************************************************************************
 *
 * Authors:    Sjors Scheres
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
//#define METADATA
#ifndef METADATA

#include <data/args.h>
#include <data/image.h>
#include <data/selfile.h>
#include <data/docfile.h>

void Usage();

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    bool            round_shifts = false;
    float           xx, yy;
    FileName        fn_img, fn_out;
    SelFile         SF;
    DocFile         DF;
    ImageXmipp      img;
    headerXmipp     head;

// Check command line options ===========================================
    try
    {

        SF.read(getParameter(argc, argv, "-i"));
        fn_out = getParameter(argc, argv, "-o");
        round_shifts = checkParameter(argc, argv, "-round_shifts");

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
    }

// Extracting information  ==================================================
    try
    {

        DF.reserve(SF.ImgNo());
        Matrix1D<double> docline;
        DF.append_comment("Headerinfo columns: rot (1) , tilt (2), psi (3), Xoff (4), Yoff (5), Weight (6), Flip (7)");

        docline.initZeros(7);
        SF.go_beginning();
        while (!SF.eof())
        {
            fn_img = SF.NextImg();
            if (fn_img=="") break;
            head.read(fn_img);
            head.get_originOffsets(xx, yy);
            if (round_shifts)
            {
                xx = (float)ROUND(xx);
                yy = (float)ROUND(yy);
            }
            docline(0) = head.Phi();
            docline(1) = head.Theta();
            docline(2) = head.Psi();
            docline(3) = xx;
            docline(4) = yy;
            docline(5) = head.Weight();
            docline(6) = head.Flip();
            DF.append_comment(fn_img);
            DF.append_data_line(docline);
        }
        DF.write(fn_out);
        std::cerr << " done!" << std::endl;
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }

}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    printf("Purpose:\n");
    printf(" Extracts the geometric transformation (angles & shifts) in the header of 2D-images.\n");
    printf("Usage:\n");
    printf("   header_extract \n");
    printf("        -i <selfile>       : input selfile\n");
    printf("        -o <docfile>       : output document file\n");
    printf("       [-round_shifts]     : Round shifts to integers \n");
    exit(1);
}
#else

#include <data/args.h>
#include <data/image.h>
#include <data/metadata.h>

void Usage();

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    bool            round_shifts = false;
    float           xx, yy;
    FileName        fn_img, fn_out,fn_in;
    ImageXmipp      img;
    headerXmipp     head;
    FileName        inputFile;
    

// Check command line options ===========================================
    try
    {
        
        fn_in  = getParameter(argc, argv, "-i"); 
	    if (checkParameter(argc, argv, "-o"))
            fn_out = getParameter(argc, argv, "-o");
	    else
	        fn_out = fn_in;   
            round_shifts = checkParameter(argc, argv, "-round_shifts");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
    }
    MetaData SF(fn_in);
    SF.removeObjects( MDL_ENABLED, -1 );
    // Extracting information  ==================================================
    try
    {

	    long int ret=SF.firstObject();
	    if(ret==MetaData::NO_OBJECTS_STORED)
	    {
	        std::cerr << "Empty inputFile File\n";
	        exit(1);
	    }
	    do
        {
	        fn_img = SF.image();
            if (fn_img=="") break;
            head.read(fn_img);
            head.get_originOffsets( xx, yy );
            if (round_shifts)
            {
                xx = (float)ROUND(xx);
                yy = (float)ROUND(yy);
            }
	        SF.setAngleRot( head.Phi());
	        SF.setAngleTilt( head.Theta());
    	    SF.setAnglePsi( head.Psi());
	        SF.setShiftX( xx );
	        SF.setShiftY( yy );
	        SF.setWeight( head.Weight());
	        SF.setFlip( head.Flip());
        }
        while (SF.nextObject()!= MetaData::NO_MORE_OBJECTS);
	
        SF.write(fn_out);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }

}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cout << " Purpose:\n";
    std::cout << " Extracts the geometric transformation (angles & shifts) in the header of 2D-images.\n";
    std::cout << " Usage:\n";
    std::cout << "    header_extract \n";
    std::cout << "         -i <selfile>       : input selfile\n";
    std::cout << (std::string)"       [-o <docfile> ]     : output metaData file, by default data\n" +
                                         "\t\t is stored in input metaData file\n";
    std::cout << "        [-round_shifts]     : Round shifts to integers \n";
    exit(1);
}

#endif
