/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres" and "Takanori Nakane"
 * MRC Laboratory of Molecular Biology
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
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <src/projector.h>
#include <src/fftw.h>
#include <src/args.h>
#include <src/euler.h>
#include <src/transformations.h>
//#include <src/symmetries.h>
//#include <src/time.h>

class rotate_map
{
private:
	Matrix2D<RFLOAT> A3D;
	MultidimArray<Complex> F2D;
	MultidimArray<RFLOAT> dummy;
	FourierTransformer transformer;

public:

	FileName fn_in, fn_out;
	int padding_factor, r_min_nn;
	RFLOAT angpix, rot, tilt, psi;

	// I/O Parser
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_in = parser.getOption("--i", "Input map to be rotated");
		fn_out = parser.getOption("--o", "Rootname for output rotated map", "rotated.mrc");
		rot = textToFloat(parser.getOption("--rot", "Rot angle (in degrees)", "0"));
		tilt = textToFloat(parser.getOption("--tilt", "Tilt angle (in degrees)", "0"));
		psi = textToFloat(parser.getOption("--psi", "Psi angle (in degrees)", "0"));
                padding_factor = textToInteger(parser.getOption("--pad", "Padding factor", "2"));

		// Hidden
		r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}


	void project()
	{
		Image<RFLOAT> vol_in;
		int orig_size, interpolator;

		interpolator = TRILINEAR;

		std::cout << " Reading map: " << fn_in << std::endl;
		vol_in.read(fn_in);
		orig_size = XSIZE(vol_in());
		std:: cout << " The input box size: " << orig_size << std::endl;
		if (orig_size % 2 != 0)
			REPORT_ERROR("The input box size must be an even number.");

		angpix = vol_in.samplingRateX();
		std::cout << " Using the pixel size in the input image header: " << angpix << " A/px" << std::endl;

		// Set up the projector
		int data_dim = 3;

		std::cout << " Rotating the map by given euler angles: ROT = " << rot << " TILT = " << tilt << " PSI = " << psi << std::endl << std::endl;

		Projector full_projector(orig_size, interpolator, padding_factor, r_min_nn, data_dim);
		Image<RFLOAT> vol_out;
		FourierTransformer final_transformer;

		full_projector.computeFourierTransformMap(vol_in(), dummy, 2 * orig_size);
		Euler_rotation3DMatrix(rot, tilt, psi, A3D);
		F2D.initZeros(orig_size, orig_size, orig_size / 2 + 1);
		vol_out().reshape(vol_in());
		full_projector.get2DFourierTransform(F2D, A3D, IS_NOT_INV);

		transformer.inverseFourierTransform(F2D, vol_out());
		CenterFFT(vol_out(), false);
		vol_out.setSamplingRateInHeader(angpix);
		vol_out.write(fn_out);
		std::cout << " The rotated map has been written to " << fn_out << std::endl;

	} // end project function
};

int main(int argc, char *argv[])
{
	rotate_map app;

	try
	{
		app.read(argc, argv);
		app.project();
	}

	catch (RelionError XE)
	{
	        //prm.usage();
        	std::cerr << XE;
	        exit(1);
	}

	return 0;
}
