/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



#include <iostream>
#include <chrono>

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>


using namespace ASTex;


int main()
{
	// All following explanation could also apply to RGBA

	ImageRGBu8 image(4,4);
	// ImageRGBu8::PixelType => itk::RGBPixel<uint8>

	// pixelAbsolute return ref or const-ref to PixelType
	ImageRGBu8::PixelType& p1 = image.pixelAbsolute(0,0);
	ImageRGBu8::PixelType& p2 = image.pixelAbsolute(1,0);
	ImageRGBu8::PixelType& p3 = image.pixelAbsolute(2,0);

	// iterator.Value() return (const) PixelType
	ImageRGBu8::ConstIterator it = image.beginConstIterator();
	const ImageRGBu8::PixelType& v=it.Value();
	std::cout << v << std::endl;

	// for_all traversal fonction lambda take PixelType parameter for reference (const)
	image.for_all_pixels([] (ImageRGBu8::PixelType& p)
	{
		p[0]=65;
	});


	// Why using DoublePixelEigen:
	// - more nice constructor
	// - more nices operators (*/)

	// Examples of constructors
	p1 = ImageRGBu8::itkPixel(134);
	p2 = ImageRGBu8::itkPixel(128,128,128);

	p3 = p1 + p2; // possible overflow !
//	p3 /= 2; not possible with pixel type
	std::cout<< p1 << " + "<< p2 << " (compute uint8) = " << p3 << std::endl;

	using DoublePix = ImageRGBu8::DoublePixelEigen;

	// easy type conversion to any compatible type:

	p3 = ImageRGBu8::itkPixel((ImageRGBu8::eigenPixel(p1)+ImageRGBu8::eigenPixel(p2))/2);

	std::cout<< p1 << " + "<< p2 << " /2 (compute uint16) = " << p3 << std::endl;


	DoublePix ep1 = image.pixelEigenAbsolute(0,0);
	DoublePix ep2 = image.pixelEigenAbsolute(1,0);
	DoublePix ep3 = image.pixelEigenAbsolute(2,0);

	image.pixelEigenAbsolute(1,1) = (ep1+ep2+ep3)/3 ;

	std::cout<< image.pixelAbsolute(0,0) << " + "<< image.pixelAbsolute(1,0) << " + "<< image.pixelAbsolute(2,0) << " /3 = " << " + "<< image.pixelAbsolute(1,1) << std::endl;


	return EXIT_SUCCESS;
}

