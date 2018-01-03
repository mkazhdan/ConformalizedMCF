/*
Copyright (c) 2018, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef GRID_INCLUDED
#define GRID_INCLUDED

#include <stdio.h>

template<class Data>
class Grid
{
protected:
	Data* values;
	int resX,resY;
public:
	Grid(void);
	~Grid(void);

	// Returns the dimension of the array
	void resolution(int& resX,int& resY) const;
	int width(void) const;
	int height(void) const;
	// Allocates memory for the array
	int resize(const int& resX,const int& resY);

	// Clears the values of the array to 0
	void clear(void);

	// Returns a reference to the indexed  array element
	Data& operator() (const int& x,const int& y);
	const Data& operator() (const int& x,const int& y) const;
	Data* operator[] (const int& x);
	// Returns the linear interpolation of the value at the spedified index
	Data operator() (const double& x,const double& y);
	Data operator() (const double& x,const double& y) const;

	// Returns the square of the L2-norm of the array elements
	template<class Real>
	Real squareNorm(void) const;

	// Reads in an array from the specified file
	int read(const char* fileName);
	int read(FILE* fp);

	// Writes out the array to the specified file
	int write(const char* fileName) const;
	int write(FILE* fp) const;

	// Returns the square of the L2-difference between two arrays
	template<class Real>
	static Real SquareDifference(const Grid& g1,const Grid& g2);

	// Returns the dot-product of two arrays
	template<class Real>
	static Real Dot(const Grid& g1,const Grid& g2);
};

#include "Grid.inl"
#endif // GRID_INCLUDED
