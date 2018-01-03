/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

#ifndef POLYNOMIAL_INCLUDED
#define POLYNOMIAL_INCLUDED

#include <vector>

template< unsigned int Degree >
class Polynomial
{
public:
	double coefficients[ Degree+1 ];

	Polynomial( void );
	Polynomial( double c );
	template< unsigned int Degree2 > Polynomial( const Polynomial< Degree2 >& P );
	double operator()( double t ) const;
	double integral( double tMin , double tMax ) const;
	double& operator[] ( unsigned int idx );
	const double& operator[] ( unsigned int idx ) const;

	int operator == (const Polynomial& p) const;
	int operator != (const Polynomial& p) const;
	int isZero( void ) const;
	void setZero( void );

	template< unsigned int Degree2 >
	Polynomial& operator  = ( const Polynomial< Degree2 > &p );
	Polynomial& operator += ( const Polynomial& p );
	Polynomial& operator -= ( const Polynomial& p );
	Polynomial  operator -  ( void ) const;
	Polynomial  operator +  ( const Polynomial& p ) const;
	Polynomial  operator -  ( const Polynomial& p ) const;

	template< unsigned int Degree2 > Polynomial< ( Degree>Degree2 ? Degree : Degree2 ) > operator + ( const Polynomial< Degree2 >& p ) const;
	template< unsigned int Degree2 > Polynomial< ( Degree>Degree2 ? Degree : Degree2 ) > operator - ( const Polynomial< Degree2 >& p ) const;
	template< unsigned int Degree2 > Polynomial< Degree+Degree2 > operator * ( const Polynomial< Degree2 >& p ) const;

	Polynomial& operator += ( double s );
	Polynomial& operator -= ( double s );
	Polynomial& operator *= ( double s );
	Polynomial& operator /= ( double s );
	Polynomial  operator +  ( double s ) const;
	Polynomial  operator -  ( double s ) const;
	Polynomial  operator *  ( double s ) const;
	Polynomial  operator /  ( double s ) const;

	Polynomial shrink( double s ) const;
	Polynomial scale ( double s ) const;
	Polynomial shift ( double t ) const;

	Polynomial< (Degree>0 ? Degree-1 : 0) > derivative( void ) const;
	Polynomial< Degree+1 > integral( void ) const;

	void printnl( char x='x' , const char* header=NULL ) const;

	Polynomial& addScaled( const Polynomial& p , double scale );

	static void Negate   ( const Polynomial& in , Polynomial& out );
	static void Subtract ( const Polynomial& p1 , const Polynomial& p2 , Polynomial& q );
	static void Scale    ( const Polynomial& p  , double w , Polynomial& q );
	static void AddScaled( const Polynomial& p1 , double w1 , const Polynomial& p2 , double w2 , Polynomial& q );
	static void AddScaled( const Polynomial& p1 , const Polynomial& p2 , double w2 , Polynomial& q );
	static void AddScaled( const Polynomial& p1 , double w1 , const Polynomial& p2 , Polynomial& q );

	void getSolutions( double c , std::vector< double >& roots , double EPS ) const;
};
//template< int Degree1 , int Degree2 >
//bool Factor( const Polynomial< Degree1 >& numerator , const Polynomial< Degree2 >& denominator , Polynomial< ( Degree1>Degree2 ? Degree1-Degree2 : 0 ) >& quotient , Polynomial< ( Degree2>0 ? Degree2-1 : 0 ) >& remainder );

template< unsigned int Degree1 , unsigned int Degree2 >
class Polynomial2
{
public:
	double coefficients[ Degree1+1 ][ Degree2+1 ];

	Polynomial2( void );
	template< unsigned int _Degree1 , unsigned int _Degree2 > Polynomial2( const Polynomial2< _Degree1 , _Degree2 >& p );
	double operator()( double s , double t ) const;
	double integral( double sMin , double sMax , double tMin , double tMax ) const;

	bool operator == ( const Polynomial2& p ) const;
	bool operator != ( const Polynomial2& p ) const;
	bool isZero( void ) const;
	void setZero( void );

	template< unsigned int _Degree1 , unsigned int _Degree2 > Polynomial2& operator = ( const Polynomial2< _Degree1 , _Degree2 > &p );
	Polynomial2& operator += ( const Polynomial2& p );
	Polynomial2& operator -= ( const Polynomial2& p );
	Polynomial2  operator -  ( void ) const;
	Polynomial2  operator +  ( const Polynomial2& p ) const;
	Polynomial2  operator -  ( const Polynomial2& p ) const;
	template< unsigned int _Degree1 , unsigned int _Degree2 >
	Polynomial2< Degree1+_Degree1 , Degree2+_Degree2 > operator * ( const Polynomial2< _Degree1 , _Degree2 >& p ) const;

	Polynomial2& operator += ( double s );
	Polynomial2& operator -= ( double s );
	Polynomial2& operator *= ( double s );
	Polynomial2& operator /= ( double s );
	Polynomial2  operator +  ( double s ) const;
	Polynomial2  operator -  ( double s ) const;
	Polynomial2  operator *  ( double s ) const;
	Polynomial2  operator /  ( double s ) const;
	
	Polynomial2< Degree1+Degree2 , Degree1+Degree2 > xform( const double* M ) const;	// This is a 2x2 matrix indexed in row-major format
	Polynomial2 scale( double s , double t ) const;
	Polynomial2 shift( double s , double t ) const;

	std::pair< Polynomial2< (Degree1>0 ? Degree1-1 : 0) , Degree2 > , Polynomial2< Degree1 , (Degree2>0 ? Degree2-1 : 0) > > derivative( void ) const;

	void printnl( char x='x' , char y='y' , const char* header=NULL ) const;
};

#include "Polynomial.inl"
#endif // POLYNOMIAL_INCLUDED
