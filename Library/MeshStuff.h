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

#ifndef MESH_STUFF_INCLUDED
#define MESH_STUFF_INCLUDED

#include <omp.h>
#include "PPolynomial.h"
#include "SparseMatrix.h"
#include "HalfEdge.h"

#define AREA_CUT_OFF 0
//#define CONFORMAL_CUT_OFF 1e-16
#define CONFORMAL_CUT_OFF 1e-15


// \int_0^1 \int_0^{1-a}(1-a-b)^l * a^m * b^n db da = l! * m! * n! / (l+m+n+2)!
template< int I > long long Factorial( void );
template< > long long Factorial< 0 >( void ) { return 1; }
template< int I > long long Factorial( void ) { return Factorial< I-1 >() * I; }
template< int I1 , int I2 , int I3 > long long MultiFactorial( void ) { return Factorial< I1 >() * Factorial< I2 >() * Factorial< I2 >(); }
long long factorial( int N  )
{
	long long f = 1;
	for( int i=1 ; i<=N ; i++ ) f *= i;
	return f;
}
template< typename Real , int Degree >
inline Real TriangleIntegral( const Real a[][3] )
{
	Real integral = Real( 0 );
	long long size = 1;
	for( int i=0 ; i<Degree ; i++ ) size *= 3;
	for( int i=0 ; i<size ; i++ )
	{
		int idx = i;
		Real product = Real( 1. );
		int count[] = { 0 , 0 , 0 };
		for( int d=0 ; d<Degree ; d++ )
		{
			int ii = idx%3;
			idx /= 3;
			product *= a[d][ii];
			count[ii]++;
		}
		integral += product * factorial( count[0] ) * factorial( count[1] ) * factorial( count[2] );
	}
	return integral / factorial( Degree+2 );
}

template< typename Real , class Vertex >
inline Real Jacobian( const Vertex& v0 , const Vertex& v1 , const Vertex& v2 )
{
	Point3D< Real > A = Point3D< Real >( v0 ) , B = Point3D< Real >( v1 ) , C = Point3D< Real >( v2 );
	Point3D< Real > N = Point3D< Real >::CrossProduct( B-A , C-A );
	return Real( Length( N ) );
}
template< typename Real , class Vertex >
inline Real Integral( const Vertex &v0 , const Vertex &v1 , const Vertex &v2 )
{
#if 1
	Real a[1][3]; // Should be a[0][3] but that's illegal
	return TriangleIntegral< Real , 0 >( a ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#else
	return 
		(
		Real
		(
		1
		) * MultiFactorial< 0 , 0 , 0 >()
		) / Real( Factorial< 2 >() ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#endif
}
template< typename Real , class Vertex >
inline Real Integral( const Vertex &v0 , const Vertex &v1 , const Vertex &v2 , Real a0 , Real a1 , Real a2 )
{
#if 1
	Real a[][3] = { { a0 , a1 , a2 } };
	return TriangleIntegral< Real , 1 >( a ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#else
	// (a0*x + a1*y + a2*z)
	return
		(
		(
		a0+a1+a2
		) * MultiFactorial< 1 , 0 , 0 >()
		) / Real( Factorial< 3 >() ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#endif
}
template< typename Real , class Vertex >
inline Real Integral( const Vertex &v0 , const Vertex &v1 , const Vertex &v2 , Real a0 , Real a1 , Real a2 , Real b0 , Real b1 , Real b2 )
{
#if 1
	Real a[][3] = { { a0 , a1 , a2 } , { b0 , b1 , b2 } };
	return TriangleIntegral< Real , 2 >( a ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#else
	// (a0*x + a1*y + a2*z) * (b0*x + b1*y + b2*z)
	return
		(
		(
		a0*b0 + a1*b1 + a2*b2
		) * MultiFactorial< 2 , 0 , 0 >() +
		(
		a0*b1 + a1*b0 + a0*b2 + a2*b0 + a1*b2 + a2*b1
		) * MultiFactorial< 1 , 1 , 0 >()
		) / Real( Factorial< 4 >() ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#endif
}
template< typename Real , class Vertex >
inline Real Integral( const Vertex &v0 , const Vertex &v1 , const Vertex &v2 , Real a0 , Real a1 , Real a2 , Real b0 , Real b1 , Real b2 , Real c0 , Real c1 , Real c2 )
{
#if 1
	Real a[][3] = { { a0 , a1 , a2 } , { b0 , b1 , b2 } , { c0 , c1 , c2 } };
	return TriangleIntegral< Real , 3 >( a ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#else
	// (a0*x + a1*y + a2*z) * (b0*x + b1*y + b2*z) * (c0*x + c1*y + c2*z)
	return
		(
		(
		a0*b0*c0 + a1*b1*c1 + a2*b2*c2
		) * MultiFactorial< 3 , 0 , 0 >() +
		(
		a0*b0*(c1+c2) + a1*b1*(c0+c2) + a2*b2*(c0+c1) +
		a0*c0*(b1+b2) + a1*c1*(b0+b2) + a2*c2*(b0+b1) +
		b0*c0*(a1+a2) + b1*c1*(a0+b2) + b2*c2*(a0+a1) 
		) * MultiFactorial< 2 , 1 , 0 >() +
		(
		a0*(b1*c2+b2*c1) + a1*(b0*c2+b2*c0) + a2*(b0*c1+b1*c0)
		) * MultiFactorial< 1 , 1 , 1 >()
		) / Real( Factorial< 5 >() ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#endif
}
template< typename Real , class Vertex >
inline Real Integral( const Vertex &v0 , const Vertex &v1 , const Vertex &v2 , Real a0 , Real a1 , Real a2 , Real b0 , Real b1 , Real b2 , Real c0 , Real c1 , Real c2 , Real d0 , Real d1 , Real d2 )
{
#if 1
	Real a[][3] = { { a0 , a1 , a2 } , { b0 , b1 , b2 } , { c0 , c1 , c2 } , { d0 , d1 , d2 } };
	return TriangleIntegral< Real , 4 >( a ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#else
	// (a0*x + a1*y + a2*z) * (b0*x + b1*y + b2*z) * (c0*x + c1*y + c2*z) * (d0*x + d1*y + d2*z)
	return
		(
		(
		a0*b0*c0*d0 + a1*b1*c1*d1 + a2*b2*c2*d2
		) * MultiFactorial< 4 , 0 , 0 >() +
		(
		a0*b0*c0*(d1+d2) + a1*b1*c1*(d0+d2) + a2*b2*c2(d0+d1) +
		a0*b0*d0*(c1+c2) + a1*b1*d1*(c0+c2) + a2*b2*d2(c0+c1) +
		a0*c0*d0*(b1+b2) + a1*c1*d1*(b0+b2) + a2*c2*d2(b0+b1) +
		b0*c0*d0*(a1+a2) + b1*c1*d1*(a0+a2) + b2*c2*d2(a0+a1)
		) * MultiFactorial< 3 , 1 , 0 >() +
		(
		a0*b0*(c1*d1+c2*d2) + a1*b1*(c0*d0+c2*d2) + a2*b2*(c0*d0+c1*d1) +
		a0*c0*(b1*d1+b2*d2) + a1*c1*(b0*d0+b2*d2) + a2*c2*(b0*d0+b1*d1) +
		a0*d0*(c1*b1+c2*b2) + a1*d1*(c0*b0+c2*b2) + a2*d2*(c0*b0+c1*b1)
		) * MultiFactorial< 2 , 2 , 0 >() +
		(
		a0*b0*(c1*d2+c2*d1) + a1*b1*(c0*d2+c2*d0) + a2*b2*(c0*d1+c1*d0) +
		a0*c0*(b1*d2+b2*d1) + a1*c1*(b0*d2+b2*d0) + a2*c2*(b0*d1+b1*d0) +
		a0*d0*(c1*b2+c2*b1) + a1*d1*(c0*b2+c2*b0) + a2*d2*(c0*b1+c1*b0)
		) * MultiFactorial< 2 , 1 , 1 >()
		) / Real( Factorial< 6 >() ) * Jacobian< Real , Vertex >( v0 , v1 , v2 );
#endif
}

// Given the triangle T=ABC with A = (0,0) , B = (1,0) , C = (0,1)
// The barycentric coordinates of a point in the triangle are ( 1-x-y , x , y )
// The integral of a function f over the triangle is:
// I[f] = \int_T f(p) dp
//		= \int_0^1 \int_0^{1-x} f(x,y) dy dx
// Setting f(x,y) = x^2 gives:
// I[f] = \int_0^1 x^2 ( \int_0^{1-x} dy ) dx
//		= \int_0^1 (x^2-x^3) dx
//		= [ 1/3 x^3 - 1/4 x^4 ]_0^1
//		= 1/12
// Setting f(x,y) = xy gives:
// I[f] = \int_0^1 x ( \int_0^{1-x} y dy ) dx
//		= \int_0^1 x ( (1-x)^2 / 2 ) dx
//		= \int_0^1 ( ( x - 2*x^2 + x^3 ) / 2 ) dx
//		= [ x^2 / 2 - 2/3 x^3 + x^4 / 4 ]_0^1 / 2
//		= [ 1/2 - 2/3 + 1/4 ]_0^1 / 2
//		= [ 6/12 - 8/3 + 3/4 ]_0^1 / 2
//		= 1 / 24
// Setting f(x,y) = (1-x-y) * x gives:
// I[f] = \int_0^1 (1-x)*x ( \int_0^{1-x} dy ) dx - 1/24
//		= \int_0^1 (1-x)^2 * x dx - 1/24
//		= 1/12 - 1/24
//		= 1 / 24
// Setting \Phi(x,y) = v_1 + (v2-v1) x + (v3-v1) y we have:
//                      | <v2-v1,v2-v1>   <v2-v1,v3-v1> |
// g = d\Phi^ T d\Phi = |                               |
//                      | <v2-v1,v3-v1>   <v3-v1,v3-v1> |
// Which gives:
//                     1               | <v3-v1,v3-v1>   -<v2-v1,v3-v1> |
// g^{-1} = -------------------------- |                                |
//          || (v2-v1) x (v3-v1 ) ||^2 |-<v2-v1,v3-v1>    <v2-v1,v2-v1> |
// And the integral of g^{-1} over T is:
//                          1             | <v3-v1,v3-v1>   -<v2-v1,v3-v1> |
// I[g^{-1}] = -------------------------- |                                |
//             2 || (v2-v1) x (v3-v1 ) || |-<v2-v1,v3-v1>    <v2-v1,v2-v1> |

struct NullData {};
typedef HalfEdgeMesh< NullData , NullData , NullData > EmptyHEMesh;



template< class Real , class Vertex >
void GetVertexNormals( const Vertex* vertices , int vCount , bool circular , Point2D< Real >* n )
{
#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ ) n[i] = Point2D< Real >( );

#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ )
	{
		if( i || circular )
		{
			int j = (i+vCount-1) % vCount;
			Point2D< Real > e = vertices[i]-vertices[j];
			n[i] += Point2D< Real >( e[1] , -e[0] );
		}
		if( i<vCount-1 || circular )
		{
			int j = (i+1) % vCount;
			Point2D< Real > e = vertices[j]-vertices[i];
			n[i] += Point2D< Real >( e[1] , -e[0] );
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ ) n[i] /= Real( sqrt( Point2D< Real >::SquareNorm( n[i] ) ) );
}
template< class Real , class Vertex >
void GetVertexNormals( const std::vector< Vertex >& vertices , bool circular , std::vector< Point2D< Real > >& n )
{
	n.resize( vertices.size() );
	GetVertexNormals( &vertices[0] , vertices.size() , circular , &n[0] );

}
template< class Real , class Vertex , class HEMesh >
void GetVertexNormals( const Vertex* vertices , const HEMesh& mesh , Point3D< Real >* n )
{
#pragma omp parallel for
	for( int i=0 ; i<mesh.vertex_size() ; i++ ) n[i] = Point3D< Real >( );

#pragma omp parallel for
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet(i);
		int v[] = { (int)mesh.index( f.halfedge().vertex() ) , (int)mesh.index( f.halfedge().next().vertex() ) , (int)mesh.index( f.halfedge().next().next().vertex() ) };
		Point3D< Real > A = Point3D< Real >( vertices[ v[0] ] );
		Point3D< Real > B = Point3D< Real >( vertices[ v[1] ] );
		Point3D< Real > C = Point3D< Real >( vertices[ v[2] ] );
		Point3D< Real > N = Point3D< Real >::CrossProduct( B-A , C-A );
		for( int j=0 ; j<3 ; j++ )
#pragma omp critical
		{
			n[ v[j] ] += N;
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<mesh.vertex_size() ; i++ ) n[i] /= Real( sqrt( Point3D< Real >::SquareNorm( n[i] ) ) );
}
template< class Real , class Vertex >
void GetSoRVertexNormals( const Vertex* vertices , Point3D< Real >* n , int sz )
{
#pragma omp parallel for
	for( int i=0 ; i<sz ; i++ ) n[i] = Point3D< Real >( );

#pragma omp parallel for
	for( int i=0 ; i<sz ; i++ )
	{
		Point2D< Real > d;
		if     ( i==0    ) d = Point2D< Real >( vertices[i+1] ) - Point2D< Real >( vertices[i  ] );
		else if( i==sz-1 ) d = Point2D< Real >( vertices[i  ] ) - Point2D< Real >( vertices[i-1] );
		else               d = Point2D< Real >( vertices[i+1] ) - Point2D< Real >( vertices[i-1] );
		n[i] = Point3D< Real >( -d[1] , d[0] , 0 );
	}
#pragma omp parallel for
	for( int i=0 ; i<sz ; i++ ) n[i] /= Real( sqrt( Point3D< Real >::SquareNorm( n[i] ) ) );
}

template< class Real , class Vertex , class HEMesh >
void GetVertexNormals( const std::vector< Vertex >& vertices , const HEMesh& mesh , std::vector< Point3D< Real > >& n )
{
	n.clear();
	n.resize( vertices.size() );
	GetVertexNormals( &vertices[0] , mesh , &n[0] );
}
template< class Real , class Vertex >
void GetSoRVertexNormals( const std::vector< Vertex >& vertices , std::vector< Point3D< Real > >& n )
{
	n.clear();
	n.resize( vertices.size() );
	GetSoRVertexNormals( &vertices[0] , &n[0] , vertices.size() );
}
template< class Real , class Vertex , class HEMesh >
void GetFaceNormals( const std::vector< Vertex >& vertices , const HEMesh& mesh , std::vector< Point3D< Real > >& n )
{
	n.clear();
	n.resize( mesh.facet_size() );

#pragma omp parallel for
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet(i);
		Point3D< Real > A = Point3D< Real >( vertices[ mesh.index( f.halfedge().vertex() ) ] );
		Point3D< Real > B = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().vertex() ) ] );
		Point3D< Real > C = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().next().vertex() ) ] );
		n[i] = Point3D< Real >::CrossProduct( B-A , C-A );
	}
	for( int i=0 ; i<n.size() ; i++ ) n[i] /= Real( sqrt( Point3D< Real >::SquareNorm( n[i] ) ) );
}
template< class Real , class Vertex >
Real SoRVolumeCenter( const std::vector< Vertex >& vertices )
{
	Real volume = Real(0) , center = Real(0);
#pragma omp parallel for reduction( + : volume , center )
	for( int i=0 ; i<vertices.size()-1 ; i++ )
	{
		Point2D< Real > p1 = Point2D< Real >( vertices[i] ) , p2 = Point2D< Real >( vertices[i+1] );
		//		Real v = fabs( p1[1]-p2[1] ) * fabs( p1[0]+p2[0] ) / 2;
		Real v = ( p1[1]-p2[1] ) * ( p1[0]+p2[0] ) / 2;
		Real c = ( p1[1]+p2[1] ) / 2 * v;
		center += c;
		volume += v;
	}
	return center / volume;
}
template< class Real , class Vertex , class HEMesh >
Point3D< Real > VolumeCenter( const std::vector< Vertex >& vertices , const HEMesh& mesh )
{
	Real volume = Real(0);
	Real centerX , centerY , centerZ;
	centerX = centerY = centerZ = 0;
	Point3D< Real > center;
#pragma omp parallel for reduction( + : volume , centerX , centerY , centerZ )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet(i);
		int idx[] = { (int)mesh.index( f.halfedge().vertex() ) , (int)mesh.index( f.halfedge().next().vertex() ) , (int)mesh.index( f.halfedge().next().next().vertex() ) };
		Point3D< Real > V[] = { Point3D< Real >( vertices[idx[0]] ) , Point3D< Real >( vertices[idx[1]] ) , Point3D< Real >( vertices[idx[2]] ) };
		XForm3x3< Real > xForm;
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) xForm(i,j) = V[i][j];
		Real v = xForm.determinant();
		Point3D< Real > c = ( V[0]+V[1]+V[2] )/4;
		c *= v;
		centerX += c[0];
		centerY += c[1];
		centerZ += c[2];
		volume += v;
	}
	return Point3D< Real >( centerX , centerY , centerZ ) / volume;
}
template< class Real , class Vertex >
Point2D< Real > AreaCenter( const std::vector< Vertex >& vertices , bool circular )
{
	Real area = Real(0);
	Real centerX , centerY;
	centerX = centerY = 0;
	Point2D< Real > center;
#pragma omp parallel for reduction( + : area , centerX , centerY )
	for( int i=0 ; i<( circular ? vertices.size() : vertices.size()-1 ) ; i++ )
	{
		int j = (i+1)%vertices.size();
		Point2D< Real > V[] = { Point2D< Real >( vertices[i] ) , Point2D< Real >( vertices[j] ) };
		XForm2x2< Real > xForm;
		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) xForm(i,j) = V[i][j];
		Real a = xForm.determinant();
		Point2D< Real > c = ( V[0]+V[1] )/3;
		c *= a;
		centerX += c[0];
		centerY += c[1];
		area += a;
	}
	return Point2D< Real >( centerX , centerY ) / area;
}
template< class Real , class Vertex >
Real Area( const Vertex* vertices , const std::vector< TriangleIndex >& triangles )
{
	Real area = Real(0);
#pragma omp parallel for reduction( + : area )
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { Point3D< Real >( vertices[triangles[i][0]] ) , Point3D< Real >( vertices[triangles[i][1]] ) , Point3D< Real >( vertices[triangles[i][2]] ) };
		Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
		Real a = sqrt( Point3D< Real >::SquareNorm( n ) );
		area += a;
	}
	return area/2;
}

template< class Real , class Vertex >
Point3D< Real > AreaCenter( const Vertex* vertices , const std::vector< TriangleIndex >& triangles )
{
	Real area = Real(0);
	Real centerX , centerY , centerZ;
	centerX = centerY = centerZ = Real(0);
#pragma omp parallel for reduction( + : area , centerX , centerY , centerZ )
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { Point3D< Real >( vertices[triangles[i][0]] ) , Point3D< Real >( vertices[triangles[i][1]] ) , Point3D< Real >( vertices[triangles[i][2]] ) };
		Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
		Real a = sqrt( Point3D< Real >::SquareNorm( n ) );
		Point3D< Real > c = (v[0]+v[1]+v[2]) * a / 3.;
		centerX += c[0];
		centerY += c[1];
		centerZ += c[2];
		area += a;
	}
	return Point3D< Real >( centerX , centerY , centerZ ) / area;
}
template< class Real , class Vertex >
Real SoRAreaCenter( const Vertex* vertices , int vCount )
{
	double area=0 , yCenter=0;
#pragma omp parallel for reduction( + : area , yCenter )
	for( int i=0 ; i<vCount-1 ; i++ )
	{
		Point2D< Real > v1 = Point2D< Real >( vertices[i+0] );
		Point2D< Real > v2 = Point2D< Real >( vertices[i+1] );
		if( v1[0]*v2[0]<0 )
		{
			Point2D< Real > v = (v1 * v2[0] - v2 * v1[0]) /( v1[0] - v2[0] );
			double l1 = Length< Real >( v1-v );
			double l2 = Length< Real >( v2-v );
			double x1 = fabs( v1[0] + v[0] ) / 2;
			double x2 = fabs( v2[0] + v[0] ) / 2;
			double y1 =     ( v1[1] + v[1] ) / 2;
			double y2 =     ( v2[1] + v[1] ) / 2;
			area += l1 * x1 + l2 * y2;
			yCenter += y1 * l1 * x1 + y2 * l2 * x2;
		}
		else
		{
			double l = Length< Real >( v2-v1 );
			double x = fabs( v1[0] + v2[0] ) / 2;
			double y =     ( v1[1] + v2[1] ) / 2;
			area += l * x ;
			yCenter += y * l * x;
		}
	}
	return yCenter/area;
}
template< class Real , class Vertex >
Point2D< Real > LengthCenter( const Vertex* vertices , int vCount , bool circular )
{
	Real length = Real(0);
	Real centerX , centerY;
	centerX = centerY = Real(0);
#pragma omp parallel for reduction( + : length , centerX , centerY )
	for( int i=0 ; i<( circular ? vCount : vCount-1 ) ; i++ )
	{
		int j = (i+1)%vCount;
		Real l = sqrt( Point2D< Real >::SquareNorm( Point2D< Real >( vertices[j] - vertices[i] ) ) );
		Point2D< Real > c = (vertices[i]+vertices[j]) / 2;
		centerX += c[0] * l;
		centerY += c[1] * l;
		length += l;
	}
	return Point2D< Real >( centerX , centerY ) / length;
}

template< class Real , class Vertex , class HEMesh >
Point3D< Real > AreaCenter( const Vertex* vertices , const HEMesh& mesh )
{
	Real area = Real(0);
	Real centerX , centerY , centerZ;
	centerX = centerY = centerZ = Real(0);
#pragma omp parallel for reduction( + : area , centerX , centerY , centerZ )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet(i);
		Point3D< Real > A = Point3D< Real >( vertices[ mesh.index( f.halfedge().vertex() ) ] );
		Point3D< Real > B = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().vertex() ) ] );
		Point3D< Real > C = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().next().vertex() ) ] );
		Point3D< Real > N = Point3D< Real >::CrossProduct( B-A , B-C );
		Real a = sqrt( Point3D< Real >::SquareNorm( N ) );
		Point3D< Real > D = (A+B+C) * a / 3.;
		centerX += D[0];
		centerY += D[1];
		centerZ += D[2];
		area += a;
	}
	return Point3D< Real >( centerX , centerY , centerZ ) / area;
}
template< class Real , class Vertex >
Real BoundingRadius( const Vertex* vertices , int vCount , Point3D< Real > center=Point3D< Real >() )
{
	Real radius2 = 0;
	for( int i=0 ; i<vCount ; i++ )
	{
		Real r2 = Point3D< Real >::SquareNorm( center-Point3D< Real >( vertices[i] ) );
		if( r2>radius2 ) radius2=r2;
	}
	return sqrt( radius2 );
}
template< class Real , class Vertex >
Real BoundingRadius2D( const Vertex* vertices , int vCount , Point2D< Real > center=Point2D< Real >() )
{
	Real radius2 = 0;
	for( int i=0 ; i<vCount ; i++ )
	{
		Real r2 = Point2D< Real >::SquareNorm( center-Point2D< Real >( vertices[i] ) );
		if( r2>radius2 ) radius2=r2;
	}
	return sqrt( radius2 );
}
template< class Real , class Vertex >
Real SoRBoundingRadius( const Vertex* vertices , int vCount , Real yCenter=0 )
{
	Point2D< Real > center( 0 , yCenter );
	Real radius = Real( 0.f );
	for( int i=0 ; i<vCount ; i++ )
	{
		Point2D< Real > v = Point2D< Real >( vertices[i] );
		Real tmp = Length< Real >( v-center );
		if( tmp > radius ) radius = tmp;
	}
	return radius;
}

template< class Real , class Vertex >
Real CurveLength( const Vertex* vertices , int vCount , bool circular )
{
	double length=0;
#pragma omp parallel for reduction( + : length )
	for( int i=0 ; i<( circular ? vCount : vCount-1 ) ; i++ )
	{
		Point2D< Real > v1 = Point2D< Real >( vertices[ i+0        ] );
		Point2D< Real > v2 = Point2D< Real >( vertices[(i+1)%vCount] );
		length += sqrt( Point2D< Real >::SquareNorm( v2-v1 ) );
	}
	return length;
}

template< class Real , class Vertex >
Real SoRArea( const Vertex* vertices , int vCount )
{
	double area=0;
#pragma omp parallel for reduction( + : area )
	for( int i=0 ; i<vCount-1 ; i++ )
	{
		Point2D< Real > v1 = Point2D< Real >( vertices[i+0] );
		Point2D< Real > v2 = Point2D< Real >( vertices[i+1] );
		if( v1[0]*v2[0]<0 )
		{
			Point2D< Real > v = ( v1 * v2[0] - v2 * v1[0] ) /( v1[0] - v2[0] );
			double l1 = Length< Real >( v1-v );
			double l2 = Length< Real >( v2-v );
			double x1 = fabs( v1[0] + v[0] ) / 2;
			double x2 = fabs( v2[0] + v[0] ) / 2;
			area += l1 * x1 + l2 * x2;
		}
		else
		{
			double l = Length< Real >( v2-v1 );
			double x = fabs( v1[0] + v2[0] ) / 2;
			area += l * x;
		}
	}
	return area * 2. * M_PI;
}

template< class Real , class Vertex , class HEMesh >
Real Area( const Vertex* vertices , const HEMesh& mesh )
{
	Real area = Real(0);
#pragma omp parallel for reduction( + : area )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet(i);
		Point3D< Real > A = Point3D< Real >( vertices[ mesh.index( f.halfedge().vertex() ) ] );
		Point3D< Real > B = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().vertex() ) ] );
		Point3D< Real > C = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().next().vertex() ) ] );
		Point3D< Real > N = Point3D< Real >::CrossProduct( B-A , B-C );
		Real a = sqrt( Point3D< Real >::SquareNorm( N ) );
		area += a;
	}
	return area/2;
}
#if 0
template< class Real , class Vertex , class HEMesh >
Real Area( const std::vector< Vertex >& vertices , const HEMesh& mesh )
{
	return Area( &vertices[0] , mesh );
}
#endif
template< class Real , class Vertex >
Real GetSoRSphericalVariation( const std::vector< Vertex >& vertices , Real& radius )
{
	// Get the vertex areas
	std::vector< Real > vAreas( vertices.size() , Real(0) );
	for( int i=0 ; i<vertices.size()-1 ; i++ )
	{
		Point2D< Real > p1 = Point2D< Real >( vertices[i] ) , p2 = Point2D< Real >( vertices[i+1] );
#if 1
		double area;
		if( p1[0]*p2[0]<0 )
		{
			Point2D< Real > p = (p1 * p2[0] - p2 * p1[0]) /( p1[0] - p2[0] );
			double l1 = Length< Real >( p1-p );
			double l2 = Length< Real >( p2-p );
			double x1 = fabs( p1[0] + p[0] ) / 2;
			double x2 = fabs( p2[0] + p[0] ) / 2;
			area = l1 * x1 + l2 * x2;
		}
		else
		{
			double l = Length< Real >( p2-p1 );
			double x = fabs( p1[0] + p2[0] ) / 2;
			area = l * x;
		}
		area *= 2. * M_PI;
#else
		double area = sqrt( Point2D< Real >::SquareNorm( p2-p1 ) ) * fabs( (p1[0]+p2[0])/2 ) * 2. * M_PI;
#endif
		if( i ) vAreas[i-1] += area;
		vAreas[i] += area;
	}

	// Get the vertex-area-weighted average
	Real areaSum = Real(0);
	Point2D< Real > center;
	for( int i=0 ; i<vertices.size() ; i++ ) center[1] += vertices[i][1] * vAreas[i] , areaSum += vAreas[i];
	center /= areaSum;

	// Get the radii and the average radius
	radius = 0;
	std::vector< Real > radii( vertices.size() , Real(0) );
	for( int i=0 ; i<vertices.size() ; i++ ) radii[i] = sqrt( Point2D< Real >::SquareNorm( Point2D< Real >( vertices[i]-center ) ) ) , radius += radii[i] * vAreas[i];
	radius /= areaSum;

	// Get the radius variation
	Real var = Real(0);
	for( int i=0 ; i<radii.size() ; i++ ) var += ( radii[i]-radius ) * ( radii[i]-radius ) * vAreas[i];
	var /= areaSum;
	return sqrt( var ) / radius;
}
template< class Real , class Vertex >
Real GetCircularVariation( const std::vector< Vertex >& vertices , bool circular , Real& radius )
{
	// Get the vertex areas
	std::vector< Real > vLengths( vertices.size() , Real(0) );
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		vLengths[i] = 0;
		if( i || circular )
		{
			int j = (i+vertices.size()-1) % vertices.size();
			vLengths[i] += sqrt( Point2D< Real >::SquareNorm( vertices[i] - vertices[j] ) );
		}
		if( i<vertices.size()-1 || circular )
		{
			int j = (i+1) % vertices.size();
			vLengths[i] += sqrt( Point2D< Real >::SquareNorm( vertices[i] - vertices[j] ) );
		}
	}

	// Get the vertex-length-weighted average
	Real lengthSum = Real(0);
	Point2D< Real > center;
	for( int i=0 ; i<vertices.size() ; i++ ) center += vertices[i] * vLengths[i] , lengthSum += vLengths[i];
	center /= lengthSum;

	// Get the radii and the average radius
	radius = 0;
	std::vector< Real > radii( vertices.size() , Real(0) );
	for( int i=0 ; i<vertices.size() ; i++ ) radii[i] = sqrt( Point2D< Real >::SquareNorm( Point2D< Real >( vertices[i]-center ) ) ) , radius += radii[i] * vLengths[i];
	radius /= lengthSum;

	// Get the radius variation
	Real var = Real(0);
	for( int i=0 ; i<radii.size() ; i++ ) var += ( radii[i]-radius ) * ( radii[i]-radius ) * vLengths[i];
	var /= lengthSum;
	return sqrt( var ) / radius;
}
template< class Real , class Vertex , class HEMesh >
Real GetSphericalVariation( const std::vector< Vertex >& vertices , const HEMesh& mesh , Real& radius )
{
	// Get the vertex areas
	std::vector< Real > vAreas( vertices.size() , Real(0) );
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet(i);
		size_t v[] = { mesh.index( f.halfedge().vertex() ) , mesh.index( f.halfedge().next().vertex() ) , mesh.index( f.halfedge().next().next().vertex() ) };
		Point3D< Real > V[] = { Point3D< Real >( vertices[v[0]] ) , Point3D< Real >( vertices[v[1]] ) , Point3D< Real >( vertices[v[2]] ) };
		double area = sqrt( Point3D< Real >::SquareNorm( Point3D< Real >::CrossProduct( V[2]-V[0] , V[1]-V[0] ) ) );
		for( int j=0 ; j<3 ; j++ ) vAreas[v[j]] += area;
	}
	// Get the vertex-area-weighted average
	Real areaSum = Real(0);
	Point3D< Real > center;
	for( int i=0 ; i<vertices.size() ; i++ ) center += vertices[i] * vAreas[i] , areaSum += vAreas[i];
	center /= areaSum;

	// Get the radii and the average radius
	radius = 0;
	std::vector< Real > radii( vertices.size() , Real(0) );
	for( int i=0 ; i<vertices.size() ; i++ ) radii[i] = sqrt( Point3D< Real >::SquareNorm( Point3D< Real >( vertices[i]-center ) ) ) , radius += radii[i] * vAreas[i];
	radius /= areaSum;

	// Get the radius variation
	Real var = Real(0);
	for( int i=0 ; i<radii.size() ; i++ ) var += ( radii[i]-radius ) * ( radii[i]-radius ) * vAreas[i];
	var /= areaSum;
	return sqrt( var ) / radius;
}
template< class Real >
struct SoRMetricData
{
	Real d00inv , d11inv , area;
	SoRMetricData( void ) { d00inv = d11inv = area = Real(0); }
	void setTraceAndDeterminant( const Point2D< Real >* v , Real& trace , Real& determinant ) const
	{
		Real d00 = ( v[0][0] + v[1][0] ) * ( v[0][0] + v[1][0] );
		Real d11 = Point2D< Real >::SquareNorm( v[0] - v[1] );
		d00 *= d00inv , d11 *= d11inv;
		trace = d00 + d11;
		determinant = d00 * d11;
	}
};
template< class Real >
struct MetricData
{
	XForm2x2< Real > dinv;
	Real area;
	MetricData( void ) { area = Real(0.); }
	void setTraceAndDeterminant( const Point3D< Real >* v , Real& trace , Real& determinant ) const
	{
		Point3D< Real > t[] = { v[1]-v[0] , v[2]-v[0] };
		XForm2x2< Real > d , d2;
		for( int  j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) d2(j,k) = Point3D< Real >::Dot( t[j] , t[k] );

		// Phi_1(x,y) -> x * t1[0] + y * t1[1]
		// Phi_2(x,y) -> x * t2[0] + y * t2[1]
		// D_1[i][j] = < t1[i] , t1[j] >
		// D_2[i][j] = < t2[i] , t2[j] >
		// Set s1[0] = D_1^{-1/2}(1,0)
		//     s1[1] = D_1^{-1/2}(0,1)
		// Then Phi_1(s1[0]) and Phi_1(s1[1]) are orthonormal
		// Want the matrix D[i][j] = < Phi_2(Phi_1(s1[i])) , Phi_2(Phi_1(s1[j])) >
		//                         = Phi_1(s1[i])^t * D_2 * Phi_1(s1[j])
		//                         = D_1^{-1/2} * D_2 * D_1{-1/2}
		// The eigenvectors of the matrix D are also the eigenvectors of the matrix
		// E = D_1^{-1/2} * D * D^{1/2} = D_1^{-1} * D_2
		d = dinv * d2;
		trace = d(0,0) + d(1,1);
		determinant = d(0,0)*d(1,1) - d(0,1)*d(1,0);
	}
};
template< class Real , class Vertex >
void InitSoRMetricData( const std::vector< Vertex >& v , std::vector< SoRMetricData< Real > >& cData )
{
	cData.resize( v.size()-1 );
	Real areaSum = Real(0);
#pragma omp parallel for reduction( + : areaSum )
	for( int i=0 ; i<cData.size() ; i++ )
	{
		Point2D< Real > p1 = Point2D< Real >( v[i] ) , p2 = Point2D< Real >( v[i+1] );
		cData[i].d00inv = ( p1[0]+p2[0] ) * ( p1[0]+p2[0] );
		cData[i].d11inv = Point2D< Real >::SquareNorm( p1-p2 );
		cData[i].area = sqrt( cData[i].d00inv * cData[i].d11inv );
		cData[i].d00inv = Real( 1./cData[i].d00inv );
		cData[i].d11inv = Real( 1./cData[i].d11inv );
		areaSum += cData[i].area;
	}
	for( int i=0 ; i<cData.size() ; i++ ) cData[i].area /= areaSum;
}
template< class Real , class Vertex , class HEMesh >
void InitMetricData( const std::vector< Vertex >& v , const HEMesh& mesh , std::vector< MetricData< Real > >& cData )
{
	cData.resize( mesh.facet_size() );
	Real areaSum = Real(0);
#pragma omp parallel for reduction( + : areaSum )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		size_t vi[] = { mesh.index( mesh.facet(i).halfedge().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().next().vertex() ) };
		Point3D< Real > _v[] = { Point3D< Real >( v[vi[0]] ), Point3D< Real >( v[vi[1]] ) , Point3D< Real >( v[vi[2]] ) };
		Point3D< Real > _t[] = { _v[1]-_v[0] , _v[2]-_v[0] };
		cData[i].area = sqrt( Point3D< Real >::SquareNorm( Point3D< Real >::CrossProduct( _t[0] , _t[1] ) ) ) / Real(2);
		areaSum += cData[i].area;
		if( cData[i].area<=Real(CONFORMAL_CUT_OFF) ) continue;
		XForm2x2< Real > d;
		for( int  j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) d(j,k) = Point3D< Real >::Dot( _t[j] , _t[k] );
		cData[i].dinv = d.inverse();
	}
#pragma omp parallel for
	for( int i=0 ; i<cData.size() ; i++ ) cData[i].area /= areaSum;
}
#if 1
template< class Real , class Vertex , class HEMesh >
Real GetInitializedAreaError( const std::vector< MetricData< Real > >& v1 , const std::vector< Vertex >& v2 , const HEMesh& mesh , bool useDeterminant=true )
{
	Real error = Real(0);
#pragma omp parallel for reduction( + : error )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		size_t v[] = { mesh.index( mesh.facet(i).halfedge().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().next().vertex() ) };
		Point3D< Real > _v2[] = { Point3D< Real >( v2[v[0]] ), Point3D< Real >( v2[v[1]] ) , Point3D< Real >( v2[v[2]] ) };
		Real det , trace;
		v1[i].setTraceAndDeterminant( _v2 , trace , det );
		Real tError = det;
		if( useDeterminant ) tError /= Real( sqrt(det) );
		else                 tError /= trace/2;
		tError *= v1[i].area;
		error += tError;
	}
	return error;
}
template< class Real , class Vertex , class HEMesh >
Real GetInitializedConformalError( const std::vector< MetricData< Real > >& v1 , const std::vector< Vertex >& v2 , const HEMesh& mesh , bool useDeterminant=true )
{
	Real error = Real(0);
#pragma omp parallel for reduction( + : error )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		size_t v[] = { mesh.index( mesh.facet(i).halfedge().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().next().vertex() ) };
		Point3D< Real > _v2[] = { Point3D< Real >( v2[v[0]] ), Point3D< Real >( v2[v[1]] ) , Point3D< Real >( v2[v[2]] ) };
		Real det , trace;
		v1[i].setTraceAndDeterminant( _v2 , trace , det );
		Real tError = trace*trace - det * 4;
		if( useDeterminant ) tError /= Real( sqrt(det) );
		else                 tError /= trace/2;
		tError *= v1[i].area;
		error += tError;
	}
	return error/2;
}
#else
template< class Real , class Vertex , class HEMesh >
Real GetInitializedAreaError( const std::vector< MetricData< Real > >& v1 , const std::vector< Vertex >& v2 , const HEMesh& mesh , bool useDeterminant=true )
{
	Real error = Real(0);
#pragma omp parallel for reduction( + : error )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		int v[] = { mesh.index( mesh.facet(i).halfedge().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().next().vertex() ) };
		Point3D< Real > _v2[] = { Point3D< Real >( v2[v[0]] ), Point3D< Real >( v2[v[1]] ) , Point3D< Real >( v2[v[2]] ) };
		Real det , trace;
		v1[i].setTraceAndDeterminant( _v2 , trace , det );
		trace /= Real(2);
		Real tError = det;
		if( useDeterminant ) tError /= Real( sqrt(det) );
		else                 tError /= trace;
		tError *= v1[i].area;
		error += tError;
	}
	return error;
}
template< class Real , class Vertex , class HEMesh >
Real GetInitializedConformalError( const std::vector< MetricData< Real > >& v1 , const std::vector< Vertex >& v2 , const HEMesh& mesh , bool useDeterminant=true )
{
	Real error = Real(0);
#pragma omp parallel for reduction( + : error )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		int v[] = { mesh.index( mesh.facet(i).halfedge().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().next().vertex() ) };
		Point3D< Real > _v2[] = { Point3D< Real >( v2[v[0]] ), Point3D< Real >( v2[v[1]] ) , Point3D< Real >( v2[v[2]] ) };
		Real det , trace;
		v1[i].setTraceAndDeterminant( _v2 , trace , det );
		trace /= Real(2);
		Real tError = trace*trace - det; // trace * trace - 4 * det
		if( useDeterminant ) tError /= Real( sqrt(det) );
		else                 tError /= trace;
		tError *= v1[i].area;
		error += tError;
	}
	return error;
}
#endif
template< class Real , class Vertex , class HEMesh >
Real GetInitializedAreaAndConformalError( const std::vector< MetricData< Real > >& v1 , const std::vector< Vertex >& v2 , const HEMesh& mesh , bool useDeterminant=true )
{
	Real error = Real(0);
#pragma omp parallel for reduction( + : error )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		int v[] = { mesh.index( mesh.facet(i).halfedge().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().next().vertex() ) };
		Point3D< Real > _v2[] = { Point3D< Real >( v2[v[0]] ), Point3D< Real >( v2[v[1]] ) , Point3D< Real >( v2[v[2]] ) };
		Real det , trace;
		v1[i].setTraceAndDeterminant( _v2 , trace , det );
		trace /= Real(2);
		Real tError = trace*trace;
		if( useDeterminant ) tError /= Real( sqrt(det) );
		else                 tError /= trace;
		tError *= v1[i].area;
		error += tError;
	}
	return error;
}

template< class Real , class Vertex , class HEMesh >
Real GetInitializedConformalRatio( const std::vector< MetricData< Real > >& v1 , const std::vector< Vertex >& v2 , const HEMesh& mesh , bool clampTarget=false )
{
	Real error = Real(0) , areaSum = Real(0);
#pragma omp parallel for reduction( + : error , areaSum )
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		size_t v[] = { mesh.index( mesh.facet(i).halfedge().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().vertex() ) , mesh.index( mesh.facet(i).halfedge().next().next().vertex() ) };
		Point3D< Real > _v2[] = { Point3D< Real >( v2[v[0]] ), Point3D< Real >( v2[v[1]] ) , Point3D< Real >( v2[v[2]] ) };
#if 1
		Real det , trace;
		v1[i].setTraceAndDeterminant( _v2 , trace , det );
		trace /= Real(2);
		if( clampTarget && det/Real(4)<=Real(CONFORMAL_CUT_OFF*CONFORMAL_CUT_OFF) ) continue;
#else
		Point3D< Real > _t2[] = { _v2[1]-_v2[0] , _v2[2]-_v2[0] };
		Real tArea2 = sqrt( Point3D< Real >::SquareNorm( Point3D< Real >::CrossProduct( _t2[0] , _t2[1] ) ) ) / Real(2);
		if( tArea2<=Real(CONFORMAL_CUT_OFF) ) continue;
		XForm2x2< Real > d , d2;
		for( int  j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) d2(j,k) = Point3D< Real >::Dot( _t2[j] , _t2[k] );

		// Phi_1(x,y) -> x * t1[0] + y * t1[1]
		// Phi_2(x,y) -> x * t2[0] + y * t2[1]
		// D_1[i][j] = < t1[i] , t1[j] >
		// D_2[i][j] = < t2[i] , t2[j] >
		// Set s1[0] = D_1^{-1/2}(1,0)
		//     s1[1] = D_1^{-1/2}(0,1)
		// Then Phi_1(s1[0]) and Phi_1(s1[1]) are orthonormal
		// Want the matrix D[i][j] = < Phi_2(Phi_1(s1[i])) , Phi_2(Phi_1(s1[j])) >
		//                         = Phi_1(s1[i])^t * D_2 * Phi_1(s1[j])
		//                         = D_1^{-1/2} * D_2 * D_1{-1/2}
		// The eigenvectors of the matrix D are also the eigenvectors of the matrix
		// E = D_1^{-1/2} * D * D^{1/2} = D_1^{-1} * D_2
		d = v1[i].dinv * d2;
		Real trace = ( d(0,0)+d(1,1) ) / Real(2);
		Real det   =   d(0,0)*d(1,1) - d(0,1)*d(1,0);
#endif
		// The eigenvectors of this matrix are the roots of
		// P(x) = [x-d(0,0)]*[x-d(1,1)] - d(1,0)*d(0,1)
		//      = x^2 - x * Tr(d) + Det(d)
		// x = Tr(d)/2 +/- sqrt( Tr^2(d)/4 - Det(d) )
		Real disc  = trace*trace-det;
		if( disc<=Real(0) ) disc = 0;
		else                disc = sqrt( disc );
		Real x1 = trace - disc;
		Real x2 = trace + disc;
		if( x1<=Real(0) ) continue;
		Real tError = sqrt( x2/x1 ) * v1[i].area;
		error += tError;
		areaSum += v1[i].area;
	}
	return error / areaSum;
}
template< class Real , class Vertex >
Real GetInitializedSoRConformalRatio( const std::vector< SoRMetricData< Real > >& v1 , const std::vector< Vertex >& v2 , bool clampTarget=false )
{
	Real error = Real(0) , areaSum = Real(0);
#pragma omp parallel for reduction( + : error , areaSum )
	for( int i=0 ; i<v2.size()-1 ; i++ )
	{
		Point2D< Real > _v2[] = { Point2D< Real >( v2[i] ) , Point2D< Real >( v2[i+1] ) };

		Real det , trace;
		v1[i].setTraceAndDeterminant( _v2 , trace , det );
		trace /= Real(2);
		if( clampTarget && det/Real(4)<=Real(CONFORMAL_CUT_OFF*CONFORMAL_CUT_OFF) ) continue;
		// The eigenvectors of this matrix are the roots of
		// P(x) = [x-d(0,0)]*[x-d(1,1)] - d(1,0)*d(0,1)
		//      = x^2 - x * Tr(d) + Det(d)
		// x = Tr(d)/2 +/- sqrt( Tr^2(d)/4 - Det(d) )
		Real disc  = trace*trace-det;
		if( disc<=Real(0) ) disc = 0;
		else                disc = sqrt( disc );
		Real x1 = trace - disc;
		Real x2 = trace + disc;
		if( x1<=Real(0) ) continue;
		Real tError = sqrt( x2/x1 ) * v1[i].area;
		error += tError;
		areaSum += v1[i].area;
	}
	return error / areaSum;
}
template< class Real , class Vertex1 , class Vertex2 , class HEMesh >
Real GetConformalRatio( const std::vector< Vertex1 >& v1 , const std::vector< Vertex2 >& v2 , const HEMesh& mesh )
{
	std::vector< MetricData< Real > > _v1;
	InitMetricData( v1 , mesh , _v1 );
	return GetInitializedConformalRatio( _v1 , v2 , mesh );
}
template< class Real >
Real GetConformalRatio( const Point3D< Real >* v1 , const Point3D< Real >* v2 )
{
	Point3D< Real > t1[] = { v1[1]-v1[0] , v1[2]-v1[0] };
	Point3D< Real > t2[] = { v2[1]-v2[0] , v2[2]-v2[0] };

	XForm2x2< Real > d , d1 , d2;
	for( int  j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
		d1(j,k) = Point3D< Real >::Dot( t1[j] , t1[k] ) , d2(j,k) = Point3D< Real >::Dot( t2[j] , t2[k] );

	// Phi_1(x,y) -> x * t1[0] + y * t1[1]
	// Phi_2(x,y) -> x * t2[0] + y * t2[1]
	// D_1[i][j] = < t1[i] , t1[j] >
	// D_2[i][j] = < t2[i] , t2[j] >
	// Set s1[0] = D_1^{-1/2}(1,0)
	//     s1[1] = D_1^{-1/2}(0,1)
	// Then Phi_1(s1[0]) and Phi_1(s1[1]) are orthonormal
	// Want the matrix D[i][j] = < Phi_2(Phi_1(s1[i])) , Phi_2(Phi_1(s1[j])) >
	//                         = Phi_1(s1[i])^t * D_2 * Phi_1(s1[j])
	//                         = D_1^{-1/2} * D_2 * D_1{-1/2}
	// The eigenvectors of the matrix D are also the eigenvectors of the matrix
	// E = D_1^{-1/2} * D * D^{1/2} = D_1^{-1} * D_2
	d = d1.inverse() * d2;
	// The eigenvectors of this matrix are the roots of
	// P(x) = [x-d(0,0)]*[x-d(1,1)] - d(1,0)*d(0,1)
	//      = x^2 - x * Tr(d) + Det(d)
	// x = Tr(d)/2 +/- sqrt( Tr^2(d)/4 - Det(d) )
	Real trace = ( d(0,0)+d(1,1) ) / Real(2);
	Real det   =   d(0,0)*d(1,1) - d(0,1)*d(1,0);
	Real disc  = trace*trace-det;
	//	if( disc<Real(0) ) fprintf( stderr , "[WARNING] negative discriminant set to zero: %g\n" , disc ) , disc = 0;
	if( disc<=Real(0) ) disc = 0;
	else          disc = sqrt( disc );
	Real x1 = trace - disc;
	Real x2 = trace + disc;
	if( x1<Real(0) ) return -1.;
	return sqrt( x2/x1 );
}

template< class Real >
void SetSoRStencil( Point3D< Real > p0 , Point3D< Real > p1 , Point3D< Real > p2 , Real dStencil[] , Real lStencil[] )
{
	Point3D< Real > v1 , v2;
	Real area;

	// Edge 0-1
	v1 = p0-p2 , v2 = p1-p2;
	area = Length( Point3D< Real >::CrossProduct( v1 , v2 ) );
	if( area>Real(AREA_CUT_OFF) )
	{
		dStencil[0] = area / Real(24.);
		lStencil[0] = -Point3D< Real >::Dot( v1 , v2 ) / area / Real(2.);
	}
	else dStencil[0] = lStencil[0] = Real(0.);

	// Edge 1-2
	v1 = p1-p0 , v2 = p2-p0;
	area = Length( Point3D< Real >::CrossProduct( v1 , v2 ) );
	if( area>Real(AREA_CUT_OFF) )
	{
		dStencil[1] = area / Real(24.);
		lStencil[1] = -Point3D< Real >::Dot( v1 , v2 ) / area / Real(2.);
	}
	else dStencil[1] = lStencil[1] = Real(0.);

	// Edge 2-0
	v1 = p2-p1 , v2 = p0-p1;
	area = Length( Point3D< Real >::CrossProduct( v1 , v2 ) );
	if( area>Real(AREA_CUT_OFF) )
	{
		dStencil[2] = area / Real(24.);
		lStencil[2] = -Point3D< Real >::Dot( v1 , v2 ) / area / Real(2.);
	}
	else dStencil[2] = lStencil[2] = Real(0.);

	dStencil[3] =   dStencil[2] + dStencil[0];
	lStencil[3] = - lStencil[2] - lStencil[0];
	dStencil[4] =   dStencil[0] + dStencil[1];
	lStencil[4] = - lStencil[0] - lStencil[1];
	dStencil[5] =   dStencil[1] + dStencil[2];
	lStencil[5] = - lStencil[1] - lStencil[2];
}
template< class Real , class CReal , class Vertex >
void GetSoRMatrices( const std::vector< Vertex > &vertices , int axialRes , SparseMatrix< Real , int >* DMatrix , SparseMatrix< Real , int >* LMatrix , SparseMatrix< Real , int >* SMatrix , bool resize )
{
	//============================ pre-allocate the sparse matrix ==============================
	if( resize )
	{
		for( int d=0 ; d<2 ; d++ )
		{
			if( DMatrix ) DMatrix[d].resize( vertices.size() );
			if( LMatrix ) LMatrix[d].resize( vertices.size() );
			if( SMatrix ) SMatrix[d].resize( vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ )
				if( i==0 || i==vertices.size()-1 )
				{
					if( DMatrix ) DMatrix[d].SetRowSize( i , 2 );
					if( LMatrix ) LMatrix[d].SetRowSize( i , 2 );
					if( SMatrix ) SMatrix[d].SetRowSize( i , 2 );
				}
				else
				{
					if( DMatrix ) DMatrix[d].SetRowSize( i , 3 );
					if( LMatrix ) LMatrix[d].SetRowSize( i , 3 );
					if( SMatrix ) SMatrix[d].SetRowSize( i , 3 );
				}
		}
	}

	//============================= construct the matrices ==============================
	double theta = 2. * M_PI / axialRes;
	double cTheta = cos(theta) , sTheta = sin(theta);

#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		CReal dStencil[4][6] , lStencil[4][6];

		Real dot0[2] , lap0[2] , dot1[2] , lap1[2] , dot2[2] , lap2[2];

		Point3D< CReal > q[3][2];
		{
			Point2D< CReal > p = Point2D< CReal >( vertices[i] );
			q[1][0] = Point3D< CReal >( p[0] , p[1] , 0 ) , q[1][1] = Point3D< CReal >( p[0]*cTheta , p[1] , p[0]*sTheta );
		}
		if( i>0 )
		{
			Point2D< CReal > p = Point2D< CReal >( vertices[i-1] );
			q[0][0] = Point3D< CReal >( p[0] , p[1] , 0 ) , q[0][1] = Point3D< CReal >( p[0]*cTheta , p[1] , p[0]*sTheta );
			SetSoRStencil( q[1][0] , q[1][1] , q[0][0] , dStencil[2] , lStencil[2] );
			SetSoRStencil( q[1][1] , q[0][1] , q[0][0] , dStencil[3] , lStencil[3] );
		}
		if( i<vertices.size()-1 )
		{
			Point2D< CReal > p = Point2D< CReal >( vertices[i+1] );
			q[2][0] = Point3D< CReal >( p[0] , p[1] , 0 ) , q[2][1] = Point3D< CReal >( p[0]*cTheta , p[1] , p[0]*sTheta );
			SetSoRStencil( q[1][0] , q[2][0] , q[2][1] , dStencil[0] , lStencil[0] );
			SetSoRStencil( q[1][0] , q[2][1] , q[1][1] , dStencil[1] , lStencil[1] );
		}

		dot0[0] = dot0[1] = lap0[0] = lap0[1] = Real( 0 );
		dot1[0] = dot1[1] = lap1[0] = lap1[1] = Real( 0 );
		dot2[0] = dot2[1] = lap2[0] = lap2[1] = Real( 0 );
		if( i==0 )
		{
			dot1[0] = Real( dStencil[0][3] * axialRes );
			lap1[0] = Real( lStencil[0][3] * axialRes );
			dot2[0] = Real( ( dStencil[0][0] + dStencil[0][2] ) * axialRes );
			lap2[0] = Real( ( lStencil[0][0] + lStencil[0][2] ) * axialRes );
		}
		else if( i==vertices.size()-1 )
		{
			dot1[0] = Real( dStencil[3][3] * axialRes );
			lap1[0] = Real( lStencil[3][3] * axialRes );
			dot0[0] = Real( ( dStencil[3][0] + dStencil[3][2] ) * axialRes );
			lap0[0] = Real( ( lStencil[3][0] + lStencil[3][2] ) * axialRes );
		}
		else
		{
			dot1[1] = Real( ( dStencil[1][2] + dStencil[2][0] ) * axialRes * 2 );
			lap1[1] = Real( ( lStencil[1][2] + lStencil[2][0] ) * axialRes * 2 );

			if( i==1 )
			{
				dot0[0] = Real( ( dStencil[2][1] + dStencil[2][2] ) * axialRes );
				lap0[0] = Real( ( lStencil[2][1] + lStencil[2][2] ) * axialRes );
				dot1[0] += Real( ( dStencil[2][3] + dStencil[2][4] ) * axialRes );
				lap1[0] += Real( ( lStencil[2][3] + lStencil[2][4] ) * axialRes );
			}
			else
			{
				dot0[0] = Real( ( dStencil[2][2] + dStencil[3][0] ) * axialRes );
				lap0[0] = Real( ( lStencil[2][2] + lStencil[3][0] ) * axialRes );
				dot0[1] = Real( ( dStencil[2][1] + dStencil[3][2] ) * axialRes );
				lap0[1] = Real( ( lStencil[2][1] + lStencil[3][2] ) * axialRes );
				dot1[0] += Real( ( dStencil[2][3] + dStencil[2][4] + dStencil[3][3] ) * axialRes );
				lap1[0] += Real( ( lStencil[2][3] + lStencil[2][4] + lStencil[3][3] ) * axialRes );
			}
			if( i==vertices.size()-2 )
			{
				dot2[0] = Real( ( dStencil[1][0] + dStencil[1][1] ) * axialRes );
				lap2[0] = Real( ( lStencil[1][0] + lStencil[1][1] ) * axialRes );
				dot1[0] += Real( ( dStencil[1][3] + dStencil[1][5] ) * axialRes );
				lap1[0] += Real( ( lStencil[1][3] + lStencil[1][5] ) * axialRes );
			}
			else
			{
				dot2[0] = Real( ( dStencil[0][0] + dStencil[1][1] ) * axialRes );
				lap2[0] = Real( ( lStencil[0][0] + lStencil[1][1] ) * axialRes );
				dot2[1] = Real( ( dStencil[0][2] + dStencil[1][0] ) * axialRes );
				lap2[1] = Real( ( lStencil[0][2] + lStencil[1][0] ) * axialRes );
				dot1[0] += Real( ( dStencil[0][3] + dStencil[1][3] + dStencil[1][5] ) * axialRes );
				lap1[0] += Real( ( lStencil[0][3] + lStencil[1][3] + lStencil[1][5] ) * axialRes );
			}
		}
		{
			int idx = 0;
			if( DMatrix )
			{
				if( i==0 || i==vertices.size()-1 ) DMatrix[0][i][0] = MatrixEntry< Real , int >( i , dot1[0] );
				else                               DMatrix[0][i][0] = MatrixEntry< Real , int >( i , (dot1[0]+dot1[1]*cTheta)/2 );
				DMatrix[1][i][0] = MatrixEntry< Real , int >( i , dot1[0] + dot1[1] );
			}
			if( LMatrix )
			{
				if( i==0 || i==vertices.size()-1 ) LMatrix[0][i][0] = MatrixEntry< Real , int >( i , lap1[0] );
				else                               LMatrix[0][i][0] = MatrixEntry< Real , int >( i , (lap1[0]+lap1[1]*cTheta)/2 );
				LMatrix[1][i][0] = MatrixEntry< Real , int >( i , lap1[0] + lap1[1] );
			}
			if( SMatrix ) SMatrix[0][i][0] = SMatrix[1][i][0] = MatrixEntry< Real , int >( i , Real(0.) );
			idx++;
			if( i )
			{
				if( DMatrix )
				{
					if( i==1 || i==vertices.size()-1 ) DMatrix[0][i][idx] = MatrixEntry< Real , int >( i-1 , 0 );
					else                               DMatrix[0][i][idx] = MatrixEntry< Real , int >( i-1 , (dot0[0]+dot0[1]*cTheta)/2 );
					DMatrix[1][i][idx] = MatrixEntry< Real , int >( i-1 , dot0[0]+dot0[1] );
				}
				if( LMatrix )
				{
					if( i==1 || i==vertices.size()-1 ) LMatrix[0][i][idx] = MatrixEntry< Real , int >( i-1 , 0 );
					else                               LMatrix[0][i][idx] = MatrixEntry< Real , int >( i-1 , (lap0[0]+lap0[1]*cTheta)/2 );
					LMatrix[1][i][idx] = MatrixEntry< Real , int >( i-1 , lap0[0]+lap0[1] );
				}
				if( SMatrix ) SMatrix[0][i][idx] = SMatrix[1][i][idx] = MatrixEntry< Real , int >( i-1 , Real(0.) );
				idx++;
			}
			if( i<vertices.size()-1 )
			{
				if( DMatrix )
				{
					if( i==0 || i==vertices.size()-2 ) DMatrix[0][i][idx] = MatrixEntry< Real , int >( i+1 , 0 );
					else                               DMatrix[0][i][idx] = MatrixEntry< Real , int >( i+1 , (dot2[0]+dot2[1]*cTheta)/2 );
					DMatrix[1][i][idx] = MatrixEntry< Real , int >( i+1 , dot2[0]+dot2[1] );
				}
				if( LMatrix )
				{
					if( i==0 || i==vertices.size()-2 ) LMatrix[0][i][idx] = MatrixEntry< Real , int >( i+1 , 0 );
					else                               LMatrix[0][i][idx] = MatrixEntry< Real , int >( i+1 , (lap2[0]+lap2[1]*cTheta)/2 );
					LMatrix[1][i][idx] = MatrixEntry< Real , int >( i+1 , lap2[0]+lap2[1] );
				}
				if( SMatrix ) SMatrix[0][i][idx] = SMatrix[1][i][idx] = MatrixEntry< Real , int >( i+1 , Real(0.) );
				idx++;
			}
		}
	}
}

template< class Real , class CReal , class Vertex , class HEMesh >
void GetMatrices( const std::vector< Vertex > &vertices , const HEMesh& mesh ,  SparseMatrix< Real , int >* DMatrix , SparseMatrix< Real , int >* LMatrix , SparseMatrix< Real ,int >* SMatrix , bool resize , bool threadSafe )
{
	//=========================== Compute the list of bad triangles ===========================
	std::vector< CReal > triangleAreas( mesh.facet_size() );
	int _threads = threadSafe ? omp_get_num_procs() : 1;


	if( DMatrix || LMatrix )
#pragma omp parallel for num_threads( _threads )
		for( int i=0 ; i<mesh.facet_size() ; i++ )
		{
			typename HEMesh::Halfedge_const_handle he = mesh.facet(i).halfedge();
			int A = (int)mesh.index( he.vertex() );
			int B = (int)mesh.index( he.next().vertex() );
			int C = (int)mesh.index( he.next().next().vertex() );
			Point3D< CReal > AB = Point3D< CReal >( vertices[A] - vertices[B] );
			Point3D< CReal > CA = Point3D< CReal >( vertices[C] - vertices[A] );
			triangleAreas[i] = sqrt( Point3D< CReal >::SquareNorm( Point3D< CReal >::CrossProduct( AB , CA ) ) ) / CReal( 2. );
		}

		//============================ pre-allocate the sparse matrix ==============================
		if( resize )
		{
			if( DMatrix ) DMatrix->resize( (int)vertices.size() );
			if( LMatrix ) LMatrix->resize( (int)vertices.size() );
			if( SMatrix ) SMatrix->resize( (int)vertices.size() );
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ )
			{
				typename HEMesh::Vertex_around_vertex_const_circulator iter = mesh.vertex_around_vertex_begin( i );
				typename HEMesh::Vertex_around_vertex_const_circulator end = iter;
				int nCount = 0;
				do
				{
					// Count the number of adjacent vertices
					nCount++;
					iter++;
				}
				while( iter!=end );

				// Set the row size
				if( DMatrix ) DMatrix->SetRowSize( i , nCount+1 );
				if( LMatrix ) LMatrix->SetRowSize( i , nCount+1 );
				if( SMatrix ) SMatrix->SetRowSize( i , nCount+1 );
			}
		}

		//============================= construct the matrices ==============================
#pragma omp parallel for
		for( int i=0 ; i<mesh.vertex_size() ; i++ )
		{
			if( DMatrix ) DMatrix->rowSizes[i] = 1;
			if( LMatrix ) LMatrix->rowSizes[i] = 1;
			if( SMatrix ) SMatrix->rowSizes[i] = 1;

			typename HEMesh::Halfedge_around_vertex_const_circulator iter = mesh.halfedge_around_vertex_begin(i);
			typename HEMesh::Halfedge_around_vertex_const_circulator end = iter;

			CReal area = 0 , weight = 0;
			do
			{
				CReal a = 0 , w = 0;
				int A = (int)mesh.index( iter.start_vertex() );
				int B = (int)mesh.index( iter.end_vertex() );

				// \int_T (1-x) * y dp = 2 * |T| \int_0^1 \int_0^x (1-x) * y dy dx
				//                     = 2 * |T| \int_0^1 (1-x) [y^2/2]_0^{1-x} dx
				//                     =     |T| \int_0^1 (1-x) * x^3 dx
				//                     =     |T| \int_0^1 (x^2-x^3) dx
				//                     =     |T| [1/3x^3-1/4x^4]_0^1
				//                     =     |T| [1/3-1/4]_0^1
				//                     =     |T| / 12

				// d_1 = v_1 - v_0 , d_2 = v_2 - v_0
				// D_{ij} = < d_i , d_j >
				// \sqrt{ \det(D) } = 2 * |T|
				//          \det(D) = 4 * |T|^2
				// \int_T < \grad (1-x) , \grad y > dp = |T| \grad^t(1-x) D^{-1} \grad(y)
				//                                     = |T| / \det(D) < d_i , d_j >
				//                                     = < d_i , d_j > / ( 4 |T| )

				typename HEMesh::Halfedge_const_handle h;

				if( DMatrix || LMatrix )
				{
					h = iter;
					if( h.facet() && triangleAreas[ mesh.index( h.facet() ) ]>CReal(AREA_CUT_OFF) )
					{
						int C = (int)mesh.index( h.next().end_vertex() );

						const CReal tArea = triangleAreas[ mesh.index( h.facet() ) ];
						if( DMatrix ) a += tArea/CReal(12.);
						if( LMatrix )
						{
							Point3D< CReal > AC = Point3D< CReal >( vertices[A] - vertices[C] );
							Point3D< CReal > BC = Point3D< CReal >( vertices[B] - vertices[C] );
							w -= Point3D< CReal >::Dot( AC , BC ) / ( tArea * CReal(4.) );
						}
					}

					h = iter.opposite();
					if( h.facet()&& triangleAreas[ mesh.index( h.facet() ) ]>CReal(AREA_CUT_OFF) )
					{
						int D =(int) mesh.index( h.next().end_vertex() );

						const CReal tArea = triangleAreas[ mesh.index( h.facet() ) ];
						if( DMatrix ) a += tArea / CReal(12.);
						if( LMatrix )
						{
							Point3D< CReal > AD = Point3D< CReal >( vertices[A] - vertices[D] );
							Point3D< CReal > BD = Point3D< CReal >( vertices[B] - vertices[D] );
							w -= Point3D< CReal >::Dot( AD , BD ) / ( tArea * CReal(4.) );
						}
					}
				}

				if( DMatrix ) (*DMatrix)[B][ DMatrix->rowSizes[B]++ ] = MatrixEntry< Real , int >( A , Real(a) ) , area += a;
				if( LMatrix ) (*LMatrix)[B][ LMatrix->rowSizes[B]++ ] = MatrixEntry< Real , int >( A , Real(w) ) , weight -= w;
				if( SMatrix ) (*SMatrix)[B][ SMatrix->rowSizes[B]++ ] = MatrixEntry< Real , int >( A , Real(0.) );

				iter++;
			}
			while( iter!=end );
			if( DMatrix ) (*DMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( area ) );
			if( LMatrix ) (*LMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( weight ) );
			if( SMatrix ) (*SMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( 0. ) );
		}
}
template< class Real , class CReal , class Vertex , class HEMesh >
void GetGraphMatrices( const std::vector< Vertex > &vertices , const HEMesh& mesh ,  SparseMatrix< Real , int >* DMatrix , SparseMatrix< Real , int >* LMatrix , SparseMatrix< Real ,int >* SMatrix , bool twoDMassMatrix , bool resize , bool threadSafe )
{
	int _threads = threadSafe ? omp_get_num_procs() : 1;



	//============================ pre-allocate the sparse matrix ==============================
	if( resize )
	{
		if( DMatrix ) DMatrix->resize( vertices.size() );
		if( LMatrix ) LMatrix->resize( vertices.size() );
		if( SMatrix ) SMatrix->resize( vertices.size() );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ )
		{
			typename HEMesh::Vertex_around_vertex_const_circulator iter = mesh.vertex_around_vertex_begin( i );
			typename HEMesh::Vertex_around_vertex_const_circulator end = iter;
			int nCount = 0;
			do
			{
				// Count the number of adjacent vertices
				nCount++;
				iter++;
			}
			while( iter!=end );

			// Set the row size
			if( DMatrix ) DMatrix->SetRowSize( i , nCount+1 );
			if( LMatrix ) LMatrix->SetRowSize( i , nCount+1 );
			if( SMatrix ) SMatrix->SetRowSize( i , nCount+1 );
		}
	}

	//============================= construct the matrices ==============================
#pragma omp parallel for num_threads( _threads )
	for( int i=0 ; i<mesh.vertex_size() ; i++ )
	{
		if( DMatrix ) DMatrix->rowSizes[i] = 1;
		if( LMatrix ) LMatrix->rowSizes[i] = 1;
		if( SMatrix ) SMatrix->rowSizes[i] = 1;

		typename HEMesh::Halfedge_around_vertex_const_circulator iter = mesh.halfedge_around_vertex_begin(i);
		typename HEMesh::Halfedge_around_vertex_const_circulator end = iter;

		int count = 0;
		CReal length = 0;
		do
		{
			CReal l = 0;
			int A = mesh.index( iter.start_vertex() );
			int B = mesh.index( iter.end_vertex() );

			l = CReal( sqrt( Point3D< CReal >::SquareNorm( vertices[A]-vertices[B] ) ) );
			// Note that though the graph matrix wants to be 1D, we make it 2D by having
			// the mass-matrix grow quadratically with scale
			if( twoDMassMatrix ) l *= l;
			if( DMatrix ) (*DMatrix)[B][ DMatrix->rowSizes[B]++ ] = MatrixEntry< Real , int >( A , Real(l/6) );
			if( LMatrix ) (*LMatrix)[B][ LMatrix->rowSizes[B]++ ] = MatrixEntry< Real , int >( A , Real(-1) );
			if( SMatrix ) (*SMatrix)[B][ SMatrix->rowSizes[B]++ ] = MatrixEntry< Real , int >( A , Real(0.) );

			length += l;
			count++;
			iter++;
		}
		while( iter!=end );
		if( DMatrix ) (*DMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( length/3 ) );
		if( LMatrix ) (*LMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( count ) );
		if( SMatrix ) (*SMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( 0. ) );
	}
}
template< class Real , class CReal , class Vertex >
void GetCurveMatrices( const std::vector< Vertex >& vertices , bool circular , SparseMatrix< Real , int >* DMatrix , SparseMatrix< Real , int >* LMatrix , SparseMatrix< Real , int >* SMatrix , bool resize , bool threadSafe )
{
	//=========================== Compute the list of bad triangles ===========================
	std::vector< CReal > edgeLengths( circular ? vertices.size() : vertices.size()-1 );
	int _threads = threadSafe ? omp_get_num_procs() : 1;


	if( DMatrix || LMatrix )
#pragma omp parallel for num_threads( _threads )
		for( int i=0 ; i<( circular ? vertices.size() : vertices.size()-1 ) ; i++ )
			edgeLengths[i] = sqrt( Point2D< CReal >::SquareNorm( vertices[i] - vertices[(i+1)%vertices.size()] ) );

	//============================ pre-allocate the sparse matrix ==============================
	if( resize )
	{
		if( DMatrix ) DMatrix->resize( vertices.size() );
		if( LMatrix ) LMatrix->resize( vertices.size() );
		if( SMatrix ) SMatrix->resize( vertices.size() );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ )
		{
			int nCount = 1;
			if( i                   || circular ) nCount++;
			if( i<vertices.size()-1 || circular ) nCount++;
			if( DMatrix ) DMatrix->SetRowSize( i , nCount );
			if( LMatrix ) LMatrix->SetRowSize( i , nCount );
			if( SMatrix ) SMatrix->SetRowSize( i , nCount );
		}
	}

	//============================= construct the matrices ==============================
#pragma omp parallel for num_threads( _threads )
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		if( DMatrix ) DMatrix->rowSizes[i] = 1;
		if( LMatrix ) LMatrix->rowSizes[i] = 1;
		if( SMatrix ) SMatrix->rowSizes[i] = 1;

		CReal length = 0 , weight = 0;
		if( i || circular )
		{
			int j = ( i+(vertices.size()-1) ) % vertices.size();
			CReal l =  edgeLengths[ j ] / 6;
			CReal w = -edgeLengths[ j ];
			length += l*2;
			weight -= w;
			if( DMatrix ) (*DMatrix)[i][ DMatrix->rowSizes[i]++ ] = MatrixEntry< Real , int >( j , Real(l) );
			if( LMatrix ) (*LMatrix)[i][ LMatrix->rowSizes[i]++ ] = MatrixEntry< Real , int >( j , Real(w) );
			if( SMatrix ) (*SMatrix)[i][ SMatrix->rowSizes[i]++ ] = MatrixEntry< Real , int >( j , Real(0) );
		}
		if( i<vertices.size()-1 || circular )
		{
			int j = ( i+1 ) % vertices.size();
			CReal l =  edgeLengths[i] / 6;
			CReal w = -edgeLengths[i];
			length += l*2;
			weight -= w;
			if( DMatrix ) (*DMatrix)[i][ DMatrix->rowSizes[i]++ ] = MatrixEntry< Real , int >( j , Real(l) );
			if( LMatrix ) (*LMatrix)[i][ LMatrix->rowSizes[i]++ ] = MatrixEntry< Real , int >( j , Real(w) );
			if( SMatrix ) (*SMatrix)[i][ SMatrix->rowSizes[i]++ ] = MatrixEntry< Real , int >( j , Real(0) );
		}
		if( DMatrix ) (*DMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( length ) );
		if( LMatrix ) (*LMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( weight ) );
		if( SMatrix ) (*SMatrix)[i][0] = MatrixEntry< Real , int >( i , Real( 0. ) );
	}
}
int SoRVertexCount( int cRes , int aRes ) { return ( cRes - 2 ) * aRes + 2; }

template< class Real >
void GetSoRProlongation( int cRes , int aRes , SparseMatrix< Real , int > P[2] )
{
	for( int d=0 ; d<2 ; d++ ) P[d].resize( cRes );
#pragma omp parallel for
	for( int i=0 ; i<cRes ; i++ )
		if( i==0 )
		{
			P[0].SetRowSize( i , 1 ) , P[1].SetRowSize( i , 1 );
			P[0][i][0].Value = P[1][i][0].Value = 1;
			P[0][i][0].N = P[1][i][0].N = 0;
		}
		else if( i==cRes-1 )
		{
			P[0].SetRowSize( i , 1 ) , P[1].SetRowSize( i , 1 );
			P[0][i][0].Value = P[1][i][0].Value = 1;
			P[0][i][0].N = P[1][i][0].N = SoRVertexCount( cRes , aRes )-1;
		}
		else
		{
			P[0].SetRowSize( i , aRes ) , P[1].SetRowSize( i , aRes );
			for( int j=0 ; j<aRes ; j++ )
			{
				double theta = ( 2. * M_PI * j ) / aRes;
				P[0][i][j].Value = cos(theta) , P[1][i][j].Value = 1;
				P[0][i][j].N = P[1][i][j].N = 1 + (i-1)*aRes + j;
			}
		}
}

template< int Bins >
class Histogram
{
	double _min , _max;
public:
	double values[Bins];
	Histogram( const Histogram& h )
	{
		_min = h._min , _max = h._max;
		memcpy( values , h.values , sizeof(double)*Bins );
	}
	Histogram( double min=0. , double max=1. )
	{
		if( min>=max ) fprintf( stderr , "Bad bounds %g<%g\n" , min , max ) , exit( 0 );
		_min = min , _max = max;
		memset( values , 0 , sizeof( double ) * Bins );
	}
	void add( double value , double weight=1.0 )
	{
		double b = (value-_min) / (_max-_min) * Bins;
		int bb = std::max< int >( 0 , std::min< int >( Bins-1 , int(b) ) );
		values[bb] += weight;
	}
	void normalize( void )
	{
		double sum = 0;
		for( int i=0 ; i<Bins ; i++ ) sum += values[i];
		if( sum ) for( int i=0 ; i<Bins ; i++ ) values[i] /= sum;
	}
	void print( FILE* fp=stdout )
	{
		for( int i=0 ; i<Bins ; i++ )
		{
			double x = (_max-_min) * (i+0.5) / Bins + _min;
			fprintf( fp , "%f %g\n" , x , values[i] );
		}
	}
};

template< class Real >
void SetConvexConcaveEdgeValueAndWeight( Point3D< Real > A , Point3D< Real > B , Point3D< Real > C , Point3D< Real > D , Real& value , Real& weight )
{
	value = weight = 0.;
	Point3D< Real > AB = A - B;

	Point3D< Real > AC = A - C;
	Point3D< Real > BC = B - C;
	Point3D< Real > AD = A - D;
	Point3D< Real > BD = B = D;

	Point3D< Real > n1 = Point3D< Real >::CrossProduct( AC , BC ); // Upward facing normal
	Point3D< Real > n2 = Point3D< Real >::CrossProduct( BD , AD ); // Upward facing normal

	n1 = Point3D< Real >::CrossProduct( n1 , AB );	// Outward facing tangent
	Real l;
	l = Point3D< Real >::SquareNorm( n2 );
	if( l )
	{
		n2 /= sqrt( l );
		l = Point3D< Real >::SquareNorm( n1 );
		if( l )
		{
			n1 /= sqrt( l );
			weight = sqrt( Point3D< Real >::SquareNorm( AB ) );
			value = Point3D< Real >::Dot( n1 , n2 );
		}
	}
}
template< class Real , class Vertex , class HEMesh , int Bins >
Histogram< Bins > GetConvexConcaveEdgeDistribution( const std::vector< Vertex > &vertices , const HEMesh& mesh , double range=1.1 )
{
	Histogram< Bins > histogram( -range , range );
	for( int i=0 ; i<mesh.halfedge_size() ; i++ )
	{
		typename HEMesh::Halfedge_const_handle h1 = mesh.halfedge(i);
		typename HEMesh::Halfedge_const_handle h2 = h1.opposite();

		if( h1.facet() && h2.facet() )
		{
			//           B
			//         / | \
			//        /  |  \
			//       C   |   D
			//        \  |  /
			//         \ | /
			//           A
			Point3D< Real > A = Point3D< Real >( vertices[ mesh.index( h1.start_vertex() ) ] );
			Point3D< Real > B = Point3D< Real >( vertices[ mesh.index( h2.start_vertex() ) ] );
			Point3D< Real > C = Point3D< Real >( vertices[ mesh.index( h1.next().end_vertex() ) ] );
			Point3D< Real > D = Point3D< Real >( vertices[ mesh.index( h2.next().end_vertex() ) ] );

			Real value , weight;
			SetConvexConcaveEdgeValueAndWeight( A , B , C , D , value , weight );
			histogram.add( value , weight );
		}
	}
	histogram.normalize();
	return histogram;
}
template< class Real , class Vertex , class HEMesh , int Bins >
std::pair< Histogram< Bins > , int > GetStarShapedTriangleDistribution( const std::vector< Vertex > &vertices , const HEMesh& mesh , Point3D< Real > center , double range=1.1 )
{
	std::pair< Histogram< Bins > , int > histogram;
	histogram.first = Histogram< Bins >( -range , range );
	histogram.second = 0;
#pragma omp parallel for
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet(i);
		Point3D< Real > A = Point3D< Real >( vertices[ mesh.index( f.halfedge().vertex() ) ] );
		Point3D< Real > B = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().vertex() ) ] );
		Point3D< Real > C = Point3D< Real >( vertices[ mesh.index( f.halfedge().next().next().vertex() ) ] );
		Point3D< Real > N1 = Point3D< Real >::CrossProduct( B-A , C-A );
		Point3D< Real > N2 = (A+B+C)/3 - center;
		Real l = sqrt( Point3D< Real >::SquareNorm( N2 ) );
		if( l )
		{
			N2 /= l;
			l = sqrt( Point3D< Real >::SquareNorm( N1 ) );
			if( l )
			{
				N1 /= l;
				Real dot = Point3D< Real >::Dot( N1 , N2 );
#pragma omp critical
				{
					histogram.first.add( dot , l );
					if( dot<0 ) histogram.second++;
				}
			}
		}
	}
	histogram.first.normalize();
	return histogram;
}
template< class Real , class Vertex , class HEMesh , int Bins >
std::pair< Histogram< Bins > , int > GetStarShapedTriangleDistribution( const std::vector< Vertex > &vertices , const HEMesh& mesh , double range=1.1 )
{
	return GetStarShapedTriangleDistribution< Real , Vertex , HEMesh , Bins >( vertices , mesh, AreaCenter< Real , Vertex , HEMesh >( &vertices[0] , mesh ) , range );
}
template< class Real , class HEMesh >
Point3D< Real > TranslateAreaCenterToOrigin( std::vector< Point3D< Real > >& vertices , HEMesh& mesh )
{
	Point3D< Real > center = -AreaCenter< Real , Point3D< Real > , HEMesh >( vertices , mesh );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] += center;
	return center;
}
template< class Real >
Real TranslateSoRVolumeCenterToOrigin( std::vector< Point2D< Real > >& vertices )
{
	Real center = -SoRVolumeCenter< Real , Point2D< Real > >( vertices );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i][1] += center;
	return center;
}
template< class Real , class HEMesh >
Point3D< Real > TranslateVolumeCenterToOrigin( std::vector< Point3D< Real > >& vertices , HEMesh& mesh )
{
	Point3D< Real > center = -VolumeCenter< Real , Point3D< Real > , HEMesh >( vertices , mesh );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] += center;
	return center;
}
template< class Real >
Point2D< Real > TranslateAreaCenterToOrigin( std::vector< Point2D< Real > >& vertices , bool circular )
{
	Point2D< Real > center = -AreaCenter< Real , Point2D< Real > >( vertices , circular );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] += center;
	return center;
}
template< class Real >
void FitVertices( std::vector< Point3D< Real > >& vertices , const Point3D< Real >& center , const Real& width , Point3D< Real >& translation , Real& scale )
{
	Point3D< Real > min , max;
	for( size_t i=0 ; i<vertices.size() ; i++ )
		for( int c=0 ; c<3 ;c++ )
		{
			if( !i || vertices[i][c]<min[c] )	min[c] = vertices[i][c];
			if( !i || vertices[i][c]>max[c] )	max[c] = vertices[i][c];
		}

		for( int c=0 ; c<3 ; c++ )
			if( !c || scale<max[c]-min[c] ) scale = max[c]-min[c];
		scale /= width;
		for( int c=0 ; c<3 ; c++ ) translation[c] = (max[c]+min[c])/Real(2) - center[c]*scale;

		for( size_t s=0 ; s<vertices.size() ; s++ ) vertices[s] = ( Point3D< Real >( vertices[s] ) - translation ) / scale;
}
template< class Real >
void FitVertices( std::vector< Point3D< Real > >& vertices , const Point3D< Real >& center , const Real& width )
{
	Point3D< Real > t;
	Real s;
	FitVertices( vertices , center , width , t , s );
}
template< class Real , class Vertex >
void BoundingBox( const std::vector< Vertex >& vertices , Point3D< Real >& min , Point3D< Real >& max )
{
	min = max = Point3D< Real >( vertices[0] );
	for( int i=1 ; i<vertices.size() ; i++ )
		for( int j=0 ; j<3 ; j++ )
		{
			if( Point3D< Real >( vertices[i] )[j]<min[j] ) min[j] = Point3D< Real >( vertices[i] )[j];
			if( Point3D< Real >( vertices[i] )[j]>max[j] ) max[j] = Point3D< Real >( vertices[i] )[j];
		}
}
template< class Real >
Real MakeUnitLength( std::vector< Point2D< Real > >& vertices , bool circular )
{
	Real length = Real(0);
	for( int i=0 ; i<( circular ? vertices.size() : vertices.size()-1 ) ; i++ )
	{
		int j = (i+1)%vertices.size();
		length += sqrt( Point2D< Real >::SquareNorm( vertices[j] - vertices[i] ) );
	}
	Real scl = Real(1.) / sqrt( length );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] *= scl;
	return Real(1.) / scl;
}
template< class Real , class HEMesh >
Real MakeUnitArea( std::vector< Point3D< Real > >& vertices , const HEMesh& mesh )
{
	Real area = Real(0);
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		typename HEMesh::Facet_const_handle f = mesh.facet( i );
		Point3D< Real > p[] = { vertices[ mesh.index( f.halfedge().vertex() ) ] , vertices[ mesh.index( f.halfedge().next().vertex() ) ] , vertices[ mesh.index( f.halfedge().next().next().vertex() ) ] };
		area += sqrt( Point3D< Real >::SquareNorm( Point3D< Real >::CrossProduct( p[1]-p[0] , p[2]-p[0] ) ) );
	}
	area /= Real(2);
	Real scl = Real(1.) / sqrt( area );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] *= scl;
	return Real(1.) / scl;
}
template< class Real >
Real SoRMakeUnitArea( std::vector< Point2D< Real > >& vertices )
{
	Real area = SoRArea< Real , Point2D< Real > >( &vertices[0] , vertices.size() );
	Real scl = Real(1.) / sqrt( area );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] *= scl;
	return Real(1.) / scl;
}
class RandomField3D
{
private:
	struct _RandomField3D
	{
		struct _Cache
		{
			bool set;
			int idx[3];
			_Cache( void ) { set=false; }
			Point3D< double > lastPoints[4][4][4];
		};
		double _a;
		int _r , _seed;
		PPolynomial< 3 > wFunction;

		_RandomField3D( void ) { wFunction = PPolynomial< 3 >::GaussianApproximation( ) , _seed = 0; }
		static long long InterleaveBits( int x , int y , int z )
		{
			long long key = 0;
			for( int i=0 ; i<32 ; i++ ) key |= ( ( x & (1<<i) )<<i ) | ( ( y & (1<<i) )<<(i+1) ) | ( ( z & (1<<i) )<<(i+2) );
			return key;
		}

		Point3D< double > _offset( int i , int j , int k )
		{
#if 1
			srand( (unsigned int)InterleaveBits( i , j , k ) );
#else
			srand(i) ,                   i = rand();
			srand(j) , rand() ,          j = rand();
			srand(k) , rand() , rand() , k = rand();
			srand( i*_r*_r + j*_r + k );
#endif
			return Point3D< double >( 2.*Random< double >()-1. , 2.*Random< double >()-1. , 2.*Random< double >()-1. ) * _a;
		}
		Point3D< double > operator() ( Point3D< double > p )
		{
			_Cache c=_Cache();
			return operator()( p , c );
		}
		Point3D< double > operator() ( Point3D< double > p , _Cache& c )
		{
			int idx[3] , *_idx = c.idx;
			double d[3];
			double weights[3][4];
			for( int i=0 ; i<3 ; i++ ) idx[i] = int( floor( p[i]*_r ) )-1 , d[i] = p[i]*_r - idx[i];
			if( !c.set || idx[0]!=_idx[0] || idx[1] !=_idx[1] || idx[2]!=_idx[2] )
				for( int j=0 ; j<4 ; j++ ) 
					for( int k=0 ; k<4 ; k++ ) 
						for( int l=0 ; l<4 ; l++ ) 
							c.lastPoints[j][k][l] = _offset( idx[0]+j , idx[1]+k , idx[2]+l );
			c.set = true;
			_idx[0] = idx[0] , _idx[1] = idx[1] , _idx[2] = idx[2];
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<4 ; j++ ) weights[i][j] = wFunction( d[i]-j );

			Point3D< double > off;
			for( int j=0 ; j<4 ; j++ ) 
				for( int k=0 ; k<4 ; k++ ) 
					for( int l=0 ; l<4 ; l++ ) 
						off += c.lastPoints[j][k][l] * weights[0][j] * weights[1][k] * weights[2][l];
			return off;
		}
	};
	std::vector< _RandomField3D > d;
public:
	struct Cache
	{
		std::vector< _RandomField3D::_Cache > caches;
	};
	RandomField3D( int depth , double amplitude=1. , double persistence=0.5 )
	{
		d.resize( depth );
		for( int i=0 ; i<d.size() ; i++ ) d[i]._r = 2<<i , d[i]._a = (amplitude*=persistence);
	}

	template< class Real >
	Point3D< Real > operator() ( Point3D< Real > p , Cache& c=Cache() )
	{
		if( !c.caches.size() ) c.caches.resize( d.size() );
		Point3D< double > q;
		for( int i=0 ; i<d.size() ; i++ ) q += d[i]( Point3D< double >( p ) , c.caches[i] );
		return Point3D< Real >( q );
	}
};
template< class Real >
struct BoundingRectangle
{
	Point2D< Real > p[2];
	BoundingRectangle( void ){;}
	BoundingRectangle( Point2D< Real > q ){ p[0] = p[1] = q; }
	BoundingRectangle( Point2D< Real > q1 , Point2D< Real > q2 )
	{
		p[0][0] = std::min< Real >( q1[0] , q2[0] );
		p[0][1] = std::min< Real >( q1[1] , q2[1] );
		p[1][0] = std::max< Real >( q1[0] , q2[0] );
		p[1][1] = std::max< Real >( q1[1] , q2[1] );
	}
	void insert( Point2D< Real > q )
	{
		p[0][0] = std::min< Real >( p[0][0] , q[0] );
		p[0][1] = std::min< Real >( p[0][1] , q[1] );
		p[1][0] = std::max< Real >( p[1][0] , q[0] );
		p[1][1] = std::max< Real >( p[1][1] , q[1] );
	}
	BoundingRectangle& operator += ( BoundingRectangle br )
	{
		insert( br.p[0] ) , insert( br.p[1] );
		return *this;
	}
	BoundingRectangle operator + ( BoundingRectangle br ) const
	{
		Point2D< Real > q1 = Point2D< Real >( std::min< Real >( p[0][0] , br.p[0][0] ) , std::min< Real >( p[0][1] , br.p[0][1] ) );
		Point2D< Real > q2 = Point2D< Real >( std::max< Real >( p[1][0] , br.p[1][0] ) , std::max< Real >( p[1][1] , br.p[1][1] ) );
		return BoundingRectangle( q1 , q2 );
	}
	void expand( Real eps )
	{
		p[0][0] -= eps , p[0][1] -= eps;
		p[1][0] += eps , p[1][1] += eps;
	}
	char isOutside( Point2D< Real > q ) const
	{
		char outside = 0;
		if( q[0]<=p[0][0] ) outside |= 1;
		if( q[0]> p[1][0] ) outside |= 2;
		if( q[1]<=p[0][1] ) outside |= 4;
		if( q[1]> p[1][1] ) outside |= 8;
		return outside;
	}

	static bool Overlap( BoundingRectangle r1 , BoundingRectangle r2 )
	{
		return !( r1.p[0][0]>r2.p[1][0] || r1.p[0][1]>r2.p[1][1] || r1.p[1][0]<r2.p[0][0] || r1.p[1][1]<r2.p[0][1] );
	}
};
template< class Real >
bool Intersect( const std::pair< Point2D< Real > , Point2D< Real > >& ls1 , const std::pair< Point2D< Real > , Point2D< Real > >& ls2 , Real eps=Real(0) )
{
	Point2D< Real > d1 = ls1.second-ls1.first , d2 = ls2.second-ls2.first;
	Point2D< Real > n1( d1[1] , -d1[0] ) , n2( d2[1] , -d2[0] );
	if( Point2D< Real >::Dot( n1 , ls2.first-ls1.first )*Point2D< Real >::Dot( n1 , ls2.second-ls1.first )>=eps ) return false;
	if( Point2D< Real >::Dot( n2 , ls1.first-ls2.first )*Point2D< Real >::Dot( n2 , ls1.second-ls2.first )>=eps ) return false;
	// Solve for t such that:
	// 0 = < ls1.first + t * d1 - ls2.first , n2 >
	//   = < ls1.first-ls2.first , n2 > + t * < d1 , n2 >
	// t = < ls2.first-ls1.first , n2 > / < d1 , n2 >
	Real t = Point2D< Real >::Dot( ls2.first-ls1.first , n2 ) / Point2D< Real >::Dot( d1 , n2 );
	if( t<-eps || t>1+eps ) return false;
	else return true;
}

int SoRVertexIndex( int cRes , int aRes , int idx )
{
	if( idx==0 ) return 0;
	else return ( (idx-1) / aRes ) + 1;
}
template< class Vertex2D , class Real >
void SetSoRVertices( const Vertex2D* v2 , Point3D< Real >* v3 , int cRes , int aRes )
{
	Point2D< Real > p;
	p = Point2D< Real >( v2[0] );
	v3[0] = Point3D< Real >( p[0] , p[1] , 0 );
#pragma omp parallel for
	for( int i=1 ; i<cRes-1 ; i++ )
	{
		p = Point2D< Real >( v2[i] );
		for( int j=0 ; j<aRes ; j++ )
		{
			double theta = ( 2. * M_PI * j ) / aRes;
			v3[ 1 + (i-1)*aRes + j ] = Point3D< Real >( p[0] * cos(theta) , p[1] , p[0] * sin(theta) );
		}
	}
	p = Point2D< Real >( v2[cRes-1] );
	v3[ SoRVertexCount( cRes , aRes ) - 1 ] = Point3D< Real >( p[0] , p[1] , 0 );
}
template< class Vertex3D , class Real >
void UnsetSoRVertices( const Vertex3D* v3 , Point2D< Real >* v2 , int cRes , int aRes )
{
	Point3D< Real > p;
	p = Point3D< Real >( v3[0] );
	v2[0] = Point2D< Real >( p[0] , p[1] );
#pragma omp parallel for
	for( int i=1 ; i<cRes-1 ; i++ )
	{
#if 1
		Point2D< Real > _p;
		for( int j=0 ; j<aRes ; j++ )
		{
			Point3D< Real > p = Point3D< Real >( v3[ 1 + (i-1)*aRes + j ] );
			double theta = ( 2. * M_PI * j ) / aRes;
			_p[0] += Real( cos(theta) * p[0] + sin(theta) * p[1] );
			_p[1] += p[1];
		}
		v2[i] = _p / aRes;
#else
		Point3D< Real > p = Point3D< Real >( v3[ 1 + (i-1)*aRes ] );
#endif
		v2[i] = Point2D< Real >( p[0] , p[1] );
	}
	p = Point3D< Real >( v3[ SoRVertexCount( cRes , aRes ) - 1 ] );
	v2[cRes-1] = Point2D< Real >( p[0] , p[1] );
}

void SetSoRTriangles( int cRes , int aRes , std::vector< TriangleIndex >& triangles )
{
	int vCount = SoRVertexCount( cRes , aRes );
	triangles.reserve( 2*aRes + 2 * ( cRes - 3 ) * aRes );
	TriangleIndex tri;
	// Add the caps
	{
		int i=0;
		int start = 0 , end = aRes;

		for( int j=start ; j<end ; j++ )
		{
			tri[0] = 0;
			tri[1] = 1 + j;
			tri[2] = 1 + ( (j+1)%aRes );
			triangles.push_back( tri );
		}
	}
	{
		int i=cRes-1;
		int start = 0 , end = aRes;

		for( int j=start ; j<end ; j++ )
		{
			tri[0] = vCount-1;
			tri[1] = 1 + (cRes-3)*aRes + ( (j+1)%aRes );
			tri[2] = 1 + (cRes-3)*aRes + j;
			triangles.push_back( tri );
		}
	}
	// Add everything else
	for( int i=1 ; i<cRes-2 ; i++ )
	{
		int i1=i-1 , i2=i;
		int start = 0 , end = aRes;
		for( int j=start ; j<end ; j++ )
		{
			int j1=j , j2=(j+1)%aRes;
			tri[0] = 1 + i1*aRes + j1;
			tri[1] = 1 + i2*aRes + j2;
			tri[2] = 1 + i1*aRes + j2;
			triangles.push_back( tri );

			tri[0] = 1 + i2*aRes + j2;
			tri[1] = 1 + i1*aRes + j1;
			tri[2] = 1 + i2*aRes + j1;
			triangles.push_back( tri );
		}
	}
}
#endif // MESH_STUFF_INCLUDED