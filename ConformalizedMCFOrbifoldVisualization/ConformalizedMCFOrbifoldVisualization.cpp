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


#undef USE_CHOLMOD
#undef USE_SUITESPARSE

#define USE_EIGEN
#undef EIGEN_USE_MKL_ALL

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <omp.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "Library/CmdLineParser.h"
#include "Library/Ply.h"
#include "Library/Timer.h"
#include "Library/Image.h"
#include "Symmetry.h"
#include "SystemData.h"
#include "VisualizationData.h"

const char DefaultSymmetry[] = "C1";

static inline void gets( char* str , int strLen ){ fgets( str , strLen , stdin ); }

cmdLineParameter< char* > In( "in" ) , Symmetry( "sym" );
cmdLineParameter< int > Threads( "threads" , omp_get_num_procs() );
cmdLineParameter< float > StepSize( "stepSize" , 1.f );
cmdLineSequence< int > Cones( "cones" );
cmdLineReadable  Verbose( "verbose" ) , Simple( "simple" );
cmdLineReadable* params[] = { &In , &StepSize , &Threads , &Symmetry , &Verbose , &Simple , &Cones , NULL };

void Usage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <symmetry type>=%s]\n" , Symmetry.name , DefaultSymmetry );
	printf( "\t\tC<N>: Cyclic symmetry of order N\n" );
	printf( "\t\tD<N>: Dihedral symmetry of order 2N\n" );
	printf( "\t\tT:    Tetrahedral symmetry\n" );
	printf( "\t\tO:     Octahedral symmetry\n" );
	printf( "\t\tI:    Icosahedral symmetry\n" );
	printf( "\t[--%s <step size>=%f]\n" , StepSize.name , StepSize.value );
	printf( "\t[--%s <source> or <start:end> or <start:source:end>\n" , Cones.name );
	printf( "\t[--%s]\n" , Simple.name );
	printf( "\t[--%s]\n" , Verbose.name );
}


void SetLoops( const std::vector< std::pair< unsigned int , unsigned int > >& edges , std::vector< std::vector< unsigned int > >& loops )
{
	std::vector< std::pair< unsigned int , unsigned int > > _edges = edges;
	while( _edges.size() )
	{
		std::vector< unsigned int > loop;
		std::pair< unsigned int , unsigned int > edge = _edges.back() ; _edges.pop_back();
		int start = edge.first , current = edge.second;
		loop.push_back( start );
		while( true )
		{
			int i;
			for( i=0 ; i<_edges.size() ; i++ ) if( _edges[i].first==current )
			{
				loop.push_back( _edges[i].first );
				current = _edges[i].second;
				break;
			}
			if( i==_edges.size() ) fprintf( stderr , "[ERROR] Could not close loop\n" ) , exit( 0 );
			_edges[i] = _edges.back() ; _edges.pop_back();
			if( current==start ) break;
		}
		loops.push_back( loop );
	}
}
void SetLoops( const std::vector< TriangleIndex >& triangles , std::vector< std::vector< unsigned int > >& loops )
{
	auto EdgeIndex = []( unsigned long long v1 , unsigned long long v2 ){ return (v1<<32) | v2; };
	auto FactorEdgeIndex = []( unsigned long long idx , unsigned int& v1 , unsigned int&v2 ){ v1 = (unsigned int)(idx>>32) , v2 =  (unsigned int)( (idx<<32)>>32 ); };
	std::unordered_set< unsigned long long > edgeMap;
	std::vector< std::pair< unsigned int , unsigned int > > edges;
	for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) edgeMap.insert( EdgeIndex( triangles[i][j] , triangles[i][(j+1)%3] ) );

	for( auto i=edgeMap.begin() ; i!=edgeMap.end() ; i++ )
	{
		unsigned int v1 , v2;
		FactorEdgeIndex( *i , v1 , v2 );
		if( edgeMap.find( EdgeIndex(v2,v1) )==edgeMap.end() ) edges.push_back( std::pair< unsigned int , unsigned int >( (int)v1 , (int)v2 ) );
	}
	SetLoops( edges , loops );
}

template< class Real >
Point3D< Real > NormalColor( Point3D< Real > n )
{
	Point3D< Real > c = ( -n + Point3D< Real >( 1 , 1 , 1 ) ) * Real(128);
	for( int d=0 ; d<3 ; d++ )
	{
		if( c[d]>Real(255) ) c[d] = Real(255);
		if( c[d]<Real(  0) ) c[d] = Real(  0);
	}
	return c;
}
template< class Real >
void SetNormalColors( const std::vector< TriangleIndex >& triangles , std::vector< PlyColorVertex< Real > >& verts )
{
	for( int i=0 ; i<verts.size() ; i++ ) verts[i].color = Point3D< Real >();
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { verts[triangles[i][0]].point , verts[triangles[i][1]].point , verts[triangles[i][2]].point };
		Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
		for( int j=0 ; j<3 ; j++ ) verts[ triangles[i][j] ].color += n;
	}
	for( int i=0 ; i<verts.size() ; i++ ) verts[i].color = NormalColor< Real >( verts[i].color/(Real)Length( verts[i].color ) );
}

struct DijkstraVertex
{
	int vertex , prev;
	double dist;
	DijkstraVertex( int v=0 , int p=-1 , double d=0 ) : vertex(v) , prev(p) , dist(d) {}
};
template< typename Vertex >
std::vector< DijkstraVertex > RunDijkstra( const std::vector< TriangleIndex >& triangles , const std::vector< Vertex >& vertices , int source )
{
	auto EdgeIndex = []( unsigned long long v1 , unsigned long long v2 ){ return (v1<<32) | v2; };
	auto FactorEdgeIndex = []( unsigned long long idx , unsigned int& v1 , unsigned int&v2 ){ v1 = (unsigned int)(idx>>32) , v2 =  (unsigned int)( (idx<<32)>>32 ); };

	std::unordered_set< unsigned long long > edges;
	for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) edges.insert( EdgeIndex( triangles[i][j] , triangles[i][(j+1)%3] ) );

	std::vector< std::vector< int > > oneRings( vertices.size() );
	for( auto i=edges.begin() ; i!=edges.end() ; i++ )
	{
		unsigned int v1 , v2;
		FactorEdgeIndex( *i , v1 , v2 );
		oneRings[v1].push_back( v2 );
		oneRings[v2].push_back( v1 );
	}

	std::vector< bool > processed( vertices.size() , false );
	std::vector< DijkstraVertex > dVertices( vertices.size() );
	struct S
	{
		int vertex , prev;
		double dist;
		S( int v , int p=-1 , double d=-1 ) : vertex(v) , prev(p) , dist(d) {}
		static bool Compare( S s1 , S s2 )
		{
			if     ( s1.dist<0 ) return false;
			else if( s2.dist<0 ) return true;
			else                 return s1.dist<s2.dist;
		}
	};
	auto SCompare = []( S s1 , S s2 )
	{
		if     ( s1.dist<0 ) return true;
		else if( s2.dist<0 ) return false;
		else                 return s1.dist>s2.dist;
	};
	std::priority_queue< S , std::vector< S > , decltype( SCompare ) > q( SCompare );
	q.push( S( source , -1 , 0 ) );

	while( q.size() )
	{
		S s = q.top() ; q.pop();
		if( !processed[ s.vertex ] )
		{
			dVertices[ s.vertex ] = DijkstraVertex( s.vertex , s.prev , s.dist );
			for( int i=0 ; i<oneRings[s.vertex].size() ; i++ ) if( !processed[ oneRings[s.vertex][i] ] )
			{
				Point3D< float > d = Point3D< float >( vertices[ s.vertex ] ) - Point3D< float >( vertices[ oneRings[s.vertex][i] ] );
				q.push( S( oneRings[s.vertex][i] , s.vertex , s.dist + Length(d) ) );
			}
			processed[ s.vertex ] = true;
		}
	}
	return dVertices;
}

int main
(
	int argc ,
	char* argv[]
)
{
	typedef PlyColorVertex< float > Vertex;
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	omp_set_num_threads( Threads.value > 1 ? Threads.value : 1 );

	SymmetryGroup *G;
	OrbifoldMesh< Vertex > oMesh;
	std::vector< unsigned int > seam;

	// Read the mesh
	{
		int fileType;
		bool propertyFlags[ Vertex::ReadComponents ];
		PlyReadTriangles( In.value , oMesh.vertices , oMesh.triangles , Vertex::ReadProperties , propertyFlags , Vertex::ReadComponents , fileType );
		oMesh.vInfo.resize( oMesh.vertices.size() );
		bool hasColor = ( propertyFlags[3] || propertyFlags[6] ) && ( propertyFlags[4] || propertyFlags[7] ) && ( propertyFlags[5] || propertyFlags[8] );
		if( !hasColor ) SetNormalColors( oMesh.triangles , oMesh.vertices );

	}

	// Center at the origin
	{
		Point3D< double > center;
		double area = 0;
		for( int i=0 ; i<oMesh.triangles.size() ; i++ )
		{
			Point3D< double > v[] = { Point3D< double >( oMesh.vertices[ oMesh.triangles[i][0] ] ) , Point3D< double >( oMesh.vertices[ oMesh.triangles[i][1] ] ) , Point3D< double >( oMesh.vertices[ oMesh.triangles[i][2] ] ) };
			double a = Length( Point3D< double >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) );
			center += ( v[0] + v[1] + v[2] ) / 3. * a;
			area += a;
		}
		center /= area;
		for( int i=0 ; i<oMesh.vertices.size() ; i++ ) oMesh.vertices[i].point -= Point3D< float >( center );
	}
	
	// Compute the boundary
	{
		std::vector< std::vector< unsigned int > > boundaryLoops;
		SetLoops( oMesh.triangles , boundaryLoops );
		if( !boundaryLoops.size() ) seam = std::vector< unsigned int >();
		else if( boundaryLoops.size()!=1 ) fprintf( stderr , "[ERROR] Expected at most one boundary loop(s): %d\n" , (int)boundaryLoops.size() ) , exit( 0 );
		else seam = boundaryLoops[0];
	}

	int symmetryOrder , symmetryType;
	{
		const char* s = Symmetry.set ?  Symmetry.value : DefaultSymmetry;
		switch( s[0] )
		{
		case 'C': case 'c':
			symmetryType = SYMMETRY_CYCLIC;
			symmetryOrder = atoi( s+1 );
			if( seam.size() ) G = new CyclicSymmetryGroup< true  >( symmetryOrder , Verbose.set );
			else              G = new CyclicSymmetryGroup< false >( symmetryOrder , Verbose.set );
			break;
		case 'D': case 'd':
			symmetryType = SYMMETRY_DIHEDRAL;
			symmetryOrder = atoi( s+1 );
			if( seam.size() ) G = new DihedralSymmetryGroup< true  >( symmetryOrder , Verbose.set );
			else              G = new DihedralSymmetryGroup< false >( symmetryOrder , Verbose.set );
			break;
		case 'T': case 't':
			symmetryType = SYMMETRY_TETRAHEDRAL;
			if( seam.size() ) G = new TetrahedralSymmetryGroup< true  >( Verbose.set );
			else              G = new TetrahedralSymmetryGroup< false >( Verbose.set );
			break;
		case 'O': case 'o':
			symmetryType = SYMMETRY_OCTAHEDRAL;
			if( seam.size() ) G = new OctahedralSymmetryGroup< true  >( Verbose.set );
			else              G = new OctahedralSymmetryGroup< false >( Verbose.set );
			break;
		case 'I': case 'i':
			symmetryType = SYMMETRY_ICOSAHEDRAL;
			if( seam.size() ) G = new IcosahedralSymmetryGroup< true  >( Verbose.set );
			else              G = new IcosahedralSymmetryGroup< false >( Verbose.set );
			break;
		default:
			fprintf( stderr , "[ERROR] Unrecognized symmetry type: %c\n" , s[0] );
			exit( 0 );
		}
	}

	VisualizationData< double >* vData = new VisualizationData< double >( *G );
	vData->stepSize = StepSize.value;

	{
		if( seam.size() )
		{
			if( symmetryType==SYMMETRY_CYCLIC )
			{
				int src1 = 0 , src2 = (int)seam.size()/2 , src3 = (int)seam.size();
				oMesh.vInfo[ seam[ src1 ] ].fixedAxis = G->fixedAxisIndex( Point3D< double >( 0 , 0 , 1 ) );
				oMesh.vInfo[ seam[ src2 ] ].fixedAxis = G->fixedAxisIndex( Point3D< double >( 0 , 0 , 1 ) );
				for( int i=src1+1 ; i<src2 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = G->fixedPlaneIndex( Point3D< double >( 1 , 0 , 0 ) );
				double theta = M_PI/symmetryOrder;
				Point3D< double > n( cos(theta) , sin(theta) , 0 );
				for( int i=src2+1 ; i<src3 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = G->fixedPlaneIndex( n );

				vData->seamVertices[0] = seam[src1] , vData->seamVertices[1] = seam[src2] , vData->seamVertices[2] = vData->seamVertices[3] = -1;
			}
			else if( symmetryType==SYMMETRY_DIHEDRAL )
			{
				int src1 = 0 , src2 = (int)seam.size()/3 , src3 = ( (int)seam.size()*2 ) / 3 , src4 = (int)seam.size();
				double theta = M_PI/symmetryOrder;
				Point3D< double > v[] = { Point3D< double >( 0 , 0 , 1 ) , Point3D< double >( 1 , 0 , 0 ) , Point3D< double >( cos(theta) , sin(theta) , 0 ) };
				oMesh.vInfo[ seam[src1] ].fixedAxis = G->fixedAxisIndex( v[0] );
				oMesh.vInfo[ seam[src2] ].fixedAxis = G->fixedAxisIndex( v[1] );
				oMesh.vInfo[ seam[src3] ].fixedAxis = G->fixedAxisIndex( v[2] );
				for( int i=src1+1 ; i<src2 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = G->fixedPlaneIndex( Point3D< double >::CrossProduct( v[0] , v[1] ) );
				for( int i=src2+1 ; i<src3 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = G->fixedPlaneIndex( Point3D< double >::CrossProduct( v[1] , v[2] ) );
				for( int i=src3+1 ; i<src4 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = G->fixedPlaneIndex( Point3D< double >::CrossProduct( v[2] , v[0] ) );

				vData->seamVertices[0] = seam[src2] , vData->seamVertices[1] = seam[src3] , vData->seamVertices[2] = seam[src1] , vData->seamVertices[3] = -1;
			}
			else if( symmetryType==SYMMETRY_TETRAHEDRAL || symmetryType==SYMMETRY_OCTAHEDRAL || symmetryType==SYMMETRY_ICOSAHEDRAL )
			{
				int src1 = 0 , src2 = (int)seam.size()/3 , src3 = ( (int)seam.size()*2 ) / 3 , src4 = (int)seam.size();
				const PlatonicSymmetryGroup< true >& PG = *(PlatonicSymmetryGroup< true >*)G;

				Point3D< double > faceCenter = PG.vertices[ PG.triangles[0][0] ] + PG.vertices[ PG.triangles[0][1] ] + PG.vertices[ PG.triangles[0][2] ];
				Point3D< double > edgeCenter = PG.vertices[ PG.triangles[0][0] ] + PG.vertices[ PG.triangles[0][1] ];
				Point3D< double >     corner = PG.vertices[ PG.triangles[0][0] ];
				oMesh.vInfo[ seam[src1] ].fixedAxis = PG.fixedAxisIndex( faceCenter );
				oMesh.vInfo[ seam[src2] ].fixedAxis = PG.fixedAxisIndex( edgeCenter );
				oMesh.vInfo[ seam[src3] ].fixedAxis = PG.fixedAxisIndex(     corner );
				for( int i=src1+1 ; i<src2 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = PG.fixedPlaneIndex( Point3D< double >::CrossProduct( faceCenter , edgeCenter ) );
				for( int i=src2+1 ; i<src3 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = PG.fixedPlaneIndex( Point3D< double >::CrossProduct( edgeCenter ,     corner ) );
				for( int i=src3+1 ; i<src4 ; i++ ) oMesh.vInfo[ seam[i] ].fixedPlane = PG.fixedPlaneIndex( Point3D< double >::CrossProduct(     corner , faceCenter ) );

				vData->seamVertices[0] = seam[src1] , vData->seamVertices[1] = seam[src2] , vData->seamVertices[2] = seam[src3] , vData->seamVertices[3] = -1;
			}
		}
		else
		{
			int source , midIndex;
			if( Cones.count==0 || Cones.count==1 || Cones.count==3 )
			{
				int target1 , target2;
				// Get the source
				switch( Cones.count )
				{
					case 0: source = rand() % oMesh.vertices.size() ; break;
					case 1: source = Cones.values[0]                ; break;
					case 3: source = Cones.values[1]                ; break;
					default: fprintf( stderr , "[ERROR] Should not be here\n" );
				}
				if( source<0 || source>=oMesh.vertices.size() ) fprintf( stderr , "[ERROR] Soure index out of bounds: %d < %d\n" , source , (int)oMesh.vertices.size() ) , exit( 0 );
				std::vector< DijkstraVertex > dVerts = RunDijkstra( oMesh.triangles , oMesh.vertices , source );

				// If the end-points are not given, compute them
				if( Cones.count!=3 )
				{
					double maxValue = 0;
					for( int i=0 ; i<dVerts.size() ; i++ ) if( dVerts[i].dist>maxValue ) maxValue = dVerts[i].dist , target1 = i;
					{
						std::vector< DijkstraVertex > _dVerts = RunDijkstra( oMesh.triangles , oMesh.vertices , target1 );
						double maxValue = 0;
						for( int i=0 ; i<dVerts.size() ; i++ ) if( _dVerts[i].dist>maxValue ) maxValue = dVerts[i].dist , target2 = i;
					}
				}

				// Generate the shortest paths from the source to the end-points and create the seam
				std::vector< unsigned int > _seam;

				seam.push_back( target1 );
				while( dVerts[ target1 ].prev!=-1 ){ target1 = dVerts[target1].prev ; seam.push_back( target1 ); }
				_seam.push_back( target2 );
				while( dVerts[ target2 ].prev!=-1 ){ target2 = dVerts[target2].prev ; _seam.push_back( target2 ); }
				while( seam[ seam.size()-2 ]==_seam[ _seam.size()-2 ] ) seam.pop_back() , _seam.pop_back();
				source = seam.back();
				midIndex = (int)seam.size()-1;
				for( int i=(int)_seam.size()-2 ; i>=0 ; i-- ) seam.push_back( _seam[i] );
			}
			else if( Cones.count==2 )
			{
				int target1 = Cones.values[0] , target2 = Cones.values[1];
				std::vector< DijkstraVertex > dVerts = RunDijkstra( oMesh.triangles , oMesh.vertices , target1 );
				seam.push_back( target2 );
				while( dVerts[ target2 ].prev!=-1 ){ target2 = dVerts[target2].prev ; seam.push_back( target2 ); }
				midIndex = (int)seam.size()/2;
				source = seam[ midIndex ];
			}
			else fprintf( stderr , "[ERROR] Expected beteween 1 and 3 cones\n" ) , exit( 0 );

			if( !symmetryType==SYMMETRY_CYCLIC || symmetryOrder>1 || !Simple.set ) oMesh.split( seam );
			if( symmetryType==SYMMETRY_CYCLIC )
			{
				if( symmetryOrder==1 && Simple.set ) vData->seamVertices[0] = vData->seamVertices[1] = vData->seamVertices[2] = vData->seamVertices[3] = -1;
				else
				{
					for( int i=1 ; i<seam.size()-1 ; i++ ) oMesh.vInfo[ oMesh.vInfo[ seam[i] ].opposite ].symmetry = 1;
					// If the symmetryOrder==1, every axis is a fixed axis
					if( symmetryOrder>1 )
					{
						oMesh.vInfo[ seam[0] ].fixedAxis = G->fixedAxisIndex( Point3D< double >( 0 , 0 , 1 ) );
						oMesh.vInfo[ seam.back() ].fixedAxis = G->fixedAxisIndex( Point3D< double >( 0 , 0 , -1 ) );
					}
					vData->seamVertices[0] = seam[0] , vData->seamVertices[1] = seam.back() , vData->seamVertices[2] = vData->seamVertices[3] = -1;
				}
			}
			else if( symmetryType==SYMMETRY_DIHEDRAL )
			{
				for( int i=1 ; i<=midIndex ; i++ )            oMesh.vInfo[ oMesh.vInfo[ seam[i] ].opposite ].symmetry = 1;
				for( int i=midIndex ; i<seam.size()-1 ; i++ ) oMesh.vInfo[ oMesh.vInfo[ seam[i] ].opposite ].symmetry = 3;
				double theta = M_PI/symmetryOrder;
				if( symmetryOrder>1 )
				{
					oMesh.vInfo[              source           ].fixedAxis = G->fixedAxisIndex( Point3D< double >( 0 , 0 ,  1 ) );
					oMesh.vInfo[  oMesh.vInfo[source].opposite ].fixedAxis = G->fixedAxisIndex( Point3D< double >( 0 , 0 , -1 ) );
				}
				oMesh.vInfo[     seam[0] ].fixedAxis = G->fixedAxisIndex( Point3D< double >( 1 , 0 , 0 ) );
				oMesh.vInfo[ seam.back() ].fixedAxis = G->fixedAxisIndex( Point3D< double >( cos(theta) , sin(theta) , 0 ) );

				vData->seamVertices[0] = seam[0] , vData->seamVertices[1] = seam.back() , vData->seamVertices[2] = source , vData->seamVertices[3] = oMesh.vInfo[source].opposite;
			}
			else if( symmetryType==SYMMETRY_TETRAHEDRAL || symmetryType==SYMMETRY_OCTAHEDRAL || symmetryType==SYMMETRY_ICOSAHEDRAL )
			{
				const PlatonicSymmetryGroup< false >& PG = *(PlatonicSymmetryGroup< false >*)G;
				unsigned int eIndex;
				Point3D< double > faceCenter = PG.vertices[ PG.triangles[0][0] ] + PG.vertices[ PG.triangles[0][1] ] + PG.vertices[ PG.triangles[0][2] ];
				Point3D< double > edgeCenter = PG.vertices[ PG.triangles[0][0] ] + PG.vertices[ PG.triangles[0][1] ];
				Point3D< double >     corner = PG.vertices[ PG.triangles[0][0] ];
				eIndex = PG.elementIndex( AxisAngleRotation( faceCenter , 2.*M_PI/3. ) );
				for( int i=1 ; i<midIndex ; i++ ) oMesh.vInfo[ oMesh.vInfo[ seam[i] ].opposite ].symmetry = eIndex;
				eIndex = PG.elementIndex( AxisAngleRotation( edgeCenter , M_PI ) );
				for( int i=midIndex ; i<seam.size()-1 ; i++ ) oMesh.vInfo[ oMesh.vInfo[ seam[i] ].opposite ].symmetry = eIndex;
				oMesh.vInfo[             source           ].fixedAxis = PG.fixedAxisIndex( corner );
				oMesh.vInfo[ oMesh.vInfo[source].opposite ].fixedAxis = PG.fixedAxisIndex( corner );
				oMesh.vInfo[             seam[0]          ].fixedAxis = PG.fixedAxisIndex( faceCenter );
				oMesh.vInfo[             seam.back()      ].fixedAxis = PG.fixedAxisIndex( edgeCenter );

				vData->seamVertices[0] = seam[0] , vData->seamVertices[1] = seam.back() , vData->seamVertices[2] = source , vData->seamVertices[3] = oMesh.vInfo[source].opposite;
			}
		}
	}
	vData->init( oMesh );

	char windowName[512];
	if     ( symmetryType==SYMMETRY_CYCLIC      ) sprintf( windowName , "Orbifold CMCF (Cyclic): %d" , symmetryOrder );
	else if( symmetryType==SYMMETRY_DIHEDRAL    ) sprintf( windowName , "Orbifold CMCF (Dihedral): 2*%d" , symmetryOrder );
	else if( symmetryType==SYMMETRY_TETRAHEDRAL ) sprintf( windowName , "Orbifold CMCF (Tetrahedral)" );
	else if( symmetryType==SYMMETRY_OCTAHEDRAL  ) sprintf( windowName , "Orbifold CMCF (Octahedral)" );
	else if( symmetryType==SYMMETRY_ICOSAHEDRAL ) sprintf( windowName , "Orbifold CMCF (Icosahedral)" );
	Visualization::Viewer::Run( vData , argc , argv , windowName );

	return EXIT_SUCCESS;
}
