/*
Copyright (c) 2012, Michael Kazhdan
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
#include <omp.h>
#include "Library/CmdLineParser.h"
#include "Library/Ply.h"
#include "Library/SparseMatrix.h"
#include "Library/HalfEdge.h"
#include "Library/LinearSolvers.h"
#include "Library/MemoryUsage.h"
#include "Library/MeshStuff.h"
#include "Library/Timer.h"

const int CONFORMAL = 1<<0;
const int AUTHALIC  = 1<<1;

enum
{
	FLOW_TRADITIONAL,
	FLOW_CONFORMAL,
	FLOW_LINEAR,
	FLOW_COUNT
};
cmdLineParameter< char* > In( "in" ) , OutHeader( "outHeader" ) , Log( "log" ) , XForm( "xForm" );
cmdLineParameter< int > Threads( "threads" , omp_get_num_procs() ) , FlowType( "flow" , FLOW_CONFORMAL+1 );
cmdLineParameter< float > StepSize( "stepSize" , 0.0001f ) , Multiplier( "multiplier" , 1.f );
cmdLineSequence< int > Steps( "steps" );
cmdLineReadable Verbose( "verbose" );

void Usage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <flow type>=%d]\n" , FlowType.name , FlowType.value );
	printf( "\t\t%d] Traditional MCF\n"   , FLOW_TRADITIONAL+1 );
	printf( "\t\t%d] Conformalized MCF\n" , FLOW_CONFORMAL  +1 );
	printf( "\t\t%d] Linearized MCF\n"    , FLOW_LINEAR     +1 );
	printf( "\t[--%s <start:increment:end> or <start:end> or <start>]\n" , Steps.name );
	printf( "\t[--%s <step size>=%f]\n" , StepSize.name , StepSize.value );
	printf( "\t[--%s <output mesh header>]\n" , OutHeader.name );
	printf( "\t[--%s]\n" , Verbose.name );
}
cmdLineReadable* params[] = { &In , &OutHeader , &FlowType , &Steps , &StepSize , &Threads , & Verbose , &Log , &Multiplier , &XForm , NULL };


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
struct LogData
{
	double sphericalError;
	double conformalRatio;
	double deformationScale;
	LogData( void ) { sphericalError = conformalRatio = deformationScale = -1.; }
};

template< class Real >
void PoseMesh( const EmptyHEMesh& mesh , std::vector< Point3D< Real > >& vertices , Point3D< Real >& center , Real& scale )
{
	center = AreaCenter< Real , Point3D< Real > , EmptyHEMesh >( &vertices[0] , mesh );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] -= center;
	scale = 0;
	for( int i=0 ; i<vertices.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) scale = std::max< Real >( scale , (Real)fabs(vertices[i][j]) );
	scale *= 2;
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] /= scale;
}
template< class Real >
void PoseMesh( const EmptyHEMesh& mesh , std::vector< Point3D< Real > >& vertices )
{
	Point3D< Real > center;
	Real scale;
	PoseMesh( mesh , vertices , center , scale );
}

template< class Real >
void OutputMesh
(
 char* fileName ,
 const EmptyHEMesh& mesh ,
 const std::vector< Point3D< Real > >& vertices ,
 const std::vector< Point3D< float > >& colors ,
 Point3D< Real > center , Real scale ,
 bool rescale
)
{
	std::vector< TriangleIndex > triangles( mesh.facet_size() );
	for( int i=0 ; i<mesh.facet_size() ; i++ )
	{
		EmptyHEMesh::Facet_const_handle f = mesh.facet(i);
		int v[] = { (int)mesh.index( f.halfedge().vertex() ) , (int)mesh.index( f.halfedge().next().vertex() ) , (int)mesh.index( f.halfedge().next().next().vertex() ) };
		for( int j=0 ; j<3 ; j++ ) triangles[i][j] = v[j];
	}
	std::vector< PlyColorVertex< float > > outVertices( vertices.size() );
	std::vector< Point3D< Real > > y = vertices;
	if( rescale ) PoseMesh( mesh , y );
	for( int i=0 ; i<y.size() ; i++ ) y[i] = y[i] * scale + center;
	for( int i=0 ; i<outVertices.size() ; i++ )
	{
		outVertices[i].point = Point3D< float >( y[i] );
		Point< float , 4 > n;
		outVertices[i].color = colors[i];
	}
	PlyWriteTriangles( fileName , outVertices , triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
}


template< class Real , class CReal >
bool Execute
(
 int flowType ,
 EmptyHEMesh mesh ,
 std::vector< PlyColorVertex< float > >& vertices ,
 std::vector< LogData >& logData ,
 const char* header ,
 const std::vector< int >& outSteps ,
 int steps = 0
)
{
#ifdef USE_CHOLMOD
	typedef CholmodSolver Solver;
#elif defined( USE_EIGEN )
	typedef EigenSolverCholeskyLLt< Real , ConstPointer( MatrixEntry< Real , int > ) > Solver;
#else // !(USE_CHOLMOD || USE_EIGEN)
#error "No solver defined"
#endif // (USE_CHOLMOD || USE_EIGEN)

	double stepSize = StepSize.value;
	std::vector< Point3D< float > > colors( vertices.size() );
	std::vector< Point3D< Real > > b( vertices.size() ) , x( vertices.size() ) , oldX( vertices.size() ) , n( vertices.size() );
	Vector< Real > _b( b.size() ) , _x( x.size() );
	std::vector< MetricData< Real > > cData;
	for( int i=0 ; i<vertices.size() ; i++ ) x[i] = Point3D< Real >( vertices[i].point ) , colors[i] = vertices[i].color;

	for( int i=0 ; i<outSteps.size() ; i++ ) if( outSteps[i]>steps ) steps = outSteps[i];

	InitMetricData( vertices , mesh , cData );

	Point3D< Real > center;
	Real scale = 1.;
	PoseMesh( mesh , x , center , scale );

	// Check if the mesh has any boundary vertices.
	// If it does, these vertices will be fixed throughout the flow.
	// Otherwise, the mesh needs to be rescaled so that it doesn't collapse down to a point.
	std::vector< bool > isBoundaryVertex( mesh.vertex_size() , false );
	bool hasBoundary = false;
	for( int i=0 ; i<mesh.vertex_size() ; i++ )
	{
		bool isBoundary = false;
		EmptyHEMesh::Halfedge_around_vertex_const_circulator iter = mesh.halfedge_around_vertex_begin( i );
		EmptyHEMesh::Halfedge_around_vertex_const_circulator end = iter;
		do
		{
			if( !iter.facet() ) isBoundary = true;
			iter++;
		}
		while( end!=iter );
		isBoundaryVertex[i] = isBoundary;
		hasBoundary |= isBoundary;
	}
	bool rescale = !hasBoundary;
	if( Verbose.set )
	{
		if( rescale ) printf( "Mesh rescaling enabled\n" );
		else          printf( "Mesh rescaling disabled\n" );
	}
	if( rescale ) scale *= MakeUnitArea( x , mesh );

	// Initialize the matrices.
	// M is the current matrix which is initialized with the sparse matrix topology.
	// D stores the mass matrix.
	// L stores the stiffness matrix.
	// If the flow is authalic, the mass matrix is unchanged throughout the flow so
	// we can initialize it once.
	// If the flow is conformal, the stiffness matrix is unchanged throughout the flow so
	// we can initialize it once.
	SparseMatrix< Real , int > D , L , M;
	GetMatrices< Real , CReal >( x , mesh , NULL , NULL , &M , true , true );
	if     ( (flowType&CONFORMAL) && (flowType&AUTHALIC) ) GetMatrices< Real , CReal >( x , mesh , &D ,   &L , NULL , true , true );
	else if( (flowType&CONFORMAL)                        ) GetMatrices< Real , CReal >( x , mesh , NULL , &L , NULL , true , true );
	else if(                         (flowType&AUTHALIC) ) GetMatrices< Real , CReal >( x , mesh , &D , NULL , NULL , true , true );

	double radius;
	logData.resize( steps+1 );
	logData[0].sphericalError = GetSphericalVariation< Real , Point3D< Real > , EmptyHEMesh >( x , mesh , radius );
	logData[0].conformalRatio = 1.0;
	logData[0].deformationScale = 0.;

	if( header )
	{
		bool output = false;
		for( int j=0 ; j<outSteps.size() ; j++ ) if( outSteps[j]==0 ) output = true;
		if( output )
		{
			char outFileName[1024];
			sprintf( outFileName , "%s.0.ply" , header );
			OutputMesh( outFileName , mesh , x , colors , center , scale , rescale );
		}
	}
	Solver* solver = NULL;
	for( int i=0 ; i<steps ; i++ )
	{
		double t;
		double tt = Timer::Time();

		t = Timer::Time();
		// Update the matrices that are changing
		if     ( !(flowType&CONFORMAL) && !(flowType&AUTHALIC) ) GetMatrices< Real , CReal >( x , mesh , &D ,   &L , NULL , i==0 , true );
		else if( !(flowType&CONFORMAL)                         ) GetMatrices< Real , CReal >( x , mesh , NULL , &L , NULL , i==0 , true );
		else if(                          !(flowType&AUTHALIC) ) GetMatrices< Real , CReal >( x , mesh , &D , NULL , NULL , i==0 , true );

		double mTime = Timer::Time()-t;
		// Set the new constraint vector: b = D * x
		
		D.Multiply( &x[0] , &b[0] );
		// Set the system matrix: M = D + t * L
#pragma omp parallel for
		for( int j=0 ; j<D.rows ; j++ ) for( int k=0 ; k<D.rowSizes[j] ; k++ ) M[j][k].Value = D[j][k].Value + L[j][k].Value * stepSize;
		if( hasBoundary )
#pragma omp parallel for
			for( int j=0 ; j<M.rows ; j++ )
				if( isBoundaryVertex[j] )
				{
					for( int k=0 ; k<M.rowSizes[j] ; k++ ) M[j][k].Value = M[j][k].N==j ? 1 : 0;
					b[j] = x[j];
				}
				else for( int k=0 ; k<M.rowSizes[j] ; k++ ) if( isBoundaryVertex[ M[j][k].N ] ) b[j] -= x[ M[j][k].N ] * M[j][k].Value , M[j][k].Value = 0;
		stepSize *= Multiplier.value;

		t = Timer::Time();
		// If this is the first solve perform both the symbolic and numerical factorization.
		if( !solver ) solver = new Solver( M );
		// Otherwise, if the system matrix has changed, update the numerical factorization.
		else if( !(flowType&CONFORMAL) || !(flowType&AUTHALIC) ) solver->update( M );

		double sTime1 = Timer::Time()-t;

		t = Timer::Time();
		// Solve for each of the x, y, z coefficients independently
		for( int d=0 ; d<3 ; d++ )
		{
#pragma omp parallel for
			for( int j=0 ; j<b.size() ; j++ ) _b[j] = b[j][d] , _x[j] = oldX[j][d] = x[j][d];
			solver->solve( &_b[0] , &_x[0] );
#pragma omp parallel for
			for( int j=0 ; j<x.size() ; j++ ) x[j][d] = _x[j];
		}

		double sTime2 = Timer::Time()-t;
		double mem = MemoryInfo::UsageMB();
		// Re-scale/center the mash
		if( rescale ) TranslateVolumeCenterToOrigin( x , mesh ) , MakeUnitArea( x , mesh );
		double sphericalError = GetSphericalVariation< Real , Point3D< Real > , EmptyHEMesh >( x , mesh , radius );
		double conformalRatio = GetInitializedConformalRatio< Real , Point3D< Real > , EmptyHEMesh >( cData , x , mesh , true );
		double differenceNorm=0 , oldNorm=0;

		for( int j=0 ; j<D.rows ; j++ )
			for( int k=0  ; k<D.rowSizes[j] ; k++ )
			{
				int jj = D[j][k].N;
				differenceNorm += Point3D< Real >::Dot( oldX[j]-x[j] , oldX[jj]-x[jj] ) * D[j][k].Value;
				oldNorm += Point3D< Real >::Dot( oldX[j] , oldX[jj] ) * D[j][k].Value;
			}
		double deformationScale = sqrt( differenceNorm / oldNorm );

		logData[i+1].sphericalError = sphericalError;
		logData[i+1].conformalRatio = conformalRatio;
		logData[i+1].deformationScale = deformationScale;

		if( Verbose.set )
		{
			tt = Timer::Time() - tt;
			printf( "\r[%d / %d] %4.2f(s): D-Norm=%6.5f / QC-Ratio=%6.5f / R-Var=%6.5f     " , i+1 , steps , tt , deformationScale , conformalRatio , sphericalError );
		}

		if( header )
		{
			bool output = false;
			for( int j=0 ; j<outSteps.size() ; j++ ) if( outSteps[j]==i+1 ) output = true;
			if( output )
			{
				char outFileName[1024];
				sprintf( outFileName , "%s.%d.ply" , header , i+1 );
				OutputMesh( outFileName , mesh , x , colors , center , scale , rescale );
			}
		}
	}
	if( Verbose.set && steps>0 ) printf( "\n" );
	if( rescale ) PoseMesh( mesh , x );

	for( int i=0 ; i<x.size() ; i++ ) x[i] = x[i] * scale + center;
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].point = Point3D< float >( x[i] );
	if( solver ) delete solver;
	return true;
}
int main
(
 int argc ,
 char* argv[]
)
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}

	omp_set_num_threads( Threads.value > 1 ? Threads.value : 1 );

	std::vector< int > outSteps;
	{
		int start=1 , end=1 , increment=1;
		if     ( Steps.count>=3 ) start = Steps.values[0] , increment = Steps.values[1] , end = Steps.values[2];
		else if( Steps.count>=2 ) start = Steps.values[0] , end = Steps.values[1];
		else if( Steps.count>=1 ) start = end = Steps.values[0];
		for( int i=start ; i<=end ; i+=increment ) outSteps.push_back( i );
	}

	int fileType;
	std::vector< TriangleIndex > triangles;
	std::vector< PlyColorVertex< float > > vertices;
	bool propertyFlags[ PlyColorVertex< float >::ReadComponents ];
	PlyReadTriangles( In.value , vertices , triangles , PlyColorVertex< float >::ReadProperties , propertyFlags , PlyColorVertex< float >::ReadComponents , fileType );
	bool hasColor = (propertyFlags[3]||propertyFlags[6]) && (propertyFlags[4]||propertyFlags[7]) && (propertyFlags[5]||propertyFlags[8]);

	EmptyHEMesh mesh;
	mesh.SetHalfEdgeData( triangles );

	if( !hasColor )
	{
		XForm3x3< float > xForm;
		xForm *= 0;
		xForm(0,0) = xForm(1,1) = xForm(2,2) = float(1);
		if( XForm.set )
		{
			printf( "Reading transform\n" );
			FILE* fp = fopen( XForm.value , "r" );
			if( !fp ) fprintf( stderr , "[WARNING] Failed to read in xform from: %s\n" , XForm.value );
			else
			{
				float temp;
				XForm4x4< float > _xForm;
				for( int j=0 ; j<4 ; j++ ) for( int i=0 ; i<4 ; i++ )
				{
					if( fscanf( fp , " %f " , &temp )!=1 ) fprintf( stderr , "[ERROR] Failed to read matrix coefficient\n" ) , exit( 0 );
					_xForm(i,j) = float(temp);
				}
				fclose( fp );
				_xForm = _xForm.inverse();
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) xForm(i,j) = _xForm(i,j);
			}
		}

		std::vector< Point3D< float > > normals;
		GetVertexNormals( vertices , mesh , normals );
		for( int i=0 ; i<normals.size() ; i++ )
		{
			// Negating the normal to comply with earlier version.
			normals[i] = -xForm * normals[i];
			normals[i] /= (float)Length( normals[i] );
			vertices[i].color = NormalColor( normals[i] );
		}
	}

	std::vector< LogData > logData;

	int flowType = 0;
	switch( FlowType.value )
	{
	case FLOW_TRADITIONAL+1: flowType; break;
	case FLOW_CONFORMAL+1:   flowType = CONFORMAL ; break;
	case FLOW_LINEAR+1:      flowType = CONFORMAL | AUTHALIC ; break;
	default:
		fprintf( stderr , "[ERROR] Unknown flow-type: %d\n" , FlowType.value );
		Usage( argv[0] );
		return EXIT_FAILURE;
	}

	bool success = Execute< double , double >( flowType , mesh , vertices , logData , OutHeader.value , outSteps );

	if( Log.set )
	{
		FILE* fp = fopen( Log.value , "w" );
		if( !fp ) fprintf( stderr , "[WARNING] Failed to open file for writing: %s\n" , Log.value );
		else
		{
			fprintf( fp , "Iters\tSpherical-Error\tConformal-Ratio\tDeformation-Scale\n" );
			for( int i=0 ; i<logData.size() ; i++ ) fprintf( fp , "%d\t%g\t%g\t%g\n" , i , logData[i].sphericalError , logData[i].conformalRatio , logData[i].deformationScale );
			fclose( fp );
		}
	}
	return EXIT_SUCCESS;
}