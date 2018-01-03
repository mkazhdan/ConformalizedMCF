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

#define LUMPED_MASS


#include "Library/SparseMatrix.h"
#include "Library/LinearSolvers.h"
#include "Library/HalfEdge.h"
#include "Library/MeshStuff.h"



template< class Real >
struct SystemData
{
	struct DeformationStats
	{
		double iterationTime;
		double deformationScale;
		double maxDeformation;
		double conformalRatio;
		double areaEnergy;
		double conformalEnergy;
		double sphericalRadius;
		double sphericalError;
		double flippedTriangleArea;
		int flippedTriangleCount;
		DeformationStats( void )
		{
			iterationTime = 0.;
			deformationScale = 0.;
			maxDeformation = 0.;
			conformalRatio = 1.;
			sphericalError = 0.;
			sphericalRadius = 0.;
			conformalEnergy = 0.;
			areaEnergy = 0.;
			flippedTriangleArea = 0.;
			flippedTriangleCount = 0;
		}
	};
protected:
	std::vector< Point3D< Real > > _b , _x;
public:
	const SymmetryGroup& G;
	double stepSize;
	EmptyHEMesh heMesh;
	OrbifoldMesh< Point3D< Real > > tileMesh;
	std::vector< Point3D< Real > >  oldPositions;
	std::vector< MetricData< Real > > cData;
	SparseMatrix< SquareMatrix< double , 3 > , int > D , L;
	SparseMatrix< double , int > M;
	bool noPose;

#ifdef USE_CHOLMOD
	typedef CholmodSolver Solver;
#elif defined( USE_EIGEN )
	typedef EigenSolverCholeskyLLt< Real , ConstPointer( MatrixEntry< Real , int > ) > Solver;
#else // !(USE_CHOLMOD || USE_EIGEN)
#error "No solver defined"
#endif // (USE_CHOLMOD || USE_EIGEN)

	Solver* solver;

	SystemData( const SymmetryGroup& g ) : G(g) { noPose=false; }
	void sphericalNormalize( void )
	{
		Real r = 0;
		for( int i=0 ; i<tileMesh.vertices.size() ; i++ ) r += Point3D< Real >::SquareNorm( tileMesh.vertices[i] );
		r = (Real)sqrt( r / tileMesh.vertices.size() );
		for( int i=0 ; i<tileMesh.vertices.size() ; i++ ) tileMesh.vertices[i] *= r / (Real)sqrt( Point3D< Real >::SquareNorm( tileMesh.vertices[i] ) );
	}
	void pose( void )
	{
		if( noPose ) return;
		Point3D< double > center = AreaCenter< double , Point3D< double > , EmptyHEMesh >( &tileMesh.vertices[0] , heMesh );
		double size;
		size = sqrt( Area< double , Point3D< double > , EmptyHEMesh >( &tileMesh.vertices[0] , heMesh ) * G.size() );
		Point3D< double > c;
		for( unsigned int g=0 ; g<G.size() ; g++ ) c += G.xform(g) * center;
		center = c / G.size();
#pragma omp parallel for
		for( int i=0 ; i<heMesh.vertex_size() ; i++ ) tileMesh.vertices[i] = ( tileMesh.vertices[i]-center ) / size;
	}
	bool flippedSphericalTriangles( int& flipCount , double& flipArea , std::vector< int >* flippedTriangles )
	{
		bool invert;
		int fCount = 0 , totalCount = 0;
		double fArea = 0. , totalArea = 0.;
#pragma omp parallel for reduction( + : fCount , fArea , totalCount , totalArea )
		for( int i=0 ; i<heMesh.facet_size() ; i++ )
		{
			int idx[] = { (int)heMesh.index( heMesh.facet(i).halfedge().vertex() ) , (int)heMesh.index( heMesh.facet(i).halfedge().next().vertex() ) , (int)heMesh.index( heMesh.facet(i).halfedge().next().next().vertex() ) };
			Point3D< double > v[] = { tileMesh.vertices[idx[0]] , tileMesh.vertices[idx[1]] , tileMesh.vertices[idx[2]] };
			Point3D< double > n = Point3D< double >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
			Point3D< double > c = v[0] + v[1] + v[2];
			double area = Length( n ) / 2.;
			totalCount++ , totalArea += area;
			if( Point3D< double >::Dot( c , n )<=0 ) fCount++ , fArea += area;
		}
		if( fArea>totalArea/2 )
		{
			flipCount = totalCount - fCount;
			flipArea = totalArea - fArea;
			invert = true;
		}
		else
		{
			flipCount = fCount;
			flipArea = fArea;
			invert = false;
		}
		if( flippedTriangles )
		{
			flippedTriangles->resize( flipCount );
			int fCount = 0;
			for( int i=0 ; i<heMesh.facet_size() ; i++ )
			{
				int idx[] = { (int)heMesh.index( heMesh.facet(i).halfedge().vertex() ) , (int)heMesh.index( heMesh.facet(i).halfedge().next().vertex() ) , (int)heMesh.index( heMesh.facet(i).halfedge().next().next().vertex() ) };
				Point3D< double > v[] = { tileMesh.vertices[idx[0]] , tileMesh.vertices[idx[1]] , tileMesh.vertices[idx[2]] };
				Point3D< double > n = Point3D< double >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
				Point3D< double > c = v[0] + v[1] + v[2];
				double area = Length( n ) / 2.;
				double dot = Point3D< double >::Dot( c , n ) * ( invert ? -1 : 1 );
				if( dot<=0 ) (*flippedTriangles)[fCount++] = i;
			}
		}
		return invert;
	}

	SquareMatrix< double , 3 > GetDoFMatrix( VertexInfo vi , bool fullRank=false ) const
	{
		if( vi.fixedAxis>=0 )
		{
			SquareMatrix< double , 3 > M;
			Point3D< double > axis = G.fixedAxis( vi.fixedAxis );
			for( int i=0 ; i<3 ; i++ ) M(0,i) = axis[i];
			if( fullRank )
			{
				std::pair< Point3D< double > , Point3D< double > > axes = GetNormalVectors( axis );
				for( int i=0 ; i<3 ; i++ ) M(1,i) = axes.first[i] , M(2,i) = axes.second[i];
			}
			return M;
		}
		else if( vi.fixedPlane>=0 )
		{
			SquareMatrix< double , 3 > M;
			std::pair< Point3D< double > , Point3D< double > > axes = G.fixedPlane( vi.fixedPlane );
			for( int i=0 ; i<3 ; i++ ) M(0,i) = axes.first[i] , M(1,i) = axes.second[i];
			if( fullRank )
			{
				Point3D< double > axis = Point3D< double >::CrossProduct( axes.first , axes.second );
				for( int i=0 ; i<3 ; i++ ) M(2,i) = axis[i];
			}
			return M;
		}
		else return SquareMatrix< double , 3 >::Identity();
	}

	void symmetrizeSystemMatrix( const SparseMatrix< Real , int >& in , SparseMatrix< SquareMatrix< Real , 3 > , int >& out ) const
	{
		int dofs;
		for( dofs=0 ; dofs<tileMesh.vInfo.size() && tileMesh.source( dofs )==dofs ; dofs++ ) ;
		out.resize( dofs );

		for( int i=0 ; i<in.rows ; i++ ) out.rowSizes[ tileMesh.source( i ) ] += in.rowSizes[i];

#pragma omp parallel for
		for( int i=0 ; i<dofs ; i++ )
		{
			size_t tmp = out.rowSizes[i];
			out.rowSizes[i] = 0;
			out.SetRowSize( i , tmp );
			out.rowSizes[i] = 0;
		}

#pragma omp parallel for
		for( int i=0 ; i<in.rows ; i++ )
		{
			unsigned int src1 = tileMesh.source( i );
			for( int _j=0 ; _j<in.rowSizes[i] ; _j++ )
			{
				int j = in[i][_j].N;
				unsigned int src2 = tileMesh.source( j );

				SquareMatrix< double , 3 > R;
				{
					SquareMatrix< double , 3 > P1 = GetDoFMatrix( tileMesh.vInfo[i] );
					SquareMatrix< double , 3 > P2 = GetDoFMatrix( tileMesh.vInfo[j] );
					SquareMatrix< double , 3 > g1 = G.xform( tileMesh.vInfo[i].symmetry ) * P1;
					SquareMatrix< double , 3 > g2 = G.xform( tileMesh.vInfo[j].symmetry ) * P2;
					R += g2.transpose() * g1 * G.size();
				}
				R *= in[i][_j].Value;
				out[src1][ out.rowSizes[src1]++ ] = MatrixEntry< SquareMatrix< Real , 3 > , int >( src2 , R );
			}
		}
#pragma omp parallel for
		for( int i=0 ; i<out.rows ; i++ )
		{
			bool needsCollapse = false;
			for( int j=0 ; j<out.rowSizes[i] ; j++ ) if( tileMesh.vInfo[ out[i][j].N ].opposite!=-1 ) needsCollapse = true;
			if( needsCollapse ) out.CollapseRow(i);
		}
	}

	void setVerticesFromDoFs( void )
	{
#pragma omp parallel for
		for( int i=0 ; i<tileMesh.vInfo.size() ; i++ ) tileMesh.vertices[i] = G.xform( tileMesh.vInfo[i].symmetry ) * GetDoFMatrix( tileMesh.vInfo[ tileMesh.source( i ) ] ) * _x[ tileMesh.source( i ) ];
	}
	void setDoFsFromVertices( void )
	{
#pragma omp parallel for
		for( int i=0 ; i<tileMesh.vInfo.size() ; i++ ) if( tileMesh.source(i)==i )
			_x[i] = ( G.xform( tileMesh.vInfo[i].symmetry ) * GetDoFMatrix( tileMesh.vInfo[i] , true ) ).inverse() * tileMesh.vertices[i];
	}

	template< class Vertex >
	void init( const OrbifoldMesh< Vertex >& tile )
	{
		// Copy the tile
		{
			tileMesh.triangles = tile.triangles;
			tileMesh.vInfo = tile.vInfo;
			tileMesh.vertices.resize( tile.vertices.size() );
			for( int i=0 ; i<tileMesh.vertices.size() ; i++ ) tileMesh.vertices[i] = Point3D< Real >( tile.vertices[i] );
		}
		heMesh.SetHalfEdgeData( tileMesh.triangles );
		{
			int dofs;
			for( dofs=0 ; dofs<tileMesh.vInfo.size() && tileMesh.source( dofs )==dofs ; dofs++ ) ;
			_x.resize( dofs );
		}
		pose();
		setDoFsFromVertices();

		solver = NULL;
		_b.resize( _x.size() ) , oldPositions.resize( tileMesh.vertices.size() );
		InitMetricData( tileMesh.vertices , heMesh , cData );
		{
			SparseMatrix< Real , int > _L;
			GetMatrices< Real , Real >( tileMesh.vertices , heMesh , NULL , &_L , NULL , true , true );
			symmetrizeSystemMatrix( _L , L );
#ifdef LUMPED_MASS
			for( int i=0 ; i<L.rows ; i++ ) for( int j=0 ; j<L.rowSizes[i] ; j++ ) if( L[i][j].N==i ) std::swap( L[i][j] , L[i][0] );
#endif // LUMPED_MASS
		}
		M.resize( L.rows*3 );
		for( int i=0 ; i<L.rows ; i++ ) for( int j=0 ; j<3 ; j++ )
		{
			M.SetRowSize( 3*i+j , L.rowSizes[i]*3 );
			for( int k=0 ; k<L.rowSizes[i] ; k++ ) for( int l=0 ; l<3 ; l++ ) M[3*i+j][3*k+l].N = 3*L[i][k].N+l;
		}
		setVerticesFromDoFs();
	}
	bool update( Real stepSize , DeformationStats& dStats , std::vector< int >* flippedTriangles=NULL )
	{
		Timer timer;

		if( stepSize>0 )
		{
			static int flowIteration = 0;
			flowIteration++;
			static bool firstTime = true;

			{
				static SparseMatrix< Real , int > _D;
#ifdef LUMPED_MASS
				if( firstTime )
				{
					_D.resize( tileMesh.vertices.size() );
#pragma omp parallel for
					for( int i=0 ; i<_D.rows ; i++ )
					{
						_D.SetRowSize( i , 1 );
						_D[i][0].N = i;
					}
				}
#pragma omp parallel for
				for( int i=0 ; i<_D.rows ; i++ ) _D[i][0].Value = 0;
#pragma omp parallel for
				for( int i=0 ; i<tileMesh.triangles.size() ; i++ )
				{
					Point3D< double > v[] = { tileMesh.vertices[ tileMesh.triangles[i][0] ] , tileMesh.vertices[ tileMesh.triangles[i][1] ] , tileMesh.vertices[ tileMesh.triangles[i][2] ] };
					double mass = Length( Point3D< double >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 6;
					for( int j=0 ; j<3 ; j++ )
#pragma omp atomic
						_D[ tileMesh.source( tileMesh.triangles[i][j] ) ][0].Value += mass;
				}
#else // !LUMPED_MASS
				GetMatrices< Real , Real >( tileMesh.vertices , heMesh , &_D , NULL , NULL , firstTime , true );
#endif // LUMPED_MASS
				symmetrizeSystemMatrix( _D , D );
			}
			firstTime = false;

			setDoFsFromVertices();
			D.Multiply( ( ConstPointer( Point3D< Real > ) )GetPointer( _x ) , GetPointer( _b ) );
#pragma omp parallel for
			for( int j=0 ; j<D.rows ; j++ )
			{
				if( tileMesh.vInfo[j].fixedAxis>=0  ) for( int k=0 ; k<D.rowSizes[j] ; k++ ) if( D[j][k].N==j ) for( int l=1 ; l<3 ; l++ ) D[j][k].Value(l,l) = D[j][k].Value(0,0);
				if( tileMesh.vInfo[j].fixedPlane>=0 ) for( int k=0 ; k<D.rowSizes[j] ; k++ ) if( D[j][k].N==j ) for( int l=2 ; l<3 ; l++ ) D[j][k].Value(l,l) = D[j][k].Value(0,0);
			}

#ifdef LUMPED_MASS
#pragma omp parallel for
			for( int j=0 ; j<L.rows ; j++ ) for( int k=0 ; k<L.rowSizes[j] ; k++ )
				if( j==L[j][k].N ) for( int l=0 ; l<3 ; l++ ) for( int m=0 ; m<3 ; m++ ) M[3*j+l][3*k+m].Value = L[j][k].Value(l,m) * stepSize + D[j][0].Value(l,m);
				else               for( int l=0 ; l<3 ; l++ ) for( int m=0 ; m<3 ; m++ ) M[3*j+l][3*k+m].Value = L[j][k].Value(l,m) * stepSize;
#else // !LUMPED_MASS
#pragma omp parallel for
			for( int j=0 ; j<D.rows ; j++ ) for( int k=0 ; k<D.rowSizes[j] ; k++ ) for( int l=0 ; l<3 ; l++ ) for( int m=0 ; m<3 ; m++ )
				M[3*j+l][3*k+m].Value = D[j][k].Value(l,m) + L[j][k].Value(l,m) * stepSize;
#endif // LUMPED_MASS


			if( !solver ) solver = new Solver( M );
			solver->update( M );
#pragma omp parallel for
			for( int j=0 ; j<tileMesh.vertices.size() ; j++ ) oldPositions[j] = tileMesh.vertices[j];

			solver->solve( ( ConstPointer( Real ) )GetPointer( _b ) , ( Pointer( Real ) )GetPointer( _x ) );
			setVerticesFromDoFs();
		}
		pose();

		dStats.sphericalError  = GetSphericalVariation< Real , Point3D< Real > , EmptyHEMesh >( tileMesh.vertices , heMesh , dStats.sphericalRadius );
		dStats.conformalRatio  = GetInitializedConformalRatio< Real , Point3D< Real > , EmptyHEMesh >( cData , tileMesh.vertices , heMesh );
		dStats.conformalEnergy = GetInitializedConformalError< Real , Point3D< Real > , EmptyHEMesh >( cData , tileMesh.vertices , heMesh , false );
		dStats.areaEnergy      = GetInitializedAreaError     < Real , Point3D< Real > , EmptyHEMesh >( cData , tileMesh.vertices , heMesh , false );
		bool invert = flippedSphericalTriangles( dStats.flippedTriangleCount , dStats.flippedTriangleArea , flippedTriangles );
		dStats.areaEnergy -= 1.;

		double differenceNorm=0 , oldNorm=0;
		dStats.maxDeformation = 0.;
		for( int i=0 ; i<tileMesh.vertices.size() ; i++ )
		{
			double tmp = Point3D< Real >::SquareNorm( tileMesh.vertices[i]-oldPositions[i] );
			if( tmp>dStats.maxDeformation ) dStats.maxDeformation = tmp;
		}
		dStats.maxDeformation = sqrt(dStats.maxDeformation) / dStats.sphericalRadius;
#pragma omp parallel for reduction( + : differenceNorm , oldNorm )
		for( int j=0 ; j<D.rows ; j++ ) for( int k=0  ; k<D.rowSizes[j] ; k++ )
			{
				int jj = D[j][k].N;
				SquareMatrix< Real , 3 > _D = D[j][k].Value;
				double dNorm = Point3D< Real >::Dot( _D * ( oldPositions[j] - tileMesh.vertices[j] ) , ( oldPositions[jj] - tileMesh.vertices[jj] ) );
				double oNorm = Point3D< Real >::Dot( _D *   oldPositions[j]                          ,   oldPositions[jj]                           );
				differenceNorm += dNorm;
				oldNorm += oNorm;
			}
		dStats.deformationScale = sqrt( differenceNorm / oldNorm );
		dStats.iterationTime = timer.elapsed();
		return invert;
	}
};

