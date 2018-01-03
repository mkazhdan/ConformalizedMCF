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

enum
{
	SYMMETRY_CYCLIC ,
	SYMMETRY_DIHEDRAL ,
	SYMMETRY_TETRAHEDRAL ,
	SYMMETRY_OCTAHEDRAL ,
	SYMMETRY_ICOSAHEDRAL ,
	SYMMETRY_COUNT
};

SquareMatrix< double , 3 > AxisProjection( Point3D< double > axis )
{
	axis /= Length( axis );
	SquareMatrix< double , 3 > P;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) P(i,j) = axis[i] * axis[j];
	return P;
}
SquareMatrix< double , 3 > AxisAngleRotation( Point3D< double > axis , double angle )
{
	axis /= Length( axis );
	SquareMatrix< double , 3 > M1=SquareMatrix< double , 3 >::Identity() , M2 , M3;
	M2(1,0) = -axis[2] , M2(2,0) = axis[1] , M2(2,1) = -axis[0];
	M2(0,1) = -M2(1,0) , M2(0,2) = -M2(2,0) , M2(1,2) = -M2(2,1);
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M3(i,j) = axis[i] * axis[j];
	return M1 * cos(angle) + M2 * sin(angle) + M3 * ( 1.-cos(angle) );
}

SquareMatrix< double , 3 > NormalReflection( Point3D< double > normal ){ return SquareMatrix< double , 3 >::Identity() - AxisProjection( normal ) * 2.; }

unsigned int Kernel( SquareMatrix< double , 3 > M , double epsilon2 , Point3D< double >& p )
{
	if( fabs( M.determinant() )>epsilon2 ) return 0;
	else
	{
		for( int i=0 ; i<3 ; i++ )
		{
			int _i[] = { (i+1)%3 , (i+2)%3 };
			SquareMatrix< double , 2 > A;
			Point2D< double > b;
			for( int k=0 ; k<2 ; k++ )
			{
				b[k] = - M( i , _i[k] );
				for( int j=0 ; j<2 ; j++ ) A(j,k) = M( _i[j] , _i[k] );
			}
			if( fabs( A.determinant() )>epsilon2 )
			{
				b = A.inverse() * b;
				p = Point3D< double >();
				Point3D< double > v;
				p[i] = 1 , p[ _i[0] ] = b[0] , p[ _i[1] ] = b[1];
				p /= Length( p );
				return 1;
			}
		}
		return 2;
	}
}
std::pair< Point3D< double > , Point3D< double > > GetNormalVectors( Point3D< double > n )
{
	std::pair< Point3D< double > , Point3D< double > > vs;
	n /= Length( n );
	for( int i=0 ; i<3 ; i++ )
	{
		vs.first = Point3D< double >();
		vs.first[i] = 1;
		if( fabs( Point3D< double >::Dot( n , vs.first )>0.9 ) ) continue;
		vs.second = Point3D< double >::CrossProduct( n , vs.first );
		vs.second /= Length( vs.second );
		vs.first = Point3D< double >::CrossProduct( vs.second , n );
		return vs;
	}
	fprintf( stderr , "[ERROR] GetNormalVectors: Could not generate normal vectors: %g %g %g\n" , n[0] , n[1] , n[2] ) , exit( 0 );
	return vs;
}
template< typename PrintNumber >
void PrintValue( double v , const char* header="" , bool newLine=true , PrintNumber Print=[]( double n ){ printf( "%7.4f" , n ); } )
{
	printf( "%s" , header );
	printf( " " ) , Print( v );
	if( newLine ) printf( "\n" );
}
void PrintValue( double v , const char* header="" , bool newLine=true ){ PrintValue( v , header , newLine , []( double n ){ printf( "%7.4f" , n ); } ); }

template< unsigned int Dim >
void PrintPoint( Point< double , Dim > p , const char* header="" , bool newLine=true ){ PrintPoint( p , header , newLine , []( double n ){ printf( "%7.4f" , n ); } ); }
template< unsigned int Dim , typename PrintNumber >
void PrintPoint( Point< double , Dim > p , const char* header="" , bool newLine=true , PrintNumber Print=[]( double n ){ printf( "%7.4f" , n ); } )
{
	printf( "%s" , header );
	for( int i=0 ; i<Dim ; i++ ) printf( " " ) , Print( p[i] );
	if( newLine ) printf( "\n" );
}

template< unsigned int Dim >
void PrintMatrix( SquareMatrix< double , Dim > M , const char* header="" , bool newLine=true ){ PrintMatrix( M , header , newLine , []( double n ){ printf( "%7.4f" , n ); } ); }
template< unsigned int Dim , typename PrintNumber >
void PrintMatrix( SquareMatrix< double , Dim > M , const char* header="" , bool newLine=true , PrintNumber Print=[]( double n ){ printf( "%7.4f" , n ); } )
{
	if( newLine )
	{
		for( int i=0 ; i<Dim ; i++ )
		{
			printf( "%s" , header );
			for( int j=0 ; j<Dim ; j++ ) printf( " " ) , Print( M(j,i) );
			printf( "\n" );
		}
	}
	else
	{
		printf( "%s" , header );
		for( int i=0 ; i<Dim ; i++ )
		{
			for( int j=0 ; j<Dim ; j++ ) printf( " " ) , Print( M(j,i) );
			printf( i==2 ? "\n" : "\t" );
		}
	}
}


struct SymmetryGroup
{
	virtual SquareMatrix< double , 3 > xform( unsigned int idx ) const = 0;
	virtual unsigned int size( void ) const = 0;
	virtual unsigned int inverse( unsigned int idx ) const { return elementIndex( xform(idx).inverse() ); }
	virtual unsigned int compose( unsigned int idx1 , unsigned int idx2 ) const { return elementIndex( xform(idx1) * xform(idx2) ); }

	Point3D< double > fixedAxis( unsigned int idx ) const
	{
		if( idx<_fixedAxes.size() ) return _fixedAxes[idx];
		fprintf( stderr , "[ERROR] SymmetryGroup::fixedAxis: index out of bounds: %d >= %lld\n" , idx , _fixedAxes.size() ) , exit( 0 );
		return Point3D< double >();
	}
	std::pair< Point3D< double > , Point3D< double > > fixedPlane( unsigned int idx ) const
	{
		if( idx<_fixedPlanes.size() ) return _fixedPlanes[idx];
		fprintf( stderr , "[ERROR] SymmetryGroup::fixedPlane: index out of bounds: %d >= %lld\n" , idx , _fixedPlanes.size() ) , exit( 0 );
		return std::pair< Point3D< double > , Point3D< double > >();
	}
	unsigned int elementIndex( SquareMatrix< double , 3 > R ) const
	{
		for( unsigned int g=0 ; g<size() ; g++ ) if( SquareMatrix< double , 3 >::SquareNorm( R - xform(g) )<_epsilon2 ) return g;
		fprintf( stderr , "[ERROR] SymmetryGroup::elementIndex: Couldn't match rotation to rotation list\n" );
		PrintMatrix( R );
		double err; int idx=-1;
		for( unsigned int g=0 ; g<size() ; g++ ) if( idx==-1 || SquareMatrix< double , 3 >::SquareNorm( R - xform(g) )<err ) idx = g , err = SquareMatrix< double , 3 >::SquareNorm( R - xform(g) );
		PrintMatrix( xform(idx) , "\t" );
		exit( 0 );
		return -1;
	}
	unsigned int fixedAxisIndex( Point3D< double > a ) const
	{
		a /= Length( a );
		for( int i=0 ; i<_fixedAxes.size() ; i++ )
			if( Point3D< double >::SquareNorm( a - _fixedAxes[i] )<=_epsilon2 || Point3D< double >::SquareNorm( a + _fixedAxes[i] )<=_epsilon2 ) return i;
		fprintf( stderr , "[ERROR] SymmetryGroup::fixedAxisIndex: Couldn't match axis to axis list\n" );
		exit( 0 );
		return -1;
	}
	unsigned int fixedPlaneIndex( Point3D< double > p ) const
	{
		p /= Length( p );
		for( int i=0 ; i<_fixedPlanes.size() ; i++ )
		{
			Point3D< double > n = Point3D< double >::CrossProduct( _fixedPlanes[i].first , _fixedPlanes[i].second );
			if( Point3D< double >::SquareNorm( p - n )<=_epsilon2 || Point3D< double >::SquareNorm( p + n )<=_epsilon2 ) return i;
		}
		fprintf( stderr , "[ERROR] SymmetryGroup::fixedPlaneIndex: Couldn't match plane to plane list: %g %g %g\n" , p[0] , p[1] , p[2] );
		exit( 0 );
		return -1;
	}
	size_t fixedAxisCount( void ) const { return _fixedAxes.size(); }
	size_t fixedPlaneCount( void ) const { return _fixedPlanes.size(); }

	SymmetryGroup( double epsilon2 ) : _epsilon2( epsilon2 ) {}

protected:
	std::vector< Point3D< double > > _fixedAxes;
	std::vector< std::pair< Point3D< double > , Point3D< double > > > _fixedPlanes;
	bool _fixedPoint;
	double _epsilon2;
	void _test( bool verbose=false ) const
	{
		if( verbose ) printf( "Testing (%d):\n" , size() );

		// Test identity
		if( verbose ) printf( "\tIdentity\n" );
		{
			int i = 0;
			SquareMatrix< double , 3 > g = xform(i);
			double d2 = SquareMatrix< double , 3 >::SquareNorm( g - SquareMatrix< double , 3 >::Identity() );
			if( d2>_epsilon2 ) printf( "Identity: %g\n" , d2 );
		}
		// Test inversion
		if( verbose ) printf( "\tInversion\n" );
		for( unsigned int i=0 ; i<size() ; i++ )
		{
			unsigned int i_inv = inverse( i );
			SquareMatrix< double , 3 > g = xform(i) , g_inv = xform( i_inv );
			double d2 = SquareMatrix< double , 3 >::SquareNorm( g.inverse() - g_inv );
			if( d2>_epsilon2 ) printf( "Inverse[%d] %g\n" , i , d2 );
		}
		// Test composition
		if( verbose ) printf( "\tComposition\n" );
		for( unsigned int i=0 ; i<size() ; i++ ) for( unsigned int j=0 ; j<size() ; j++ )
		{
			unsigned int ij = compose(i,j);
			SquareMatrix< double , 3 > g1 = xform(i) , g2 = xform(j) , g12 = xform(ij);
			double d2 = SquareMatrix< double , 3 >::SquareNorm( g12 - g1*g2 );
			if( d2>_epsilon2 ) printf( "Composition[%d %d] %g\n" , i , j , d2 );
		}
		// Test distinctiveness
		if( verbose ) printf( "\tDistinctiveness\n" );
		for( unsigned int i=0 ; i<size() ; i++ ) for( unsigned int j=0 ; j<i ; j++ )
		{
			SquareMatrix< double , 3 > g1 = xform(i) , g2 = xform(j);
			double d2 = SquareMatrix< double , 3 >::SquareNorm( g1 - g2 );
			if( d2<=_epsilon2 ) printf( "Distinctiveness[%d %d] %g\n" , i , j , d2 );
		}
	}
	void _init( bool verbose=false )
	{
		if( verbose ) printf( "Initializing (%d):\n" , size() );
		_fixedPoint = false;
		for( unsigned int g=1 ; g<size() ; g++ )
		{
			SquareMatrix< double , 3 > M = xform( g );
			if( M.determinant()>0 )
			{
				Point3D< double > v;
				int dim = Kernel( M - SquareMatrix< double , 3 >::Identity() , _epsilon2 , v );
				if( dim!=1 )
				{
					fprintf( stderr , "[ERROR] Could not find fixed axis: %d\n" , dim );
					PrintMatrix( M );
					exit( 0 );
				}
				v /= Length( v );
				bool found = false;
				for( int i=0 ; i<_fixedAxes.size() ; i++ )
					if( Point3D< double >::SquareNorm( v - _fixedAxes[i] )<=_epsilon2 || Point3D< double >::SquareNorm( v +  _fixedAxes[i] )<=_epsilon2 ) found = true;
				if( !found ) _fixedAxes.push_back( v );
			}
			else
			{
				Point3D< double > v;
				SquareMatrix< double , 3 > _M = M + SquareMatrix< double , 3 >::Identity();
				if( SquareMatrix< double , 3 >::SquareNorm( _M )<_epsilon2 ) _fixedPoint = true;
				else
				{
					int dim = Kernel( _M , _epsilon2 , v );
					if( dim!=1 )
					{
						fprintf( stderr , "[ERROR] Could not find fixed plane: %d\n" , dim );
						PrintMatrix( M );
						exit( 0 );
					}
					v /= Length( v );
					bool found = false;
					for( int i=0 ; i<_fixedPlanes.size() ; i++ )
					{
						Point3D< double > n = Point3D< double >::CrossProduct( _fixedPlanes[i].first , _fixedPlanes[i].second );
						if( Point3D< double >::SquareNorm( v - n )<=_epsilon2 || Point3D< double >::SquareNorm( v +  n )<=_epsilon2 ) found = true;
					}
					if( !found ) _fixedPlanes.push_back( GetNormalVectors(v) );
				}
			}
		}
		if( verbose ) printf( "Fixed axes / planes: %d %d\n" , (int)_fixedAxes.size()  , (int)_fixedPlanes.size() );
	}
};
template< bool Reflect >
struct PlatonicSymmetryGroup : public SymmetryGroup
{
	std::vector< TriangleIndex > triangles;
	std::vector< Point3D< double > > vertices;
	unsigned int size( void ) const { return (int)triangles.size() * ( Reflect ? 6 : 3 ); }

	PlatonicSymmetryGroup( double epsilon2 ) : SymmetryGroup( epsilon2 ) {}
	SquareMatrix< double , 3 > xform( unsigned int idx ) const
	{
		auto SetMatrix = [&] ( int t )
		{
			Point3D< double > v[] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
			SquareMatrix< double , 3 > M;
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M(i,j) = v[i][j];
			return M;
		};
		Point3D< double > c = vertices[ triangles[0][0] ] + vertices[ triangles[0][1] ] + vertices[ triangles[0][2] ];
		SquareMatrix< double , 3 > R1 = SetMatrix( idx/(Reflect ? 6 :3 ) ) * SetMatrix(0).inverse();
		SquareMatrix< double , 3 > R2 = AxisAngleRotation( c , ( 2. * M_PI * ( Reflect ? ( (idx%6)>>1 ) : ( idx%3 ) ) ) / 3. );
		SquareMatrix< double , 3 > Ref = NormalReflection( Point3D< double >::CrossProduct( c , vertices[ triangles[0][0] ] ) );
		if( Reflect && idx&1 ) return R1*R2*Ref;
		else                   return R1*R2;
	}
};
template< bool Reflect >
struct CyclicSymmetryGroup : public SymmetryGroup
{
	CyclicSymmetryGroup( unsigned int N , bool verbose=false , double epsilon2=1e-8 ) : SymmetryGroup( epsilon2 ) , _N(N) { _test(verbose) , _init(verbose); }
	SquareMatrix< double , 3 > xform( unsigned int n ) const
	{
		if( !Reflect ) n <<= 1;
		SquareMatrix< double , 3 > Rot = AxisAngleRotation( Point3D< double >( 0 , 0 , 1 ) , ( 2. * M_PI * (n>>1) ) / _N );
		SquareMatrix< double , 3 > Ref = (n&1) ? NormalReflection( Point3D< double >( 1 , 0 , 0 ) ) : SquareMatrix< double , 3 >::Identity();
		return Rot * Ref;
	}
	unsigned int size( void ) const { return Reflect ? 2*_N : _N; }
protected:
	unsigned int _N;
};
template< bool Reflect >
struct DihedralSymmetryGroup : public SymmetryGroup
{
	DihedralSymmetryGroup( unsigned int N , bool verbose=false , double epsilon2=1e-8 ) : SymmetryGroup( epsilon2 ) , _N(N) { _test(verbose) , _init(verbose); }
	SquareMatrix< double , 3 > xform( unsigned int n ) const
	{
		if( !Reflect ) n <<= 1;
		SquareMatrix< double , 3 > RotZ = AxisAngleRotation( Point3D< double >( 0 , 0 , 1 ) , ( 2. * M_PI * (n>>2) ) / _N );
		SquareMatrix< double , 3 > RotX = ( n&2 ) ? AxisAngleRotation( Point3D< double >( 1 , 0 , 0 ) , M_PI ) : SquareMatrix< double , 3 >::Identity();
		SquareMatrix< double , 3 > RefZ = ( n&1 ) ? NormalReflection( Point3D< double >( 0 , 0 , 1 ) ) : SquareMatrix< double , 3 >::Identity();
		return RotZ * RotX * RefZ;
	};
	unsigned int size( void ) const { return Reflect ? 4*_N : 2*_N; }
protected:
	unsigned int _N;
};

template< bool Reflect >
struct TetrahedralSymmetryGroup : public PlatonicSymmetryGroup< Reflect >
{
	TetrahedralSymmetryGroup( bool verbose , double epsilon2=1e-8 ) : PlatonicSymmetryGroup( epsilon2 )
	{
		vertices.resize( 4 ) , triangles.resize( 4 );
		vertices[0] = Point3D< double >( sqrt(8./9) , 0 , -1./3 );
		vertices[1] = Point3D< double >( -sqrt(2./9) , sqrt(2./3) , -1./3 );
		vertices[2] = Point3D< double >( -sqrt(2./9) , -sqrt(2./3) , -1./3 );
		vertices[3] = Point3D< double >( 0 , 0 , 1 );
		triangles[0] = TriangleIndex( 0 , 2 , 1 );
		triangles[1] = TriangleIndex( 0 , 1 , 3 );
		triangles[2] = TriangleIndex( 1 , 2 , 3 );
		triangles[3] = TriangleIndex( 2 , 0 , 3 );

		_test(verbose) , _init(verbose);
	}
};
template< bool Reflect >
struct OctahedralSymmetryGroup : public PlatonicSymmetryGroup< Reflect >
{
	OctahedralSymmetryGroup( bool verbose=false , double epsilon2=1e-8 ) : PlatonicSymmetryGroup( epsilon2 )
	{
		vertices.resize( 6 ) , triangles.resize( 8 );
		vertices[0] = Point3D< double >(  1 ,  0 ,  0 );
		vertices[1] = Point3D< double >(  0 ,  1 ,  0 );
		vertices[2] = Point3D< double >( -1 ,  0 ,  0 );
		vertices[3] = Point3D< double >(  0 , -1 ,  0 );
		vertices[4] = Point3D< double >(  0 ,  0 ,  1 );
		vertices[5] = Point3D< double >(  0 ,  0 , -1 );
		triangles[0] = TriangleIndex( 0 , 1 , 4 );
		triangles[1] = TriangleIndex( 1 , 2 , 4 );
		triangles[2] = TriangleIndex( 2 , 3 , 4 );
		triangles[3] = TriangleIndex( 3 , 0 , 4 );
		triangles[4] = TriangleIndex( 1 , 0 , 5 );
		triangles[5] = TriangleIndex( 2 , 1 , 5 );
		triangles[6] = TriangleIndex( 3 , 2 , 5 );
		triangles[7] = TriangleIndex( 0 , 3 , 5 );

		_test(verbose) , _init(verbose);
	}
};
template< bool Reflect >
struct IcosahedralSymmetryGroup : public PlatonicSymmetryGroup< Reflect >
{
	IcosahedralSymmetryGroup( bool verbose=false , double epsilon2=1e-8 ) : PlatonicSymmetryGroup( epsilon2 )
	{
		double x = 2./sqrt(5.) , y = 1./sqrt(5.) , theta = 2.*M_PI/10.;
		vertices.push_back( Point3D< double >( 0 , 0 ,  1 ) );
		for( int i=0 ; i<5 ; i++ ) vertices.push_back( Point3D< double >( cos((2*i+0)*theta)*x , sin((2*i+0)*theta)*x ,  y ) );
		for( int i=0 ; i<5 ; i++ ) vertices.push_back( Point3D< double >( cos((2*i+1)*theta)*x , sin((2*i+1)*theta)*x , -y ) );
		vertices.push_back( Point3D< double >( 0 , 0 , -1 ) );

		for( int i=0 ; i<5 ; i++ ) triangles.push_back( TriangleIndex( 0 , (i+0)%5 + 1 , (i+1)%5 + 1 ) );
		for( int i=0 ; i<5 ; i++ ) triangles.push_back( TriangleIndex( (i+1)%5 + 1 , (i+0)%5 + 1 , (i+0)%5 + 6 ) );
		for( int i=0 ; i<5 ; i++ ) triangles.push_back( TriangleIndex( (i+0)%5 +6 , (i+1)%5 + 6 , (i+1)%5 + 1 ) );
		for( int i=0 ; i<5 ; i++ ) triangles.push_back( TriangleIndex( 11 , (i+1)%5 + 6 , (i+0)%5 +6 ) );

		_test(verbose) , _init(verbose);
	}
};


struct VertexInfo
{
	VertexInfo( int opp=-1 , int sym=0 , int fa=-1 , int fp=-1 ) : opposite(opp) , symmetry(sym) , fixedAxis(fa) , fixedPlane(fp) {}
	int opposite , symmetry , fixedAxis , fixedPlane;
};
template< typename Vertex >
struct OrbifoldMesh
{
	std::vector< TriangleIndex > triangles;
	std::vector< VertexInfo > vInfo;
	std::vector< Vertex > vertices;

	OrbifoldMesh& operator += ( const OrbifoldMesh& tMesh );
	void split( const std::vector< unsigned int >& seam );
	void merge( unsigned int edge1 , unsigned int edge2 );
	void print( void ) const;
	void sourceFirst( void );

	unsigned int source( unsigned int v ) const
	{
		if( vInfo[v].opposite<0 || (unsigned int)vInfo[v].opposite>v ) return v;
		else return vInfo[v].opposite;
	}
};

template< typename Vertex >
void OrbifoldMesh< Vertex >::print( void ) const
{
	printf( "Triangles:\n" );
	for( int i=0 ; i<triangles.size() ; i++ ) printf( "\t%d] %d %d %d\n" , i , triangles[i][0] , triangles[i][1] , triangles[i][2] );
	printf( "Vertices:\n" );
	for( int i=0 ; i<vertices.size() ; i++ ) printf( "\t%d] opposite = %d , symmetry = %d\n" , i , vInfo[i].opposite , vInfo[i].symmetry );
}

template< typename Vertex >
void OrbifoldMesh< Vertex >::sourceFirst( void )
{
	int vMap( vertices.size() );
	int count = 0;
	for( int i=0 ; i<vertices.size() ; i++ ) if( i==source(i) ) vMap[count++] = i;
	for( int i=0 ; i<vertices.size() ; i++ ) if( i!=source(i) ) vMap[count++] = i;
	{
		std::vector< VertexInfo > _vInfo( vInfo.size() );
		for( int i=0 ; i<vInfo.size() ; i++ ) _vInfo[ vMap[i] ] = vInfo[i];
		vInfo = _vInfo;
	}
	{
		std::vector< Vertex > _vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) _vertices[ vMap[i] ] = vertices[i];
		vertices = _vertices;
	}
	for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) triangles[i][j] = vMap[ triangles[i][j] ];
	for( int i=0 ; i<vInfo.size() ; i++ ) if( vInfo[i].opposite>=0 ) vInfo[i].opposite = vMap[ vInfo[i].opposite ];
}

struct CornerIndex
{
	unsigned int t , v;
	CornerIndex( unsigned int _t=0 , unsigned int _v=0 ) : t(_t) , v(_v) {}
};
bool SetOneRing( const std::vector< TriangleIndex >& triangles , unsigned int v , std::vector< CornerIndex >& triangleOneRing )
{
	triangleOneRing.resize( 0 );

	auto Index = [&]( CornerIndex c , int d ){ return triangles[ c.t ][ (c.v+d+3)%3 ]; };

	// Get the boundary edges
	std::vector< CornerIndex > adjacent;
	for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) if( triangles[i][j]==v ) adjacent.push_back( CornerIndex( i , j ) );


	std::vector< CornerIndex > preOneRing , postOneRing;
	postOneRing.push_back( adjacent.back() );
	preOneRing.push_back( adjacent.back() );
	adjacent.pop_back();

	while( adjacent.size() )
	{
		unsigned int next = Index( postOneRing.back() , 1 );
		bool found = false;
		for( int i=0 ; i<adjacent.size() && !found ; i++ )
		{
			if( Index( adjacent[i] , -1 )==next )
			{
				postOneRing.push_back( adjacent[i] );
				adjacent[i] = adjacent.back();
				adjacent.pop_back();
				found = true;
			}
		}
		if( !found ) break;
	}
	while( adjacent.size() )
	{
		unsigned int next = Index( preOneRing.back() , -11 );
		bool found = false;
		for( int i=0 ; i<adjacent.size() && !found ; i++ )
		{
			if( Index( adjacent[i] , 1 )==next )
			{
				preOneRing.push_back( adjacent[i] );
				adjacent[i] = adjacent.back();
				adjacent.pop_back();
				found = true;
			}
		}
		if( !found ) break;
	}
	if( adjacent.size() ) fprintf( stderr , "[ERROR] SetOneRing: Could not find one ring\n" ) , exit( 0 );
	if( preOneRing.size()==1 )
	{
		triangleOneRing = postOneRing;
		return true;
	}
	else
	{
		for( int i=(int)preOneRing.size()-1 ; i>=1 ; i++ ) triangleOneRing.push_back( preOneRing[i] );
		for( int i=0 ; i<postOneRing.size() ; i++ ) triangleOneRing.push_back( postOneRing[i] );
		return false;
	}
}

template< class Vertex >
void OrbifoldMesh< Vertex >::split( const std::vector< unsigned int >& seam )
{
	auto EdgeIndex = []( unsigned long long v1 , unsigned long long v2 ){ return (v1<<32) | v2; };
	auto FactorEdgeIndex = []( unsigned long long idx , unsigned int& v1 , unsigned int&v2 ){ v1 = (unsigned int)(idx>>32) , v2 =  (unsigned int)( (idx<<32)>>32 ); };
	auto Index = [&]( CornerIndex c , int d ){ return triangles[ c.t ][ (c.v+d+3)%3 ]; };

	std::vector< std::pair< CornerIndex , int > > newIndices;

	if( seam.size()<3 ) fprintf( stderr , "[ERROR] OrbifoldMesh::split: Seams must be at least three vertices long: %d\n" , (int)seam.size() ) , exit( 0 );

	for( int i=1 ; i<(int)seam.size()-1 ; i++ )
	{
		vInfo[ seam[i] ].opposite = (unsigned int)vInfo.size();

		size_t sz = vertices.size();
		vertices.push_back( vertices[ seam[i] ] );
		vInfo.push_back( VertexInfo( seam[i] ) );

		std::vector< CornerIndex > oneRing;
		if( !SetOneRing( triangles , seam[i] , oneRing ) ) fprintf( stderr , "[ERROR] OrbifoldMesh::Split: Expected seam vertex to be non-boundary\n" ) , exit( 0 );
		int start=-1 , end=-1;
		{
			int sz = (int)oneRing.size();
			for( int j=0 ; j<sz ; j++ ) if( Index( oneRing[j] , -1 )==seam[i-1] ){ start = j ; break; }
			if( start==-1 ) fprintf( stderr , "[ERROR] OrbifoldMesh::Split: Could not find seam-start in one-ring\n" ) , exit( 0 );
			for( int j=start ; j<start+sz ; j++ ) if( Index( oneRing[j%sz] , 1 )==seam[i+1] ){ end = j+1 ; break; }
			if( end==-1 ) fprintf( stderr , "[ERROR] OrbifoldMesh::Split: Could not find seam-end in one-ring\n" ) , exit( 0 );
		}
		for( int j=start ; j<end ; j++ ) newIndices.push_back( std::pair< CornerIndex , int >( oneRing[j%oneRing.size()] , (int)sz ) );
	}
	for( int i=0 ; i<newIndices.size() ; i++ ) triangles[ newIndices[i].first.t ][ newIndices[i].first.v ] = newIndices[i].second;
}
