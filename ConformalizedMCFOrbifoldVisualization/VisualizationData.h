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

#include "Library/Camera.h"
#include "Library/Visualization.h"

template< typename Real >
struct VisualizationData : public Visualization::Viewable
{
	enum
	{
		LINE_STEP_SIZE ,
		LINE_ITERATION_TIME ,
		LINE_DEFORMATION ,
		LINE_ENERGY ,
		LINE_CONFORMAL_RATIO ,
		LINE_SPHERICAL_ERROR ,
		LINE_FLIPPED_TRIANGLES ,
		LINE_COUNT
	};
	int lineStart;
	int steps = 0;
	float zNear = 0.0001f;
	float zFar  = 10.0f;

	int rotating = 0 , panning = 0 , scaling = 0;
	int oldX , oldY , newX , newY;
	Camera camera;
	bool lightOn = true;
	bool edgeMode = false;
	bool boundaryMode = false;
	bool cornerMode = false;
	double stepSize;
	bool faceMode = true;
	int iteration;
	SystemData< double >::DeformationStats dStats;

	GLfloat lightAmbient [4];
	GLfloat lightDiffuse [4];
	GLfloat lightSpecular[4];
	GLfloat shapeSpecular[4];
	GLfloat shapeSpecularShininess = 128.;

	char videoHeader[1024];
	int videoSteps = 0 , videoIter = -1 , videoFrame = 0;

	const SymmetryGroup& G;
	SystemData< Real > sData;
	int seamVertices[4];
	double scale;
	int vCount;
	char* vDataBuffer;
	Point3D< float >* vertices;
	Point3D< float >* colors;
	Point3D< float >* normals;
	std::vector< TriangleIndex > triangles;
	std::vector< int > flippedTriangles;
	bool isInverted , flipOn , toggleFlipped , sourceOnly;
	struct EdgeData
	{
		int v1 , v2;
		int t1 , t2;
	};
	std::vector< EdgeData > edges;
	GLuint vbo=0 , ebo=0;

	static void ToggleLightCallBack       ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->lightOn       = !( (VisualizationData*)v)->lightOn;       }
	static void ToggleEdgesCallBack       ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->edgeMode      = !( (VisualizationData*)v)->edgeMode;      }
	static void ToggleBoundaryCallBack    ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->boundaryMode  = !( (VisualizationData*)v)->boundaryMode;  }
	static void ToggleCornersCallBack     ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->cornerMode    = !( (VisualizationData*)v)->cornerMode;    }
	static void ToggleFacesCallBack       ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->faceMode      = !( (VisualizationData*)v)->faceMode;      }
	static void ToggleSourceTileCallBack  ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->sourceOnly    = !( (VisualizationData*)v)->sourceOnly;    }
	static void ToggleFlipCallBack        ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->toggleFlipped = !( (VisualizationData*)v)->toggleFlipped; }
	static void TogglePoseCallBack        ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->sData.noPose  = !( (VisualizationData*)v)->sData.noPose;  }
	static void SphericalNormalizeCallBack( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->sphericalNormalize(); }
	static void SetStepSizeCallBack       ( Visualization::Viewable* v , const char* prompt ){ if( prompt ) ( (VisualizationData*)v)->stepSize = atof( prompt ); }
	static void ToggleFlowCallBack        ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->steps = ( (VisualizationData*)v)->steps!=0 ? 0 : -1; }
	static void AdvanceOneCallBack        ( Visualization::Viewable* v , const char* ){ ( (VisualizationData*)v)->steps = 1; }
	static void VideoCallBack             ( Visualization::Viewable* v , const char* prompt )
	{
		VisualizationData* vd = (VisualizationData*)v;
		if( prompt )
		{
			int videoCount;
			if( sscanf( prompt , "%s %d %d" , vd->videoHeader , &videoCount , &vd->videoSteps )==3 ) vd->videoIter = videoCount * vd->videoSteps;
			else fprintf( stderr , "[WARNING] Could not parse video data: %s\n" , prompt );
		}
	}
	static void OutputMeshCallBack        ( Visualization::Viewable* v , const char* prompt )
	{
		VisualizationData* vd = (VisualizationData*)v;
		if( prompt )
		{
			if( vd->sourceOnly )
			{
				std::vector< PlyColorVertex< float > > vertices( vd->vCount );
				for( int i=0 ; i<vd->vCount ; i++ ) vertices[i].point = Point3D< float >( vd->vertices[i] ) , vertices[i].color = vd->colors[i] * 255.f;
				PlyWriteTriangles( prompt , vertices , vd->triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
			}
			else
			{
				std::vector< TriangleIndex > triangles( vd->triangles.size() * vd->G.size() );
				std::vector< PlyColorVertex< float > > vertices( vd->vCount * vd->G.size() );
				unsigned int gSize = vd->G.size();
				for( unsigned int g=0 ; g<gSize ; g++ )
				{
					SquareMatrix< double , 3 > R = vd->G.xform(g);
					for( int i=0 ; i<vd->vCount ; i++ ) vertices[i*gSize+g].point = Point3D< float >( R * Point3D< float >( vd->vertices[i] ) ) , vertices[i*gSize+g].color = vd->colors[i] * 255.f;
					for( int i=0 ; i<vd->triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) triangles[i*gSize+g][j] = vd->triangles[i][j]*gSize+g;
				}
				PlyWriteTriangles( prompt , vertices , triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
			}
		}
	}

	VisualizationData( const SymmetryGroup& _G ) : sData(_G) , G( _G )
	{
		showFPS = false;
		camera.position = - Point3D< double >( 0 , 0 , 1. ) * 2;
		lightAmbient [0] = lightAmbient [1] = lightAmbient [2] = 0.25f , lightAmbient [3] = 1.f;
		lightDiffuse [0] = lightDiffuse [1] = lightDiffuse [2] = 0.70f , lightDiffuse [3] = 1.f;
		lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 1.00f , lightSpecular[3] = 1.f;
		shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 1.00f , shapeSpecular[3] = 1.f;
		iteration = 0;
		stepSize = 0;
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'l' , "toggle light"                               , ToggleLightCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'e' , "toggle edges"                               , ToggleEdgesCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'f' , "toggle faces"                               , ToggleFacesCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'b' , "toggle boundary"                            , ToggleBoundaryCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'c' , "toggle corners"                             , ToggleCornersCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 't' , "toggle source tile"                         , ToggleSourceTileCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'T' , "toggle flip"                                , ToggleFlipCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'n' , "spherical normalize"                        , SphericalNormalizeCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , ' ' , "toggle flow"                                , ToggleFlowCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'p' , "toggle pose"                                , TogglePoseCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , '+' , "advance single"                             , AdvanceOneCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'u' , "set step size"       , "Step Size"          , SetStepSizeCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'v' , "get video"           , "Header Count Steps" , VideoCallBack ) );
		callBacks.push_back( Visualization::KeyboardCallBack( this , 'o' , "output mesh"         , "File name"          , OutputMeshCallBack ) );
	}
	void pose( void )
	{
		scale = BoundingRadius< Real , Point3D< Real > >( GetPointer( sData.tileMesh.vertices ) , (int)sData.tileMesh.vertices.size() , Point3D< Real >() );
#pragma omp parallel for
		for( int i=0 ; i<(int)sData.tileMesh.vertices.size() ; i++ ) vertices[i] = Point3D< float >( Point3D< Real >( sData.tileMesh.vertices[i] ) / scale );

		for( int i=0 ; i<(int)sData.tileMesh.vertices.size() ; i++ ) normals[i] = Point3D< float >();
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< float > v[] = { vertices[triangles[i][0]] , vertices[triangles[i][1]] , vertices[triangles[i][2]] };
			Point3D< float > n = Point3D< float >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
			for( int j=0 ; j<3 ; j++ ) normals[ triangles[i][j] ] += n;
		}
		for( int i=0 ; i<(int)sData.tileMesh.vertices.size() ; i++ ) normals[i] /= (float)Length( normals[i] );
	}

	void init( const OrbifoldMesh< PlyColorVertex< float > >& tile )
	{
		isInverted = false;
		toggleFlipped = true;
		flipOn = false;
		sourceOnly = false;

		vCount = (int)tile.vertices.size();
		triangles = tile.triangles;

		vDataBuffer = (char*)malloc( 3 * sizeof( Point3D< float > ) * vCount  );
		vertices = ( Point3D< float >* )( vDataBuffer + 0 * vCount * sizeof( Point3D< float > ) );
		normals  = ( Point3D< float >* )( vDataBuffer + 1 * vCount * sizeof( Point3D< float > ) );
		colors   = ( Point3D< float >* )( vDataBuffer + 2 * vCount * sizeof( Point3D< float > ) );

		for( int i=0 ; i<tile.vertices.size() ; i++ ) colors[i] = tile.vertices[i].color / 255.f;

		std::unordered_map< long long , EdgeData > eMap;
		for( int i=0 ; i<triangles.size() ; i++ )
			for( int j=0 ; j<3 ; j++ )
			{
				int v1 = triangles[i][j] , v2 = triangles[i][(j+1)%3];
				long long key = EdgeKey( v1 , v2 );
				std::unordered_map< long long , EdgeData >::iterator iter = eMap.find(key);
				if( iter==eMap.end() )
				{
					EdgeData eData;
					eData.v1=v1 , eData.v2=v2 , eData.t1=i , eData.t2=-1;
					eMap[key] = eData;
				}
				else iter->second.t2 = i;
			}
		edges.reserve( eMap.size() );
		for( std::unordered_map< long long , EdgeData >::iterator iter=eMap.begin() ; iter!=eMap.end() ; iter++ ) edges.push_back( iter->second );
		sData.init( tile );
		pose();
		char* str = new char[1024];
		sprintf( str , "Vertices/Triangles = %d / %d" , (int)vCount , (int)triangles.size() );
		info.push_back( str );
		lineStart = (int)info.size();
		for( int i=0 ; i<LINE_COUNT ; i++ )
		{
			str = new char[1024];
			str[0] = 0;
			info.push_back( str );
		}
	}

	void invert( bool inv )
	{
		if( inv!=isInverted )
		{
			for( int i=0 ; i<triangles.size() ; i++ ) std::swap( triangles[i][0] , triangles[i][2] );
			glDeleteBuffers( 1 , &ebo );
			glDeleteBuffers( 1,  &vbo );
			ebo = vbo = 0;
			setVBO();
			for( int i=0 ; i<vCount ; i++ ) normals[i] = -normals[i];

			glBindBuffer( GL_ARRAY_BUFFER , vbo );
			glBufferSubData( GL_ARRAY_BUFFER , 0 , 2 * sizeof( Point3D< float > ) * vCount , vDataBuffer );
			glBindBuffer( GL_ARRAY_BUFFER , 0 );
			isInverted = !isInverted;
		}
	}
	void invert( void ){ return invert( !isInverted ); }

	void setVBO( void )
	{
		if( !vbo )
		{
			glGenBuffers( 1 , &vbo );
			glBindBuffer( GL_ARRAY_BUFFER , vbo );
			glBufferData( GL_ARRAY_BUFFER , vCount * 3 * sizeof( Point3D< float > ) , vDataBuffer , GL_DYNAMIC_DRAW );
			glBindBuffer( GL_ARRAY_BUFFER , 0 );
		}

		if( !ebo )
		{
			glGenBuffers( 1 , &ebo );
			glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , ebo );
			glBufferData( GL_ELEMENT_ARRAY_BUFFER , triangles.size() * sizeof( int ) * 3 , &triangles[0] , GL_STATIC_DRAW );
			glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
		}
	}
	void sphericalNormalize( void )
	{
		sData.sphericalNormalize();
		update( false );
	}
	void update( bool updateFlow )
	{
		invert( sData.update( updateFlow ? stepSize : 0 , dStats , &flippedTriangles ) );
		pose();

		glBindBuffer( GL_ARRAY_BUFFER , vbo );
		glBufferSubData( GL_ARRAY_BUFFER , 0 , 2 * sizeof( Point3D< float > ) * vCount , vDataBuffer );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
		if( updateFlow ) iteration++;
	}
	void drawVBO( bool drawColor )
	{
		glBindBuffer( GL_ARRAY_BUFFER , vbo );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , ebo );
		glEnableClientState( GL_VERTEX_ARRAY );
		glVertexPointer( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + 0 * sizeof( Point3D< float > ) * vCount );
		glNormalPointer(     GL_FLOAT , 0 , (GLubyte*)NULL + 1 * sizeof( Point3D< float > ) * vCount );
		glColorPointer ( 3 , GL_FLOAT , 0 , (GLubyte*)NULL + 2 * sizeof( Point3D< float > ) * vCount );
		if( drawColor ) glEnableClientState ( GL_COLOR_ARRAY ); 
		else            glDisableClientState( GL_COLOR_ARRAY );
		glEnableClientState ( GL_NORMAL_ARRAY );
		if( sourceOnly ) glDrawElements( GL_TRIANGLES , (int)triangles.size()*3 , GL_UNSIGNED_INT , NULL );
		else
		{
			for( unsigned int g=0 ; g<G.size() ; g++ )
			{
				glPushMatrix();
				SquareMatrix< double , 3 > M = G.xform(g).transpose();
				double m[4][4];
				memset( m , 0 , sizeof(m) );
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) m[i][j] = M(i,j);
				m[3][3] = 1;
				glMultMatrixd( &m[0][0] );
				glDrawElements( GL_TRIANGLES , (int)triangles.size()*3 , GL_UNSIGNED_INT , NULL );
				glPopMatrix();
			}
		}
		glDisableClientState( GL_COLOR_ARRAY );
		glDisableClientState( GL_NORMAL_ARRAY );
		glBindBuffer( GL_ARRAY_BUFFER , 0 );
		glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
	}
	void drawSeamVertices( void ) const
	{
		for( int i=0 ; i<4 ; i++ ) if( seamVertices[i]>=0 )
		{
			switch( i )
			{
			case 0:  glColor3f( 1 , 0 ,  0 ) ; break;
			case 1:  glColor3f( 0 , 0 ,  1 ) ; break;
			default: glColor3f( 0 , 1 ,  0 ) ; break;
			}
			for( unsigned int g=0 ; g<( sourceOnly ? 1 : G.size() ) ; g++ )
			{
				SquareMatrix< double , 3 > R = G.xform(g);
				Point3D< float > p = Point3D< float >( R( Point3D< double >( vertices[ seamVertices[i] ] ) ) );
				glPushMatrix();
				glTranslatef( p[0] , p[1] , p[2] );
				glutSolidSphere( 0.02 , 20 , 10 );
				glPopMatrix();
			}
		}
	}
	void drawBoundary( void ) const
	{
		glDisable( GL_LIGHTING );
		glColor3f( 0 , 0 , 0 );
		glLineWidth( 2.5 );
		glEnable( GL_BLEND );
		glEnable( GL_LINE_SMOOTH );
		glBegin( GL_LINES );
		for( unsigned int g=0 ; g<( sourceOnly ? 1 : G.size() ) ; g++ )
		{
			SquareMatrix< double , 3 > R = G.xform(g);
			for( int i=0 ; i<edges.size() ; i++ )
				if( edges[i].t2==-1 )
				{
					Point3D< float > v1 = Point3D< float >( R( Point3D< double >( vertices[edges[i].v1] ) ) );
					Point3D< float > v2 = Point3D< float >( R( Point3D< double >( vertices[edges[i].v2] ) ) );
					glVertex3f( v1[0] , v1[1] , v1[2] ) , glVertex3f( v2[0] , v2[1] , v2[2] );
				}
		}
		glEnd();
	}
	void drawFlippedTriangles( void )
	{
		if( flipOn )
		{
			glDisable( GL_POLYGON_OFFSET_FILL );
			glDisable( GL_LIGHTING );
			glColor3f( 0.f , 0.f , 0.f );
			glBegin( GL_TRIANGLES );
			for( unsigned int g=0 ; g<( sourceOnly ? 1 : G.size() ) ; g++ )
			{
				SquareMatrix< double , 3 > R = G.xform(g);
				for( int i=0 ; i<flippedTriangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ )
				{
					int v = triangles[ flippedTriangles[i] ][j];
					Point3D< float > p = Point3D< float >( R( Point3D< double >( vertices[v] ) ) );
					glVertex3f( p[0] , p[1] , p[2] );
				}
			}
			glEnd();
		}
		flipOn = !flipOn;
	}
	float sphericalScale( void ) { return (float)fabs( Point3D< float >::Dot( camera.position , camera.forward ) ) - 1.f; }
	float sphericalTangentialMotion( float t ) { return sphericalScale() * t *2; }
	float sphericalNormalMotion( float  t )
	{
		// Want a function F(t) s.t.:
		// 1] F(0) = l
		// 2] F'(s) = a * ( F(s)-b )
		// [2] => F(s) = e^{as+c} + b
		// [1] => l-b = e^c
		//     => F(s) = e^{as}*(l-b) + b
		//	float speed = Length( camera.position ) - (1.+zNear);
		float speed = (float)Length( camera.position ) - (1.f+zNear);
		return exp( t ) * speed - speed;
	}
	void idle( void )
	{
		if( steps || videoIter>=0 )
		{
			update( true );
			if( steps ) steps--;
			glutPostRedisplay();
		}
		else if( toggleFlipped && flippedTriangles.size() ) glutPostRedisplay();
	}

	void mouseFunc(int button, int state, int x, int y)
	{
		newX = x; newY = y;

		if( button==GLUT_LEFT_BUTTON ) rotating = !state;
		else if( button==GLUT_RIGHT_BUTTON ) scaling = !state;
		else printf( "Button: %d\n" , button );
	}
	void mouseWheelFunc( int button , int dir , int x , int y )
	{
		float stepSize = 10.f / ( screenWidth + screenHeight );
		float scaleSize = stepSize*2;

		scaleSize *= sphericalScale();

		if( dir==1 ) camera.translate(  camera.forward*scaleSize );
		else if( dir==-1 ) camera.translate( -camera.forward*scaleSize );
		glutPostRedisplay();
	}
	void motionFunc( int x , int y )
	{
		oldX = newX; oldY = newY;
		newX = x;    newY = y;

		int imageSize = std::min< int >( screenWidth , screenHeight );
		float rel_x = (newX - oldX) / (float)imageSize;
		float rel_y = (newY - oldY) / (float)imageSize;
		float pRight = -rel_x , pUp = rel_y;
		float sForward = rel_y;
		float rRight = rel_y , rUp = rel_x;

		{
			// Rotating by an angle of t at a distance of r from the origin
			// corresponds to moving a distance of t*r in the tangential direction (for small values of t)
			rRight /= (float)Length( camera.position ) , rUp /= (float)Length( camera.position );
			pRight   = sphericalTangentialMotion( pRight );
			pUp      = sphericalTangentialMotion( pUp );
			rRight   = sphericalTangentialMotion( rRight );
			rUp      = sphericalTangentialMotion( rUp );
			sForward *= sphericalScale();
		}
		if( rotating ) camera.rotateUp( rUp ) , camera.rotateRight( rRight );
		else if( scaling ) camera.translate( camera.forward * sForward );
		else if( panning ) camera.translate( camera.right * pRight + camera.up * pUp );
		glutPostRedisplay();
	}
	void keyboardFunc( unsigned char key , int x , int y )
	{
		float speed = 1.f;
		float scale = (float)M_PI / 180.f;
		speed = sphericalTangentialMotion( speed );
		switch( key )
		{
		case 'A': camera.rotateRight  ( -22.5 *scale ) ; break;
		case 'a': camera.rotateRight  ( -speed*scale ) ; break;
		case 'Z': camera.rotateRight  (  22.5 *scale ) ; break;
		case 'z': camera.rotateRight  (  speed*scale ) ; break;
		case 'Q': camera.rotateUp     ( -22.5 *scale ) ; break;
		case 'q': camera.rotateUp     ( -speed*scale ) ; break;
		case 'W': camera.rotateUp     (  22.5 *scale ) ; break;
		case 'w': camera.rotateUp     (  speed*scale ) ; break;
		case 'S': camera.rotateForward( -22.5 *scale ) ; break;
		case 's': camera.rotateForward( - 1.  *scale ) ; break;
		case 'X': camera.rotateForward(  22.5 *scale ) ; break;
		case 'x': camera.rotateForward(   1.  *scale ) ; break;
		}
	}
	void specialFunc( int key, int x, int y )
	{
		float stepSize = 10.f / ( screenWidth + screenHeight );
		if( glutGetModifiers()&GLUT_ACTIVE_CTRL ) stepSize /= 16;
		float panSize = stepSize , scaleSize = stepSize*2;

		{
			panSize = sphericalTangentialMotion( panSize );
			scaleSize *= sphericalScale();
		}
		switch( key )
		{
		case Visualization::KEY_UPARROW:    camera.translate(  camera.forward*scaleSize ) ; break;
		case Visualization::KEY_DOWNARROW:  camera.translate( -camera.forward*scaleSize ) ; break;
		}
		glutPostRedisplay();
	}

	void display( void )
	{
		static bool firstTime = true;
		setVBO();
		firstTime = false;

		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		{
			double s = sphericalScale();
			double ar = double( screenWidth ) / screenHeight;
			double bx , by;
			if( ar<1 ) bx = s , by = s/ar;
			else       bx = s*ar , by = s;
			glOrtho( -bx , bx , -by , by , s*0.9 , s + 2.*1.1  );
		}
		glMatrixMode( GL_MODELVIEW );

		glLoadIdentity();
		camera.draw();

		GLfloat lPosition[4];

		{
			Point3D< float > d = camera.up + camera.right - camera.forward * 4;
			lPosition[0] = d[0] , lPosition[1] = d[1] , lPosition[2] = d[2];
		}
		lPosition[3] = 0.0;
		glLightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER , GL_TRUE );
		glLightModeli( GL_LIGHT_MODEL_TWO_SIDE , GL_FALSE );
		glLightfv( GL_LIGHT0 , GL_AMBIENT , lightAmbient );
		glLightfv( GL_LIGHT0 , GL_DIFFUSE , lightDiffuse );
		glLightfv( GL_LIGHT0 , GL_SPECULAR , lightSpecular );
		glLightfv( GL_LIGHT0 , GL_POSITION , lPosition );
		glEnable( GL_LIGHT0 );

		///////////////
		// Draw Here //
		///////////////
		glClearColor( 1 , 1 , 1 , 1 );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

		glDisable( GL_CULL_FACE );
		glEnable( GL_DEPTH_TEST );

		glColorMaterial( GL_FRONT_AND_BACK , GL_AMBIENT_AND_DIFFUSE );
		glEnable( GL_COLOR_MATERIAL );

		if( boundaryMode ) drawBoundary();

		glMaterialfv( GL_FRONT_AND_BACK , GL_SPECULAR  , shapeSpecular );
		glMaterialf ( GL_FRONT_AND_BACK , GL_SHININESS , shapeSpecularShininess );

		if( faceMode )
		{
			glPolygonOffset( 3.f , 3.f );
			glEnable( GL_POLYGON_OFFSET_FILL );
			if( lightOn ) glEnable( GL_LIGHTING );
			else          glDisable( GL_LIGHTING );
			glColor3f( 0.8f , 0.8f , 0.8f );

			drawVBO( true );

			if( toggleFlipped ) drawFlippedTriangles();

			glDisable( GL_POLYGON_OFFSET_FILL );
		}
		if( edgeMode )
		{
			double s;
			s = sphericalScale();
			Point3D< double > t = camera.forward * s / 256;
			glEnable( GL_POLYGON_OFFSET_FILL );
			if( lightOn ) glEnable( GL_LIGHTING );
			else          glDisable( GL_LIGHTING );
			GLint src , dst;
			glGetIntegerv( GL_BLEND_SRC , &src );
			glGetIntegerv( GL_BLEND_DST , &dst );
			glTranslatef( (float)-t[0] , (float)-t[1] , (float)-t[2] );
			if( !faceMode )
			{
				if( lightOn ) glEnable( GL_LIGHTING );
				else          glDisable( GL_LIGHTING );
			}
			glColor3f( 0.125 , 0.125 , 0.125 );
			glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
			glEnable( GL_BLEND );
			glEnable( GL_LINE_SMOOTH );
			glLineWidth( 0.25f );
			glPolygonMode( GL_FRONT_AND_BACK , GL_LINE );
			drawVBO( !faceMode );
			glPolygonMode( GL_FRONT_AND_BACK , GL_FILL );
			glDisable( GL_LINE_SMOOTH );
			glDisable( GL_BLEND );
			glBlendFunc( src , dst );
		}

		if( cornerMode )
		{
			glEnable( GL_POLYGON_OFFSET_FILL );
			if( lightOn ) glEnable( GL_LIGHTING );
			else          glDisable( GL_LIGHTING );
			drawSeamVertices();
		}

		sprintf( info [lineStart+LINE_STEP_SIZE ] , "Step-size = %.1e" , stepSize );
		sprintf( info [lineStart+LINE_FLIPPED_TRIANGLES ] , "Flipped Triangles = %d / %.1e" , dStats.flippedTriangleCount , dStats.flippedTriangleArea );
		sprintf( info [lineStart+LINE_SPHERICAL_ERROR ] , "Spherical Error = %.1e * %.1e" , dStats.sphericalError , dStats.sphericalRadius );
		sprintf( info [lineStart+LINE_CONFORMAL_RATIO ] , "Conformal Ratio = %.4f" , dStats.conformalRatio );
		sprintf( info [lineStart+LINE_ENERGY ] , "Energy = %.4e" , dStats.areaEnergy + dStats.conformalEnergy );
		sprintf( info [lineStart+LINE_DEFORMATION ] , "Deformation Avg/Max = %.1e / %.1e" , dStats.deformationScale , dStats.maxDeformation );
		sprintf( info [lineStart+LINE_ITERATION_TIME ] , "Iteration Time[%d] = %.2f (s)" , iteration , dStats.iterationTime );

		if( videoIter>=0 )
		{
			if( (videoIter%videoSteps)==0 )
			{
				char fileName[1024];
				sprintf( fileName , "%s.%d.jpg" , videoHeader , videoFrame );
				saveFrameBuffer( fileName );
				videoFrame++;
			}
			videoIter--;
			if( videoIter==-1 ) exit( 0 );
		}
	}
};
