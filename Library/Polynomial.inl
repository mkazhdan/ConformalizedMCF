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

#include <float.h>
#include <math.h>
#include <algorithm>

////////////////
// Polynomial //
////////////////
template< unsigned int Degree > Polynomial< Degree >::Polynomial( void ) { memset( coefficients , 0 , sizeof(double)*(Degree+1) ); }
template< unsigned int Degree > Polynomial< Degree >::Polynomial( double c ) { memset( coefficients , 0 , sizeof(double)*(Degree+1) ) , coefficients[0] = c; }

template< unsigned int Degree >
template< unsigned int Degree2 >
Polynomial<Degree>::Polynomial(const Polynomial<Degree2>& P){
	memset(coefficients,0,sizeof(double)*(Degree+1));
	for(int i=0;i<=Degree && i<=Degree2;i++){coefficients[i]=P.coefficients[i];}
}
template< unsigned int Degree >       double& Polynomial< Degree >::operator[]( unsigned int idx )       { return coefficients[idx]; }
template< unsigned int Degree > const double& Polynomial< Degree >::operator[]( unsigned int idx ) const { return coefficients[idx]; }

template< unsigned int Degree >
template< unsigned int Degree2 >
Polynomial<Degree>& Polynomial<Degree>::operator  = (const Polynomial<Degree2> &p){
	int d=Degree<Degree2?Degree:Degree2;
	memset(coefficients,0,sizeof(double)*(Degree+1));
	memcpy(coefficients,p.coefficients,sizeof(double)*(d+1));
	return *this;
}

template< unsigned int Degree >
Polynomial< (Degree>0 ? Degree-1 : 0) > Polynomial<Degree>::derivative(void) const
{
	Polynomial< (Degree>0 ? Degree-1 : 0) > p;
	for( int i=0 ; i<Degree ; i++ ) p.coefficients[i]=coefficients[i+1]*(i+1);
	return p;
}

template< unsigned int Degree >
Polynomial<Degree+1> Polynomial<Degree>::integral(void) const{
	Polynomial<Degree+1> p;
	p.coefficients[0]=0;
	for(int i=0;i<=Degree;i++){p.coefficients[i+1]=coefficients[i]/(i+1);}
	return p;
}
#if 1
template< > double Polynomial< 0 >::operator() ( double t ) const { return coefficients[0]; }
template< > double Polynomial< 1 >::operator() ( double t ) const { return coefficients[0] + coefficients[1]*t; }
template< > double Polynomial< 2 >::operator() ( double t ) const { return coefficients[0] + (coefficients[1] + coefficients[2]*t) * t; }
template< unsigned int Degree >
double Polynomial<Degree>::operator() ( double t ) const
{
	double temp = coefficients[ Degree ];
	for( int i=Degree ; i>0 ; i-- )
	{
		temp *= t;
		temp += coefficients[i-1];
	}
	return temp;
}
#else
template< unsigned int Degree >
double Polynomial<Degree>::operator() ( double t) const{
	double temp=1;
	double v=0;
	for(int i=0;i<=Degree;i++){
		v+=temp*coefficients[i];
		temp*=t;
	}
	return v;
}
#endif
template< unsigned int Degree >
double Polynomial< Degree >::integral( double tMin , double tMax ) const
{
	double v = 0 , t1=tMin , t2=tMax;
	for( int i=0 ; i<=Degree ; i++ )
	{
		v += coefficients[i]*(t2-t1)/(i+1);
		if( t1!=-DBL_MAX && t1!=DBL_MAX ) t1*=tMin;
		if( t2!=-DBL_MAX && t2!=DBL_MAX ) t2*=tMax;
	}
	return v;
}
template< unsigned int Degree >
int Polynomial<Degree>::operator == (const Polynomial& p) const{
	for(int i=0;i<=Degree;i++){if(coefficients[i]!=p.coefficients[i]){return 0;}}
	return 1;
}
template< unsigned int Degree >
int Polynomial<Degree>::operator != (const Polynomial& p) const{
	for(int i=0;i<=Degree;i++){if(coefficients[i]==p.coefficients[i]){return 0;}}
	return 1;
}
template< unsigned int Degree >
int Polynomial<Degree>::isZero(void) const{
	for(int i=0;i<=Degree;i++){if(coefficients[i]!=0){return 0;}}
	return 1;
}
template< unsigned int Degree >
void Polynomial<Degree>::setZero(void){memset(coefficients,0,sizeof(double)*(Degree+1));}

template< unsigned int Degree >
Polynomial<Degree>& Polynomial<Degree>::addScaled(const Polynomial& p , double s){
	for(int i=0;i<=Degree;i++){coefficients[i]+=p.coefficients[i]*s;}
	return *this;
}
template< unsigned int Degree >
Polynomial<Degree>& Polynomial<Degree>::operator += (const Polynomial<Degree>& p){
	for(int i=0;i<=Degree;i++){coefficients[i]+=p.coefficients[i];}
	return *this;
}
template< unsigned int Degree >
Polynomial<Degree>& Polynomial<Degree>::operator -= (const Polynomial<Degree>& p){
	for(int i=0;i<=Degree;i++){coefficients[i]-=p.coefficients[i];}
	return *this;
}
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::operator + (const Polynomial<Degree>& p) const{
	Polynomial q;
	for(int i=0;i<=Degree;i++){q.coefficients[i]=(coefficients[i]+p.coefficients[i]);}
	return q;
}
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::operator - (const Polynomial<Degree>& p) const{
	Polynomial q;
	for(int i=0;i<=Degree;i++)	{q.coefficients[i]=coefficients[i]-p.coefficients[i];}
	return q;
}
template< unsigned int Degree >
void Polynomial<Degree>::Scale(const Polynomial& p,double w,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p.coefficients[i]*w;}
}
template< unsigned int Degree >
void Polynomial<Degree>::AddScaled(const Polynomial& p1,double w1,const Polynomial& p2,double w2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]*w1+p2.coefficients[i]*w2;}
}
template< unsigned int Degree >
void Polynomial<Degree>::AddScaled(const Polynomial& p1,double w1,const Polynomial& p2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]*w1+p2.coefficients[i];}
}
template< unsigned int Degree >
void Polynomial<Degree>::AddScaled(const Polynomial& p1,const Polynomial& p2,double w2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]+p2.coefficients[i]*w2;}
}

template< unsigned int Degree >
void Polynomial<Degree>::Subtract(const Polynomial &p1,const Polynomial& p2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]-p2.coefficients[i];}
}
template< unsigned int Degree >
void Polynomial<Degree>::Negate(const Polynomial& in,Polynomial& out){
	out=in;
	for(int i=0;i<=Degree;i++){out.coefficients[i]=-out.coefficients[i];}
}

template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::operator - (void) const{
	Polynomial q=*this;
	for(int i=0;i<=Degree;i++){q.coefficients[i]=-q.coefficients[i];}
	return q;
}
template< unsigned int Degree >
template< unsigned int Degree2 >
Polynomial< ( Degree>Degree2 ? Degree : Degree2 ) > Polynomial<Degree>::operator + ( const Polynomial< Degree2 >& p ) const
{
	Polynomial< ( Degree>Degree2 ? Degree : Degree2 ) > q;
	for( int i=0 ; i<=Degree  ; i++ ) q.coefficients[i]  =   coefficients[i];
	for( int i=0 ; i<=Degree2 ; i++ ) q.coefficients[i] += p.coefficients[i];
	return q;
}
template< unsigned int Degree >
template< unsigned int Degree2 >
Polynomial< ( Degree>Degree2 ? Degree : Degree2 ) > Polynomial<Degree>::operator - ( const Polynomial< Degree2 >& p ) const
{
	Polynomial< ( Degree>Degree2 ? Degree : Degree2 ) > q;
	for( int i=0 ; i<=Degree  ; i++ ) q.coefficients[i]  =   coefficients[i];
	for( int i=0 ; i<=Degree2 ; i++ ) q.coefficients[i] -= p.coefficients[i];
	return q;
}
template< unsigned int Degree >
template< unsigned int Degree2 >
Polynomial<Degree+Degree2> Polynomial<Degree>::operator * (const Polynomial<Degree2>& p) const{
	Polynomial<Degree+Degree2> q;
	for(int i=0;i<=Degree;i++){for(int j=0;j<=Degree2;j++){q.coefficients[i+j]+=coefficients[i]*p.coefficients[j];}}
	return q;
}

template< unsigned int Degree >
Polynomial<Degree>& Polynomial<Degree>::operator += (double s){
	coefficients[0]+=s;
	return *this;
}
template< unsigned int Degree >
Polynomial<Degree>& Polynomial<Degree>::operator -= (double s){
	coefficients[0]-=s;
	return *this;
}
template< unsigned int Degree >
Polynomial<Degree>& Polynomial<Degree>::operator *= (double s){
	for(int i=0;i<=Degree;i++){coefficients[i]*=s;}
	return *this;
}
template< unsigned int Degree >
Polynomial<Degree>& Polynomial<Degree>::operator /= (double s){
	for(int i=0;i<=Degree;i++){coefficients[i]/=s;}
	return *this;
}
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::operator + (double s) const{
	Polynomial<Degree> q=*this;
	q.coefficients[0]+=s;
	return q;
}
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::operator - (double s) const{
	Polynomial q=*this;
	q.coefficients[0]-=s;
	return q;
}
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::operator * (double s) const{
	Polynomial q;
	for(int i=0;i<=Degree;i++){q.coefficients[i]=coefficients[i]*s;}
	return q;
}
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::operator / (double s) const{
	Polynomial q(this->degree());
	for(int i=0;i<=Degree;i++){q.coefficients[i]=coefficients[i]/s;}
	return q;
}
#if 1
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::shrink(double s) const{
	Polynomial q=*this;
	double s2=1.0;
	for(int i=0;i<=Degree;i++){
		q.coefficients[i]*=s2;
		s2 *= s;
	}
	return q;
}
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::scale(double s) const { return shrink( 1.0/s ); }
#else
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::scale(double s) const{
	Polynomial q=*this;
	double s2=1.0;
	for(int i=0;i<=Degree;i++){
		q.coefficients[i]*=s2;
		s2/=s;
	}
	return q;
}
#endif
template< unsigned int Degree >
Polynomial<Degree> Polynomial<Degree>::shift(double t) const{
	Polynomial<Degree> q;
	for(int i=0;i<=Degree;i++){
		double temp=1;
		for(int j=i;j>=0;j--){
			q.coefficients[j]+=coefficients[i]*temp;
			temp*=-t*j;
			temp/=(i-j+1);
		}
	}
	return q;
}
template< unsigned int Degree >
void Polynomial< Degree >::printnl( char x , const char* header ) const
{
	bool first = true;
	if( header ) printf( "%s" , header );
	for( int j=0 ; j<=Degree ; j++ )
	{
		char sign = coefficients[j]<0 ? '-' : ( first ? ' ' : '+' );
		double c = fabs( coefficients[j] );
		if( !c );
		else if( j==0 )
		{
			if( c==1 )             printf( " %c 1" , sign );
			else if( floor(c)==c ) printf( " %c %d" , sign , (int)c );
#if 1
			else                   printf( " %c %g" , sign , c );
#else
			else                   printf( " %c %8.6f" , sign , c );
#endif
		}
		else if( j==1 )
		{
			if( c==1 )             printf( " %c %c" , sign , x );
			else if( floor(c)==c ) printf( " %c %d %c" , sign , (int)c , x );
#if 1
			else                   printf( " %c %g %c" , sign , c , x );
#else
			else                   printf( " %c %8.6f %c" , sign , c , x );
#endif
		}
		else
		{
			if( c==1 )             printf( " %c %c^%d" , sign , x , j );
			else if( floor(c)==c ) printf( " %c %d %c^%d" , sign , (int)c , x , j );
#if 1
			else                   printf( " %c %g %c^%d" , sign , c , x , j );
#else
			else                   printf( " %c %8.6f %c^%d" , sign , c , x , j );
#endif
		}
		if( c ) first = false;
	}
	if( first ) printf( " 0" );
	printf( "\n" );
}
template< unsigned int Degree >
void Polynomial<Degree>::getSolutions(double c,std::vector<double>& roots,double EPS) const {
	double r[4][2];
	int rCount=0;
	roots.clear();
	switch(Degree){
	case 1:
		rCount=Factor(coefficients[1],coefficients[0]-c,r,EPS);
		break;
	case 2:
		rCount=Factor(coefficients[2],coefficients[1],coefficients[0]-c,r,EPS);
		break;
	case 3:
		rCount=Factor(coefficients[3],coefficients[2],coefficients[1],coefficients[0]-c,r,EPS);
		break;
//	case 4:
//		rCount=Factor(coefficients[4],coefficients[3],coefficients[2],coefficients[1],coefficients[0]-c,r,EPS);
//		break;
	default:
		printf("Can't solve polynomial of degree: %d\n",Degree);
	}
	for(int i=0;i<rCount;i++){
		if(fabs(r[i][1])<=EPS){
			roots.push_back(r[i][0]);
//printf("%d] %f\t%f\n",i,r[i][0],(*this)(r[i][0])-c);
		}
	}
}
#if 0
template< unsigned int Degree >
bool Factor( const Polynomial< Degree >& numerator , const Polynomial< Degree >& denominator , Polynomial< 0 >& quotient , Polynomial< Degree-1 >& remainder )
{
	if( denominator.coefficients[Degree] )
	{
		quotient.coefficient[0] = numerator.coefficient[0] / denominator.coefficient[0];
		remainder = numerator - denominator * quotient.coefficient[0];
		return true;
	}
	else return false;
}
template< int Degree1 , int Degree2 >
bool Factor( const Polynomial< Degree1 >& numerator , const Polynomial< Degree2 >& denominator , Polynomial< ( Degree1>Degree2 ? Degree1-Degree2 : 0 ) >& quotient , Polynomial< ( Degree2>0 ? Degree2-1 : 0 ) >& remainder )
{
	if( !denominator.coefficients[Degree2] ) return false;
	if( Degree2>Degree1 )
	{
		remainder = numerator;
		quotient.coefficients[0] = 0;
	}
	else if( Degree2==Degree1 )
	{
		quotient.coefficient[0] = numerator.coefficient[0] / denominator.coefficient[0];
		remainder = numerator - denominator * quotient.coefficient[0];
	}
	else
	{
		Polyonimal< Degree2-1 > _remainder;
		Polynomial< ( Degree1-1>Degree2 ? Degree1-1-Degree2 : 0 ) > _quotient;

		quotient.setZero();
		quotient.coefficient[Degree1-Degree2] = numerator.coefficient[Degree1] / denominator.coefficient[2];
		_remainder = numerator - denominator * quotient;
		Factor( _remainder , numerator , _quotient , remainder );
		quotient += _quotient;
	}
	return true;
}
#endif

/////////////////
// Polynomial2 //
/////////////////
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 >::Polynomial2( void ){ memset( coefficients , 0 , sizeof(double)*(Degree1+1)*(Degree2+1) ); }

template< unsigned int  Degree1 , unsigned int  Degree2 >
template< unsigned int _Degree1 , unsigned int _Degree2 >
Polynomial2< Degree1 , Degree2 >::Polynomial2( const Polynomial2< _Degree1 , _Degree2 >& p )
{
	memset( coefficients , 0 , sizeof(double)*(Degree1+1)*(Degree2+1) );
	for( int i=0 ; i<=Degree1 && i<=_Degree1 ; i++ ) for( int j=0 ; j<=Degree2 && j<=_Degree2 ; j++ ) coefficients[i][j] = p.coefficients[i][j];
}


template< unsigned int  Degree1 , unsigned int  Degree2 >
template< unsigned int _Degree1 , unsigned int _Degree2 >
Polynomial2< Degree1 , Degree2 >& Polynomial2< Degree1 , Degree2 >::operator = ( const Polynomial2< _Degree1 , _Degree2 > &p )
{
	memset( coefficients , 0 , sizeof(double)*(Degree1+1)*(Degree2+1) );
	for( int i=0 ; i<=Degree1 && i<=_Degree1 ; i++ ) for( int j=0 ; j<=Degree2 && j<=_Degree2 ; j++ ) coefficients[i][j] = p.coefficients[i][j];
	return *this;
}

template< unsigned int Degree1 , unsigned int Degree2 >
std::pair< Polynomial2< (Degree1>0 ? Degree1-1 : 0) , Degree2 > , Polynomial2< Degree1 , (Degree2>0 ? Degree2-1 : 0) > > Polynomial2< Degree1 , Degree2 >::derivative(void) const
{
	std::pair< Polynomial2< (Degree1>0 ? Degree1-1 : 0) , Degree2 > , Polynomial2< Degree1 , (Degree2>0 ? Degree2-1 : 0) > > p;
	for( int i=0 ; i< Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) p.first.coefficients[i][j] = coefficients[i+1][j] * (i+1);
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j< Degree2 ; j++ ) p.second.coefficients[i][j] = coefficients[i][j+1] * (j+1);
	return p;
}

template< > double Polynomial2< 0 , 0 >::operator() ( double s , double t ) const { return coefficients[0][0]; }
template< > double Polynomial2< 1 , 0 >::operator() ( double s , double t ) const { return coefficients[0][0] + coefficients[1][0]*s; }
template< > double Polynomial2< 0 , 1 >::operator() ( double s , double t ) const { return coefficients[0][0] + coefficients[0][1]*t; }
template< > double Polynomial2< 1 , 1 >::operator() ( double s , double t ) const { return coefficients[0][0] + coefficients[1][0]*s + coefficients[0][1]*t + coefficients[1][1]*s*t; }
template< unsigned int Degree1 , unsigned int Degree2 >
double Polynomial2< Degree1 , Degree2 >::operator() ( double s , double t ) const
{
	double value = 0 , _t = 1;
	for( int j=0 ; j<=Degree2 ; j++ )
	{
		double temp = coefficients[ Degree1 ][j];
		for( int i=Degree1 ; i>0 ; i-- )
		{
			temp *= s;
			temp += coefficients[i-1][j];
		}
		value += temp * _t;
		_t *= t;
	}
	return value;
}
template< unsigned int Degree1 , unsigned int Degree2 >
double Polynomial2< Degree1 , Degree2 >::integral( double sMin , double sMax , double tMin , double tMax ) const
{
	double v = 0 , s1 = sMin , s2 = sMax;
	for( int i=0 ; i<=Degree1 ; i++ )
	{
		double t1=tMin , t2=tMax;
		for( int j=0 ; j<=Degree2 ; j++ )
		{
			v += coefficients[i][j]*(t2-t1)*(s2-s1)/( (i+1) * (j+1) );
			t1 *= tMin , t2 *= tMax;
		}
		s1 *= sMin , s2 *= sMax;
	}
	return v;
}
template< unsigned int Degree1 , unsigned int Degree2 >
bool Polynomial2< Degree1 , Degree2 >::operator == ( const Polynomial2& p ) const
{
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) if( coefficients[i][j]!=p.coefficients[i][j] ) return false;
	return true;
}
template< unsigned int Degree1 , unsigned int Degree2 >
bool Polynomial2< Degree1 , Degree2 >::operator != ( const Polynomial2& p ) const
{
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) if( coefficients[i][j]!=p.coefficients[i] ) return true;
	return false;
}
template< unsigned int Degree1 , unsigned int Degree2 >
bool Polynomial2< Degree1 , Degree2 >::isZero(void) const{
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) if( coefficients[i] ) return false;
	return true;
}
template< unsigned int Degree1 , unsigned int Degree2 >
void Polynomial2< Degree1 , Degree2 >::setZero( void ){ memset( coefficients , 0 , sizeof(double)*(Degree1+1)*(Degree2+1) ); }

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 >& Polynomial2< Degree1 , Degree2 >::operator += ( const Polynomial2< Degree1 , Degree2 >& p )
{
	for( int i=0 ; i<=Degree1 ; i++ )for( int j=0 ; j<=Degree2 ; j++ ) coefficients[i][j] += p.coefficients[i][j];
	return *this;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 >& Polynomial2< Degree1 , Degree2 >::operator -= ( const Polynomial2< Degree1 , Degree2 >& p )
{
	for( int i=0 ; i<=Degree1 ; i++ )for( int j=0 ; j<=Degree2 ; j++ ) coefficients[i][j] -= p.coefficients[i][j];
	return *this;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::operator + ( const Polynomial2< Degree1 , Degree2 >& p ) const
{
	Polynomial2 q;
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) q.coefficients[i][j] = coefficients[i][j] + p.coefficients[i][j];
	return q;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::operator - ( const Polynomial2< Degree1 , Degree2 >& p ) const
{
	Polynomial2 q;
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) q.coefficients[i][j] = coefficients[i][j] - p.coefficients[i][j];
	return q;
}

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::operator - ( void ) const
{
	Polynomial2 q = *this;
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) q.coefficients[i][j] = -q.coefficients[i][j];
	return q;
}
template< unsigned int  Degree1 , unsigned int  Degree2 >
template< unsigned int _Degree1 , unsigned int _Degree2 >
Polynomial2< Degree1+_Degree1 , Degree2+_Degree2 > Polynomial2< Degree1 , Degree2 >::operator * ( const Polynomial2< _Degree1 , _Degree2 >& p ) const
{
	Polynomial2< Degree1+_Degree1 , Degree2+_Degree2> q;
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ )
		for( int ii=0 ; ii<=_Degree1 ; ii++ ) for( int jj=0 ; jj<=_Degree2 ; jj++ )
			q.coefficients[i+ii][j+jj] += coefficients[i][j] * p.coefficients[ii][jj];
	return q;
}

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 >& Polynomial2< Degree1 , Degree2 >::operator += ( double s )
{
	coefficients[0][0] += s;
	return *this;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 >& Polynomial2< Degree1 , Degree2 >::operator -= ( double s )
{
	coefficients[0][0] -= s;
	return *this;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 >& Polynomial2< Degree1 , Degree2 >::operator *= ( double s )
{
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) coefficients[i][j] *= s;
	return *this;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 >& Polynomial2< Degree1 , Degree2 >::operator /= ( double s )
{
	double _s = 1./s;
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) coefficients[i][j] *= _s;
	return *this;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::operator + ( double s ) const
{
	Polynomial2< Degree1 , Degree2 > q = *this;
	q.coefficients[0][0] += s;
	return q;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::operator - ( double s ) const
{
	Polynomial2< Degree1 , Degree2 > q = *this;
	q.coefficients[0][0] -= s;
	return q;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::operator * ( double s ) const
{
	Polynomial2< Degree1 , Degree2 > q;
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) q.coefficients[i][j] = coefficients[i][j] * s;
	return q;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::operator / ( double s ) const
{
	double _s = 1./s;
	Polynomial2< Degree1 , Degree2 > q;
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ ) q.coefficients[i][j] = coefficients[i][j] * _s;
	return q;
}

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1+Degree2 , Degree1+Degree2 > Polynomial2< Degree1 , Degree2 >::xform( const double* M ) const
{
	Polynomial2< Degree1+Degree2 , Degree1+Degree2 > q;
	double det = M[0]*M[3] - M[1]*M[2];
	double M_inv[] = { M[3]/det , -M[1]/det , -M[2]/det , M[0]/det };
	Polynomial2< 1 , 1 > X , Y;
	X.coefficients[1][0] = M_inv[0] , X.coefficients[0][1] = M_inv[1];
	Y.coefficients[1][0] = M_inv[2] , Y.coefficients[0][1] = M_inv[3];

	for( int d1=0 ; d1<=Degree1 ; d1++ )
	{
		Polynomial2< Degree1+Degree2 , Degree1+Degree2 > p1 ; p1.coefficients[0][0] = 1;
		for( int d2=0 ; d2<=Degree2 ; d2++ )
		{
			Polynomial2< Degree1+Degree2 , Degree1+Degree2 > p2 ; p2.coefficients[0][0] = 1;
			q += p1 * p2 * coefficients[d1][d2];
			p2 = p2 * Y;
		}
		p1 = p1 * X;
	}
	return q;
}
template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::scale( double s , double t ) const
{
	Polynomial2 q = *this;
	double _s = 1.0;
	for( int i=0 ; i<=Degree1 ; i++ )
	{
		double _t = _s;
		for( int j=0 ; j<=Degree2 ; j++ )
		{
			q.coefficients[i][j] *= _t;
			_t /= t;
		}
		_s /= s;
	}
	return q;
}

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial2< Degree1 , Degree2 > Polynomial2< Degree1 , Degree2 >::shift( double s , double t ) const
{
	Polynomial2 q;
	for( int j=0 ; j<=Degree2 ; j++ ) for( int i=0 ; i<=Degree1 ; i++ )
	{
		double temp=1;
		for( int k=i ; k>=0 ; k-- )
		{
			q.coefficients[k][j] += coefficients[i][j]*temp;
			temp *= -s*k;
			temp /= (i-k+1);
		}
	}
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ )
	{
		double temp=1;
		for( int k=j ; k>=0 ; k-- )
		{
			q.coefficients[i][k] += coefficients[i][j]*temp;
			temp *= -t*k;
			temp /= (j-k+1);
		}
	}
	return q;
}

template< unsigned int Degree1 , unsigned int Degree2 >
void Polynomial2< Degree1 , Degree2 >::printnl( char x , char y , const char* header ) const
{
	bool first = true;
	if( header ) printf( "%s" , header );
	for( int i=0 ; i<=Degree1 ; i++ ) for( int j=0 ; j<=Degree2 ; j++ )
	{
		char sign = coefficients[i][j]<0 ? '-' : ( first ? ' ' : '+' );
		double c = fabs( coefficients[i][j] );
		if( !c );
		else
		{
			printf( " %c" , sign );
			if( i==0 && j==0 )
			{
				if( c==1 )             printf( " 1" );
				else if( floor(c)==c ) printf( " %d" , (int)c );
				else                   printf( " %g" , c );
			}
			else
			{
				if( c==1 );
				else if( floor(c)==c ) printf( " %d" , (int)c );
				else                   printf( " %g" , c );

				if( i==0 );
				else if( i==1 ) printf( " %c" , x );
				else            printf( " %c^%d" , x , i );

				if( j==0 );
				else if( j==1 ) printf( " %c" , y );
				else            printf( " %c^%d" , y , j );
			}
		}
		if( c ) first = false;
	}
	if( first ) printf( " 0" );
	printf( "\n" );
}
