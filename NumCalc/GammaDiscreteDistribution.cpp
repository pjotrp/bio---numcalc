//
// File: GammaDiscreteDistribution.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Oct 26 20:36:12 2003
//

#include "GammaDiscreteDistribution.h"

// From Utils:
#include <Utils/MapTools.h>

// From the STL:
#include <cmath>

const double GammaDiscreteDistribution::VERYBIG = static_cast<double>(1.7E+23);
const int    GammaDiscreteDistribution::ITMAX = 100;
const double GammaDiscreteDistribution::EPS = static_cast<double>(0.0000003);
const double GammaDiscreteDistribution::FPMIN = static_cast<double>(1.0e-30);
const double GammaDiscreteDistribution::ERR_FOR_GAMMA_CALC = static_cast<double>(0.00001);

/** Constructor: **************************************************************/

GammaDiscreteDistribution::GammaDiscreteDistribution(unsigned int n, double alpha) : AbstractDiscreteDistribution()
{
	// We use a lower bound of 0.05 for alpha to prohibe errors due to computer
	// floating precision: if alpha is quite low (gamma -> constant), some classes
	// may have the same category value, leading to a classe number lower than expected.
	// NB: if this is the case, then a warning is shown. This may happen in optimization
	// algorithms.
	_alphaConstraint = new IncludingPositiveReal(0.05);
	_parameters.addParameter(Parameter("alpha", alpha, _alphaConstraint));
	applyParameters(n);
}


GammaDiscreteDistribution::~GammaDiscreteDistribution() {
	delete _alphaConstraint;
}

/******************************************************************************/

void GammaDiscreteDistribution::fireParameterChanged() {
	//cout << "Parameter changed, alpha = " << getParameter("alpha") << endl;
	applyParameters(getNumberOfCategories());	
}

/******************************************************************************/

Domain GammaDiscreteDistribution::getDomain() const
{
    return Domain(_bonderi, MapTools::getKeys<double, double, AbstractDiscreteDistribution::Order>(_distribution));
}

/******************************************************************************/

void GammaDiscreteDistribution::applyParameters(unsigned int numberOfCategories) {
	if(numberOfCategories <= 0)
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in GammaDiscreteDistribution::applyParameters()." << endl;
	_distribution.clear();
	_bonderi.clear();
	_bonderi.resize(numberOfCategories + 1);
	//cout<<"number of categories is: "<<categories()<<endl;
	//cout<<"alpha is: "<<alpha<<endl;
	if (numberOfCategories == 1) {
		_distribution[0] = 1.0;
		return;
	} else if (numberOfCategories > 1) {
		//fillMedian(numberOfCategories);
		fillMean(numberOfCategories);
		if(getNumberOfCategories() != numberOfCategories) {
			cout << "WARNING!!! Couldn't create " << numberOfCategories << " distinct categories." << endl;
			cout << "WARNING!!! This may occure if you specified a too low alpha parameter." << endl;
		}
		return ;
	} else {
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in GammaDiscreteDistribution::applyParameters()." << endl;
	}
}

/******************************************************************************/

void GammaDiscreteDistribution::fillMean(unsigned int numberOfCategories) {
	if(numberOfCategories <= 0)
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in GammaDiscreteDistribution::applyParameters()." << endl;
	fillBonderi(numberOfCategories);
	int i;
	int K = numberOfCategories;
	double alpha = getParameter("alpha");
	for (i = 0; i < K; ++i) {
		double rate = the_average_r_in_category_between_a_and_b(_bonderi[i], _bonderi[i + 1], alpha, alpha, K);
		_distribution[rate] = 1. / K;
	}
}

/******************************************************************************/

//Added by Julien Dutheil on 30/04/03;
void GammaDiscreteDistribution::fillMedian(unsigned int numberOfCategories) {
	if(numberOfCategories <= 0)
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in GammaDiscreteDistribution::applyParameters()." << endl;
	fillBonderi(numberOfCategories);
	int i;
	double t, lnga1;
	int median = 0;
	double alpha = getParameter("alpha");
	double beta = alpha;
	int K = numberOfCategories;
	double factor = K;
	Vdouble freqK(K);
	vector<double> rates;
	rates.resize(K);
	if (median) {
    	for(i = 0       ; i < K; i++ ) rates[i] = PointGamma((i * 2. + 1) / (2. * K), alpha, beta);
      	for(i = 0, t = 0; i < K; i++) t += rates[i];
      	for(i = 0       ; i < K; i++) rates[i] *= factor / t;
   	} else {
   		lnga1 = LnGamma(alpha + 1);
      	for (i = 0; i < K - 1; i++)	freqK[i] = PointGamma((i + 1.0) / K, alpha, beta);
      	for (i = 0; i < K - 1; i++)	freqK[i] = IncompleteGamma(freqK[i] * beta, alpha + 1, lnga1);
      	rates[0]     = freqK[0] * factor;
      	rates[K - 1] = (1 - freqK[K - 2]) * factor;
      	for (i = 1; i < K - 1; i++)  rates[i] = (freqK[i] - freqK[i - 1]) * factor;
   	}
	for(i = 0; i < K; i ++) _distribution[rates[i]] = 1. / K;
}

/******************************************************************************/

void GammaDiscreteDistribution::fillBonderi(unsigned int numberOfCategories) {
	if(numberOfCategories <= 0)
		cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in GammaDiscreteDistribution::applyParameters()." << endl;
	double alpha = getParameter("alpha");
	unsigned int i;
	for (i = 1; i < numberOfCategories; ++i) {
		_bonderi[i] = search_for_z_in_dis_with_any_beta(alpha, alpha, (double)i / numberOfCategories);
	}
	_bonderi[0] = 0;
	_bonderi[i] = VERYBIG;
}

/******************************************************************************/

double GammaDiscreteDistribution::the_average_r_in_category_between_a_and_b(
	double left,
	double right,
	double alpha,
	double beta,
	int k)
{
	// and and b are the border of ahoson k
	double tmp;
	tmp = gammp(alpha + 1, right * beta) - gammp(alpha + 1, left * beta);
	tmp = (tmp * alpha / beta) * k;
	return tmp;
}

/******************************************************************************/

//From Yang's PAML package:

/******************************************************************************/

double GammaDiscreteDistribution::PointNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
       points of the normal distribution.  26: 118-121.

*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) z=999;
   else {
      y = sqrt (log(1 / (p1 * p1)));
      z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   }
   return (p<0.5 ? -z : z);
}

/******************************************************************************/

double GammaDiscreteDistribution::CDFNormal (double x)
{
/* Hill ID  (1973)  The normal integral.  Applied Statistics, 22:424-427.
   Algorithm AS 66.   (error < ?)
   adapted by Z. Yang, March 1994.  Hill's routine does not look good, and I
   haven't consulted
      Adams AG  (1969)  Algorithm 39.  Areas under the normal curve.
      Computer J. 12: 197-198.
*/
    int invers=0;
    double p, limit=10, t=1.28, y=x*x/2;

    if (x<0) {  invers=1;  x=-x; }
    if (x>limit)  return (invers?0:1);
    if (x<t)
       p = .5 - x * (    .398942280444 - .399903438504 * y
                   /(y + 5.75885480458 - 29.8213557808
                   /(y + 2.62433121679 + 48.6959930692
                   /(y + 5.92885724438))));
    else
       p = 0.398942280385 * exp(-y) /
           (x - 3.8052e-8 + 1.00000615302 /
           (x + 3.98064794e-4 + 1.98615381364 /
           (x - 0.151679116635 + 5.29330324926 /
           (x + 4.8385912808 - 15.1508972451 /
           (x + 0.742380924027 + 30.789933034 /
           (x + 3.99019417011))))));
    return (invers ? p : 1-p);
}

/******************************************************************************/

double GammaDiscreteDistribution::LnGamma (double x) throw (BadNumberException)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.

   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z;

   if(x <= 0) {
      throw BadNumberException("GammaDiscreteDistribution::LnGamma. lnGamma not implemented for x<0", x);
      //if((int)x-x==0) { cout << "lnGamma undefined" << endl; return(-1); }
      //for (fneg=1; x<0; x++) fneg/=x;
      //if(fneg<0) errorMsg::reportError("strange!! check lngamma");
      //fneg=log(fneg);
   }
   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  fneg+ f + (x-0.5)*log(x) - x + .918938533204673
          + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
               +.083333333333333)/x;
}

/******************************************************************************/

double GammaDiscreteDistribution::DFGamma(double x, double alpha, double beta) throw (BadNumberException)
{
/* mean=alpha/beta; var=alpha/beta^2
*/
   if (alpha <= 0) throw BadNumberException("GammaDiscreteDistribution::DFGamma. alpha <= 0.", alpha);
   if (beta  <= 0) throw BadNumberException("GammaDiscreteDistribution::DFGamma. beta <= 0." , beta );
   if (alpha > 100) throw BadNumberException("GammaDiscreteDistribution::DFGamma. large alpha in DFGamma.", alpha);
   return pow(beta * x, alpha) / x * exp(-beta * x - LnGamma(alpha));
}

/******************************************************************************/

double GammaDiscreteDistribution::IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   /* double accurate=1e-8, overflow=1e30; */
   double accurate=1e-10, overflow=1e60;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=exp(p*log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}

/******************************************************************************/

double GammaDiscreteDistribution::PointChi2 (double prob, double v) throw (Exception)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g, small=1e-6;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<small)   return(0);
   if (p>1-small) return(9999);
   if (v<=0)      return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0)
      throw Exception("GammaDiscreteDistribution::PointChi2. IncompleteGamma.");
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}

double GammaDiscreteDistribution::PointGamma(double prob, double alpha, double beta) throw (Exception) {
	return PointChi2(prob, 2.0 * alpha) / (2.0 * beta);
}

/******************************************************************************/

double GammaDiscreteDistribution::search_for_z_in_dis_with_any_beta(double alpha,double beta, double ahoson)
{
	return (search_for_z_in_dis_with_beta_1(alpha,ahoson)/beta);
}

/******************************************************************************/

double GammaDiscreteDistribution::search_for_z_in_dis_with_beta_1(double alpha, double ahoson) throw (BadNumberException)
{
	if ( ahoson > 1 || ahoson < 0 ) throw BadNumberException("GammaDiscreteDistribution::search_for_z_in_dis_with_beta_1.", ahoson);
	double left    = 0;
	double right   = 99999.0;
	double tmp     = 5000.0;
	double results = 0.0;

	for (int i = 0; i < 100000000 ; i++){
		results = gammp(alpha, tmp);
		if (results > ahoson){
			right = tmp;
		} else left = tmp;
		tmp = (right + left) / 2;
		if (fabs(ahoson - results) < ERR_FOR_GAMMA_CALC) return tmp;
	}
	throw Exception("GammaDiscreteDistribution::search_for_z_in_dis_with_beta_1. First bonderi is 0");
	exit(-1);// also quit the program
	return 0;
}

/******************************************************************************/

void GammaDiscreteDistribution::gser(double *gamser, double a, double x, double *gln) throw (BadNumberException)
{
	//double gammln(double xx);

	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) throw BadNumberException("GammaDiscreteDistribution::gser. x less than 0.", x);
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		throw BadNumberException("GammaDiscreteDistribution::gser. a too large, ITMAX too small.", a);
		return;
	}
}

/******************************************************************************/

void GammaDiscreteDistribution::gcf(double *gammcf, double a, double x, double *gln) throw (BadNumberException)
{
	//double gammln(double xx);
	int i;
	double an, b, c, d, del, h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) throw BadNumberException("GammaDiscreteDistribution::gcf. a too large, ITMAX too small.", a);
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

#undef ITMAX
#undef EPS
#undef FPMIN

/******************************************************************************/

double GammaDiscreteDistribution::gammln(double xx)
{
	double x, y, tmp, ser;
	static double cof[6]={
		static_cast<double>(76.18009172947146),
		static_cast<double>(-86.50532032941677),
		static_cast<double>(24.01409824083091),
		static_cast<double>(-1.231739572450155),
		static_cast<double>(0.1208650973866179e-2),
		static_cast<double>(-0.5395239384953e-5)
	};
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015f;
	for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005 * ser / x);
}

/******************************************************************************/

double GammaDiscreteDistribution::gammp(double a, double x) throw (BadNumberException)
{
	//void gcf(double *gammcf, double a, double x, double *gln);
	//void gser(double *gamser, double a, double x, double *gln);
	double gamser, gammcf, gln;

	if (x <  0.0) throw BadNumberException("GammaDiscreteDistribution::gammp. x < 0.", x);
	if (a <= 0.0) throw BadNumberException("GammaDiscreteDistribution::gammp. a <= 0.", a);
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

/******************************************************************************/
