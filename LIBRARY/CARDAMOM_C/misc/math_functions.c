
#include <math.h>




double randn()
{

    double pi=3.141592653589793;
    double r1=(double)random()/RAND_MAX;
    double r2=(double)random()/RAND_MAX;

    double rn=sqrt(-2*log(r1)) * cos(2*pi*r2);

    return rn;
}

double std(double a[], int n) 
{
    if(n == 0){return 0.0;}
    double sum = 0;
	int i;
    for(i = 0; i < n; i++)
       sum += a[i];
    double mean = sum / n;
    double sq_diff_sum = 0;
    for(i=0;i<n;i++) {
       double diff = a[i] - mean;
       sq_diff_sum = sq_diff_sum + diff * diff;
    }
    double variance = sq_diff_sum / (n-1);
    return sqrt(variance);
}

double par2nor(double p, double mn, double mx)
{
    return log(p/mn)/log(mx/mn);
}

double nor2par(double p, double mn, double mx)
{
    return mn*pow(mx/mn,p);
}

