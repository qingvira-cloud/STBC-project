#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <complex> 
#include <vector>
#include <cstdlib>
#include <iomanip>





#define  OUTPUTFILE  "coherent_uncoded_bpsk_awgn.txt"
#define  OPENSTRING  "#Eb/No(dB)        BitErrRate          BitErrSample    TotalSample \n"
#define  DATASTRING  "%5.2f             %18.16f             %7d             %11u \n"
#define  DB0		0
#define  DBSTEP		5.0
#define  POINTNO	4
#define  ERRNOSTEP  100
#define  MAXERRNO	300
#define  SAMPLE_NUM 500000
#define  NORMALIZE  0.70710678118654752440084436210485


using namespace std;



#define  M3  0x00000004L
#define  M31 0x40000000L
#define  eulen(r1,i1,r2,i2) ((r1)-(r2))*((r1)-(r2))+((i1)-(i2))*((i1)-(i2))



FILE   *fp;
static char  filename[80]=OUTPUTFILE;
int	 point,biterrno,errlevel,u0,u1,t,de_info,de_info1,dsignal,samp;
long    pnstate,pntemp;
double   snrdb,snr,deviate;

double   drand48()
{
	double  w;
	/* RAND_MAX = 0xffff;  */
	w=(double)rand()/RAND_MAX;

	return  w;
}

void initial()
{
	if (fopen_s(&fp, filename, "a") != 0 || fp == nullptr) {
		printf("\nOpen file error!\n");
		exit(1);
	}
	else {
		fputs(OPENSTRING, fp);
		fclose(fp);
		printf("\nProgram Start...\n\n");
	}
}

double normal()
{
  static int iset=0;
  static double gset;
  double fac,r,v1,v2;
  if(iset==0)
  { do
    {
      v1=2.0*drand48()-1.0;
      v2=2.0*drand48()-1.0;
      r=v1*v1+v2*v2;
    }
    while(r>=1.0||r==0.0);
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return(v2*fac);
  }  else
  {
    iset=0;
    return gset;
  }
}

void record()
{
	FILE* fp = fopen(filename, "a");
	if (!fp)
		return;

	fprintf(fp, DATASTRING, snrdb, (double)biterrno / dsignal, biterrno, dsignal);
	fclose(fp);
}


void receive()
{
	double sx,sy;

	complex<double> r0,r1,s0,s1;
	complex<double> h0 (NORMALIZE*normal(),NORMALIZE*normal());
	complex<double> n0 (deviate*normal(),deviate*normal());
	complex<double> h1 (NORMALIZE*normal(),NORMALIZE*normal());
	complex<double> n1 (deviate*normal(),deviate*normal());

	complex<double> qpsk_map[2] = {complex<double>(NORMALIZE,0),complex<double>(-NORMALIZE,0)};
	
	u0=rand()%2;
	u1=rand()%2;


	r0= h0*qpsk_map[u0]+h1*qpsk_map[u1]+n0;
	r1=-h0*conj(qpsk_map[u1])+h1*conj(qpsk_map[u0])+n1;
	
	s0=conj(h0)*r0+h1*conj(r1);
	s1=conj(h1)*r0-h0*conj(r1);


	 
	double a = real((s0 - qpsk_map[0]) * (conj(s0) - conj(qpsk_map[0])));
	double b = real((s0 - qpsk_map[1]) * (conj(s0) - conj(qpsk_map[1])));


	if     (a<b)
	{ 
		de_info=0;
	}
	else
	{ 
		de_info=1;
	}
	int bit_er=de_info-u0;
///////////////////////////////////////////////////
	double a1 = real((s1 - qpsk_map[0]) * (conj(s1) - conj(qpsk_map[0]))); ;
	double b1 = real((s1 - qpsk_map[1]) * (conj(s1) - conj(qpsk_map[1])));
	if      (a1<b1)
	{
		de_info1=0;
	}

	else
	{ 
		de_info1=1;
	}
	int bit_er1=de_info1-u1;
///////////////////////////////////////////////////


	//printf("(bit_er=%d)\n",bit_er);

	if(bit_er==0)
	{
		dsignal=dsignal+1;
	}
	else 
	{
		biterrno=biterrno+1;
		dsignal=dsignal+1;
	}
/////////////////////////////////////////////////////
	if(bit_er1==0)
	{
		dsignal=dsignal+1;
	}
	else 
	{
		biterrno=biterrno+1;
		dsignal=dsignal+1;
		
	}
	//	cout << biterrno << endl;
/////////////////////////////////////////////////////	
}


int main()
{
	srand((unsigned)time(NULL));
	initial();
	for (point = 0; point < POINTNO; point++)
	{
		pnstate = 0xaaaaaaaaL;
		snrdb = DB0 + point * DBSTEP;
		snr = pow(10.0, 0.1 * snrdb);
		deviate = sqrt(0.5 / snr);
		biterrno = 0;
		errlevel = ERRNOSTEP;
		dsignal = 0;

		while (biterrno < MAXERRNO && dsignal < 5000000)
		{
			if (biterrno > errlevel)
			{
				errlevel += ERRNOSTEP;
			}
			receive();
		}
		record();

		printf("%5.2f   %18.16f   %7d   %11u \n", (snrdb), (double)biterrno / dsignal, biterrno, dsignal);
	}
	return 0;
}