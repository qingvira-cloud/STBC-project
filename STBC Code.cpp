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
#define  OPENSTRING  "#Eb/No(dB)        BitErrRate                    BitErrSample           TotalSample \n"
#define  DATASTRING  "%5.2f             %18.16f             %7d             %11u \n",\
                     (snrdb),(double)biterrno/dsignal,biterrno,dsignal,deviate
#define  DB0		0
#define  DBSTEP		5.0
#define  POINTNO	7
#define  ERRNOSTEP  100
#define  MAXERRNO	300
#define  SAMPLE_NUM 5000000
#define  NORMALIZE  0.70710678118654752440084436210485


using namespace std;



#define  M3  0x00000004L
#define  M31 0x40000000L
#define  eulen(r1,i1,r2,i2) ((r1)-(r2))*((r1)-(r2))+((i1)-(i2))*((i1)-(i2))



FILE* fp;
static char  filename[80] = OUTPUTFILE;
int	 point, biterrno, errlevel, u0, u1, t, de_info, de_info1, dsignal, samp;
long    pnstate, pntemp;
double   snrdb, snr, deviate;

static double   drand48()
{
	double  w;
	/* RAND_MAX = 0xffff;  */
	w = (double)rand() / RAND_MAX;

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
	static int iset = 0;
	static double gset;
	double fac, r, v1, v2;
	if (iset == 0)
	{
		do
		{
			v1 = 2.0 * drand48() - 1.0;
			v2 = 2.0 * drand48() - 1.0;
			r = v1 * v1 + v2 * v2;
		} while (r >= 1.0 || r == 0.0);
		fac = sqrt(-2.0 * log(r) / r);
		gset = v1 * fac;
		iset = 1;
		return(v2 * fac);
	}
	else
	{
		iset = 0;
		return gset;
	}
}

void record()
{
	FILE* fp = nullptr;
	fopen_s(&fp, filename, "a");      // 或用 fopen 並定義 _CRT_SECURE_NO_WARNINGS
	if (!fp) return;

	fprintf(fp, DATASTRING,
		snrdb, (double)biterrno / dsignal, biterrno, dsignal);

	fclose(fp);
}


void receive_2t1r()
{
	complex<double> r0, r1, s0, s1;
	complex<double> h0(NORMALIZE * normal(), NORMALIZE * normal());
	complex<double> n0(deviate * normal(), deviate * normal());
	complex<double> h1(NORMALIZE * normal(), NORMALIZE * normal());
	complex<double> n1(deviate * normal(), deviate * normal());

	complex<double> bpsk_map[2] = { 
		complex<double>(NORMALIZE,0),
		complex<double>(-NORMALIZE,0) 
	};

	u0 = rand() % 2;
	u1 = rand() % 2;

	r0 = h0 * bpsk_map[u0] + h1 * bpsk_map[u1] + n0;
	r1 = -h0 * conj(bpsk_map[u1]) + h1 * conj(bpsk_map[u0]) + n1;

	s0 = conj(h0) * r0 + h1 * conj(r1);
	s1 = conj(h1) * r0 - h0 * conj(r1);

	double a = real((s0 - bpsk_map[0]) * (conj(s0) - conj(bpsk_map[0]))); ;
	double b = real((s0 - bpsk_map[1]) * (conj(s0) - conj(bpsk_map[1])));

	if (a < b)
	{
		de_info = 0;
	}

	else
	{
		de_info = 1;
	}
	int bit_er = de_info - u0;
	///////////////////////////////////////////////////
	double a1 = real((s1 - bpsk_map[0]) * (conj(s1) - conj(bpsk_map[0]))); ;
	double b1 = real((s1 - bpsk_map[1]) * (conj(s1) - conj(bpsk_map[1])));
	if (a1 < b1)
	{
		de_info1 = 0;
	}

	else
	{
		de_info1 = 1;
	}
	int bit_er1 = de_info1 - u1;
	///////////////////////////////////////////////////


		//printf("(bit_er=%d)\n",bit_er);

	if (bit_er == 0)
	{
		dsignal = dsignal + 1;
	}
	else
	{

		biterrno = biterrno + 1;
		dsignal = dsignal + 1;

	}
	/////////////////////////////////////////////////////
	if (bit_er1 == 0)
	{
		dsignal = dsignal + 1;
	}
	else
	{

		biterrno = biterrno + 1;
		dsignal = dsignal + 1;

	}
	//	cout << biterrno << endl;
/////////////////////////////////////////////////////	
}

void receive_2t2r()
{
	complex<double> r0, r1, s0, s1, s2, s3, r2, r3;
	complex<double> h0(NORMALIZE * normal(), NORMALIZE * normal());
	complex<double> n0(deviate * normal(), deviate * normal());
	complex<double> h1(NORMALIZE * normal(), NORMALIZE * normal());
	complex<double> n1(deviate * normal(), deviate * normal());
	complex<double> h2(NORMALIZE * normal(), NORMALIZE * normal());
	complex<double> n2(deviate * normal(), deviate * normal());
	complex<double> h3(NORMALIZE * normal(), NORMALIZE * normal());
	complex<double> n3(deviate * normal(), deviate * normal());

	complex<double> bpsk_map[2] = {
		complex<double>(NORMALIZE,0),
		complex<double>(-NORMALIZE,0)
	};

	u0 = rand() % 2;
	u1 = rand() % 2;
	complex<double> s0sym = bpsk_map[u0];
	complex<double> s1sym = bpsk_map[u1];

	// --- 接收兩個時間槽 ---
	// Rx1：
	complex<double> r0 = h0 * s0sym + h1 * s1sym + n0;                      // r0 = h0*s0 + h1*s1 + n0
	complex<double> r1 = -h0 * conj(s1sym) + h1 * conj(s0sym) + n1;         // r1 = -h0*s1* + h1*s0* + n1 
	// Rx2：
	complex<double> r2 = h2 * s0sym + h3 * s1sym + n2;
	complex<double> r3 = -h2 * conj(s1sym) + h3 * conj(s0sym) + n3;         // ★ 修正這行

	// --- 各 Rx 做 Alamouti 合併 ---
	// ŝ0^(m) = h0^* r0 + h1 r1^*
	// ŝ1^(m) = h1^* r0 - h0 r1^*
	complex<double> s0_hat_1 = conj(h0) * r0 + h1 * conj(r1);
	complex<double> s1_hat_1 = conj(h1) * r0 - h0 * conj(r1);

	complex<double> s0_hat_2 = conj(h2) * r2 + h3 * conj(r3);
	complex<double> s1_hat_2 = conj(h3) * r2 - h2 * conj(r3);

	// --- 兩路 Rx 做 MRC 相加 ---
	complex<double> s0_hat = s0_hat_1 + s0_hat_2;
	complex<double> s1_hat = s1_hat_1 + s1_hat_2;

	// （可選）做總增益正規化：G = |h0|^2 + |h1|^2 + |h2|^2 + |h3|^2
	double G = norm(h0) + norm(h1) + norm(h2) + norm(h3);
	s0_hat /= G;
	s1_hat /= G;

	// --- BPSK 判決：最近點（也可改成看 real(ŝ) 的正負）---
	double d0p = real((s0_hat - bpsk_map[0]) * (conj(s0_hat) - conj(bpsk_map[0])));
	double d0m = real((s0_hat - bpsk_map[1]) * (conj(s0_hat) - conj(bpsk_map[1])));
	de_info = (d0p < d0m) ? 0 : 1;

	double d1p = real((s1_hat - bpsk_map[0]) * (conj(s1_hat) - conj(bpsk_map[0])));
	double d1m = real((s1_hat - bpsk_map[1]) * (conj(s1_hat) - conj(bpsk_map[1])));
	de_info1 = (d1p < d1m) ? 0 : 1;

	// --- 計數：一次處理兩個 bit ---
	int bit_er = de_info - u0;
	int bit_er1 = de_info1 - u1;

	if (bit_er == 0)  dsignal += 1;
	else { ++biterrno; dsignal += 1; }

	if (bit_er1 == 0) dsignal += 1;
	else { ++biterrno; dsignal += 1; }
}

int main()
{
	srand((unsigned)time(NULL));
	initial();

	for (point = 0; point < POINTNO; point++)
	{
		snrdb = DB0 + point * DBSTEP;
		snr = pow(10.0, 0.1 * snrdb);
		deviate = sqrt(0.5 / snr);

		// ==== 2T1R 測試（假設你原本的 receive() 是 2T1R）====
		biterrno = 0;
		dsignal = 0;
		errlevel = ERRNOSTEP;

		while (biterrno < MAXERRNO)
		{
			if (biterrno > errlevel) errlevel += ERRNOSTEP;
			receive_2t1r();                    // <-- 這裡用你原本的 2T1R 版本
		}
		record();
		printf("2T1R  %5.2f   %18.16f   %7d   %11u \n",
			snrdb, (double)biterrno / dsignal, biterrno, dsignal);

	}
	return 0;
}