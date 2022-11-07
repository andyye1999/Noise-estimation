#ifndef SPEECH_ENHANCE_H_
#define SPEECH_ENHANCE_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "expint.h"

typedef short Word16;
typedef long  Word32;
typedef unsigned short UWord16;
typedef unsigned long  UWord32;

#define FS    8000
#define FRM_LEN     160
#define FRM_LEN2    160       
#define	TRUE			1
#define	FALSE			0
#define	PI		   3.1415926535897932384626433832795
#define NOISE_THD 40
#define DOWN_FAC 0.316227766f
#define STAGES     8
#define ANA_LEN   (1<<STAGES)
#define HALF_ANAL_BLOCKL  ((1<<(STAGES-1))+1)
#define SIMULT              3
#define HIST_PAR_EST            1000

#define	NUM_STAGE  7
#define	FFT_LEN	   (2<<NUM_STAGE)
#define	DELAY			24
#define	NUM_CHAN		16
#define	MID_CHAN		5
#define	LO_CHAN			0
#define	HI_CHAN			15
#define	UPDATE_THLD		35
#define	METRIC_THLD		45
#define	INDEX_THLD		12
#define	SETBACK_THLD	12
#define	SNR_THLD		6
#define	INDEX_CNT_THLD	5
#define	UPDATE_CNT_THLD	50
#define	NORM_ENRG		(1.0f)
#define	NOISE_FLOOR		(1.0f / NORM_ENRG)
#define	MIN_CHAN_ENRG	(0.0625f/ NORM_ENRG)
#define	INE		     	(16.0f / NORM_ENRG)
#define	MIN_GAIN		(-13.0f)
#define	GAIN_SLOPE		0.39
#define	CNE_SM_FAC		0.1f
#define	CEE_SM_FAC		(0.55f)
#define	EMP_FAC	    	0.8
#define	HYSTER_CNT_THLD	6
#define	HIGH_TCE_DB		(50.0f)
#define	LOW_TCE_DB		(30.0f)
#define	TCE_RANGE		(HIGH_TCE_DB - LOW_TCE_DB)
#define	HIGH_ALPHA		0.99f
#define	LOW_ALPHA		0.50f
#define	ALPHA_RANGE		(HIGH_ALPHA - LOW_ALPHA)
#define	DEV_THLD		28.0
#define	FAST_STEP_SIZE	   0.05f
#define	SLOW_STEP_SIZE	   0.005f
#define UPDATE_NUM 15
#define WIN_FREQ_NUM 3


typedef struct
{
	short send_frame_cnt;
	short send_update_cnt;
	short send_last_update_cnt;
	short send_hyster_cnt;
	float send_pre_emp_mem;
	float send_in_overlap[DELAY];
	float send_ch_enrg[NUM_CHAN];
	float send_ch_noise[NUM_CHAN];
	float send_ch_enrg_long_db[NUM_CHAN];
	float send_de_emp_mem;
	float send_out_overlap[FFT_LEN - FRM_LEN];

	float ref_pre_emp_mem;
	float ref_in_overlap[DELAY];

	float ener_ref_aver[FFT_LEN / 2];
	float rab[FFT_LEN];
	float raa[FFT_LEN / 2];
	short anc_flag;
	short sle_flag;
	short loop;
	//float Smin[UPDATE_NUM][FFT_LEN / 2];
	//float Smin_sw[UPDATE_NUM][FFT_LEN / 2];
	short framenum;
	float lambda_dav[FFT_LEN / 2];
	float lambda_d[FFT_LEN / 2];
	float gamma[FFT_LEN / 2];
	float Smin[FFT_LEN / 2];
	float S[FFT_LEN / 2];
	float St[FFT_LEN / 2];
	float GH1[FFT_LEN / 2];
	float Smint[FFT_LEN / 2];
	float Smin_sw[FFT_LEN / 2];
	float Smint_sw[FFT_LEN / 2];
	float eta_2term[FFT_LEN / 2];
}ANC_STRUCT;

void ANC_init(ANC_STRUCT *st);
void ANC_run(ANC_STRUCT *st, float *send_in, float *send_out);


void c_fft(float *dat_vec, short isign);
void r_fft(float *dat_vec, short isign);

#endif
