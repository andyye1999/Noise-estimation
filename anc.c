#include "speech_enhance.h"

short ch_tbl[NUM_CHAN][2] = { {2, 6}, {7, 10}, {11, 14}, {15, 18}, {19, 22}, {23, 26}, {27, 32}, {33, 38}, {39, 44}, {45, 52}, {53, 60}, {61, 70}, {71, 82}, {83, 96}, {97, 110}, {111, 127} };
short vm_tbl[90] = { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, 17, 17, 18, 19, 20, 20, 21, 22, 23, 24, 24, 25, 26, 27, 28, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50 };
float window[DELAY + FRM_LEN];

double phs_tbl[FFT_LEN];
float alpha_eta = 0.92;
float alpha_s = 0.9;
float alpha_d = 0.85;
float beta = 2;
float eta_min = 0.0158;
float GH0 = 0.1256980508997653471;
float gama0 = 4.6;
float gama1 = 3;
float zeta0 = 1.67;
float Bmin = 1.66;
float l_mod_lswitch = 0;
float Vwin = 15;
float Nwin = 8;

void c_fft(float *dat_vec, short isign)
{
	short i, j, k, ii, jj, kk, ji, kj;
	float ftmp, re7, im7;

	for (i = j = 0; i < FFT_LEN - 2; i += 2)
	{
		if (j > i)
		{
			ftmp = dat_vec[i];
			dat_vec[i] = dat_vec[j];
			dat_vec[j] = ftmp;

			ftmp = dat_vec[i + 1];
			dat_vec[i + 1] = dat_vec[j + 1];
			dat_vec[j + 1] = ftmp;
		}

		k = FFT_LEN / 2;
		while (j >= k)
		{
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	if (isign == 1)
	{
		for (i = 0; i < NUM_STAGE; i++)
		{
			jj = 2 << i;
			kk = 4 << i;
			ii = (short)(FFT_LEN >> 1 >> i);

			for (j = 0; j < jj; j += 2)
			{
				ji = j * ii;
				for (k = j; k < FFT_LEN; k += kk)
				{
					kj = k + jj;
					re7 = dat_vec[kj] * (float)(phs_tbl[ji]) - dat_vec[kj + 1] * (float)(phs_tbl[ji + 1]);
					im7 = dat_vec[kj + 1] * (float)(phs_tbl[ji]) + dat_vec[kj] * (float)(phs_tbl[ji + 1]);
					dat_vec[kj] = (dat_vec[k] - re7) / 2;
					dat_vec[kj + 1] = (dat_vec[k + 1] - im7) / 2;
					dat_vec[k] = (dat_vec[k] + re7) / 2;
					dat_vec[k + 1] = (dat_vec[k + 1] + im7) / 2;
				}
			}
		}
	}
	else
	{
		for (i = 0; i < NUM_STAGE; i++)
		{
			jj = 2 << i;
			kk = 4 << i; 
			ii = (short)(FFT_LEN >> 1 >> i); 

			for (j = 0; j < jj; j += 2)
			{
				ji = j * ii;	
				for (k = j; k < FFT_LEN; k += kk) 
				{
					kj = k + jj; 

					re7 = dat_vec[kj] * (float)phs_tbl[ji] + dat_vec[kj + 1] * (float)phs_tbl[ji + 1];
					im7 = dat_vec[kj + 1] * (float)phs_tbl[ji] - dat_vec[kj] * (float)phs_tbl[ji + 1];
					dat_vec[kj] = dat_vec[k] - re7;
					dat_vec[kj + 1] = dat_vec[k + 1] - im7;
					dat_vec[k] = dat_vec[k] + re7;
					dat_vec[k + 1] = dat_vec[k + 1] + im7;
				}
			}
		}
	}
}

void r_fft(float *dat_vec, short isign)
{
	short i, j;
	float re8, im8, re9, im9;

	if (isign == 1)
	{
		c_fft(dat_vec, isign);

		re8 = dat_vec[0];
		re9 = dat_vec[1];
		dat_vec[0] = re8 + re9;
		dat_vec[1] = re8 - re9;

		for (j = FFT_LEN - 2, i = 2; i <= FFT_LEN / 2; i += 2, j -= 2)
		{
			re8 = dat_vec[i] + dat_vec[j];
			im8 = dat_vec[i + 1] - dat_vec[j + 1];
			re9 = dat_vec[i + 1] + dat_vec[j + 1];
			im9 = dat_vec[j] - dat_vec[i];
			dat_vec[i] = (float)(re8 + phs_tbl[i] * re9 - phs_tbl[i + 1] * im9) / 2;
			dat_vec[i + 1] = (float)(im8 + phs_tbl[i] * im9 + phs_tbl[i + 1] * re9) / 2;
			dat_vec[j] = (float)(re8 + phs_tbl[j] * re9 + phs_tbl[j + 1] * im9) / 2;
			dat_vec[j + 1] = (float)(-im8 - phs_tbl[j] * im9 + phs_tbl[j + 1] * re9) / 2;
		}
	}
	else
	{
		re8 = dat_vec[0];
		re9 = dat_vec[1];
		dat_vec[0] = (re8 + re9) / 2;
		dat_vec[1] = (re8 - re9) / 2;

		for (j = FFT_LEN - 2, i = 2; i <= FFT_LEN / 2; i += 2, j -= 2)
		{
			re8 = dat_vec[i] + dat_vec[j];
			im8 = dat_vec[i + 1] - dat_vec[j + 1];
			re9 = -(*(dat_vec + i + 1) + *(dat_vec + j + 1));
			im9 = -(*(dat_vec + j) - *(dat_vec + i));
			dat_vec[i] = (float)(re8 + phs_tbl[i] * re9 + phs_tbl[i + 1] * im9) / 2;
			dat_vec[i + 1] = (float)(im8 + phs_tbl[i] * im9 - phs_tbl[i + 1] * re9) / 2;
			dat_vec[j] = (float)(re8 + phs_tbl[j] * re9 - phs_tbl[j + 1] * im9) / 2;
			dat_vec[j + 1] = (float)(-im8 - phs_tbl[j] * im9 - phs_tbl[j + 1] * re9) / 2;
		}
		c_fft(dat_vec, isign);
	}
}

void ANC_init(ANC_STRUCT *st)
{
	short i;
	float arg;

	arg = (float)PI;
	for (i = 0; i < DELAY; i++) window[FRM_LEN + DELAY - 1 - i] = window[i] = (float)pow(sin((2 * i + 1)*arg / DELAY / 4), 2);
	for (i = DELAY; i < FRM_LEN; i++) window[i] = 1;

	for (i = 0; i < FFT_LEN / 2; i++)
	{
		phs_tbl[2 * i] = cos(-i * 2 * PI / FFT_LEN);
		phs_tbl[2 * i + 1] = sin(-i * 2 * PI / FFT_LEN);
	}

	st->send_hyster_cnt = st->send_last_update_cnt = st->send_update_cnt = st->send_frame_cnt = 0;
	for (st->send_pre_emp_mem = 0, i = 0; i < DELAY; i++) st->send_in_overlap[i] = 0;
	for (i = 0; i < NUM_CHAN; i++) st->send_ch_enrg_long_db[i] = st->send_ch_noise[i] = st->send_ch_enrg[i] = 0;
	for (st->send_de_emp_mem = 0, i = 0; i < FFT_LEN - FRM_LEN; i++) st->send_out_overlap[i] = 0;

	for (st->ref_pre_emp_mem = 0, i = 0; i < DELAY; i++) st->ref_in_overlap[i] = 0;

	for (i = 0; i < FFT_LEN / 2; i++) st->ener_ref_aver[i] = 0;

	for (i = 0; i < FFT_LEN; i++) st->rab[i] = 0;
	for (i = 0; i < FFT_LEN / 2; i++) st->raa[i] = 0;

	st->sle_flag = 0;
	st->anc_flag = 0;
	st->loop = 1;
	st->framenum = 0;
}

void ANC_run(ANC_STRUCT *st, float *send_in, float *send_out)
{
	short i, j, send_update_flag, index_cnt, vm_sum, ch_snr[NUM_CHAN], freq_num, hist[100];
	float send_dat_buf[FFT_LEN], enrg, tne, tce, ftmp1, ch_enrg_dev, ch_enrg_db[NUM_CHAN], alpha, vv, ch_gain[FFT_LEN / 2];
	float ref_dat_buf[FFT_LEN], ener_ref[FFT_LEN / 2], ener_send[FFT_LEN / 2], corr_ref_ref, corr_send_send, corr_ref_send, aver_level;
	float err[FFT_LEN], ref_ener[FFT_LEN / 2];
	float Ya2[FFT_LEN / 2], lambda_dav[FFT_LEN / 2], lambda_d[FFT_LEN / 2], gamma[FFT_LEN / 2], Smin[FFT_LEN / 2], S[FFT_LEN / 2], St[FFT_LEN / 2], GH1[FFT_LEN / 2], Smint[FFT_LEN / 2], Smin_sw[FFT_LEN / 2], Smint_sw[FFT_LEN / 2], eta_2term[FFT_LEN / 2], Sf[FFT_LEN / 2];
	float eta[FFT_LEN / 2], v[FFT_LEN / 2], gama_min[FFT_LEN / 2], zeta[FFT_LEN / 2], I_f[FFT_LEN / 2], conv_I[FFT_LEN / 2], Sft[FFT_LEN / 2], conv_Y[FFT_LEN / 2];
	float win_freq[WIN_FREQ_NUM] = { 0.5,1.0,0.5 };
	float Sf_conv[FFT_LEN / 2 + WIN_FREQ_NUM - 1], I_f_conv[FFT_LEN / 2 + WIN_FREQ_NUM - 1], conv_Y_conv[FFT_LEN / 2 + WIN_FREQ_NUM - 1];
	float gamma_mint[FFT_LEN / 2], zetat[FFT_LEN / 2], qhat[FFT_LEN / 2], phat[FFT_LEN / 2];
	float alpha_dt[FFT_LEN / 2], G[FFT_LEN / 2];



	st->framenum++;
	send_dat_buf[DELAY] = (float)(send_in[0] - EMP_FAC * st->send_pre_emp_mem);
	for (i = 1; i < FRM_LEN; i++) send_dat_buf[DELAY + i] = (float)(send_in[i] - EMP_FAC * send_in[i - 1]);
	st->send_pre_emp_mem = send_in[FRM_LEN - 1];
	for (i = 0; i < DELAY; i++) send_dat_buf[i] = st->send_in_overlap[i];
	for (i = 0; i < DELAY; i++) st->send_in_overlap[i] = send_dat_buf[FRM_LEN + i];
	for (i = 0; i < FRM_LEN + DELAY; i++) send_dat_buf[i] *= window[i];
	for (i = DELAY + FRM_LEN; i < FFT_LEN; i++) send_dat_buf[i] = 0;
	r_fft(send_dat_buf, 1);


	for (i = 1; i < FFT_LEN / 2; i++)  Ya2[i] = send_dat_buf[2 * i] * send_dat_buf[2 * i] + send_dat_buf[2 * i + 1] * send_dat_buf[2 * i + 1];
	conv(win_freq, Ya2 + 1, 3, FFT_LEN / 2 - 1, Sf_conv);
	for (i = 1; i < FFT_LEN / 2; i++)  Sf[i] = Sf_conv[i - 1];
	if (st->loop == 1)
	{
		for (i = 1; i < FFT_LEN / 2; i++)
		{
			lambda_dav[i] = Ya2[i];
			lambda_d[i] = Ya2[i];
			gamma[i] = 1;
			Smin[i] = Sf[i];
			S[i] = Sf[i];
			St[i] = Sf[i];
			GH1[i] = 1;
			Smint[i] = Sf[i];
			Smin_sw[i] = Sf[i];
			Smint_sw[i] = Sf[i];
			eta_2term[i] = GH1[i] * GH1[i] * gamma[i];
			
		}
		st->loop = 0;
	}
	else
	{
		for (i = 1; i < FFT_LEN / 2; i++)
		{
			lambda_dav[i] = st->lambda_dav[i];
			lambda_d[i] = st->lambda_d[i];
			gamma[i] = st->gamma[i];
			Smin[i] = st->Smin[i];
			S[i] = st->S[i];
			St[i] = st->St[i];
			GH1[i] = st->GH1[i];
			Smint[i] = st->Smint[i];
			Smin_sw[i] = st->Smin_sw[i];
			Smint_sw[i] = st->Smint_sw[i];
			eta_2term[i] = st->eta_2term[i];
		}
		
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		gamma[i] = Ya2[i] / max(lambda_d[i], 0.00001);
	}
	
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		eta[i] = alpha_eta * eta_2term[i] + (1 - alpha_eta) * max(gamma[i] - 1, 0);
		eta[i] = max(eta[i], eta_min);
		v[i] = gamma[i] * eta[i] / (1 + eta[i]);
		GH1[i] = eta[i] / (1 + eta[i]) * exp(0.5 * expint_E1(v[i], 0));
	}
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		S[i] = alpha_s * S[i] + (1 - alpha_s) * Sf[i];
	}


	if (st->framenum < UPDATE_NUM)
	{
		for (i = 1; i < FFT_LEN / 2; i++)
		{
			Smin[i] = min(Smin[i], S[i]);
			Smin_sw[i] = min(Smin_sw[i], S[i]);
		}
	}
	else if (st->framenum == UPDATE_NUM)
	{
		for (i = 1; i < FFT_LEN / 2; i++)
		{
			Smin[i] = min(Smin_sw[i], S[i]);
			Smin_sw[i] = S[i];
		}
	}
	
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		gama_min[i] = Ya2[i] / Bmin / Smin[i];
		zeta[i] = S[i] / Bmin / Smin[i];
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		if (gama_min[i] < gama0 && zeta[i] < zeta0)
		{
			I_f[i] = 1;
		}
		else
		{
			I_f[i] = 0;
		}
	}
	conv(win_freq, I_f + 1, 3, FFT_LEN / 2 - 1, I_f_conv);
	for (i = 1; i < FFT_LEN / 2; i++)  conv_I[i] = I_f_conv[i - 1];
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		Sft[i] = St[i];
	}
	for (i = 1; i < FFT_LEN / 2; i++) conv_Y[i] = conv_I[i] * Ya2[i];
	conv(win_freq, conv_Y + 1, 3, FFT_LEN / 2 - 1, conv_Y_conv);
	for (i = 1; i < FFT_LEN / 2; i++)  conv_Y[i] = conv_Y_conv[i - 1];
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		if (conv_I[i] != 0)
		{
			Sft[i] = conv_Y[i] / conv_I[i];
		}
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		St[i] = alpha_s * St[i] + (1 - alpha_s) * Sft[i];
	}

	if (st->framenum < UPDATE_NUM)
	{
		for (i = 1; i < FFT_LEN / 2; i++)
		{
			Smint[i] = min(Smin[i], St[i]);
			Smint_sw[i] = min(Smint_sw[i], St[i]);
		}
	}
	else if (st->framenum == UPDATE_NUM)
	{
		for (i = 1; i < FFT_LEN / 2; i++)
		{
			Smint[i] = min(Smint_sw[i], St[i]);
			Smint_sw[i] = St[i];
			st->framenum = 0;
		}
	}
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		gamma_mint[i] = Ya2[i] / Bmin / Smint[i];
		zetat[i] = S[i] / Bmin / Smint[i];
	}
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		qhat[i] = 1;
		phat[i] = 0;
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		if (gamma_mint[i] > 1 && gamma_mint[i] < gama1 && zetat[i] < zeta0)
		{
			qhat[i] = (gama1 - gamma_mint[i]) / (gama1 - 1);
		}
		else if (gamma_mint[i] >= gama1 || zetat[i] >= zeta0)
		{
			qhat[i] = 0;
		}
	}
	for (i = 1; i < FFT_LEN / 2; i++)
	{
		phat[i] = 1.0 / (1.0 + qhat[i] * (1.0 + eta[i]) * exp(-v[i]) / (1.0 - qhat[i] + 0.00001) );
		if (gamma_mint[i] >= gama1 || zetat[i] >= zeta0) phat[i] = 1;
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		alpha_dt[i] = alpha_d + (1 - alpha_d) * phat[i];
		lambda_dav[i] = alpha_dt[i] * lambda_dav[i] + (1 - alpha_dt[i]) * Ya2[i];
		lambda_d[i] = lambda_dav[i] * beta;
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		gamma[i] = Ya2[i] / max(lambda_d[i], 0.00001);
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		eta[i] = alpha_eta * eta_2term[i] + (1 - alpha_eta) * max(gamma[i] - 1, 0);
		eta[i] = max(eta[i], eta_min);
		v[i] = gamma[i] * eta[i] / (1 + eta[i]);
		GH1[i] = eta[i] / (1 + eta[i]) * exp(0.5 * expint_E1(v[i], 0));
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		G[i] = pow(GH1[i], phat[i]) * pow(GH0, (1 - phat[i]));
		eta_2term[i] = pow(GH1[i],2) * gamma[i];
	}

	for (i = 1; i < FFT_LEN / 2; i++)
	{
		st->lambda_dav[i] = lambda_dav[i];
		st->lambda_d[i] = lambda_d[i];
		st->gamma[i] = gamma[i];
		st->Smin[i] = Smin[i];
		st->S[i] = S[i];
		st->St[i] = St[i];
		st->GH1[i] = GH1[i];
		st->Smint[i] = Smint[i];
		st->Smin_sw[i] = Smin_sw[i];
		st->Smint_sw[i] = Smint_sw[i];
		st->eta_2term[i] = eta_2term[i];
	}

	for (i = 1; i < FFT_LEN / 2; i++) send_dat_buf[2 * i] *= G[i];
	for (i = 1; i < FFT_LEN / 2; i++) send_dat_buf[2 * i + 1] *= G[i];

	r_fft(send_dat_buf, -1);

	for (i = 0; i < FFT_LEN - FRM_LEN; i++) send_dat_buf[i] += st->send_out_overlap[i];
	for (i = 0; i < FFT_LEN - FRM_LEN; i++) st->send_out_overlap[i] = send_dat_buf[i + FRM_LEN];
	for (vv = st->send_de_emp_mem, i = 0; i < FRM_LEN; i++)
	{
		vv = (float)(send_dat_buf[i] + EMP_FAC * vv);
		send_out[i] = max(-32768.0f, min(32767.0f, vv));
	}
	st->send_de_emp_mem = vv;
}
