#define _CRT_SECURE_NO_WARNINGS
#include "speech_enhance.h"

Word32 frame;
void main(void)
{
	Word16 in_send[FRM_LEN], in_air[FRM_LEN], in_receive[FRM_LEN], out_send[FRM_LEN], out_receive[FRM_LEN], ns_flag = 0;
	float in_airf[FRM_LEN], in_receivef[FRM_LEN], out_receivef[FRM_LEN], in_f32[FRM_LEN], out_f32[FRM_LEN];
	float eng_in, eng_out, eng_ratio = 1;
	Word32 i;
	FILE *fq, *fs;

	
	ANC_STRUCT anc_st;
	ANC_init(&anc_st);
	
	fq = fopen("NOISY_fileid_0.pcm", "rb");	//主麦输入
	fs = fopen("NOISY_fileid_0out.pcm", "wb");			//上行输出
	if ((fq == NULL)) { printf("Cann't open file!!!\n");  exit(0); }

	for (frame = 0; ; frame++)
	{
		if (fread(in_send, sizeof(short), FRM_LEN, fq) != FRM_LEN) break;
		for (i = 0; i < FRM_LEN2; i++) in_f32[i] = (float)in_send[i];
		ANC_run(&anc_st, in_f32, out_f32);	
		for (i = 0; i < FRM_LEN2; i++) out_send[i] = (short)out_f32[i];
		fwrite(out_send, sizeof(short), FRM_LEN, fs);				// 上行输出

	}
	fclose(fq);
	fclose(fs);

}
