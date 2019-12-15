#include <iostream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
using namespace std;

const int IMG_WIDTH = 1920;
const int IMG_HEIGHT = 1080;
const int IMG_SIZE = 1920 * 1080;
const char* YUV_FILE1 = "demo/dem1.yuv";
const char* YUV_FILE2 = "demo/dem2.yuv";
const char* OUT_FILE1 = "demo/nosimd_fade.yuv";
const char* OUT_FILE2 = "demo/nosimd_overlay.yuv";

unsigned char yuv_y[2][IMG_SIZE];
unsigned char yuv_u[2][IMG_SIZE / 4];
unsigned char yuv_v[2][IMG_SIZE / 4];
short new_yuv_y[IMG_SIZE];
short new_yuv_u[IMG_SIZE / 4];
short new_yuv_v[IMG_SIZE / 4];
short rgb_r[2][IMG_SIZE];
short rgb_g[2][IMG_SIZE];
short rgb_b[2][IMG_SIZE];
short new_rgb_r[IMG_SIZE];
short new_rgb_g[IMG_SIZE];
short new_rgb_b[IMG_SIZE];
unsigned char raw_yuv_y[IMG_SIZE];
unsigned char raw_yuv_u[IMG_SIZE / 4];
unsigned char raw_yuv_v[IMG_SIZE / 4];


void read_yuv(const char* filename, int file_id)
{
	FILE* file = fopen(filename, "rb");
	int index = file_id - 1;  // index in the array  = 0 / 1
	fread(yuv_y[index], 1, IMG_SIZE, file);
	fread(yuv_u[index], 1, IMG_SIZE / 4, file);
	fread(yuv_v[index], 1, IMG_SIZE / 4, file);
	fclose(file);
}


void write_yuv(FILE* file)
{
	for(int i = 0; i < IMG_SIZE; i++)
		raw_yuv_y[i] = new_yuv_y[i];

	for(int i = 0; i < IMG_SIZE / 4; i++) {
		raw_yuv_u[i] = new_yuv_u[i];
		raw_yuv_v[i] = new_yuv_v[i];
	}

	fwrite(raw_yuv_y, 1, IMG_SIZE, file);
	fwrite(raw_yuv_u, 1, IMG_SIZE / 4, file);
	fwrite(raw_yuv_v, 1, IMG_SIZE / 4, file);

	return;
}

void alpha_blending(int alpha)
{
	for(int i = 0; i < IMG_WIDTH; i++)
		for(int j = 0; j < IMG_HEIGHT; j++) {
			int row_y = i * IMG_HEIGHT + j;
			// R' = A * R / 256
			new_rgb_r[row_y] = alpha * rgb_r[0][row_y] / 256;
			// G' = A * G / 256			
			new_rgb_g[row_y] = alpha * rgb_g[0][row_y] / 256;
			// B' = A * B / 256
			new_rgb_b[row_y] = alpha * rgb_b[0][row_y] / 256;
		}
}

void image_overlay(int alpha)
{
	for(int i = 0; i < IMG_WIDTH; i++)
		for(int j = 0; j < IMG_HEIGHT; j++) {
			int row_y = i * IMG_HEIGHT + j;
			// R' = (A * R1 + (256 - A) * R2) / 256
			new_rgb_r[row_y] = (alpha * rgb_r[0][row_y] + (256 - alpha) * rgb_r[1][row_y]) / 256;
			// G' = (A * G1 + (256 - A) * G2) / 256
			new_rgb_g[row_y] = (alpha * rgb_g[0][row_y] + (256 - alpha) * rgb_g[1][row_y]) / 256;			
			// B' = (A * B1 + (256 - A) * B2) / 256
			new_rgb_b[row_y] = (alpha * rgb_b[0][row_y] + (256 - alpha) * rgb_b[1][row_y]) / 256;
		}
}


void yuv2rgb(int index)
{
	for(int i = 0; i < IMG_WIDTH; i++)
		for(int j = 0; j < IMG_HEIGHT; j++) {

			int row_y = i * IMG_HEIGHT + j;
			int row_uv = (i / 2) * (IMG_HEIGHT / 2) + j / 2;
			// R = Y + 1.140 * V
			rgb_r[index][row_y] = yuv_y[index][row_y] + 1.140 * (yuv_v[index][row_uv] - 128);
			// G = Y - 0.394 * U - 0.581 * V
			rgb_g[index][row_y] = yuv_y[index][row_y] - 0.394 * (yuv_u[index][row_uv] - 128) 
									- 0.581 * (yuv_v[index][row_uv] - 128);
			// B = Y + 2.032 * U
			rgb_b[index][row_y] = yuv_y[index][row_y] + 2.032 * (yuv_u[index][row_uv] - 128);
		}
}


void rgb2yuv()
{
	for(int i = 0; i < IMG_WIDTH; i++)
		for(int j = 0; j < IMG_HEIGHT; j++) {

			int row_y = i * IMG_HEIGHT + j;
			// Y = 0.299 * R + 0.587 * G + 0.114 * B
			new_yuv_y[row_y] = 0.299 * new_rgb_r[row_y] + 0.587 * new_rgb_g[row_y] + 0.114 * new_rgb_b[row_y];

			if(i % 2 == 1 && j % 2 == 1) {
				int row_uv = (IMG_HEIGHT / 2) * (i / 2) + (j / 2);
				// U = -0.147 * R - 0.289 * G + 0.436 * B
				new_yuv_u[row_uv] = - 0.147 * new_rgb_r[row_y] - 0.289 * new_rgb_g[row_y] + 0.436 * new_rgb_b[row_y] + 128;
				// V = 0.615 * R - 0.515 * G - 0.100 * B
				new_yuv_v[row_uv] = 0.615 * new_rgb_r[row_y] - 0.515 * new_rgb_g[row_y] - 0.100 * new_rgb_b[row_y] + 128;
			}
		}
}


void process_picture(int task)
{
	clock_t start_clock = clock();
	FILE* fout = NULL;

	if (task == 1) {
		fout = fopen(OUT_FILE1, "wb");
		yuv2rgb(0);
		for(int alpha = 1; alpha < 255; alpha += 3) {
			alpha_blending(alpha);
			rgb2yuv();
			write_yuv(fout);
		}
	}
	else if (task == 2) {
		fout = fopen(OUT_FILE2, "wb");
		yuv2rgb(0);
		yuv2rgb(1);
		for(int alpha = 1; alpha < 255; alpha += 3) {
			image_overlay(alpha);
			rgb2yuv();
			write_yuv(fout);
		}
	}

	fclose(fout);
	clock_t total_time = clock() - start_clock;

	if(task == 1)
		cout << "Task1: Image Fading" << endl;
	else if (task == 2)
		cout << "Task2: Image Overlay" << endl;

	cout << "Time Elapse: " << fixed << setprecision(3) << (double)total_time / (double)CLOCKS_PER_SEC << "s" << endl;
}


int main()
{	
	read_yuv(YUV_FILE1, 1);
	read_yuv(YUV_FILE2, 2);
	process_picture(1);
	process_picture(2);
	return 0;
}
