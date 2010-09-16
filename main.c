//      main.c
//      
//      Copyright 2010 Ilya <ilya@laptop>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include <stdio.h>
#include <ao/ao.h>
#include <math.h>
#include <inttypes.h>

#include <fftw3.h>

/* Notes
 * C - 0 C# - 1
 * D - 2 D# - 3
 * E - 4 E# - 5
 * F - 6 G - 7 G# - 8
 * A - 9 A# - 10
 * H - 11
 */

/*
 *	Octaves:
 * 	0 - subcontr
 * 	1 - contr
 * 	2 - large
 * 	3 - small
 * 	4 - first
 *  5 - second
 *  ...
 */

double compute_freq(int note, int octave) {
	return 27.5 * pow(2, (octave * 12 + note) / 12.0L);
}

struct Harmonic {
	double freq;
	double phase;
	double ampl;
};

int main(int argc, char** argv)
{
	
	struct Harmonic * harmonics;
	struct Harmonic spectre[10];
	int current, pos, ch, count;
	double constant;
	double * phases;
	
	char * sample_buffer;
	char * temp_buffers[5];
	int sample_rate;
	int sample_bits;
	int sample_channels;
	int sample_buffer_size;
	int sample;
	int samples[5];
	int sample_msec;
	int sample_size;
	float freq;
	double ampl, phase, frequency;
	int i, j, k;
	int frames;  // FIXME! 4 byte!
	char tmp_buf[2];
	
	FILE * fd;
	
	double dpi = 6.2831852;
	
	ao_device *device;
	ao_sample_format format;
	int default_driver;
	
    fftw_complex *in, *out;
    fftw_plan p;	
	

	
	ao_initialize();
	default_driver = ao_default_driver_id();
	format.bits = 16;
	format.channels = 2;
	format.rate = 44100;
	format.byte_format = AO_FMT_LITTLE;
	format.matrix = "L,R";
	device = ao_open_live(default_driver, &format, NULL);	
	

    fd = fopen("test.wav", "rb");
	fseek(fd, 40, SEEK_SET);
	fread(&frames, sizeof(int), 1, fd);
	printf("%d", frames);
	

	sample_rate = 44100; // Hz
	sample_bits = 16;
	sample_channels = 2;
	sample_msec = 100; // For 10 Hz
	sample_size = (int)(sample_rate * (sample_msec / 1000.0F));
	sample_buffer_size =  sample_size * sample_bits/8 * sample_channels;
	sample_buffer = (char *) calloc(sample_buffer_size, sizeof(char));
	
	harmonics = (struct Harmonic *) malloc(sample_size * sizeof(struct Harmonic) / 2);
	phases = (double *) malloc(5 * pos * sizeof(double));
	
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sample_size);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sample_size);
    p = fftw_plan_dft_1d(sample_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	
	for (i = 0; i < (sample_rate / sample_size) * (frames / (sample_rate * sample_channels * (sample_bits / 8))) ; ++i) {
		
		for (j = 0; j < sample_size; ++j) {
			fread(&tmp_buf, sizeof(char), 2, fd);
			fread(&tmp_buf, sizeof(char), 2, fd);
			in[j][0] = 0xff * tmp_buf[1] + tmp_buf[0] * ((tmp_buf[1] < 0) ? -1 : 1); // Little endian :)
			in[j][1] = 0;	
		}
		
		fftw_execute(p);
		
		for (j = 0; j < sample_size / 2; ++j) {
			harmonics[j].ampl = sqrt(out[j][0] * out[j][0] + out[j][1] * out[j][1]) / sample_size;
			harmonics[j].phase = atan2(out[j][1], out[j][0]);
			harmonics[j].freq = j / ((double)sample_msec / 1000);
		}
		
		spectre[0].ampl = harmonics[0].ampl;
		for (j = 1; j < sample_size / 2; ++j) {
			if (harmonics[j].ampl > spectre[0].ampl) {
				current = j;
				spectre[0].ampl = harmonics[j].ampl;
				spectre[0].phase = harmonics[j].phase;
				spectre[0].freq = harmonics[j].freq;
			}
		}
		
		current *= 2;
		pos = 1;
		while (pos < 10) {
			
			spectre[pos].ampl = harmonics[current - 10].ampl;
			spectre[pos].phase = harmonics[current - 10].phase;
			spectre[pos].freq = harmonics[current - 10].freq;
			ch = current - 10;
			
			for (j = current - 9; j < current + 10; ++j) {
				if (harmonics[j].ampl > spectre[pos].ampl) {
					ch = j;
					spectre[pos].ampl = harmonics[j].ampl;
					spectre[pos].phase = harmonics[j].phase;
					spectre[pos].freq = harmonics[j].freq;
				}
			}
			current = ch * 2;
			
			if (current + 10 > sample_size / 2)
				break;
			
			++pos;
		}
		
		ampl = 0;
		for (j = 0; j < pos; ++j)
			ampl += spectre[j].ampl;
			
		printf("Sample # %d:\n", i);
		for (j = 0; j < pos; ++j)
			printf("\t%d harmonic: f = %.4lf Hz, a = %.4lf \%\n", j, spectre[j].freq, spectre[j].ampl * 100 / ampl);
		
		for (j = 0; j < pos; ++j) {
			phases[j * 3] = dpi / (float)sample_rate * compute_freq(0, 3) * (j + 1);
			phases[j * 3 + 1] = dpi / (float)sample_rate * compute_freq(4, 3) * (j + 1);
			phases[j * 3 + 2] = dpi / (float)sample_rate * compute_freq(7, 3) * (j + 1);
		}
		
		// Generate and upmix chord (C major)
		
		for (j = 0; j < sample_size; ++j) {
			samples[0] = samples[1] = samples[2] = 0;
			for (k = 0; k < pos; ++k) {
				samples[0] += (int)(spectre[k].ampl * sin((i * sample_size + j) * phases[k * 3]) / 3);
				samples[1] += (int)(spectre[k].ampl * sin((i * sample_size + j) * phases[k * 3 + 1]) / 3);
				samples[2] += (int)(spectre[k].ampl * sin((i * sample_size + j) * phases[k * 3 + 2]) / 3);
			}
			sample = samples[0] + samples[1] + samples[2];	
			
			sample_buffer[4 * j] = sample_buffer[4 * j + 2] = sample & 0xff;
			sample_buffer[4 * j + 1] = sample_buffer[4 * j + 3] = (sample >> 8) & 0xff;				
		}
		
		ao_play(device, sample_buffer, sample_buffer_size);
	}
	
	
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
	
	ao_close(device);
	ao_shutdown();
	
	return 0;
}
