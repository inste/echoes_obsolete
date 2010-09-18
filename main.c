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
#include <math.h>
#include <stdlib.h>
#include <ao/ao.h>
#include <inttypes.h>
#include <string.h>

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

struct Note {
	int note;
	int octave;
};

struct Chord {
	struct Note note[3];
};

void read_data_into_fftw(FILE * fd, fftw_complex * in, int sample_size) {
	int j;
	char tmp_buf[2];
		
	for (j = 0; j < sample_size; ++j) {
		fread(&tmp_buf, sizeof(char), 2, fd);
		fread(&tmp_buf, sizeof(char), 2, fd);
		in[j][0] = 0xff * tmp_buf[1] + tmp_buf[0] * ((tmp_buf[1] < 0) ? -1 : 1); // Little endian :)
		in[j][1] = 0;	
	}
}

void process_fftw_out_into_harmonics(fftw_complex * out, struct Harmonic * harmonics, int sample_size, int sample_msec) {
	int j;
	
	for (j = 0; j < sample_size / 2; ++j) {
		harmonics[j].ampl = sqrt(out[j][0] * out[j][0] + out[j][1] * out[j][1]) / sample_size;
		harmonics[j].phase = atan2(out[j][1], out[j][0]);
		harmonics[j].freq = j / ((double)sample_msec / 1000);
	}
}

int find_main_harmonic(struct Harmonic * harmonics, struct Harmonic * dest, int sample_size) {
	int j, current = 0;
	
	memcpy(dest, &harmonics[0], sizeof(struct Harmonic));
	for (j = 1; j < sample_size / 2; ++j)
		if ( (harmonics[j].ampl > dest->ampl) && (current = j) ) 
			memcpy(dest, &harmonics[j], sizeof(struct Harmonic));
	
	return current;
}

int lookup_for_harmonics(struct Harmonic * harmonics, struct Harmonic * spectre, int current, int sample_size) {
	int j, ch, pos = 1;	
	
	current *= 2;
	while (pos < 10) {
		memcpy(&spectre[pos], &harmonics[current - 10], sizeof(struct Harmonic));
		ch = current - 10;
		// Looking window contains 20 frequencies around n * f (needed for accurate syntheses)
		for (j = current - 9; j < current + 10; ++j)
			if ( (harmonics[j].ampl > spectre[pos].ampl) && (ch = j) ) 
				memcpy(&spectre[pos], &harmonics[j], sizeof(struct Harmonic));
		current = ch * 2;
		// Skip harmonics beyond FFT's resolution
		if (current + 10 > sample_size / 2)
			break;
		++pos;
	}
	return pos;
}

void generate_chord(struct Harmonic * spectre, double * out_buffer, struct Chord chord, int sample_rate, int sample_size, int pos, int offset) {
	int samples[5];
	int j, k, mp;
	double phases[ 3 * 10 ];
	double phase, amplitude, sample;
	const double dpi = 6.2831852;
	
	for (j = 0; j < pos; ++j) {
		mp = j * 3;
		phases[mp] = dpi / (float)sample_rate * compute_freq(chord.note[0].note, chord.note[0].octave) * (j + 1);
		phases[mp + 1] = dpi / (float)sample_rate * compute_freq(chord.note[1].note, chord.note[1].octave) * (j + 1);
		phases[mp + 2] = dpi / (float)sample_rate * compute_freq(chord.note[2].note, chord.note[2].octave) * (j + 1);
	}
	// Generate and upmix chord directly into out buffer
	for (j = 0; j < sample_size; ++j) {
		sample = 0;
		phase = offset + j;
		for (k = 0; k < pos; ++k) {
			mp = k * 3;
			sample += (spectre[k].ampl / 3) * (sin(phase * phases[mp]) + sin(phase * phases[mp + 1]) + sin(phase * phases[mp + 2]));	
		}
		out_buffer[2 * j] = out_buffer[2 * j + 1] = sample;
	}
}

void normalize_out_buffer(double * out_buffer, char * sample_buffer, int sample_size) {
	int i, sample;
	double t, max = 0;
	double multiplier;
	
	for (i = 0; i < sample_size; ++i) {
		t = fabs(out_buffer[2 * i]);
		if (t > max)
			max = t;
	}
	multiplier = (0xffff / 2 - 0xff) / max;
	for (i = 0; i < sample_size; ++i) {
		out_buffer[2 * i] *= multiplier;
		out_buffer[2 * i + 1] *= multiplier;
		sample = (int)out_buffer[2 * i];
		sample_buffer[4 * i] = sample_buffer[4 * i + 2] = sample & 0xff;
		sample_buffer[4 * i + 1] = sample_buffer[4 * i + 3] = (sample >> 8) & 0xff;		
	}
}

int main(int argc, char** argv)
{
	
	struct Harmonic * harmonics;
	struct Harmonic spectre[10];
	int current, pos;
	
	char * sample_buffer;
	double * out_buffer;
	int sample_rate;
	int sample_bits;
	int sample_channels;
	int sample_buffer_size;
	int sample;
	int sample_msec;
	int sample_size;
	float freq;
	double ampl, phase, frequency;
	int i, j, k;
	uint64_t frames;  // FIXME! 4 byte!
	
	FILE * fd;
	
	struct Chord chord;

	
	// Sound output through libao
	
	ao_device * device;
	ao_sample_format format;
	int	default_driver;

	
	// Structures for FFT
	
	fftw_complex *in, *out;
	fftw_plan p;	
	
	// Initialiazing audio output
	

	ao_initialize();
	default_driver = ao_default_driver_id();
	format.bits = 16;
	format.channels = 2;
	format.rate = 44100;
	format.byte_format = AO_FMT_LITTLE;
	format.matrix = "L,R";
	device = ao_open_live(default_driver, &format, NULL);	

	
	// Opening sound file for reading source data

	fd = fopen("test.wav", "rb");
	fseek(fd, 40, SEEK_SET);
	fread(&frames, sizeof(uint64_t), 1, fd);
	printf("%d", frames);
	
	// Hardcoded parameters for DSP

	sample_rate = 44100; // Desc. frequency, Hz
	sample_bits = 16; // Bits per sample
	sample_channels = 2; // Stereo
	sample_msec = 100; // Size of sample, 100 ms for 10 Hz FFT freq. resolution
	
	sample_size = (int)(sample_rate * (sample_msec / 1000.0F));
	sample_buffer_size =  sample_size * sample_bits/8 * sample_channels;
	sample_buffer = (char *) calloc(sample_buffer_size, sizeof(char));
	out_buffer = (double *) calloc(sample_buffer_size / (sample_bits/8), sizeof(double));
	
	// Structures for stroing harmonics of input signal
	
	harmonics = (struct Harmonic *) malloc(sample_size * sizeof(struct Harmonic) / 2);

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sample_size);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sample_size);
	p = fftw_plan_dft_1d(sample_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// D minor
	chord.note[0].note = 2; chord.note[0].octave = 3; 
	chord.note[1].note = 6; chord.note[1].octave = 3; 
	chord.note[2].note = 9; chord.note[2].octave = 3; 	
	
	// Reading util end of file (frames - number of frames in WAV)
	
	for (i = 0; i < (sample_rate / sample_size) * (frames / (sample_rate * sample_channels * (sample_bits / 8))) ; ++i) {
		
		// Reading sample from file and copying to in[][] for FFT
		read_data_into_fftw(fd, in, sample_size);
		// FFT 
		fftw_execute(p);
		// Transform FFT's out[][] into array of harmonics
		process_fftw_out_into_harmonics(out, harmonics, sample_size, sample_msec);
		// Searching main harmonic
		current = find_main_harmonic(harmonics, &spectre[0], sample_size);
		// Searching other harmonics as multiplies of main harmonic (2f, 3f, 4f, ...)
		pos = lookup_for_harmonics(harmonics, spectre, current, sample_size);
		// Debug output
		ampl = 0;
		for (j = 0; j < pos; ++j)
			ampl += spectre[j].ampl;
			
		printf("Sample # %d:\n", i);
		for (j = 0; j < pos; ++j)
			printf("\t%d harmonic: f = %.4lf Hz, a = %.4lf \%\n", j, spectre[j].freq, spectre[j].ampl * 100 / ampl);
		generate_chord(spectre, out_buffer, chord, sample_rate, sample_size, pos, i * sample_size);
		normalize_out_buffer(out_buffer, sample_buffer, sample_size);
		ao_play(device, sample_buffer, sample_buffer_size);
	}
	
	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
	
	ao_close(device);
	ao_shutdown();
	
	return 0;
}
