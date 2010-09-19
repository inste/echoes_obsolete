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
#include <inttypes.h>
#include <string.h>

#include <fftw3.h>

#include "dsp.h"
#include "types.h"
#include "wavfile.h"
#include "ao.h"

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

	struct AOutput * device;

	// Structures for FFT

	fftw_complex *in, *out;
	fftw_plan p;	

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

	device = aout_init();

	// Reading util end of file (frames - number of frames in WAV)

	for (i = 0; i < (sample_rate / sample_size) * (frames / (sample_rate * sample_channels * (sample_bits / 8))) ; ++i) {

		// Reading sample from file and copying to in[][] for FFT
		read_sample_into_buffer(fd, sample_buffer, sample_size);
		
		push_raw_le_data_into_fftw(sample_buffer, in, sample_size);
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
		aout_play(device, sample_buffer, sample_buffer_size);
	}

	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
	
	aout_close(device);

	return 0;
}
