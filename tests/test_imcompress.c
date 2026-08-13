/*
 * Tests for imcompress.c - image compression/decompression
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "test_macros.h"

#define test_path "test_imcompress.fits"
#define test_path2 "test_imcompress2.fits"

/*
 * Helper to create a simple test image
 */
static void
create_test_image(fitsfile **fptr, int bitpix, long nx, long ny, int *status)
{
	long naxes[2] = { nx, ny };

	fits_create_file(fptr, "!" test_path, status);
	fail_if(*status != 0);

	fits_create_img(*fptr, bitpix, 2, naxes, status);
	fail_if(*status != 0);
}

/*
 * Test setting and getting compression type
 */
static void
test_compression_type(void)
{
	fitsfile *fptr;
	int status = 0;
	int ctype;

	create_test_image(&fptr, SHORT_IMG, 64, 64, &status);

	/* Test RICE compression */
	fits_set_compression_type(fptr, RICE_1, &status);
	fail_if(status != 0);

	fits_get_compression_type(fptr, &ctype, &status);
	fail_if(status != 0);
	fail_if(ctype != RICE_1);

	/* Test GZIP compression */
	fits_set_compression_type(fptr, GZIP_1, &status);
	fail_if(status != 0);

	fits_get_compression_type(fptr, &ctype, &status);
	fail_if(status != 0);
	fail_if(ctype != GZIP_1);

	/* Test GZIP_2 compression */
	fits_set_compression_type(fptr, GZIP_2, &status);
	fail_if(status != 0);

	fits_get_compression_type(fptr, &ctype, &status);
	fail_if(status != 0);
	fail_if(ctype != GZIP_2);

	/* Test PLIO compression */
	fits_set_compression_type(fptr, PLIO_1, &status);
	fail_if(status != 0);

	fits_get_compression_type(fptr, &ctype, &status);
	fail_if(status != 0);
	fail_if(ctype != PLIO_1);

	/* Test HCOMPRESS */
	fits_set_compression_type(fptr, HCOMPRESS_1, &status);
	fail_if(status != 0);

	fits_get_compression_type(fptr, &ctype, &status);
	fail_if(status != 0);
	fail_if(ctype != HCOMPRESS_1);

	fits_close_file(fptr, &status);
	fail_if(status != 0);
}

/*
 * Test setting and getting tile dimensions
 */
static void
test_tile_dimensions(void)
{
	fitsfile *fptr;
	int status = 0;
	long dims_in[2] = { 32, 32 };
	long dims_out[2] = { 0, 0 };

	create_test_image(&fptr, SHORT_IMG, 64, 64, &status);

	fits_set_tile_dim(fptr, 2, dims_in, &status);
	fail_if(status != 0);

	fits_get_tile_dim(fptr, 2, dims_out, &status);
	fail_if(status != 0);
	fail_if(dims_out[0] != 32);
	fail_if(dims_out[1] != 32);

	fits_close_file(fptr, &status);
	fail_if(status != 0);
}

/*
 * Test setting and getting quantize level
 */
static void
test_quantize_level(void)
{
	fitsfile *fptr;
	int status = 0;
	float qlevel;

	create_test_image(&fptr, SHORT_IMG, 64, 64, &status);

	fits_set_quantize_level(fptr, 16.0f, &status);
	fail_if(status != 0);

	fits_get_quantize_level(fptr, &qlevel, &status);
	fail_if(status != 0);
	fail_if(fabs(qlevel - 16.0f) > 0.001f);

	fits_close_file(fptr, &status);
	fail_if(status != 0);
}

/*
 * Test setting and getting noise bits
 */
static void
test_noise_bits(void)
{
	fitsfile *fptr;
	int status = 0;
	int noisebits;

	create_test_image(&fptr, SHORT_IMG, 64, 64, &status);

	fits_set_noise_bits(fptr, 4, &status);
	fail_if(status != 0);

	fits_get_noise_bits(fptr, &noisebits, &status);
	fail_if(status != 0);
	fail_if(noisebits != 4);

	fits_close_file(fptr, &status);
	fail_if(status != 0);
}

/*
 * Test RICE compression roundtrip with short image
 */
static void
test_rice_compress_short(void)
{
	fitsfile *infptr, *outfptr;
	int status = 0;
	long naxes[2] = { 64, 64 };
	short *original, *decompressed;
	int i, is_compressed;

	original = malloc(64 * 64 * sizeof *original);
	decompressed = malloc(64 * 64 * sizeof *decompressed);
	fail_if(original == NULL || decompressed == NULL);

	/* Create gradient data */
	for (i = 0; i < 64 * 64; i += 1) {
		original[i] = (short)(i % 1000);
	}

	/* Create input image */
	fits_create_file(&infptr, "!" test_path, &status);
	fail_if(status != 0);
	fits_create_img(infptr, SHORT_IMG, 2, naxes, &status);
	fail_if(status != 0);
	fits_write_img(infptr, TSHORT, 1, 64 * 64, original, &status);
	fail_if(status != 0);
	fits_close_file(infptr, &status);
	fail_if(status != 0);

	/* Open for reading and compress to new file */
	fits_open_file(&infptr, test_path, READONLY, &status);
	fail_if(status != 0);

	fits_create_file(&outfptr, "!" test_path2, &status);
	fail_if(status != 0);

	fits_set_compression_type(outfptr, RICE_1, &status);
	fail_if(status != 0);

	fits_img_compress(infptr, outfptr, &status);
	fail_if(status != 0);

	fits_close_file(infptr, &status);
	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	/* Verify compressed file - compressed image is in HDU 2 */
	fits_open_file(&outfptr, test_path2, READONLY, &status);
	fail_if(status != 0);

	fits_movabs_hdu(outfptr, 2, NULL, &status);
	fail_if(status != 0);

	is_compressed = fits_is_compressed_image(outfptr, &status);
	fail_if(status != 0);
	fail_if(is_compressed != 1);

	/* Read back decompressed data */
	fits_read_img(outfptr, TSHORT, 1, 64 * 64, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);

	/* Verify lossless roundtrip */
	for (i = 0; i < 64 * 64; i += 1) {
		fail_if(decompressed[i] != original[i]);
	}

	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	free(original);
	free(decompressed);
}

/*
 * Test GZIP compression roundtrip
 */
static void
test_gzip_compress(void)
{
	fitsfile *infptr, *outfptr;
	int status = 0;
	long naxes[2] = { 32, 32 };
	int *original, *decompressed;
	int i;

	original = malloc(32 * 32 * sizeof *original);
	decompressed = malloc(32 * 32 * sizeof *decompressed);
	fail_if(original == NULL || decompressed == NULL);

	/* Create test data */
	for (i = 0; i < 32 * 32; i += 1) {
		original[i] = i * 100;
	}

	/* Create input image */
	fits_create_file(&infptr, "!" test_path, &status);
	fail_if(status != 0);
	fits_create_img(infptr, LONG_IMG, 2, naxes, &status);
	fail_if(status != 0);
	fits_write_img(infptr, TINT, 1, 32 * 32, original, &status);
	fail_if(status != 0);
	fits_close_file(infptr, &status);
	fail_if(status != 0);

	/* Compress with GZIP */
	fits_open_file(&infptr, test_path, READONLY, &status);
	fail_if(status != 0);

	fits_create_file(&outfptr, "!" test_path2, &status);
	fail_if(status != 0);

	fits_set_compression_type(outfptr, GZIP_1, &status);
	fail_if(status != 0);

	fits_img_compress(infptr, outfptr, &status);
	fail_if(status != 0);

	fits_close_file(infptr, &status);
	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	/* Read back - compressed image is in HDU 2 */
	fits_open_file(&outfptr, test_path2, READONLY, &status);
	fail_if(status != 0);

	fits_movabs_hdu(outfptr, 2, NULL, &status);
	fail_if(status != 0);

	fits_read_img(outfptr, TINT, 1, 32 * 32, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);

	/* Verify lossless */
	for (i = 0; i < 32 * 32; i += 1) {
		fail_if(decompressed[i] != original[i]);
	}

	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	free(original);
	free(decompressed);
}

/*
 * Test PLIO compression (good for mask images)
 */
static void
test_plio_compress(void)
{
	fitsfile *infptr, *outfptr;
	int status = 0;
	long naxes[2] = { 64, 64 };
	short *original, *decompressed;
	int i;

	original = malloc(64 * 64 * sizeof *original);
	decompressed = malloc(64 * 64 * sizeof *decompressed);
	fail_if(original == NULL || decompressed == NULL);

	/* Create mask-like data (mostly zeros with some ones) */
	for (i = 0; i < 64 * 64; i += 1) {
		original[i] = (i % 10 == 0) ? 1 : 0;
	}

	/* Create input image */
	fits_create_file(&infptr, "!" test_path, &status);
	fail_if(status != 0);
	fits_create_img(infptr, SHORT_IMG, 2, naxes, &status);
	fail_if(status != 0);
	fits_write_img(infptr, TSHORT, 1, 64 * 64, original, &status);
	fail_if(status != 0);
	fits_close_file(infptr, &status);
	fail_if(status != 0);

	/* Compress with PLIO */
	fits_open_file(&infptr, test_path, READONLY, &status);
	fail_if(status != 0);

	fits_create_file(&outfptr, "!" test_path2, &status);
	fail_if(status != 0);

	fits_set_compression_type(outfptr, PLIO_1, &status);
	fail_if(status != 0);

	fits_img_compress(infptr, outfptr, &status);
	fail_if(status != 0);

	fits_close_file(infptr, &status);
	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	/* Read back - compressed image is in HDU 2 */
	fits_open_file(&outfptr, test_path2, READONLY, &status);
	fail_if(status != 0);

	fits_movabs_hdu(outfptr, 2, NULL, &status);
	fail_if(status != 0);

	fits_read_img(outfptr, TSHORT, 1, 64 * 64, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);

	/* Verify lossless */
	for (i = 0; i < 64 * 64; i += 1) {
		fail_if(decompressed[i] != original[i]);
	}

	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	free(original);
	free(decompressed);
}

/*
 * Test HCOMPRESS compression
 */
static void
test_hcompress_compress(void)
{
	fitsfile *infptr, *outfptr;
	int status = 0;
	long naxes[2] = { 64, 64 };
	int *original, *decompressed;
	int i;

	original = malloc(64 * 64 * sizeof *original);
	decompressed = malloc(64 * 64 * sizeof *decompressed);
	fail_if(original == NULL || decompressed == NULL);

	/* Create gradient data */
	for (i = 0; i < 64 * 64; i += 1) {
		original[i] = (i % 64) + (i / 64) * 100;
	}

	/* Create input image */
	fits_create_file(&infptr, "!" test_path, &status);
	fail_if(status != 0);
	fits_create_img(infptr, LONG_IMG, 2, naxes, &status);
	fail_if(status != 0);
	fits_write_img(infptr, TINT, 1, 64 * 64, original, &status);
	fail_if(status != 0);
	fits_close_file(infptr, &status);
	fail_if(status != 0);

	/* Compress with HCOMPRESS */
	fits_open_file(&infptr, test_path, READONLY, &status);
	fail_if(status != 0);

	fits_create_file(&outfptr, "!" test_path2, &status);
	fail_if(status != 0);

	fits_set_compression_type(outfptr, HCOMPRESS_1, &status);
	fail_if(status != 0);

	fits_img_compress(infptr, outfptr, &status);
	fail_if(status != 0);

	fits_close_file(infptr, &status);
	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	/* Read back - compressed image is in HDU 2 */
	fits_open_file(&outfptr, test_path2, READONLY, &status);
	fail_if(status != 0);

	fits_movabs_hdu(outfptr, 2, NULL, &status);
	fail_if(status != 0);

	fits_read_img(outfptr, TINT, 1, 64 * 64, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);

	/* Verify lossless (with default settings) */
	for (i = 0; i < 64 * 64; i += 1) {
		fail_if(decompressed[i] != original[i]);
	}

	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	free(original);
	free(decompressed);
}

/*
 * Test fits_is_compressed_image on uncompressed image
 */
static void
test_is_compressed_uncompressed(void)
{
	fitsfile *fptr;
	int status = 0;
	int is_compressed;

	create_test_image(&fptr, SHORT_IMG, 32, 32, &status);

	is_compressed = fits_is_compressed_image(fptr, &status);
	fail_if(status != 0);
	fail_if(is_compressed != 0);

	fits_close_file(fptr, &status);
	fail_if(status != 0);
}

/*
 * Test byte image compression
 */
static void
test_compress_byte_image(void)
{
	fitsfile *infptr, *outfptr;
	int status = 0;
	long naxes[2] = { 32, 32 };
	unsigned char *original, *decompressed;
	int i;

	original = malloc(32 * 32 * sizeof *original);
	decompressed = malloc(32 * 32 * sizeof *decompressed);
	fail_if(original == NULL || decompressed == NULL);

	/* Create test data */
	for (i = 0; i < 32 * 32; i += 1) {
		original[i] = (unsigned char)(i % 256);
	}

	/* Create input image */
	fits_create_file(&infptr, "!" test_path, &status);
	fail_if(status != 0);
	fits_create_img(infptr, BYTE_IMG, 2, naxes, &status);
	fail_if(status != 0);
	fits_write_img(infptr, TBYTE, 1, 32 * 32, original, &status);
	fail_if(status != 0);
	fits_close_file(infptr, &status);
	fail_if(status != 0);

	/* Compress */
	fits_open_file(&infptr, test_path, READONLY, &status);
	fail_if(status != 0);

	fits_create_file(&outfptr, "!" test_path2, &status);
	fail_if(status != 0);

	fits_set_compression_type(outfptr, RICE_1, &status);
	fail_if(status != 0);

	fits_img_compress(infptr, outfptr, &status);
	fail_if(status != 0);

	fits_close_file(infptr, &status);
	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	/* Read back - compressed image is in HDU 2 */
	fits_open_file(&outfptr, test_path2, READONLY, &status);
	fail_if(status != 0);

	fits_movabs_hdu(outfptr, 2, NULL, &status);
	fail_if(status != 0);

	fits_read_img(outfptr, TBYTE, 1, 32 * 32, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);

	/* Verify */
	for (i = 0; i < 32 * 32; i += 1) {
		fail_if(decompressed[i] != original[i]);
	}

	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	free(original);
	free(decompressed);
}

/*
 * Compress an n-dimensional image with fits_img_compress (what fpack does),
 * then read it back one pixel range at a time.  Images with more than 3
 * dimensions used to fail with DATA_DECOMPRESSION_ERR ("only 1D, 2D, or 3D
 * images are currently supported") - see heasarc/cfitsio issue #171.
 */
static void
compress_and_read_ndim(int naxis, long *naxes)
{
	fitsfile *infptr, *outfptr;
	int status = 0;
	long npix = 1;
	short *original, *decompressed;
	long i, firstelem, nelem;

	for (i = 0; i < naxis; i += 1) {
		npix *= naxes[i];
	}

	original = malloc(npix * sizeof *original);
	decompressed = malloc(npix * sizeof *decompressed);
	fail_if(original == NULL || decompressed == NULL);

	for (i = 0; i < npix; i += 1) {
		original[i] = (short)(i * 3 + 1);
	}

	/* Create the uncompressed input image */
	fits_create_file(&infptr, "!" test_path, &status);
	fail_if(status != 0);
	fits_create_img(infptr, SHORT_IMG, naxis, naxes, &status);
	fail_if(status != 0);
	fits_write_img(infptr, TSHORT, 1, npix, original, &status);
	fail_if(status != 0);
	fits_close_file(infptr, &status);
	fail_if(status != 0);

	/* Compress it, as "fpack -g" would */
	fits_open_file(&infptr, test_path, READONLY, &status);
	fail_if(status != 0);
	fits_create_file(&outfptr, "!" test_path2, &status);
	fail_if(status != 0);
	fits_set_compression_type(outfptr, GZIP_1, &status);
	fail_if(status != 0);
	fits_img_compress(infptr, outfptr, &status);
	fail_if(status != 0);
	fits_close_file(infptr, &status);
	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	/* Read the whole compressed image back */
	fits_open_file(&outfptr, test_path2, READONLY, &status);
	fail_if(status != 0);
	fits_movabs_hdu(outfptr, 2, NULL, &status);
	fail_if(status != 0);
	fail_if(fits_is_compressed_image(outfptr, &status) == 0);
	fail_if(status != 0);

	memset(decompressed, 0, npix * sizeof *decompressed);
	fits_read_img(outfptr, TSHORT, 1, npix, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);
	for (i = 0; i < npix; i += 1) {
		fail_if(decompressed[i] != original[i]);
	}

	/*
	 * Read a range that starts and ends part way through a row, so that
	 * the read is split up into partial rows and whole planes instead of
	 * being satisfied by a single section read.
	 */
	firstelem = 8;
	nelem = npix - 13;
	fail_if(nelem < 1);
	memset(decompressed, 0, npix * sizeof *decompressed);
	fits_read_img(outfptr, TSHORT, firstelem, nelem, NULL, decompressed,
		NULL, &status);
	fail_if(status != 0);
	for (i = 0; i < nelem; i += 1) {
		fail_if(decompressed[i] != original[firstelem - 1 + i]);
	}

	/* A range wholly inside one row */
	memset(decompressed, 0, npix * sizeof *decompressed);
	fits_read_img(outfptr, TSHORT, 2, 3, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);
	for (i = 0; i < 3; i += 1) {
		fail_if(decompressed[i] != original[1 + i]);
	}

	fits_close_file(outfptr, &status);
	fail_if(status != 0);

	free(original);
	free(decompressed);
}

/*
 * Write an n-dimensional image directly into a compressed HDU, using linear
 * pixel ranges, then read it back.  This exercises the write counterpart of
 * the issue #171 code path.
 */
static void
write_and_read_ndim(int naxis, long *naxes)
{
	fitsfile *fptr;
	int status = 0;
	long npix = 1;
	short *original, *decompressed;
	long i, firstelem, nelem;

	for (i = 0; i < naxis; i += 1) {
		npix *= naxes[i];
	}

	original = malloc(npix * sizeof *original);
	decompressed = malloc(npix * sizeof *decompressed);
	fail_if(original == NULL || decompressed == NULL);

	for (i = 0; i < npix; i += 1) {
		original[i] = (short)(i * 7 + 5);
	}

	fits_create_file(&fptr, "!" test_path, &status);
	fail_if(status != 0);

	/* A compressed image is a binary table, so it needs a primary HDU */
	fits_create_img(fptr, SHORT_IMG, 0, NULL, &status);
	fail_if(status != 0);

	fits_set_compression_type(fptr, GZIP_1, &status);
	fail_if(status != 0);
	fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
	fail_if(status != 0);
	fail_if(fits_is_compressed_image(fptr, &status) == 0);
	fail_if(status != 0);

	/* Write the whole image, then rewrite a range of it */
	fits_write_img(fptr, TSHORT, 1, npix, original, &status);
	fail_if(status != 0);

	firstelem = 6;
	nelem = npix - 11;
	fail_if(nelem < 1);
	for (i = 0; i < nelem; i += 1) {
		original[firstelem - 1 + i] = (short)(-i - 1);
	}
	fits_write_img(fptr, TSHORT, firstelem, nelem,
		original + firstelem - 1, &status);
	fail_if(status != 0);

	fits_close_file(fptr, &status);
	fail_if(status != 0);

	fits_open_file(&fptr, test_path, READONLY, &status);
	fail_if(status != 0);
	fits_movabs_hdu(fptr, 2, NULL, &status);
	fail_if(status != 0);

	memset(decompressed, 0, npix * sizeof *decompressed);
	fits_read_img(fptr, TSHORT, 1, npix, NULL, decompressed, NULL,
		&status);
	fail_if(status != 0);
	for (i = 0; i < npix; i += 1) {
		fail_if(decompressed[i] != original[i]);
	}

	fits_close_file(fptr, &status);
	fail_if(status != 0);

	free(original);
	free(decompressed);
}

/*
 * Test compressed images with more than 3 dimensions (issue #171)
 */
static void
test_compress_ndim_image(void)
{
	long naxes2[2] = { 7, 5 };
	long naxes3[3] = { 7, 5, 3 };
	long naxes4[4] = { 5, 4, 3, 2 };
	long naxes5[5] = { 4, 3, 2, 2, 2 };

	compress_and_read_ndim(2, naxes2);
	compress_and_read_ndim(3, naxes3);
	compress_and_read_ndim(4, naxes4);
	compress_and_read_ndim(5, naxes5);

	write_and_read_ndim(2, naxes2);
	write_and_read_ndim(3, naxes3);
	write_and_read_ndim(4, naxes4);
	write_and_read_ndim(5, naxes5);
}

/*
 * Test dither seed setting
 */
static void
test_dither_seed(void)
{
	fitsfile *fptr;
	int status = 0;
	int seed;

	create_test_image(&fptr, SHORT_IMG, 32, 32, &status);

	fits_set_dither_seed(fptr, 1234, &status);
	fail_if(status != 0);

	fits_get_dither_seed(fptr, &seed, &status);
	fail_if(status != 0);
	fail_if(seed != 1234);

	fits_close_file(fptr, &status);
	fail_if(status != 0);
}

int
main(void)
{
	test_compression_type();
	test_tile_dimensions();
	test_quantize_level();
	test_noise_bits();
	test_rice_compress_short();
	test_gzip_compress();
	test_plio_compress();
	test_hcompress_compress();
	test_is_compressed_uncompressed();
	test_compress_byte_image();
	test_compress_ndim_image();
	test_dither_seed();

	remove(test_path);
	remove(test_path2);

	return 0;
}
