#! /usr/bin/env python

# Author: Victor Terron (c) 2015
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

description = """
Extract an extension from a multi-extension FITS file and write it as a
standalone FITS file. This new file inherits the keywords from the primary
header, regardless of the presence or not of the INHERIT keyword. Thanks to
Javier Blasco for the suggestion and help writing this script.
"""

import argparse
import astropy.io.fits
import os

def extract_extension(input_path, output_path, n):

    with astropy.io.fits.open(input_path, mode='readonly') as hdu:
        primary_header   = hdu[0].header
        extension_data   = hdu[n].data
        extension_header = hdu[n].header
        extension_header += primary_header

    astropy.io.fits.writeto(output_path, extension_data, extension_header,
                            output_verify='fix')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_path', metavar='INPUT_PATH', type=str,
                        help="the multi-extension FITS file")
    parser.add_argument('output_path', metavar='OUTPUT_PATH', type=str,
                        help="the file to which to write the extension")
    parser.add_argument('nextension', metavar='N', type=int,
                        help="the extension to extract (1-based)")
    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help="overwrite OUTPUT_PATH if it already exists")

    args = parser.parse_args()

    if os.path.exists(args.output_path) and args.overwrite:
        os.unlink(args.output_path)

    extract_extension(args.input_path, args.output_path, args.nextension)
