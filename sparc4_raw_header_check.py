#!/bin/python3

import re
import csv
import argparse
from math import isclose
from functools import partial
from astropy.io import fits
from dataclasses import dataclass
from typing import List, Type


__version__ = '2023.04.28'


_chagelog = {
    '2023.04.28': 'First stable version. Compilant with SPARC4 keys table as '
                  'of 2023.04.28',
}


class OPDCoords:
    """Coordinates of OPD"""
    lat = -22.53444444444445
    long = -45.5825
    alt = 1864.0


@dataclass
class Field:
    """Single header card definitions."""
    dtype: Type  # data type of the keyword.: str, int, float, bool
    desc: str  # description of the keyword
    allowed_values: List = None  # list of allowed values
    range: List = None  # range of allowed values (min, max)
    re: str = None  # regular expression to check the value for str values
    re_format: str = None  # explanation of re format

    def check(self, value, comment):
        """Value checking function. Returns the error message."""
        if comment != self.desc:
            return f'Description does not match: found "{comment}" '\
                   f'expected "{self.desc}"'

        # check type
        if self.dtype is not None and \
           not isinstance(value, self.dtype):
            return f'{type(value)} do not comply {self.dtype} format.'

        # check allowed values
        if self.allowed_values is not None and \
           value not in self.allowed_values:
            return 'value is not in {self.allowed_values} allowed values.'

        # re for str values
        if self.re is not None and not re.match(self.re, value):
            return f'Value does not match {self.re_format} format.'

        # range for int and float values
        if self.range is not None:
            if not self.range[0] <= value <= self.range[1]:
                return f'Value not in range {list(self.range)}'

        return self._check(value)

    def _check(self, value):
        """Custom checking function."""
        return None


_number = (int, float)


_defs = {
    # FITS general and image
    'SIMPLE': Field(dtype=bool,
                    desc='Conforms to FITS standard'),
    'BITPIX': Field(dtype=int,
                    desc='Bits per data value',
                    allowed_values=[-64, -32, 8, 16, 32]),
    'NAXIS': Field(dtype=int,
                   desc='Number of array dimensions',
                   allowed_values=[2]),
    'NAXIS1': Field(dtype=int, desc='Number of columns'),
    'NAXIS2': Field(dtype=int, desc='Number of rows'),
    'BSCALE': Field(dtype=(int, float),
                    desc='Linear factor in scaling equation'),
    'BZERO': Field(dtype=(int, float),
                   desc='Zero point in scaling equation'),

    # Image Rect
    'VBIN': Field(dtype=int, desc='Vertical binning (pix)',
                  range=[1, 1024]),
    'HBIN': Field(dtype=int, desc='Horizontal binning (pix)',
                  range=[1, 1024]),
    'INITLIN': Field(dtype=int, desc='Initial line (pix)',
                     range=[1, 1024]),
    'INITCOL': Field(dtype=int, desc='Initial column (pix)',
                     range=[1, 1024]),
    'FINALLIN': Field(dtype=int, desc='Final line (pix)',
                      range=[1, 1024]),
    'FINALCOL': Field(dtype=int, desc='Final column (pix)',
                      range=[1, 1024]),

    # Observatory
    'OBSLAT': Field(dtype=float,
                    desc='Observatory North Latitude (DEG, -90 to 90)',
                    range=[-90, 90]),
    'OBSLONG': Field(dtype=float,
                     desc='Observatory East Longitude (DEG, -180 to 180)',
                     range=[-180, 180]),
    'OBSALT': Field(dtype=float,
                    desc='Observatory elevation above sea level (m)'),

    # Guiding Camera
    'GUIDERA': Field(dtype=str, desc='RA of guider object',
                     re=r'\d{2}:\d{2}:\d{2}',
                     re_format='HH:MM:SS'),
    'GUIDEDEC': Field(dtype=str, desc='DEC of guider object',
                      re=r'\d{2}:\d{2}:\d{2}',
                      re_format='+-DD:MM:SS'),
    'GMIR': Field(dtype=_number,
                  desc='Rotation angle of the guiding mirror (deg)',
                  range=[0, 360]),
    'GFOC': Field(dtype=_number,
                  desc='Guiding focus position (mm)'),

    # Observation
    'OBSERVER': Field(dtype=str,
                      desc='Name(s) of observer(s)'),
    'OBJECT': Field(dtype=str,
                    desc='Object name'),
    'OBSTYPE': Field(dtype=str,
                     desc='Image type: OBJECT, ZERO, FLAT, DARK, FOCUS',
                     allowed_values=['OBJECT', 'ZERO', 'FLAT',
                                     'DARK', 'FOCUS']),
    'INSTRUME': Field(dtype=str, desc='Instrument name',
                      allowed_values=['SPARC4']),
    'FILENAME': Field(dtype=str, desc='File name'),
    'DATE-OBS': Field(dtype=str,
                      desc='UTC at start of observation (isot)',
                      re=r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}',
                      re_format='YYYY-MM-DDTHH:MM:SS.SSS'),
    'UTTIME': Field(dtype=str,
                    desc='UTC time at start of observation',
                    re=r'\d{2}:\d{2}:\d{2}\.\d{3}',
                    re_format='HH:MM:SS.SSS'),
    'UTDATE': Field(dtype=str,
                    desc='UTC date at start of observation',
                    re=r'\d{4}-\d{2}-\d{2}',
                    re_format='YYYY-MM-DD'),
    'RA': Field(dtype=str, desc='Right ascension: HH:MM:SS',
                re=r'\d{2}:\d{2}:\d{2}',
                re_format='HH:MM:SS'),
    'DEC': Field(dtype=str, desc='Declination: +- DD:MM:SS',
                 re=r'[+-]?\d{2}:\d{2}:\d{2}',
                 re_format='+-DD:MM:SS'),
    'EQUINOX': Field(dtype=float, desc='Equinox  of coordinates'),
    'EXPTIME': Field(dtype=_number, desc='Exposure time (s)'),

    # SPARC4 Cycles, sequences and frames
    'NCYCLES': Field(dtype=int, desc='Number of cycles'),
    'CYCLIND': Field(dtype=int, desc='Cycle index'),
    'NFRAMES': Field(dtype=int, desc='Number of frames in sequence'),
    'FRAMEIND': Field(dtype=int, desc='Frame index in sequence'),
    'NSEQ': Field(dtype=int, desc='Number of sequences in one cycle'),
    'SEQINDEX': Field(dtype=int, desc='Sequence index in cycle'),

    # SPARC4 Instrument
    'CHANNEL': Field(dtype=int,
                     desc='Instrument channel: 1 (g), 2 (r), 3 (i), 4 (z)',
                     allowed_values=[1, 2, 3, 4]),
    'FILTER': Field(dtype=str, desc='Filter'),
    'CCDSERN': Field(dtype=int,
                     desc='CCD Serial number'),
    'INSTMODE': Field(dtype=str, desc='Instrument mode: PHOT or POLAR',
                      allowed_values=['PHOT', 'POLAR']),
    'WPSEL': Field(dtype=str,
                   desc='Waveplate: half-wave, quarter-wave or None',
                   allowed_values=['L2', 'L4', 'NONE']),
    'WPPOS': Field(dtype=int,
                   desc='Waveplate index position from 1 to 16, 0 if none.',
                   range=list(range(17))),
    'CALW': Field(dtype=str,
                  desc='Selected element in calibration wheel',
                  allowed_values=['POLARIZER', 'DEPOLARIZER', 'NONE',
                                  'PINHOLE', 'POS5']),
    'ASEL': Field(dtype=bool, desc='Savart prism analyzer selected'),

    # Detector
    'ACQMODE': Field(dtype=str, desc='Acquisition mode',
                     allowed_values=['Kinetic', 'Single Scan']),
    'PREAMP': Field(dtype=str, desc='Pre-amplifier gain',
                    allowed_values=['Gain 1', 'Gain 2']),
    'READRATE': Field(dtype=_number, desc='Readout rate (MHz)'),
    'VSHIFT': Field(dtype=float, desc='Vertical shift speed (ms)'),
    'TRIGGER': Field(dtype=str, desc='Trigger mode',
                     allowed_values=['Internal', 'External']),
    'EMMODE': Field(dtype=str, desc='Output amplifier mode',
                    allowed_values=['Conventional', 'EM']),
    'EMGAIN': Field(dtype=int, desc='Electron multiplier gain'),
    'SHUTTER': Field(dtype=str, desc='Shutter mode',
                     allowed_values=['Open', 'Closed']),
    'COOLER': Field(dtype=bool, desc='CCD cooler: T or F'),
    'CCDTEMP': Field(dtype=_number, desc='CCD temperature (deg C)'),
    'TGTEMP': Field(dtype=_number, desc='CCD target temperature (deg C)'),
    'TEMPST': Field(dtype=str, desc='CCD temperature status',
                    allowed_values=['TEMPERATURE_OFF',
                                    'TEMPERATURE_NOT_REACHED',
                                    'TEMPERATURE_NOT_STABILIZED',
                                    'TEMPERATURE_STABILIZED']),
    'FRAMETRF': Field(dtype=bool, desc='Frame transfer: T or F'),
    'VCLKAMP': Field(dtype=str, desc='Clock amplitude: Normal, +1, +2, +3, +4',
                     allowed_values=['Normal', '+1', '+2', '+3', '+4']),
    'GAIN': Field(dtype=_number, desc='Gain (e-/ADU)'),
    'RDNOISE': Field(dtype=float, desc='Read noise (e-)'),

    # Telescope and TCS
    'TELFOCUS': Field(dtype=int,
                      desc='Telescope focus position (arbitrary units)'),
    'TCSHA': Field(dtype=str, desc='TCS hour angle: HH:MM:SS',
                   re=r'\d{2}:\d{2}:\d{2}',
                   re_format='HH:MM:SS'),
    'TCSDATE': Field(dtype=str, desc='TCS UT date (isot)',
                     re=r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}',
                     re_format='YYYY-MM-DDTHH:MM:SS.SSS'),
    'INSTROT': Field(dtype=_number,
                     desc='Instrument rotator angle in degrees',
                     range=[0, 360]),

    # Weather
    'EXTTEMP': Field(dtype=_number,
                     desc='Temperature (deg C), weather tower'),
    'AIRMASS': Field(dtype=float,
                     desc='Airmass at start of observation',
                     range=[1, float('inf')]),
    'PRESSURE': Field(dtype=_number,
                      desc='Barometric pressure (mb), weather tower'),
    'HUMIDITY': Field(dtype=_number,
                      desc='Relative humidity (%), weather tower'),

    # Software
    'ACSVRSN': Field(dtype=str, desc='Software version of the ACS',
                     re=r'v[\d.]+',
                     re_format='v#.#.#'),
    'CTRLINTE': Field(dtype=str, desc='Graphical control interface of ACS',
                      allowed_values=['S4GUI', 'GEI']),
    'ACSMODE': Field(dtype=bool,
                     desc='ACS in real (T) or simulated (F) mode'),
    'ICSMODE': Field(dtype=bool,
                     desc='ICS in real (T) or simulated (F) mode'),
    'TCSMODE': Field(dtype=bool,
                     desc='TCS in  real (T) or simulated (F) mode')
}


class CustomFunc:
    @staticmethod
    def _is_positive(value):
        if not value >= 0:
            return f'value {value} is not positive'

    @staticmethod
    def _float_equal(value, expect, message, tolerance=1e-5):
        if not isclose(value, expect, abs_tol=tolerance):
            return message.format(value=value, expect=expect)


# all this keywords must be positive
for i in ['NAXIS', 'NAXIS1', 'NAXIS2', 'OBSALT', 'EXPTIME',
          'NCYCLES', 'CYCLIND', 'NFRAMES', 'NSEQ', 'SEQINDEX',
          'READRATE', 'VSHIFT', 'GAIN', 'RDNOISE', 'HUMIDITY']:
    _defs[i]._check = CustomFunc._is_positive

# float equal for values that must be fixed
_defs['OBSLAT']._check = partial(CustomFunc._float_equal, expect=OPDCoords.lat,
                                 message='{value} coordinate do not match '
                                         'OPD {expect}')
_defs['OBSLONG']._check = partial(CustomFunc._float_equal,
                                  expect=OPDCoords.long,
                                  message='{value} coordinate do not match '
                                          'OPD {expect}')
_defs['OBSALT']._check = partial(CustomFunc._float_equal, expect=OPDCoords.alt,
                                 message='{value} altitude do not match '
                                         'OPD {expect}')
_defs['EQUINOX']._check = partial(CustomFunc._float_equal, expect=2000.0,
                                  message='{value} expected to be {expect}'
                                          ' equinox')


class Printer:
    """Printer class to print to stdout or file."""
    file: str = None
    fmt: str = None  # txt or csv
    c = None  # csv writer
    pad = '  -> '

    def __init__(self, file=None, fmt='txt'):
        self.file = file
        if file is not None:
            self.file = open(file, 'w')
        self.fmt = fmt

        # only execute write when the two are set.
        if self.fmt in ['txt', 'csv'] and file is not None:
            self._do_write = True
        else:
            self._do_write = False

        if self.fmt == 'csv':
            self.c = csv.writer(self.file)

    def _get_str(self, keyword, value, error):
        """Get the string for printing and writing to file."""
        if isinstance(value, str):
            value = f'"{value}"'
        return self.pad + f'{keyword}:{value} ==> {error}'

    def print(self, filename, keyword, value, error):
        """Print to stdout and save to report."""
        s = self._get_str(keyword, value, error)
        print(s)
        if self._do_write and self.fmt == 'txt':
            self.file.write(s+'\n')
        if self._do_write and self.fmt == 'csv':
            self.c.writerow([filename, keyword, value,
                             type(value).__name__, error])

    def init_report(self, filename):
        """Print the report initialization."""
        print(filename + ' report:')
        if self._do_write and self.fmt == 'txt':
            self.file.write(filename + 'report:\n')


def main(files, printer):
    for file in files:
        printer.init_report(file)
        with fits.open(file) as hdul:
            # SPARC4 Raw images must have only one HDU
            if len(hdul) != 1:
                printer.print(file, None, None,
                              f'Only one HDU expected. {len(hdul)} found.')
                continue

            # Check header keys
            hdu = hdul[0]
            # non-std keys found
            for k in hdu.header.keys():
                if k not in _defs:
                    h_v = hdu.header[k]
                    if k not in ['COMMENT', 'HISTORY']:
                        printer.print(file, k, h_v,
                                      f'{k} not in the standard.')
            # check all std keys
            for k, v in _defs.items():
                if k not in hdu.header:
                    printer.print(file, k, None, f'{k} not found in header.')
                    continue

                h_v = hdu.header[k]
                # check description
                comm = hdu.header.comments[k]
                error = v.check(h_v, comm)
                if error is not None:
                    printer.print(file, k, h_v, error)
                    continue

            # special check for WPPOS and WPSEL
            wppos = hdu.header['WPPOS']
            wpsel = hdu.header['WPSEL']
            if wpsel == 'NONE' and wppos != 0:
                printer.print(file, 'WPPOS', wppos,
                              'WPPOS must be zero if WPSEL is NONE')
            if wpsel != 'NONE' and wppos == 0:
                printer.print(file, 'WPPOS', wppos,
                              'WPPOS must not be zero is WPSEL is not NONE')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='parc4_raw_header_checker',
        description='Check SPARC4 raw image headers standards.')
    parser.add_argument('-o', '--output', type=str,
                        help='Output file name to save the report.')
    parser.add_argument('-f', '--format', type=str,
                        choices=['txt', 'csv'], default='txt',
                        help='Format of the data to save the report. '
                        'txt: plain text with the printed messages here.'
                        ' csv: comma separated values table.')
    parser.add_argument('--changelog', action='store_true')
    parser.add_argument('files', type=str, nargs='*')

    args = parser.parse_args()
    if args.changelog:
        for k, v in _chagelog.items():
            print(k, ':', v)
    if args.files is not None:
        printer = Printer(args.output, args.format)
        main(args.files, printer)
