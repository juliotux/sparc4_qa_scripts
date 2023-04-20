#!/bin/python3

import re
import csv
import argparse
from astropy.io import fits
from dataclasses import dataclass
from typing import List, Type


@dataclass
class Field:
    """Single header card definitions."""
    dtype: Type  # data type of the keyword.: str, int, float, bool
    desc: str  # description of the keyword
    allowed_values: List = None  # list of allowed values
    range: List = None  # range of allowed values (min, max)
    re: str = None  # regular expression to check the value for str values
    re_format: str = None  # explanation of re format


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
    'VBIN': Field(dtype=int, desc='Vertical binning (pix)'),
    'HBIN': Field(dtype=int, desc='Horizontal binning (pix)'),
    'INITLIN': Field(dtype=int, desc='Initial line (pix)'),
    'INITCOL': Field(dtype=int, desc='Initial column (pix)'),
    'FINALLIN': Field(dtype=int, desc='Final line (pix)'),
    'FINALCOL': Field(dtype=int, desc='Final column (pix)'),

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
    'GFOCUS': Field(dtype=None, desc='Guider focus position (mm)'),  # TODO: not defined yet
    'GUIDERA': Field(dtype=None, desc='RA of guider object'),  # TODO: not defined yet
    'GUIDEDEC': Field(dtype=None, desc='DEC of guider object'),  # TODO: not defined yet
    'GMIR': Field(dtype=None, desc='Rotation angle of the guiding mirror (deg)'),  # TODO: not defined yet
    'GFOC': Field(dtype=None, desc='Guiding focus position (mm)'),  # TODO: not defined yet

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
                 re=r'[+-]\d{2}:\d{2}:\d{2}',
                 re_format='+-DD:MM:SS'),
    'EQUINOX': Field(dtype=float, desc='Equinox  of coordinates',
                     allowed_values=[2000.0]),
    'EXPTIME': Field(dtype=(int, float), desc='Exposure time (s)',
                     range=[0, float('inf')]),

    # SPARC4 Cycles, sequences and frames
    'NCYCLES': Field(dtype=int, desc='Number of cycles'),
    'CYCLIND': Field(dtype=int, desc='Cycle index'),
    'NFRAMES': Field(dtype=int, desc='Total number of frames in sequence'),
    'FRAMEIND': Field(dtype=int, desc='Frame index in sequence'),
    'NSEQ': Field(dtype=int, desc='Total number of sequences in cycle'),
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
                   allowed_values=['L2', 'L4', 'None']),
    'WPPOS': Field(dtype=int, desc='Waveplate index position'),
    'CALW': Field(dtype=str,
                  desc='Calibration wheel: POLARIZER, DEPOLARIZER, None',
                  allowed_values=['POLARIZER', 'DEPOLARIZER', 'None']),
    'ASEL': Field(dtype=bool, desc='Savart prism analyzer selected'),

    # Detector
    'ACQMODE': Field(dtype=str, desc='Acquisition mode',
                     allowed_values=['Kinetic', 'Single Scan']),
    'PREAMP': Field(dtype=str, desc='Pre-amplifier gain',
                    allowed_values=['Gain 1', 'Gain 2']),
    'READRATE': Field(dtype=(int, float), desc='Readout rate (MHz)'),
    'VSHIFT': Field(dtype=float, desc='Vertical shift speed (ms)'),
    'TRIGGER': Field(dtype=str, desc='Trigger mode',
                     allowed_values=['Internal', 'External']),
    'EMMODE': Field(dtype=str, desc='Output amplifier mode',
                    allowed_values=['Conventional', 'EM']),
    'EMGAIN': Field(dtype=int, desc='Electron multiplier gain'),
    'SHUTTER': Field(dtype=str, desc='Shutter mode',
                     allowed_values=['Open', 'Closed']),
    'COOLER': Field(dtype=bool, desc='CCD cooler: T or F'),
    'CCDTEMP': Field(dtype=(int, float), desc='CCD temperature (deg C)'),
    'TGTEMP': Field(dtype=(int, float), desc='CCD target temperature (deg C)'),
    'TEMPST': Field(dtype=str, desc='CCD temperature status',
                    allowed_values=['TEMPERATURE_OFF',
                                    'TEMPERATURE_NOT_REACHED',
                                    'TEMPERATURE_NOT_STABILIZED',
                                    'TEMPERATURE_STABILIZED']),
    'FRAMETRF': Field(dtype=bool, desc='Frame transfer: T or F'),
    'VCLKAMP': Field(dtype=str, desc='Clock amplitude: Normal, +1, +2, +3, +4',
                     allowed_values=['Normal', '+1', '+2', '+3', '+4']),
    'GAIN': Field(dtype=(int, float), desc='Gain (e-/ADU)'),
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
    'INSTROT': Field(dtype=(int, float),
                     desc='Instrument rotator angle in degrees',
                     range=[0, 360]),

    # Weather
    'EXTTEMP': Field(dtype=(int, float),
                     desc='Temperature (deg C), weather tower'),
    'AIRMASS': Field(dtype=float,
                     desc='Airmass at start of observation'),
    'PRESSURE': Field(dtype=(int, float),
                      desc='Barometric pressure (mb), weather tower'),
    'HUMIDITY': Field(dtype=(int, float),
                      desc='Relative humidity (%), weather tower'),

    # Software
    'ACSVRSN': Field(dtype=str, desc='Software version of the ACS',
                     re=r'v[\d.]+',
                     re_format='v#.#.#'),
    'CTRLINTE': Field(dtype=str, desc='Graphical control interface of ACS',
                      allowed_values=['S4GUI']),
    'ACSMODE': Field(dtype=str, desc='ACS in simulated mode',
                     allowed_values=['Simulated', 'Real']),
    'ICSMODE': Field(dtype=str, desc='ICS in simulated mode',
                     allowed_values=['Simulated', 'Real']),
    'TCSMODE': Field(dtype=str, desc='TCS in simulated mode',
                     allowed_values=['Simulated', 'Real'])
}


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
        return self.pad + f'{keyword}:{value} {error}'

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
                if v.desc != comm:
                    printer.print(file, k, h_v,
                                  f'description does not match. '
                                  f'Expected: "{v.desc}". '
                                  f'Found "{comm}".')
                # check type
                if v.dtype is not None and \
                   not isinstance(h_v, v.dtype):
                    printer.print(file, k, h_v, f'{type(h_v)} '
                                  f'do not comply {v.dtype} format.')
                    continue
                # check allowed values
                if v.allowed_values is not None and \
                   h_v not in v.allowed_values:
                    printer.print(file, k, h_v,
                                  'value is not in '
                                  f'{v.allowed_values} allowed values.')
                # re for str values
                if v.re is not None and not re.match(v.re, h_v):
                    printer.print(file, k, h_v,
                                  f'Value does not match {v.re_format}.')
                # range for int and float values
                if v.range is not None and not v.range[0] <= h_v <= v.range[1]:
                    printer.print(file, k, h_v,
                                  f'Value not in range {list(v.range)}')


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
    parser.add_argument('files', type=str, nargs='+')

    args = parser.parse_args()
    printer = Printer(args.output, args.format)
    main(args.files, printer)
