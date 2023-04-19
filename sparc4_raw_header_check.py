#!/bin/python3

import sys
import re
from astropy.io import fits
from dataclasses import dataclass
from typing import List, Type


pad = '  -> '


# TODO: use argparse
# TODO: output format, text, csv
# TODO: output file


@dataclass
class Field:
    dtype: Type  # data type of the keyword. Fits only accepts str, int, float, bool
    desc: str  # description of the keyword
    allowed_values: List = None  # list of allowed values
    range: List = None  # range of allowed values (min, max)
    range_error: str = None  # error message if the value is out of range
    re: str = None  # regular expression to check the value for str values
    re_error: str = None  # error message if the value does not match the regular expression


_defs = {
    # FITS general and image
    'SIMPLE': Field(dtype=bool,
                    desc='conforms to FITS standard'),
    'BITPIX': Field(dtype=int,
                    desc='array data type',
                    allowed_values=[8, 16, 32]),
    'NAXIS': Field(dtype=int,
                   desc='number of array dimensions',
                   allowed_values=[2]),
    'NAXIS1': Field(dtype=int, desc=''),
    'NAXIS2': Field(dtype=int, desc=''),
    'VBIN': Field(dtype=int, desc='Vertical binning (pix)'),
    'HBIN': Field(dtype=int, desc='Horizontal binning (pix)'),
    'INITLIN': Field(dtype=int, desc='Initial line (pix)'),
    'INITCOL': Field(dtype=int, desc='Initial column (pix)'),
    'FINALLIN': Field(dtype=int, desc='Final line (pix)'),
    'FINALCOL': Field(dtype=int, desc='Final column (pix)'),
    'BSCALE': Field(dtype=(int, float), desc=''),  # type depends on BITPIX
    'BZERO': Field(dtype=(int, float), desc=''),  # type depends on BITPIX

    # Observatory
    'OBSLAT': Field(dtype=float,
                    desc='Observatory North Latitude (DEG, -90 to 90)',
                    range=[-90, 90],
                    range_error='Latitude must be between -90 and 90'),
    'OBSLONG': Field(dtype=float,
                     desc='Observatory East Longitude (DEG, 0 to 360)',
                     range=[0, 360],
                     range_error='Longitude must be between 0 and 360'),

    # Guiding Camera
    'GFOCUS': Field(dtype=None, desc='Guider focus position (mm)'),  # TODO: not defined yet
    'GUIDERA': Field(dtype=None, desc='RA of guider object'),  # TODO: not defined yet
    'GUIDEDEC': Field(dtype=None, desc='DEC of guider object'),  # TODO: not defined yet
    'GMIR': Field(dtype=None, desc='Rotation angle of the guiding mirror  (deg)'),  # TODO: not defined yet
    'GFOC': Field(dtype=None, desc='Guiding focus position (mm)'),  # TODO: not defined yet

    # Observation
    'OBSERVER': Field(dtype=str,
                      desc='Name(s) of observer(s)'),
    'OBJECT': Field(dtype=str,
                    desc='Object name'),
    'OBSTYPE': Field(dtype=str,
                     desc='Image type: OBJECT, BIAS, FLAT, DARK, FOCUS',
                     allowed_values=['OBJECT', 'BIAS', 'FLAT', 'DARK', 'FOCUS']),
    'FILENAME': Field(dtype=str, desc='File name'),  # TODO: include RE?
    'DATE-OBS': Field(dtype=str,
                      desc='UTC at start of observation (isot)',
                      re=r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}',
                      re_error='Date format must be YYYY-MM-DDTHH:MM:SS.SSS'),
    'UTTIME': Field(dtype=str,
                    desc='UTC time at start of observation',
                    re=r'\d{2}:\d{2}:\d{2}\.\d{3}',
                    re_error='Time format must be HH:MM:SS.SSS'),
    'UTDATE': Field(dtype=str,
                    desc='UTC date at start of observation',
                    re=r'\d{4}-\d{2}-\d{2}',
                    re_error='Date format must be YYYY-MM-DD'),
    'RA': Field(dtype=str, desc='Right ascension: HH:MM:SS.SSS',
                re=r'\d{2}:\d{2}:\d{2}\.\d{3}',
                re_error='RA format must be HH:MM:SS.SSS'),
    'DEC': Field(dtype=str, desc='Declination: +- DD:MM:SS.SSS',
                 re=r'[+-]\d{2}:\d{2}:\d{2}\.\d{3}',
                 re_error='DEC format must be +-DD:MM:SS.SSS'),
    'EQUINOX': Field(dtype=float, desc='Equinox  of coordinates',
                     allowed_values=[2000.0]),
    'EXPTIME': Field(dtype=(int, float), desc='Exposure time (s)'),

    # SPARC4 Specific
    'CHANNEL': Field(dtype=int,
                     desc='Instrument channel: 1 (g), 2 (r), 3 (i), 4 (z)',
                     allowed_values=[1, 2, 3, 4]),
    'CCDSERN': Field(dtype=int,
                     desc='CCD Serial number'),
    'NCYCLES': Field(dtype=int, desc='Number of cycles'),
    'CYCLIND': Field(dtype=int, desc='Cycle index in sequence'),
    'NFRAMES': Field(dtype=int, desc='Number of frames in the cycle'),
    'FRAMEIND': Field(dtype=int, desc='Frame index in the cycle'),
    'CYCLTEXP': Field(dtype=float, desc='Total cycle exposure time (s)'),
    'NSEQ': Field(dtype=int, desc='Total number of exposures in sequence'),  # special check
    'SEQINDEX': Field(dtype=int, desc='Exposure index in sequence'),
    'INSTMODE': Field(dtype=str, desc='Instrument mode: PHOT or POLAR',
                      allowed_values=['PHOT', 'POLAR']),
    'WPSEL': Field(dtype=str,
                   desc='Waveplate: L2 (​​half-wave), L4 (quarter-wave), or None',
                   allowed_values=['L2', 'L4', 'None']),
    'WPPOS': Field(dtype=int, desc='Waveplate index position'),
    'CALW': Field(dtype=str,
                  desc='Calibration wheel: POLARIZER, DEPOLARIZER, None',
                  allowed_values=['POLARIZER', 'DEPOLARIZER', 'None']),
    'ASEL': Field(dtype=bool, desc='Savart prism analyzer selected'),

    # Detector
    'ACQMODE': Field(dtype=str, desc='Acquisition mode',
                     allowed_values=['Kinetic', '']),  # TODO: all passible values
    'PREAMP': Field(dtype=str, desc='Pre-amplifier gain'),  # TODO: include format or RE?
    'READRATE': Field(dtype=str, desc='Readout rate (MHz)'),  # TODO: include format or RE?
    'VSHIFT': Field(dtype=float, desc='Vertical shift speed (ms)'),
    'TRIGGER': Field(dtype=str, desc='Trigger mode',
                     allowed_values=['Internal', 'External']),  # TODO: all passible values
    'EMMODE': Field(dtype=str, desc='Output amplifier mode',
                    allowed_values=['Conventional']),  # TODO: all passible values
    'EMGAIN': Field(dtype=int, desc='Electron multiplier gain'),
    'SHUTTER': Field(dtype=str, desc='Shutter mode',
                     allowed_values=['Open', 'Closed']),  # TODO: all passible values
    'COOLER': Field(dtype=bool, desc='CCD cooler: T or F'),
    'CCDTEMP': Field(dtype=(int, float), desc='CCD temperature (deg C)'),
    'TGTEMP': Field(dtype=(int, float), desc='CCD Target temperature (deg C)'),
    'TEMPST': Field(dtype=str, desc='Temperature status',
                    allowed_values=['TEMPERATURE_STABILIZED']),  # TODO: all passible values
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
                   re_error='TCSHA format must be HH:MM:SS'),
    'TCSDATE': Field(dtype=str, desc='TCS UT date: DD/MM/YY HH:MM:SS',
                     re=r'\d{2}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}',
                     re_error='TCSDATE format must be DD/MM/YY HH:MM:SS'),
    'INSTROT': Field(dtype=(int, float),
                     desc='Instrument rotator angle in degrees',
                     range=[0, 360],
                     range_error='INSTROT must be between 0 and 360'),

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
                     re_error='ACSVRSN format must be v#.#.#'),
    'CTRLINTE': Field(dtype=str, desc='Graphical control interface of the ACS',
                      allowed_values=['S4GUI']),
    'ACSMODE': Field(dtype=bool, desc='ACS in simulated mode'),
    'ICSMODE': Field(dtype=bool, desc='ICS in simulated mode'),
    'TCSMODE': Field(dtype=bool, desc='TCS in simulated mode')
}


def main():
    for file in sys.argv[1:]:
        print(file + ' report:')
        with fits.open(file) as hdul:
            # SPARC4 Raw images must have only one HDU
            if len(hdul) != 1:
                print(pad + f'Only one HDU expected. {len(hdul)} found.')
                continue

            # Check header keys
            hdu = hdul[0]
            # non-std keys found
            for k in hdu.header.keys():
                if k not in _defs:
                    if k not in ['COMMENT', 'HISTORY']:
                        print(pad + f'{k} not in the standard.')
            # check all std keys
            for k, v in _defs.items():
                if k not in hdu.header:
                    print(pad + f'{k} not found in header.')
                    continue

                h_v = hdu.header[k]
                # check description
                if v.desc != hdu.header.comments[k]:
                    print(pad + f'{k} description does not match. '
                          f'Expected: "{v.desc}". '
                          f'Found "{hdu.header.comments[k]}".')
                # check type
                if v.dtype is not None and \
                   not isinstance(h_v, v.dtype):
                    print(pad + f'{k}: {h_v} {type(h_v)} '
                          f'do not comply {v.dtype} format.')
                    continue
                # check allowed values
                if v.allowed_values is not None and \
                   h_v not in v.allowed_values:
                    print(pad + f'{k} value {h_v} is not in '
                          f'{v.allowed_values} allowed values.')
                # re for str values
                if v.re is not None and not re.match(v.re, h_v):
                    print(pad + f'{k}: {h_v} -> {v.re_error}.')
                # range for int and float values
                if v.range is not None and not v.range[0] <= h_v <= v.range[1]:
                    print(pad + f'{k}: {h_v} -> {v.range_error}.')


if __name__ == '__main__':
    main()
