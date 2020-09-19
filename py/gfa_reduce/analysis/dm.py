import numpy as np
from desimeter.time import mjd2lst
from desimeter.fieldmodel import FieldModel
import copy
from astropy.table import Table

def fit_dm_fieldmodel(header, ccds, _catalog):
    fm = FieldModel()

    catalog = copy.deepcopy(_catalog)

    catalog = catalog[catalog['ANG_SEP_DEG'] < 2.0/3600.0]

    if len(np.unique(catalog['EXTNAME'])) < 2:
        print('not enough sources with Gaia counterparts for a FieldModel')
        return None
        
    good_ccds = ccds[ccds['contrast'] > 2]['EXTNAME']
    
    # restrict to guide cameras with (contrast > 2)
    # if < 2 such guide cameras, then give up

    # restrict catalog to include only sources from
    # (contrast > 2) guide cameras

    keep = [(_extname in good_ccds) for _extname in catalog['extname']]
    
    # give up if fewer than some minimum number of stars ??

    catalog = catalog[keep]

    if len(np.unique(catalog['EXTNAME'])) < 2:
        print('not enough sources with Gaia counterparts for a FieldModel')
        return None

    fm.ra  = header['TARGTRA'] # other option is SKYRA
    fm.dec = header['TARGTDEC'] # other option is SKYDEC
    fm.expid = header['EXPID']
    fm.hexrot_deg = float(header['FOCUS'][5])/3600.0 # change this?

    fm.adc1 = header['ADC1PHI']
    fm.adc2 = header['ADC2PHI']

    fm.mjd = np.mean(catalog['mjd_obs'])
    fm.lst = mjd2lst(fm.mjd)

    # need to have catalog trimmed to matches w/ in 2 asec of
    # a Gaia counterpart at some point

    catalog = Table(catalog)

    fm.fit_tancorr(catalog)

    # any other metadata that I want to record?
    fm.n_cameras = len(np.unique(catalog['extname']))

    return fm
