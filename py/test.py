import astropy.io.fits as fits
import gfa_reduce.gfa_red as gfa_red

#gfa_targets = fits.getdata('/exposures/desi/20200124/00043860/fiberassign-063523.fits', extname='GFA_TARGETS')
#hdul = fits.open('gfa-43860.3.fits')
# module use /software/datasystems/desiconda/20191002/modulefiles/
# module load desiconda
# module load desi/master
# module load nightwatch/master
# export PYTHONPATH=${PYTHONPATH}:/n/home/datasystems/users/ameisner/latest/gfa_reduce/py
# export GFA_REDUCE_META=/n/home/datasystems/users/ameisner/gfa_reduce_etc
# export GAIA_CAT_DIR=/n/home/datasystems/users/ameisner/chunks-gaia-dr2-astrom
# export PS_CAT_DIR=/n/home/datasystems/users/ameisner/chunks-fitscat_qz_trim
# export PYTHONPATH=${PYTHONPATH}:/n/home/datasystems/users/ameisner/desimeter/py
# export PATH=$PATH:/n/home/datasystems/users/ameisner/desimeter/bin
# export DESIMETER_DATA=/n/home/datasystems/users/ameisner/desimeter/py/desimeter/data

def main():
    base = '/Users/dstn/data/'
    #gfa_targets = fits.getdata(base+'fiberassign-063523.fits', extname='GFA_TARGETS')
    #hdulist = fits.open(base+'gfa-43860.3.fits')
    import os
    os.environ['GFA_REDUCE_META'] = base+'gfa_reduce_etc'

    from collections import OrderedDict

    gfa_images = OrderedDict()
    
    import fitsio
    F = fitsio.FITS(base+'20201215/00067956/gfa-00067956.fits.fz')
    #gfa_targets = fitsio.read(base+'20201215/00067956/fiberassign-080605.fits',
    #                          extname='GFA_TARGETS')
    #len(gfa_targets)

    gfa_targets = fits.getdata(base+'20201215/00067956/fiberassign-080605.fits',
                               extname='GFA_TARGETS')
    
    devices = []
    gfa_data = [None]
    ccdtemp = [None]
    for hdu in F:
        info = hdu.get_info()
        ext = info.get('extname','')
        if not ext.startswith('GUIDE'):
            continue
        print('EXTNAME', ext)
        devices.append(ext)
        gfa_data.append(hdu.read())
        hdr = hdu.read_header()
        ccdtemp.append(hdr['GCCDTEMP'])
    hdr = F['GFA'].read_header()
    for item in 'EXPID EXPFRAME TILEID SEQID REQRA REQDEC FLAVOR'.split():
        if item == 'SEQID':
            gfa_images[item] = 'xxx'
        else:
            gfa_images[item] = hdr[item]
    for item in 'ST MJD-OBS DATE-OBS TIME-OBS NIGHT EXPTIME OPENSHUT REQTIME'.split():
        if item == 'TIME-OBS':
            v = hdr['DATE-OBS'].split('T')[1]
            gfa_images[item] = v
        elif item == 'REQTIME':
            gfa_images[item] = hdr['EXPTIME']
        else:
            gfa_images[item] = hdr[item]

    # added
    for item in 'FOCUS ADC1PHI ADC2PHI'.split():
        gfa_images[item] = hdr[item]

    ## per-CCD?!
    #hdr = F['GUIDE0'].read_header()
    #for item in 'GCCDTEMP'.split():
    #    gfa_images[item] = hdr[item]
    
    gfa_images['data'] = gfa_data
    gfa_images['ext_name'] = devices
    gfa_images['ccdtemp'] = ccdtemp

    # We get gfa_images

    from astropy.io.fits import HDUList, PrimaryHDU, ImageHDU, Header
    
    extnames   = gfa_images.pop('ext_name')
    image_data = gfa_images.pop('data')
    ccdtemps   = gfa_images.pop('ccdtemp')
    # drop the None off the front of the image_data, ccdtemps
    image_data.pop(0)
    ccdtemps.pop(0)
    
    hdus = []

    phdr = Header()
    for k,v in gfa_images.items():
        phdr[k] = v

    # HACK --
    #phdr['SKYRA']  = rcd['header']['SKYRA']
    #phdr['SKYDEC'] = rcd['header']['SKYDEC']
    phdr['TARGTRA' ] = phdr['REQRA']
    phdr['TARGTDEC'] = phdr['REQDEC']
    phdr['SKYRA' ] = phdr['REQRA']
    phdr['SKYDEC'] = phdr['REQDEC']
        
    hdus.append(PrimaryHDU(None, header=phdr))

    for extname,data,ccdtemp in zip(extnames, image_data, ccdtemps):
        hdr = Header()
        hdr['EXTNAME'] = extname
        hdr['GCCDTEMP'] = ccdtemp
        hdus.append(ImageHDU(data, hdr))
        print(extname)

    from collections import Counter
    hdulist = HDUList(hdus)
    #print('HDUlist:', hdulist)
    #print('GFA targets epoch:', gfa_targets['REF_EPOCH'])
    #print('GFA ref_cat:', Counter(gfa_targets['REF_CAT']))
    
    fm = gfa_red.acquire_field(gfa_targets=gfa_targets, exp_data=hdulist)

    import numpy as np
    gfa_center = OrderedDict([('ra', fm.ra),
                              ('dec', fm.dec),
                              ('hexrot', fm.hexrot_deg),
                              ('hexrate', 0.),])
    print('gfa_center:', gfa_center)

    wcsvals = []
    #for hdu in hdulist:
    for hdu in hdus:
        hdr = hdu.header
        #print('Header:', hdr)
        # Strip "GUIDE" off "GUIDE2" to return 2.
        ext = hdr.get('EXTNAME')
        print('EXTNAME', ext)
        if ext is None:
            continue
        if not ext.startswith('GUIDE'):
            continue

        print('EXT', ext)
        print('Header:', hdr)

        try:
            guide_loc = int(ext[5:], 10)
            wcsvals.append(tuple([guide_loc] + [hdr[k] for k in [
                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']]))
        except KeyError:
            print('Failed to find WCS headers for', ext)
            import traceback
            traceback.print_exc()

    print('wcsvals:', wcsvals)
    guider_wcs = np.rec.array(wcsvals,
         dtype=[('GFA_LOC', '<i8'), ('CRVAL1', '<f8'), ('CRVAL2', '<f8'),
                ('CRPIX1', '<f8'), ('CRPIX2', '<f8'),
                ('CD1_1', '<f8'), ('CD1_2', '<f8'),
                ('CD2_1', '<f8'), ('CD2_2', '<f8')])
    print('guider_wcs:', guider_wcs)


if __name__ == '__main__':
    main()
