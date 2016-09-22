#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
manifest.py
-----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
import glob
import numpy as np
import os

# Get the names of all the fits files
files = glob.glob(os.path.join('..', 'HARPS', '*.fits'))

# Get the fits headers
headers = [pyfits.getheader(f) for f in files]

# Get file info
object = [h['OBJECT'] for h in headers]
prog_id = [h['PROG_ID'] for h in headers]
pi_coi = [h['PI-COI'] for h in headers]
date_obs = [h['DATE-OBS'] for h in headers]
exptime = [h['EXPTIME'] for h in headers]
origfile = [h['ORIGFILE'] for h in headers]

# Get the total exposure time
seconds = np.sum([float(e) for e in exptime])

# Print the manifest to file
with open('../MANIFEST.txt', 'w') as file:
  print("EXOAURORA PROXIMA CENTAURI B 5577A OI LINE SEARCH", file = file)
  print("=================================================", file = file)
  print("", file = file)
  print("ORIGIN:            http://archive.eso.org/cms.html", file = file)
  print("NUMBER OF SPECTRA: %d" % len(files), file = file)
  print("TOTAL EXP TIME:    %.2f hours" % (seconds / 3600.), file = file)
  print("PROGRAM IDs:       %s" % ", ".join(sorted(list(set(prog_id)))), file = file)
  print("", file = file)
  print("DIRECT LINK 1:     http://archive.eso.org/wdb/wdb/adp/phase3_main/query?wdbo=html%2fdisplay&max_rows_returned=200&cltarget=gl551&resolver=none&wdb_input_file=&coord_sys=eq&coord1=&coord2=&box=02%2009%2000&tab_ra=on&tab_dec=on&tab_filter=on&filter=Any&tab_wavelength=on&wavelength=Any&tab_dataproduct_type=on&dataproduct_type=Any&tel_id=Any&tab_ins_id=on&ins_id=HARPS&obstech=Any&tab_date_obs=on&date_obs=&mjd_obs=&tab_exptime=on&exptime=&multi_ob=%25&tab_collection_name=on&tab_prog_id=on&prog_id=&username=&p3orig=%25&tab_origfile=on&origfile=&tab_dp_id=on&dp_id=&rel_date=&tab_referenc=on&referenc=&batch_id=&publication_date=&wdb_input_file_raw=&order_main=dummy&", file = file)
  print("DIRECT LINK 2:     http://archive.eso.org/wdb/wdb/adp/phase3_main/query?wdbo=html%2fdisplay&max_rows_returned=800&target=proxima&resolver=none&wdb_input_file=&coord_sys=eq&coord1=&coord2=&box=02%2009%2000&tab_ra=on&tab_dec=on&tab_filter=on&filter=Any&tab_wavelength=on&wavelength=Any&tab_dataproduct_type=on&dataproduct_type=Any&tel_id=Any&tab_ins_id=on&ins_id=HARPS&obstech=Any&tab_date_obs=on&date_obs=&mjd_obs=&tab_exptime=on&exptime=&multi_ob=%25&tab_collection_name=on&tab_prog_id=on&prog_id=&username=&p3orig=%25&tab_origfile=on&origfile=&tab_dp_id=on&dp_id=&rel_date=&tab_referenc=on&referenc=&batch_id=&publication_date=&wdb_input_file_raw=&order_main=dummy&", file = file)
  print("", file = file)
  print("OBJECT             PROG_ID            PI-COI             DATE-OBS                  EXPTIME           ORIGFILE", file = file)
  print("------             -------            ------             --------                  -------           --------", file = file)
  for a, b, c, d, e, f in zip(object, prog_id, pi_coi, date_obs, exptime, origfile):
    print("{:<19s}{:<19s}{:<19s}{:<26s}{:<19f}{:<s}".format(a,b,c,d,e,f), file = file)