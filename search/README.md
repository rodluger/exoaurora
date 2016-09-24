# search

In order to reproduce figures 2-6 in the paper, you'll have to download the HARPS data. These are all publicly accessible at [the ESO Archive](http://archive.eso.org/wdb/wdb/adp/phase3_main/form). Direct links to the spectra we used are
[here](http://archive.eso.org/wdb/wdb/adp/phase3_main/query?wdbo=html%2fdisplay&max_rows_returned=200&target=gl551&resolver=none&wdb_input_file=&coord_sys=eq&coord1=&coord2=&box=02%2009%2000&tab_ra=on&tab_dec=on&tab_filter=on&filter=Any&tab_wavelength=on&wavelength=Any&tab_dataproduct_type=on&dataproduct_type=Any&tel_id=Any&tab_ins_id=on&ins_id=HARPS&obstech=Any&tab_date_obs=on&date_obs=&mjd_obs=&tab_exptime=on&exptime=&multi_ob=%25&tab_collection_name=on&tab_prog_id=on&prog_id=&username=&p3orig=%25&tab_origfile=on&origfile=&tab_dp_id=on&dp_id=&rel_date=&tab_referenc=on&referenc=&batch_id=&publication_date=&wdb_input_file_raw=&order_main=dummy&)
and
[here](http://archive.eso.org/wdb/wdb/adp/phase3_main/query?wdbo=html%2fdisplay&max_rows_returned=800&target=proxima&resolver=none&wdb_input_file=&coord_sys=eq&coord1=&coord2=&box=02%2009%2000&tab_ra=on&tab_dec=on&tab_filter=on&filter=Any&tab_wavelength=on&wavelength=Any&tab_dataproduct_type=on&dataproduct_type=Any&tel_id=Any&tab_ins_id=on&ins_id=HARPS&obstech=Any&tab_date_obs=on&date_obs=&mjd_obs=&tab_exptime=on&exptime=&multi_ob=%25&tab_collection_name=on&tab_prog_id=on&prog_id=&username=&p3orig=%25&tab_origfile=on&origfile=&tab_dp_id=on&dp_id=&rel_date=&tab_referenc=on&referenc=&batch_id=&publication_date=&wdb_input_file_raw=&order_main=dummy&). Alternatively, you can check out [MANIFEST](MANIFEST.txt) for a complete manifest of the spectra we used in the paper. Once you've downloaded all of them, place the `fits` files in the `search/HARPS/` directory. Then simply run 

```
python search.py
```
