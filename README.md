# exoaurora

_There are strange things done in the midnight sun<br>
By the men who moil for gold; <br>
The Arctic trails have their secret tales <br>
That would make your blood run cold; <br>
The Northern Lights have seen queer sights, <br>
But the queerest they ever did see <br>
Was that night on the marge of Lake Lebarge <br>
I cremated Sam McGee._<br>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;__excerpt from *The Cremation of Sam McGee* by R. W. Service (1907)__

Here you will find the Python scripts used to generate the figures in Luger et al. (2016), as well as general utilities for predicting/observing exo-aurorae.

In order to reproduce Figure 2, which illustrates the search in the HARPS data, you'll need to download the HARPS fits files. These are all publicly accessible at [the ESO Archive](http://archive.eso.org/wdb/wdb/adp/phase3_main/form). Specify `Proxima Centauri` as the **Target name**, and choose `HARPS` as the **Instrument**. Request all of the science spectra, download them, and place the fits files in the directory `HARPS/`. Then simply run the Python script `/scripts/figure2.py`.

The direct links to the spectra we used are
[here](http://archive.eso.org/wdb/wdb/adp/phase3_main/query?wdbo=html%2fdisplay&max_rows_returned=200&target=gl551&resolver=none&wdb_input_file=&coord_sys=eq&coord1=&coord2=&box=02%2009%2000&tab_ra=on&tab_dec=on&tab_filter=on&filter=Any&tab_wavelength=on&wavelength=Any&tab_dataproduct_type=on&dataproduct_type=Any&tel_id=Any&tab_ins_id=on&ins_id=HARPS&obstech=Any&tab_date_obs=on&date_obs=&mjd_obs=&tab_exptime=on&exptime=&multi_ob=%25&tab_collection_name=on&tab_prog_id=on&prog_id=&username=&p3orig=%25&tab_origfile=on&origfile=&tab_dp_id=on&dp_id=&rel_date=&tab_referenc=on&referenc=&batch_id=&publication_date=&wdb_input_file_raw=&order_main=dummy&)
and
[here](http://archive.eso.org/wdb/wdb/adp/phase3_main/query?wdbo=html%2fdisplay&max_rows_returned=800&target=proxima&resolver=none&wdb_input_file=&coord_sys=eq&coord1=&coord2=&box=02%2009%2000&tab_ra=on&tab_dec=on&tab_filter=on&filter=Any&tab_wavelength=on&wavelength=Any&tab_dataproduct_type=on&dataproduct_type=Any&tel_id=Any&tab_ins_id=on&ins_id=HARPS&obstech=Any&tab_date_obs=on&date_obs=&mjd_obs=&tab_exptime=on&exptime=&multi_ob=%25&tab_collection_name=on&tab_prog_id=on&prog_id=&username=&p3orig=%25&tab_origfile=on&origfile=&tab_dp_id=on&dp_id=&rel_date=&tab_referenc=on&referenc=&batch_id=&publication_date=&wdb_input_file_raw=&order_main=dummy&).
