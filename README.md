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

In order to reproduce Figure 2, which illustrates the search in the HARPS data, you'll need to download the HARPS fits files. These are all publicly accessible at [the ESO Archive](http://archive.eso.org/wdb/wdb/adp/phase3_main/form). Specify `Proxima Centauri` as the **Target name**, and choose `HARPS` as the **Instrument**. Request all of the science spectra, download them, and place the fits files in the directory `/search/fits/HARPS/`. Then simply run the Python script `/scripts/figure2.py`.