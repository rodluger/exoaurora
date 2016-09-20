# Detectability Calculations

Here you will find the Python scripts used to generate the detectability and integration times, Table 3 in Luger et al. (2016).

In order to reproduce Table 3, which illustrates the integration times required to detect auroral signals over a range of powers,u'll need to download the HARPS fits files. These are all publicly accessible at [the ESO Archive](http://archive.eso.org/wdb/wdb/adp/phase3_main/form). Specify `Proxima Centauri` as the **Target name**, and choose `HARPS` as the **Instrument**. Request all of the science spectra, download them, and place the fits files in the directory `HARPS/`. Then simply run the Python script `/scripts/figure2.py`.

The link to the Meadows et al. 2016 spectra we used are
[here](http://vpl.astro.washington.edu/spectra/stellar/proxcen.htm)
