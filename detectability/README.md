# Detectability Calculations

Here you will find the Python scripts used to generate the detectability and integration times, Table 3 in Luger et al. (2016).

In order to reproduce Table 3, which illustrates the integration times required to detect auroral signals over a range of powers, run the following:

```
python int_time.py
```

Note: The default coronagraph design ratio is 1.0e-10, the value adopted for both the LUVOIR and the HabEX columns in Table 3.  To compute the integration times for ground-based coronagraphic telescope configurations, the coronagraph design ratio parameter can be changed in the file ```coronagraph_inputs.py```.

Aditional calculations are the TRAPPIST-1 and Proxima Centauri visible contrast ratio comparisons and the blackbody and Meadows et al. 2016 spectrum comparison.  These calculations can be reproduced by running

```
python trappist_est.py
```

and

```
blackbody_real.py
```

respectively.

The Proxima Centauri stellar spectrum used in this work is from Meadows et al. 2016 and can be found [here](http://vpl.astro.washington.edu/spectra/stellar/proxcen.htm).
