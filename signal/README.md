# auroral signal calculations

here are the scripts to calculate the expected auroral power for the 5577 line for method 1 (gnuplot script) and method 2 (python scripts) in luger et al. (2016).

fig. 1, corresponding to method 1, can be produced by running the gnuplot script below (the calculation is done in the plotting script)

```
gnuplot power_m1.plt
```

note: the gnuplot script above was using gnuplot 5.0p3. If you get errors, it is most likely due to the line style argument 'dt', which is not available in older version (4.2, 4.4, etc.). You can just remove the 'dt #' from each line and it will run in older versions. Alternatively, you can just define your own, preferred line types.

using method 2 (following wang et al (2014), steele & mcewen (1990)), the numbers in tab. 2 of luger et al. (2016) can be calculated invoking the script below:

```
python power_m2.py
```

which uses constants and functions defined in

```
auroral_signal.py
```

note: the numbers calculated by the script power\_m2.py for the neptune-like dipole are ~20% larger than those in luger et al. (2016). these are the correct values, and the paper will be corrected in revision.


scripts relevant to the estimated auroral power section of the paper.
