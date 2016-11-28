#
# This is a simple gnuplot script to
# generate Figure 1, auroral power output
# for the 5577 OI line vs. planetary magnetic
# dipole momeny
# 
# @author matilley [Matt A. Tilley, University of Washington]
# @email mtilley (at) uw (dot) edu
#
reset
set terminal postscript eps color enhanced dash size 10cm,10cm
set output "plot_swpower.eps"

#
# Set figure margins
#
set lmargin at screen 0.18
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.95

#
# Set axis, tics, key parameters
#
set key Right at screen 0.64,0.9 font "Times, 22"
set xlabel "log_{10} Magnetic dipole moment (M_{/Symbol \305})" font "Times,24" offset 0,-1
set ylabel "log_{10} 5577 @^{/Symbol \ \260}A Power Emitted (TW)" font "Times,24" offset -2
set logscale x
set xrange[1.e-1:3.e1]
set yrange[-4:0.]
set xtics font "Times,22"
set ytics font "Times,22"

#
# physical parameters
#
pi=3.1415926
m_p=1.67e-27
k_B=1.38e-23
mu_0=4*pi*1.e-7

# 
# Set stellar wind parameters for each case
#

# -- sub-Alfvenic @ Prox Cen
nsub=433.e6
vsub=6.35e5
Tsub=3.42e5
Bsub=8.25e-7
# -- super-Alfvenic @ Prox Cen
nsup=12895.e6
vsup=2.10e5
Tsup=4.77e5
Bsup=2.48e-7
# Earth @ 1 AU @ Sun
nsun=4.e6
vsun=3.8e5
Tsun=3.1e4
Bsun=5.e-9

#
# Planetary SW->aurora efficiency
#
epsE=0.01

#
# Assumed 5577 A OI fraction of auroral output
#
effOI=0.02

# Power calculation
P(n,v,T,B,eps,x)=1.e-12*effOI*eps*m_p*n*v**3.*pi*( (x*8.e15)**2./(2*mu_0*(n*m_p*v**2.) ) )**0.333333

#
# line styles
#
set style line 1 lt 1 lw 4 linecolor rgb "red"
set style line 2 lt 1 dt 3 lw 4 linecolor rgb "red"
set style line 3 lt 1 dt 4 lw 4 linecolor rgb "black"

#
# Lines denoting Earth/Neptune's dipole moments
#
set style line 11 lt 8 dt 2 lw 4 linecolor rgb "black"
set arrow from 1.,-4 to 1.,-0.9 nohead ls 11
set arrow from 1.,-0.2 to 1.,0. nohead ls 11

#set label at 1,8.84,10 ""  point pointtype 7 pointsize 3.5 lc rgb "#006600" front
set label at 0.785,-3.17,10 "{/Symbol \305}"  font "Arial,48" textcolor rgb "black" front


#
# Plot each case
#
plot log10(P(nsub,vsub,Tsub,Bsub,epsE,x)) title "ProxCen Sub-Alf" with lines ls 1, \
log10(P(nsup,vsup,Tsup,Bsup,epsE,x)) title "ProxCen Sup Alf" with lines ls 2, \
log10(P(nsun,vsun,Tsun,Bsun,epsE,x)) title "Earth/Sun" with lines ls 3
