# SYST. COMPUTION TOOL

## 0.Overview

This is a helpful toolset to compute systematic errors for ATMPD in SK.
It is dedicated to the reconstruction-related uncertainties that needs to be estimated with the "scale-and-shift" method.
Atmospheric neutrino samples are used to assess these uncertainties because no other control sample spans the same energy range and event topologies.

## 1.Scale-and-Shift Method
In the next, I will introduce the scale-and-shift method with a walk-through example of estimating the ring separation uncertainty for sub-GeV $e$-like samples with momentum $p_e$ below 400 MeV/$c$.
The MC istribution of a reconstruction variable $x$, used to separate signal and background in the event classification, is fitted to data.
The MC events are first labeled as true signal or background, and their respective $x$ distributions are fitted to data independently by applying a linear transformation,
```math
\begin{equation}
  x' = a x + b
\end{equation}
```
Four parameters will be fitted for one sample, including two scale parameters $a_1$ and $a_2$, and two shift parameters $b_1$ and $b_2$, for signal and background respectively.

In this case we are discussing, the ring counting likelihood distributions will be fitted.
The signal is composed of true single-ring events, and the background is true multi-ring events.
... TO BE CONTINUED




Step 1: Start a jupyter lab server at port 9527 (or what ever you like)
```
source start_jupyter_server.csh
```
You should see an output like
```
http://localhost:9527/lab?token=XXXXXXXXXXX
```
---
Step 2: Make a new connection to sukap, check out port 9527. For example, if you use VS code, make a new connection:
```
Ctrl+Alt+P
```
then
```
Remote-SSH: Connect to Host...
```
If you're doint this for the first time, you need to create SSH port forwarding from your local machine:
```
ssh -L 9527:localhost:9527 username@sukap02
```
Otherwise:
```
ssh username@sukap02
```
---
Step 3: From ports, find 
```
localost:9527
```
from Forwarded Address.
Then select 
```
Open in browser
```
Enter the token from step 1. You should be able to open a jupyter lab on VS code browser.
