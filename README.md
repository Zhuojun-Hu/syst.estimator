# SYST. COMPUTION TOOL

## 0.Overview

This is a helpful toolset to compute systematic errors for ATMPD in SK.
It is dedicated to the reconstruction-related uncertainties that needs to be estimated with the "scale-and-shift" method.
Atmospheric neutrino samples are used to assess these uncertainties because no other control sample spans the same energy range and event topologies.

## 1.Scale-and-Shift Method
In the next, I will introduce the scale-and-shift method with a walk-through example of estimating the ring separation uncertainty for sub-GeV $e$-like samples with momentum $p_e$ below 400 MeV/c.
The MC distribution of a reconstruction likelihood $x$, used to separate signal and background in the event classification, is fitted to data.
The MC events are first labeled as true signal or background, and their respective $x$ distributions are fitted to data independently by applying a linear transformation,
```math
\begin{equation}
  x' = a x + b
\end{equation}
```
Four parameters will be fitted for one sample, including two scaling parameters $a_1$ and $a_2$, and two shifting parameters $b_1$ and $b_2$, for signal and background respectively.

In this case we are discussing, the ring counting likelihood distributions will be fitted.
The signal is composed of true single-ring events, and the background is true multi-ring events.
The signal and background distributions is first normalized with a factor of data count over MC count.
Then, they are used as templated to fit the data distribution as described above.

After obtaining the best-fit values, we'd like to generate toy MC distributions around these best-fit values.
We want to cover the 5$\sigma$ intervals for both sides of the best-fit values.
So, first we need to find out the 5$\sigma$ intervals.

There is an interactive fitter 'interactiveFitter.ipynb' written in python on jupyter notebook for this purpose.
The reason we need an interactive fitter is that the minimizer we are using does not produce an accurate 1$\sigma$ error.
Instead of using the minimizer to determine the 5$sigma$ interval, we adopt the results from the minimizer as a seed and do grid scans.
For each sample, a grid is built wiht all four parameters, and a $chi^2$ value is computed for each grid point.
We may need several iteration to finally find the 5$sigma$ intervals.

There is an application to generate distributions of number of positive and negative events.
In the previous step, we've found the 5$sigma$ intervals for the four parameters.
We then define a grid that covers the 5$sigma$ intervals, such that the parameters are varied.
These variations shift events across the threshold used for selection, producing distributions of positive and negative counts.
These distributions are weighted by probability density functions (PDFs), constructed based on the scanned $Delta$$chi^2$ distributions following Wilks’ theorem.
The weighted distributions are further fitted by Gaussian functions to obtain means $N_{gaus}$ and standard variations $sigma_{gaus}$. 
The corresponding systematic uncertainty $epsilon$ for sub-GeV $e$-like $p<400$ 1 ring event selection is given by,
```math
\begin{equation}
  epsilon = |N_{nominal}-N_{gaus}|+sigma_{gaus} / N_{nominal},
\end{equation}
```
where $N_{nominal}$ is the nominal 1 ring event count in sub-GeV $e$-like samples with momentum $p_e$ below 400 MeV/$c$.
Similarly, the systematic uncertainty for multi ring event selection is found using the Gaussian fit results from the negative distribution.

## 2.Usage of Interactive Fitter

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

Step 4: Run the interactive fitter
Double click "interactiveFitter.ipynb" on the left panel.
Run cell [1] to load modules and functions.
Modify the input files and histogram names as required in cell [2] and run it.
Run cell [3] to launch the minimizer and obtain seeds for grid scan.
Run cell [4] for the first grid scan, and cell [5] for iterative scan.
Finally, run cell [6] to obtain the $chi^2$ distributions.

## 3.Usage of genSystErr

This app varies $a$ and $b$ and generates distributions for positive and negative counts given a cutoff.
It then weights the distributions with their corresponding $chi^2$ PDFs, before doing Gaussian fits.
Finally it computes systematic uncertainties for each sample, positive and negative, based on the Gaussian fit results.

Step 1: Define Card File
Use 
```
~/Cards/sk1.fqRC.card
```
and a template and draft your own card file.

Step 2: Compile
Simply do make.

Step 3: Run genSystErr
There are two stages, stage 1:
```
./genSystErr ../Cards/yourcardfile 1
```
generates the positive and negative counts and store them in a Tree, written to a .root file.

Stage 2:
```
./genSystErr ../Cards/yourcardfile 2
```
weights the distributions with $chi^2$ PDFs and fits them with Gaussian functions.
The weighted distributions will be stored back to the same .root file.
The final results will be printed to a text file.
