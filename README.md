## [BEAST 2](http://beast2.org) package for the PhyloPotts model

## Install through BEAUti

# Installing the package

PhyloPotts is a [BEAST2](http://beast2.org) package that requires BEAST 2 v2.7.
If you have not already done so, you can get BEAST 2 from [here](http://beast2.org).

To install PhyloPotts, it is easiest to start BEAUti (a program that is part of BEAST), and select the menu `File/Manage packages`. A package manager dialog pops up, that looks something like this:

![Package Manager](https://github.com/CompEvol/CCD/raw/master/doc/package_repos.png)

If the PhyloPotts package is listed, just click on it to select it, and hit the `Install/Upgrade` button.

If the PhyloPotts package is not listed, you may need to add a package repository by clicking the `Package repositories` button. A window pops up where you can click `Add URL` and add `https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml` in the entry. After clicking OK, the dialog should look something like this:

![Package Repositories](https://github.com/CompEvol/CCD/raw/master/doc/package_repos0.png)

Click OK and now PhyloPotts should be listed in the package manager (as in the first dialog above). Select and click Install/Upgrade to install.


## Installing by hand

Install [BEAST 2](http://beast2.org).

Download [PhyloPotts.v0.0.1.zip](https://github.com/rbouckaert/potts/releases/download/v0.0.1/PhyloPotts.v0.0.1.zip). 
Then, create a BREATH subdirectory:

```
for Windows in Users\<YourName>\BEAST\2.7\PhyloPotts
for Mac in /Users/<YourName>\/Library/Application Support/BEAST/2.7/PhyloPotts
for Linux /home/<YourName>/.beast/2.7/PhyloPotts
```

Here `<YourName>` is the username you use.
Unzip the file [PhyloPotts.v0.0.1.zip](https://github.com/rbouckaert/potts/releases/download/v0.0.1/PhyloPotts.v0.0.7.zip) in the newly created directory.


## Build from code

An alternative way to install is to build from the source code. 
First, get code for beast2, BeastFX and this repository. Then run

```bash
ant install
```

to install the package.


## Reference

Remco Bouckaert.
A first order approximation of a joint Potts and phylogenetic model.
Submitted to [BioRxiv](https://www.biorxiv.org/), 2026.

Presented at the <a href="https://biosig.lab.uq.edu.au/strphy26/">SMBE Australasian Protein Structural Phylogenetics Meeting</a>, 2026 in Brisbane
<a href="http://rbouckaert.github.io/PhyloPotts/">Slides</a>