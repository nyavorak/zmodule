# zmodule
Photometric redshift estimation with COLIBRI filter bands (but work for many others)<br />
- Use of pyGRBz (https://pygrbz.readthedocs.io/en/latest/) and pyGRBaglow (https://pygrbaglow.readthedocs.io/en/latest/)
- Tested and works on Linux, bugs on Mac M2
# Installation
- Just run `python setup.py` from the directory
# How does it work ?
- Just run `python photoz.py` from the directory
- Input name of the GRB or list of GRBs, i.e. correct format is GRBXXXXXX
- Is your input Multitargets, light curve or SED ?
- Calculation of the flux and correction for Galactic extinction if necessary
- Extraction of the SED:
    * Model of the SED is either with a single power law (model = "SPL") or a broken power law (model = "BPL")
    * At a given time (method='fixed',time_SED = ...) or time at which the flux is maximum in the reddest band (method='ReddestBand')
- MCMC sampling with the extinction laws: 'smc', 'lmc', 'mw', 'nodust'
- Comparison of all the fits statistically in term of $\chi^2$ and BIC (default threshold is 2)
* Also working Jupyter notebooks in `pyGRBz/pyGRBz/notebooks`
