# zmodule
photoz estimation <br />
- Install pyGRBz (https://pygrbz.readthedocs.io/en/latest/) and pyGRBaglow (https://pygrbaglow.readthedocs.io/en/latest/)
- Tested and works on Linux, bugs on Mac M2
# Installation
- Just run `python setup.py` from the directory
# How does it work ?
- Input name of the GRB or list of GRBs, i.e. correct format is GRBXXXXXX
- Multitargets, light curve or SED given
- Calculating flux and correcting for Galactic extinction if necessary
- Extract the SED:
    * Model of the SED either with a single power law (model = "SPL") or a broken power law (model = "BPL")
    * At a given time (method='fixed',time_SED = ...) or time at which the flux is maximum in the reddest band (method='ReddestBand')
- Put your priors
- MCMC sampling with the extinction laws: 'smc', 'lmc', 'mw', 'nodust'
- Compare all fits statistically in term of $\chi$ and BIC
