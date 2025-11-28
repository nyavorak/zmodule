# zmodule
Photometric redshift estimation with COLIBRI filter bands (but work for many others, just update the transmissions)<br />
- Use of pyGRBz (https://pygrbz.readthedocs.io/en/latest/) and pyGRBaglow (https://pygrbaglow.readthedocs.io/en/latest/)
- Tested and works on Ubuntu, Fedora, MacOs (For MacOs you may need to install `pdm` with brew: `brew install pdm`)
# Installation
- Create a python3.12 environment for the module `conda create -n yourenv python=3.12`
- Clone and go inside the source folder
- Install the requirements with `pip install -r requirements.txt`
- Just run `python setup.py` to install pyGRBz and pyGRBaglow
- To run a notebook, install your conda ipykernel with `python -m ipykernel install --user --name yourenv --display-name "yourenv"`
# How does it work ?
- Open and modify the `photoz.py` file as you want then run `python photoz.py`
- Input name of the GRB or list of GRBs, i.e. correct format is GRBXXXXXX
- Is your input Multitargets, a light curve or an SED ? Some examples for an input template are in `data/lc` and `data/sed`
- Conversion of magnitudes to flux and automatic correction from Galactic extinction if necessary
- Extraction of the SED:
    * Model of the SED is either with a single power law (model = "SPL") or a broken power law (model = "BPL")
    * At a given time (method='fixed',time_SED = ...) or time at which the flux is maximum in the reddest band (method='ReddestBand')
- MCMC sampling with the extinction laws. Accepted laws are ['smc', 'lmc', 'mw', 'nodust', 'sne']
- Comparison of all the fits statistically in term of $\chi^2$ and $\Delta$ BIC (default threshold is 2)
- Some notebooks are also available
