[![GitHub issues](https://img.shields.io/github/issues/acoiman/wildfire_modeling)](https://github.com/acoiman/wildfire_modeling/issues)  [![GitHub forks](https://img.shields.io/github/forks/acoiman/wildfire_modeling)](https://github.com/acoiman/wildfire_modeling/network)  [![GitHub stars](https://img.shields.io/github/stars/acoiman/wildfire_modeling)](https://github.com/acoiman/wildfire_modeling/stargazers)
<a href="https://buymeacoffee.com/acoiman" title="Donate to this project using Buy Me A Coffee"><img src="https://img.shields.io/badge/%20buy%20me%20a%20cofee-donate-green" alt="Buy Me A Coffee donate button"></a>



# ![fire](https://github.githubassets.com/images/icons/emoji/unicode/1f525.png) ​Wildfire Modeling in Yosemite National Park

In this tutorial, we will show you how to model wildfire events using the [r.ros](https://grass.osgeo.org/grass78/manuals/r.ros.html) and [r.spread](https://grass.osgeo.org/grass78/manuals/r.spread.html) modules of  [GRASS GIS](https://grass.osgeo.org/). We will perform fire simulations in three areas in the [Yosemite National Park](https://www.nps.gov/yose/index.htm) , California USA  during the summer season in 2020 (`from 2020-06-20 through 2020-09-22`). 

## How to start

We assume you have a running installation of  [GRASS GIS](https://grass.osgeo.org/)  and  [Anaconda](https://www.anaconda.com/) including the latest version of [Jupyter Notebook](https://jupyter.org/) on Ubuntu 20.04. You also need an active [Google Earth Engine (GEE)](https://earthengine.google.com/) account. 

Use the following commands to create a conda environment, activate it, and install packages required for this tutorial.

```
conda create -n grass_env python=3
conda activate grass_env
conda install wxpython
conda install numpy
```

You also need to install **Geopandas, GEE Python API, eeconvert, geemap** packages. In the list below I provide you some links that explain how to install them.

- Geopandas: https://geopandas.org/getting_started/install.html
- GEE Python API: https://developers.google.com/earth-engine/guides/python_install
- eeconvert: https://github.com/gee-community/eeconvert/blob/master/README.md
- geemap: https://geemap.org/installation/

Once activated the environment and installed all packages, clone this repository and move inside the directory of  `wildfire_modeling` tutorial and enter `jupyter notebook` command.

Enjoy the tutorial!!



## Author

- Abraham Coiman



## License

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.

<a href="https://www.paypal.com/paypalme/acoiman?locale.x=en_XC"
  <img height="32" src=https://raw.githubusercontent.com/everdrone/coolbadge/master/badges/Paypal/Coffee/Dark/Short.png />
  alt="Buy Me A Coffee PayPal donation button" </a>
