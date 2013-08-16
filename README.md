pyCBT
=====

Python library for Coulomb Blockade Thermometer (CBT) data fitting. It uses master equation according to *Pekola et. al. PRL*
**73** *(1994)*.

###Requirements

* Python 2.7.x http://www.python.org/
* numpy http://www.numpy.org/
* scipy http://www.scipy.org/
* matplotlib http://matplotlib.org/
* pyx http://pyx.sourceforge.net/


Packages except pyx have installers. For pyx installation, run in shell 

```
python setup.py install
```

### Usage 
Edit and run example [script](src/example_fit_and_plot.py) 
```
python example_fit_and_plot.py
```


that fits a measured CBT curve.

### Cython
There exist a cython version of file ``CBT_lib.py`` named ``CBT_lib.pyx`` that can be compiled to C-code according instructions
from [cython webpage](http://cython.org/)



#### Disclaimer

The Software is provided "as is" without warranty of any kind, either express or implied, including without limitation any implied warranties of condition, uninterrupted use, merchantability, fitness for a particular purpose, or non-infringement

