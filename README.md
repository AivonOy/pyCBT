pyCBT
=====

Python library for Coulomb Blockade Thermometer (CBT) data fitting. It uses master equation for conductance calculations according to *Pekola et. al. PRL*
**73** *(1994)*.

See 

* [Tutorial](http://nbviewer.ipython.org/urls/raw.github.com/AivonOy/pyCBT/master/pyNotebooks/pyCBT%2520tutorial.ipynb)

* [Multifit quickstart](http://nbviewer.ipython.org/urls/raw.github.com/AivonOy/pyCBT/master/pyNotebooks/pyCBT%2520multifit%2520quickstart.ipynb)


* [Multifit with explanations](http://nbviewer.ipython.org/urls/raw.github.com/AivonOy/pyCBT/master/pyNotebooks/pyCBT%2520MultiFit%2520tutorial.ipynb)

* [Capability calculations](http://nbviewer.ipython.org/github/AivonOy/pyCBT/blob/master/pyNotebooks/On%20CBT%20measurement.ipynb?create=1)



###Requirements

* Python 2.7.x http://www.python.org/
* numpy http://www.numpy.org/
* scipy http://www.scipy.org/
* matplotlib http://matplotlib.org/
* pyx http://pyx.sourceforge.net/


Packages except pyx have installers. For pyx installation, run pyx installation file ``setup.py`` in shell 

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

In Windows compiling using mingw at least works. In OSX Xcode, command line tools have been tested to work.



#### Disclaimer

The Software is provided "as is" without warranty of any kind, either express or implied, including without limitation any implied warranties of condition, uninterrupted use, merchantability, fitness for a particular purpose, or non-infringement

