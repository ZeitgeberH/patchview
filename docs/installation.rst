.. highlight:: shell

=====================
Installation & Usage
=====================
There are three ways to use Patchview

Stand-along App
-----------------------

If your operating system is Windows 10, PatchView has a prepackaged app using `Pyinstaller`_.  
You can download the zip file from the release page. After unzipped
it into a folder, you can directly double click the excutable file to start the app. No
python installation is needed. 

To remove the app, just delete the whole folder.

Note: since Pyinstaller is cross-platform, it is possible to package Patchview for other platform, Mac or Linux.
We have not try these options yet. Contributions are welcome!

Pip install
--------------

To install PatchView via `pip`_, We recommand creating a virtual python enviroment first.
If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

If you use conda, you can do this by:

.. code-block:: console
    
    $  conda create -n patchview python=3.8

Note: PatchView require Python>=3.8

After activating your virtual enviroment, run this command in your terminal:

.. code-block:: console

    $ pip install patchview

This is the preferred method to install PatchView, as it installs script that can launch
the GUI program at a command line by type

.. code-block:: console

    (patchview) $ patchview

or in a python enviroment:

.. code-block:: python

    >>> import patchview

To launch the GUI, type:

.. code-block:: python

    >>> patchview.pvGUI()

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

Note: PatchView use PyQt5 as its front-end. Pip installation may have compatible issue which may cause the following errors when launch
the GUI::

    qt.qpa.plugin: Could not find the Qt platform plugin "windows" in "". This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

If this happens, please uninstall PyQt5, PyQt5-sip, pyqtwebengine by:

.. code-block:: console

    $ pip uninstall pyqt5 pyqt5-qt5 pyqtwebengine-qt5 pyqt5-sip

then reinstall them via:

.. code-block:: console

    $ pip install pyqt5 pyqtwebengine pyqt5-sip

If this does not resolve the issue, please open a PatchView's github issue.

From sources
------------

The sources for PatchView can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/zeitgeberH/patchview

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/zeitgeberH/patchview/tarball/master

Once you have a copy of the source, use Conda to create an virtual enviroment:

.. code-block:: console

    $ conda env create -f environment.yml

then activate the enviroment and run:

.. code-block:: console

    $ conda activate patchviewPy3
    
    $ python setup.py install
    
.. _Pyinstaller: https://pyinstaller.org/en/stable/   
.. _Github repo: https://github.com/zeitgeberH/patchview
.. _tarball: https://github.com/zeitgeberH/patchview/tarball/master
