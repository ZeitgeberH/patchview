.. highlight:: shell

==============
Get started
==============
There are three ways to use Patchview

Stand-along App
-----------------------

If your operating system is Windows 10, PatchView has a prepackaged app using `Pyinstaller`_.  
You can download the zip file from its Github `releases`_ page. After unzipped
it into a folder, you can directly double click the executable file, **PatchView.exe**, to start the app. No
python installation is needed. 

.. _releases: https://github.com/ZeitgeberH/patchview/releases

To remove the app, just delete the whole folder.

**Note**

* It takes about 15 ~ 30 seconds when you run the programme the first time, as compressed libraries needs to be unzipped. 
  It should be much faster for following runs. 

* since Pyinstaller is cross-platform, it is possible to package PatchView for other platform, Mac or Linux. 
  We have not tried these options yet. Contributions are welcome!

Pip install
--------------

To install PatchView via `pip`_, We recommend creating a virtual python environment first.
If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.



If you use conda, you can do this by:

.. code-block:: console
    
    $  conda create -n patchview python=3.8

Note: PatchView require Python>=3.8

After activating your virtual environment, run this command in your terminal:

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

then activate the environment and run:

.. code-block:: console

    $ conda activate patchviewPy3
    
    $ python setup.py install
    
.. _Pyinstaller: https://pyinstaller.org/en/stable/   
.. _Github repo: https://github.com/zeitgeberH/patchview
.. _tarball: https://github.com/zeitgeberH/patchview/tarball/master


Configuration
----------------

PatchView use a Yaml file for its basic configuration.

* If you use the app version, it is located in the **patchview\\Data\\patchview.yaml** of the app folder.
* If you install it via Pip or from source, it is located in **PATH-OF-YOUR-ENVIROMENT\\Lib\\site-packages\\Patchview-xx-py3.x.egg\\patchview\\Data\\patchview.yaml**.

Open **patchview.yaml** with any text editor.  

* **RootDir** is the root node for PatchView to search files. Leave it empty (**''**) if you want to include all drives in computer. 
* **Protocols** are labels you used for series during recording. There are four default categories that PatchView uses to sort recorded series. Add any labels in the corresponding category if it is not already in the list.

After saving your changes, close the app and restart it.