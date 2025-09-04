# ILWISPy
Integrated Land and Water Information System (ILWIS) is a remote sensing and GIS software.

ILWISPy is a Geo-Processing Tool which can be used as a site-package under Python, and brings the computational functionality of ILWIS to Python. The main objective of ILWISPy is to 'achieve more with less coding' and 'increase flexibility'. ILWISPy has a large library of commonly used GIS and RS operations which can be easily coded through single line statements. The core is based on C++ and is fast. This capability, together with other Python sitepackages, offers great flexibility in data retrieval, (pre-) processing and analysis, far exceeding those offered through Graphical User Interface (GUI) Geo-processing tools.

ILWISPy can be directly imported in Python, but can also be used within a Jupyter Notebook, offering the capability of markdown text with explanations and subsequent coding fields to execute certain operations, which for educational purposes offers nice opportunities.

## Installation
The extension is available at https://pypi.org/project/ilwis/, making it possible to perform the installation through `pip install`, as follows:

`python3 -m pip install ilwis`

Add `--upgrade` to this statement, to update a previously installed version.

## Known installation issues under Linux

### When using a version of python that is installed under a conda environment

#### 1. `libpython3.11.so.1.0`
When `import ilwis` within python gives following error:

`ImportError: libpython3.11.so.1.0: cannot open shared object file: No such file or directory`

(the exact message depends on your conda-activated version of python; this message is for python 3.11)

This happens because when conda activates an environment, it doesn't automatically add the python libraries to the library-search-path.

Solution: Exit python, then (by typing commands in the terminal) add the folder of the libraries to the `LD_LIBRARY_PATH` search-path. Replace the conda environment name with the one that you use (here it is `p311`):

`export LD_LIBRARY_PATH=~/miniconda3/envs/p311/lib`

Start python in the activated environment, and try the import again.

#### 2. `libQt5Gui.so.5`
When `import ilwis` within python gives the following error:

`ImportError: libQt5Gui.so.5: cannot open shared object file: No such file or directory`

This error message comes when the `Qt5` libraries are not in the library-path or are not installed at all on the system.
ilwispy is not shipped with `Qt5`, but does depend on it. On systems that aren't shipped with `Qt5` (e.g. because they were upgraded to use `Qt6`), the safest solution is to install `Qt5` as a non-admin user, with `python-pip`.

Solution: Exit python. In the activated conda environment perform the following:

`python3 -m pip install pyqt5`

This installs the `pyqt5` extension (which is not used by ilwispy), but also downloads the `Qt5` binary libraries that are needed.

Then add the folder of the libraries to the `LD_LIBRARY_PATH` search-path. Replace the python version (here it is python 3.11) with your version, and the conda environment name (here it is p311) with the one that you use:

`export LD_LIBRARY_PATH=~/miniconda3/envs/p311/lib:~/miniconda3/envs/p311/lib/python3.11/site-packages/PyQt5/Qt5/lib`

If you are unsure if that is the correct folder, you can search for it:

`find ~ -name libQt5Gui.so.5`

Start python in the activated environment, and try `import ilwis` again.

#### 3. `ILWIS couldn't be initialized`

##### 3.1 `libgdal.so fails to load due to libldap.so`

When `import ilwis` within python gives the following error:

`ImportError: ILWIS couldn't be initialized!`

`Could not load library in module:gdalconnector, name: libgdal.so,error :Cannot load library libgdal.so: (/lib64/libldap.so.2: undefined symbol: EVP_md2, version OPENSSL_3.0.0)`

This happens because of a version-discrepancy between a system-library and a library installed in the conda-environment, preventing `libgdal.so` to load.

Solution: Exit python, and execute the following to ensure availability of the correct library version:

`conda install python-ldap --force-reinstall`

Start python in the activated environment, and try `import ilwis` again.

##### 3.2 `libgdal.so fails to load due to libstdc++.so`

When `import ilwis` within python gives the following error:

`ImportError: ILWIS couldn't be initialized!`

`Could not load library in module:gdalconnector, name: libgdal.so,error :Cannot load library libgdal.so: (/home/<username>/miniconda3/envs/p311/lib/libstdc++.so.6: version 'GLIBCXX_3.4.30' not found (required by /lib/libgdal.so))`

This happens because of a version-discrepancy between a system-library and a library installed in the conda-environment, preventing `libgdal.so` to load.

Solution: Exit python, and execute the following to ensure availability of the correct library version:

`conda install -c conda-forge gcc`

Start python in the activated environment, and try `import ilwis` again.

### When using the system-preinstalled python (`python3`)

#### 1. `libQt5Gui.so.5`
When `import ilwis` within python3 gives the following error:

`ImportError: libQt5Gui.so.5: cannot open shared object file: No such file or directory`

This error message comes when the `Qt5` libraries are not in the library-path or are not installed at all on the system.
ilwispy is not shipped with `Qt5`, but does depend on it. On systems that aren't shipped with `Qt5` (e.g. because they were upgraded to use `Qt6`), the safest solution is to install `Qt5` as a non-admin user, with `python-pip`.

Solution: Exit python, and perform the following (by typing commands in the terminal):

`python3 -m pip install pyqt5`

This installs the `pyqt5` extension (which is not used by ilwispy), but also downloads the `Qt5` binary libraries that are needed.

Then add the folder of the libraries to the `LD_LIBRARY_PATH` search-path. Replace the python version (here it is python 3.9) with your version:

`export LD_LIBRARY_PATH=~/.local/lib/python3.9/site-packages/PyQt5/Qt5/lib`

If you are unsure if that is the correct folder, you can search for it:

`find ~ -name libQt5Gui.so.5`

Start python3, and try `import ilwis` again.

### Making `LD_LIBRARY_PATH` persistent
If a change in `LD_LIBRARY_PATH` is required to make the ilwis extension work correctly, there are several ways to make this persistent.

An safe way is to make it part of the python-script. For example the `libQt5Gui.so.5` statement above would become:

`import os`

`os.environ['LD_LIBRARY_PATH'] = os.path.expanduser('~') + '/.local/lib/python3.9/site-packages/PyQt5/Qt5/lib'`

`import ilwis`

## WebSite

https://www.itc.nl/about-itc/scientific-departments/water-resources/software-tools-models/ilwis3-and-toolbox-plugins/ilwispy-getting-started/

## What's new

https://filetransfer.itc.nl/pub/52n/ilwis_py/whatsnew.txt

