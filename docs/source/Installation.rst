#################
III. Installation
#################

BASE-9 is written in C++ and designed to run on a variety of UNIX- and Linux-based operating systems. It is currently tested on

* Mac OS X 10.7 through 10.9
* Ubuntu 10.04 through 12.04
* RHEL 5 and 6
* Gentoo 13.0
* FreeBSD 9 and 10

To compile the code you will need gcc 4.7+ or clang 3.2+ (C/C++ language compilers), gsl (the gnu science library), cmake (a cross-platform build system), and Boost (a peer-reviewed, portable C++ library). To install these software packages, you may need help from your system administrator, though we provide some guidance here.

The best place to put all of this code is in the /usr/local/bin directory. If you don’t have that directory on your machine already, you can create it as follows:

::

	> sudo mkdir /usr/local
	> sudo mkdir /usr/local/bin

Note that the sudo command gives you super-user or root permission for that one command (after you enter your password at the prompt), assuming that your account has been allowed to invoke the command.

*******************************************************************
A. Installing gcc, gsl, and cmake on a Mac running OS X 10.7 – 10.9
*******************************************************************

Download the compiler. One way to do that is via downloading Xcode 4.6 or later from http://connect.apple.com. This requires that you have a developer account, but you can register for that for free. Also, it will give you 1+ GB of code and tools, most of which you will only need if you intend to develop for iPhones, Mac OS, etc. If you do follow this route, after installing Xcode, you will need to specifically install the command line tools with a window that will look similar to the one below. Click on the install button to the right of “command line tools” and it will appear as follows when done.

This will install clang and (for Xcode 4.6) GCC 4.2.

Download gsl from http://www.gnu.org/software/gsl/. Use the ftp site to obtain the source code, then

::

	mac> cd ~/Downloads
	mac> tar xzf gsl-1.15.tar.gz
	mac> cd gsl-1.15/
	mac>./configure
	mac> make
	mac> sudo make install

Download cmake from http://www.cmake.org/cmake/resources/software.html. Choose the .dmg version of the code (for the correct operating system) if you want to let Mac installation guide you through the process. We suggest placing the cmake build into ``/usr/local/bin`` by choosing that directory when prompted.

****************************************************
B. Installing gcc, gsl, and cmake on a Linux machine
****************************************************

The simplest way to install on Ubuntu is via the apt-get tool.

::

	linux> sudo apt-get install gcc
	linux> sudo apt-get install cmake cmake-curses-gui
	linux> sudo apt-get install libgsl0-dev libgsl0ldbl

A similar process will work with the yum tool on Fedora / RHEL* :

::

	linux> sudo yum install gcc gcc-c++ cmake git
	linux> sudo yum install gsl gsl-devel boost boost-devel

.. note:: RHEL 5 & 6 repositories have an old gcc version. The devtoolset package will install an alternate, up to date build environment at /opt/<distro>/devtoolset-2/

::

	linux> sudo yum install devtoolset-2-toolchain

The scl utility can create a shell referencing these alternate build tools where BASE-9 can be built:

::

	linux> scl enable  devtoolset-2 'bash'

*******************
C. Unpacking BASE-9
*******************

Create a directory where you wish to install and run the software, then download the newest code release from https://github.com/argiopetech/base/tags and the newest stellar evolution files from https://github.com/argiopetech/base-models/tags, and extract them to the appropriate directory, e.g.

::

	> tar xzf base-9.4.2.tar.gz
	> cd base-9.4.2/

Note that your computer may uncompress the .gz file for you on download, in which case the above command would instead be

::

	> tar xvf base-9.4.2.tar
	> cd base-9.4.2/

*******************
D. Installing Boost
*******************

BASE-9 has an included script to install Boost. The install location can be changed by modifying the CMAKE_INSTALL_PREFIX variable.

::

	linux> cd contrib
	linux> cmake –DCMAKE_INSTALL_PREFIX=”/usr/local” .
	linux> sudo make

Ubuntu users can save some time on this step by running

::

	linux> sudo apt-get install libboost-dev

********************
E. Installing BASE-9
********************

Once you have all of the above software in place, you are ready to install BASE-9. The following instructions should work identically for all platforms.

Change directories into the BASE-9 source directory and simply run build.sh

::

	> sudo ./build.sh


This will (if you have properly installed all libraries) build and install the BASE-9 executables and install them in the default location (generally ``/usr/local/bin``).

Alternatively, if you do not have the ability to run `sudo` on your machine, you may use

::

	> ./build_local.sh

to build and install the executables locally. The executables will be installed in the BUILD/bin
directory.