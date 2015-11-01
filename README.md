HOLA
====

HOLA ("Human-like Orthogonal Layout Algorithm") is a new graph layout algorithm
developed by [Steve Kieffer][sk], [Tim Dwyer][td], [Kim Marriott][km], and
[Michael Wybrow][mw], in the
[MArVL: Monash Adaptive Visualisation Lab][marvl] at
[Monash University][monash] in Melbourne, Australia.
The Python library was written by Steve Kieffer.

In order to use this library you need
[SWIG][swig], as well as
[Adapaca][adapaca], which is an experimental fork of the
[Adaptagrams][adaptagrams] layout library.

Currently the HOLA library provides read/write support only for the GML
file format. It comes bundled with [SPARK][spark], which it uses to parse GML.


How to cite?
------------

If you refer to HOLA in your own research or publications, please cite 
the following reference:

 *  Steve Kieffer, Tim Dwyer, Kim Marriott, and Michael Wybrow.  
    HOLA: Human-like Orthogonal Network Layout.
    In Visualization and Computer Graphics, IEEE Transactions on,
    Volume 22, Issue 1, pages 349 - 358. IEEE, 2016. [DOI][doi]


Installing
----------

The following steps can be taken on Linux or Mac; Windows users should be able
to do something similar.

 1. Clone the HOLA repository into a directory of your choice,
    `HOLADIR`.

 2. Clone [Adapaca][adapaca] into a directory of your choice,
    `ADAPDIR`.

 3. Follow the instructions at the [Adaptagrams repository page][adaptagrams]
    to build Adaptagrams.

 4. As mentioned at that page, Python bindings for Adaptagrams can be
    generated using SWIG. After building Adaptagrams you should be able
    to create the SWIG bindings with the following commands:

        $ cd ADAPDIR/cola
        $ make -f Makefile-swig-python

 5. Finally, you need to create a link from the HOLA library to Adaptagrams:

        $ cd ADAPDIR/cola
        $ touch __init__.py
        $ cd HOLADIR/main
        $ ln -s ADAPDIR/cola adaptagrams


Using HOLA
----------

For the sake of running a quick test with basic configuration, HOLA can be used
from the command line as follows:

    $ HOLADIR/hola_basic < input.gml > output.gml

For more serious use of the library you may add `HOLADIR` to your `PYTHONPATH`.
You can then import the main algorithm with

    from hola.hola import hola

and the configuration class with

    from hola.hola import HolaConfig

For further help please see the docstrings in the `hola` and other
modules of the library.


[adaptagrams]: https://github.com/mjwybrow/adaptagrams/
[adapaca]: https://github.com/skieffer/adapaca
[marvl]: http://marvl.infotech.monash.edu.au/
[monash]: http://www.infotech.monash.edu.au/about/schools/caulfield/
[td]: http://marvl.infotech.monash.edu/~dwyer/
[km]: http://www.csse.monash.edu.au/~marriott/
[mw]: http://www.csse.monash.edu.au/~mwybrow/
[sk]: http://skieffer.info
[spark]: http://pages.cpsc.ucalgary.ca/~aycock/spark/
[swig]: http://www.swig.org/
[doi]:http://dx.doi.org/10.1109/TVCG.2015.2467451

