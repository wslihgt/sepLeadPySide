==============================================================
A Graphical User Interface For Musical Audio Source Separation
==============================================================

 :Authors: Jean-Louis Durrieu 
 :Organization: LTS5, EPFL
 :Date: 2011-2012
 :Contact: http://www.durrieu.ch/research

PRELIMINARY NOTE:
=================
Check also README_v2.pdf for updated documentation.

REQUIREMENTS:
=============
 * Qt            (http://qt.nokia.com, included in PyQt4 binaries, see below)
 * Python        (use EPD: http://www.enthought.com/products/epd.php)
 * Numpy         (included in EPD, http://numpy.scipy.org)
 * Scipy         (in EPD, http://www.scipy.org)
 * Matplotlib    (in EPD, http://matplotlib.sourceforge.net)
 * PySide        (in EPD, http://www.pyside.org)

USAGE:
======

::

    python separateLeadGUI2.py

USING THE INTERFACE:
====================
Typical run:

 1) Open a file (through Menu or "Open File" button or drag and drop
    in the corresponding field)

    If desired, please precise a directory to which the program will write
    the resulting files, in the "Output Directory Suffix" line. 
    Note that for now, the program will simply concatenate the given suffix
    with the full-path of the audio file (i.e. without the filename of that 
    file). 

 2) Select desired parameters for the algorithm:

    * window length: the window size, for the analysis frames
      0.04644s is a good choice for singing voice, @44100Hz.
      Consider longer windows for bass extraction.
    * min F0: the lowest fundamental frequency F0 candidate
    * max F0: the highest fundamental frequency F0 candidate

 3) Launch first decomposition, by pushing "Load File" button.

 4) Wait for the computations to be computed (on Linux, one might
    see the evolution of the display, not available for MacOSX or Win yet)

 5) Once the program is responsive again, the main window shows the matrix
    of amplitudes, for each F0 candidate: it is a time-frequency representation
    of the signal, with x-axis as time axis, and y-axis as F0 frequency axis.
    For now, the time is indicated as the frame number (current step size
    between frames is 12.5% of window length), while the F0 scale is shown 
    on a log2 scale, proportional to the Western Muscial scale. The different
    A notes (A2: 110Hz, A3: 220Hz, A4: 440Hz,...) are also displayed, for
    ease of use. For more details on the representation, see [Durrieu2011]_.
    
    The user should then select the time-frequency zones corresponding to
    the instrument she desires to separate. The provided matrix of amplitude
    should help deciding and identifying these regions. When the user clicks,
    and moves the mouse while left button is pressed, the selected region
    corresponds to the line described by the mouse, and all the points
    between that line plus and minus half a semitone. These regions are 
    shown as red delimited contours.
    
    The user can zoom in the picture and move around thanks to the top
    toolbar. Note however that each time after using these tools, she
    needs to deactivate them again (by clicking again on the corresponding
    button) in order to be able to proceed with the source selection.
    She can also modify the display normalization, between 3 modes: 
    no-normalization (displays the actual values in dB scale), normalizing
    each frame by the maximum value, or by the sum of all values. The minimum
    and maximum values to set the color scale can also be modified.
    
    Note also that the possibility of going back in the annotation is under
    development. For now, the user can already choose the button "Delete",
    and going over the previously selected zones will "deselect" them.
    
    If the application is separateLeadGUIControlsPlay.py, then the user can 
    also have an audio feedback corresponding to the displayed time-range,
    by clicking on the image with the mouse right-button.

 6) Once satisfied with the selected regions, the user should click the
    "Separate" button, which will launch the estimation of the separated
    sources, given the user input. 
    
    Alternatively, the user could also push the "Separate (Melody)" button,
    which discards the user selected areas and, instead, automatically detects
    the most prominent melody line, as published in our previous works
    [Durrieu2010]_.


FRIENDLY FEEDBACK:
==================
The author would be greatful if the testers could provide some feedback.
More specifically, some ergonomy issues encountered would be welcome, for
instance:

    * audio feedback required (for now, we advice the use of Audacity).
    * better aspect for selected (and even more de-selected) areas.
    * adequation of the frequency scale for musicians/non-musicians?
    * readability of the representation?

KNOWN BUGS:
===========
31.8.2011:
----------

    * The Phonon module has been added, but the playback on Linux is 
      not smooth. 
    * The "right click to play" is somehow conflicting with the pan/zoom 
      tool (right click + move allows the zooming tool on the image, and 
      will start or stop the playback)

30.8.2011:
----------

    * The interface may be a bit slow, sometimes, because, under the hood, 
      each time the user changes the colormap range, or selects new regions,
      quite some computation is going on (matplotlib ``contours`` and ``set_clim``
      stuff).

REFERENCES:
===========

.. [Durrieu2010] J.-L. Durrieu, G. Richard, B. David and C. Févotte, Source/Filter Model for Main Melody Extraction From Polyphonic Audio Signals, IEEE Transactions on Audio, Speech and Language Processing, special issue on Signal Models and Representations of Musical and Environmental Sounds, March 2010, vol. 18 (3), pp. 564 -- 575.

.. [Durrieu2011] J.-L. Durrieu, G. Richard and B. David, A Musically Motivated Representation For Pitch Estimation And Musical Source Separation, IEEE Journal of Selected Topics on Signal Processing, October 2011, Vol. 5 (6), pp. 1180 - 1191.

.. [Durrieu2012] J.-L. Durrieu and J.-Ph. Thiran, Musical Audio Source Separation Based on User-Selected F0 Track, International Conference on Latent Variable Analysis and Signal Separation (LVA/ICA), March 12-15, 2012, Tel-Aviv, Israel. http://www.durrieu.ch/research/lvaica2012.html