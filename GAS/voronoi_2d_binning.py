"""
#####################################################################

Copyright (C) 2001-2014, Michele Cappellari
E-mail: michele.cappellari_at_physics.ox.ac.uk

Updated versions of the software are available from my web page
http://purl.org/cappellari/software

If you have found this software useful for your
research, we would appreciate an acknowledgment to use of
`the Voronoi binning method by Cappellari & Copin (2003)'.

This software is provided as is without any warranty whatsoever.
Permission to use, for non-commercial purposes is granted.
Permission to modify for personal or internal use is granted,
provided this copyright and disclaimer are included unchanged
at the beginning of the file. All other rights are reserved.

Sections of code authored by Erik Rosolowsky are released to the public
domain.
  

#####################################################################

NAME:
    VORONOI_2D_BINNING

AUTHOR:
      Michele Cappellari, University of Oxford
      michele.cappellari_at_physics.ox.ac.uk

PURPOSE:
      Perform adaptive spatial binning of Integral-Field Spectroscopic
      (IFS) data to reach a chosen constant signal-to-noise ratio per bin.
      This method is required for the proper analysis of IFS
      observations, but can also be used for standard photometric
      imagery or any other two-dimensional data.
      This program precisely implements the algorithm described in
      section 5.1 of the reference below.

EXPLANATION:
      Further information on VORONOI_2D_BINNING algorithm can be found in
      Cappellari M., Copin Y., 2003, MNRAS, 342, 345
      http://adsabs.harvard.edu/abs/2003MNRAS.342..345C

CALLING SEQUENCE:

    binNum, xBin, yBin, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
                x, y, signal, noise, targetSN, plot=True, quiet=False, 
                wvt=False, cvt=True, pixelsize=None)

    The function _sn_func() below returns the S/N of a bin and it can be
    changed by the user if needed.

INPUTS:
           X: Vector containing the X coordinate of the pixels to bin.
              Arbitrary units can be used (e.g. arcsec or pixels).
              In what follows the term "pixel" refers to a given
              spatial element of the dataset (sometimes called "spaxel" in
              the IFS community): it can be an actual pixel of a CCD
              image, or a spectrum position along the slit of a long-slit
              spectrograph or in the field of view of an IFS
              (e.g. a lenslet or a fiber).
              It is assumed here that pixels are arranged in a regular
              grid, so that the pixel size is a well defined quantity.
              The pixel grid however can contain holes (some pixels can be
              excluded from the binning) and can have an irregular boundary.
              See the above reference for an example and details.
           Y: Vector (same size as X) containing the Y coordinate
              of the pixels to bin.
      SIGNAL: Vector (same size as X) containing the signal
              associated with each pixel, having coordinates (X,Y).
              If the `pixels' are actually the apertures of an
              integral-field spectrograph, then the signal can be
              defined as the average flux in the spectral range under
              study, for each aperture.
              If pixels are the actual pixels of the CCD in a galaxy
              image, the signal will be simply the counts in each pixel.
       NOISE: Vector (same size as X) containing the corresponding
              noise (1 sigma error) associated with each pixel.
    TARGETSN: The desired signal-to-noise ratio in the final
              2D-binned data. E.g. a S/N~50 per pixel may be a
              reasonable value to extract stellar kinematics
              information from galaxy spectra.

KEYWORDS:
         CVT: Set this keyword to skip the Centroidal Voronoi Tessellation
              (CVT) step (vii) of the algorithm in Section 5.1 of
              Cappellari & Copin (2003).
              This may be useful if the noise is strongly non Poissonian,
              the pixels are not optimally weighted, and the CVT step
              appears to introduces significant gradients in the S/N.
              A similar alternative consists of using the /WVT keyword below.
        PLOT: Set this keyword to produce a plot of the two-dimensional
              bins and of the corresponding S/N at the end of the
              computation.
     PIXSIZE: Optional pixel scale of the input data.
              This can be the size of a pixel of an image or the size
              of a spaxel or lenslet in an integral-field spectrograph.
            - The value is computed automatically by the program, but
              this can take a long times when (X, Y) have many elements.
              In those cases the PIXSIZE keyword should be given.
       QUIET: by default the program shows the progress while accreting
              pixels and then while iterating the CVT. Set this keyword
              to avoid printing progress results.
         WVT: When this keyword is set, the routine bin2d_cvt_equal_mass is
              modified as proposed by Diehl & Statler (2006, MNRAS, 368, 497).
              In this case the final step of the algorithm, after the bin-accretion
              stage, is not a modified Centroidal Voronoi Tessellation, but it uses
              a Weighted Voronoi Tessellation.
              This may be useful if the noise is strongly non Poissonian,
              the pixels are not optimally weighted, and the CVT step
              appears to introduces significant gradients in the S/N.
              A similar alternative consists of using the /NO_CVT keyword above.
              If you use the /WVT keyword you should also include a reference to
              `the WVT modification proposed by Diehl & Statler (2006).'

OUTPUTS:
   BINNUMBER: Vector (same size as X) containing the bin number assigned
              to each input pixel. The index goes from zero to Nbins-1.
              This vector alone is enough to make *any* subsequent
              computation on the binned data. Everything else is optional!
        XBIN: Vector (size Nbins) of the X coordinates of the bin generators.
              These generators uniquely define the Voronoi tessellation.
        YBIN: Vector (size Nbins) of Y coordinates of the bin generators.
        XBAR: Vector (size Nbins) of X coordinates of the bins luminosity
              weighted centroids. Useful for plotting interpolated data.
        YBAR: Vector (size Nbins) of Y coordinates of the bins luminosity
              weighted centroids.
          SN: Vector (size Nbins) with the final SN of each bin.
     NPIXELS: Vector (size Nbins) with the number of pixels of each bin.
       SCALE: Vector (size Nbins) with the scale length of the Weighted
              Voronoi Tessellation, when the /WVT keyword is set.
              In that case SCALE is *needed* together with the coordinates
              XBIN and YBIN of the generators, to compute the tessellation
              (but one can also simply use the BINNUMBER vector).

PROCEDURES USED:
      The following procedures are contained in the main VORONOI_2D_BINNING program.
          SN_FUNC           -- Example routine to calculate the S/N of a bin.
          WEIGHTED_CENTROID -- computes weighted centroid of one bin
          BIN_ROUNDNESS     -- equation (5) of Cappellari & Copin (2003)
          BIN_ACCRETION     -- steps (i)-(v) in section 5.1
          REASSIGN_BAD_BINS -- steps (vi)-(vii) in section 5.1
          CVT_EQUAL_MASS    -- the modified Lloyd algorithm in section 4.1
          COMPUTE_USEFUL_BIN_QUANTITIES -- self explanatory
          DISPLAY_PIXELS    -- plotting of colored pixels

MODIFICATION HISTORY:
      V1.0.0: First implementation. Michele Cappellari, Leiden, June 2001
      V2.0.0: Major revisions. Stable version. MC, Leiden, 11 September 2001
      V2.1.0: First released version. Written documentation.
          MC, Vicenza, 13 February 2003
      V2.2.0: Added computation of useful bin quantities in output. Deleted some
          safety checks for zero size bins in CVT. Minor polishing of the code.
          MC, Leiden, 11 March 2003
      V2.3.0: Unified the three tests to stop the accretion of one bin.
          This can improve some bins at the border. MC, Leiden, 9 April 2003
      V2.3.1: Do *not* assume the first bin is made of one single pixel.
          Added computation of S/N scatter and plotting of 1-pixel bins.
          MC, Leiden, 13 April 2003
      V2.4.0: Addedd basic error checking of input S/N. Reintroduced the
          treatment for zero-size bins in CVT, which was deleted in V2.2.
          Thanks to Robert Sharp and Kambiz Fathi for reporting problems.
          MC, Leiden, 10 December 2003.
      V2.4.1: Added /QUIET keyword and verbose output during the computation.
          After suggestion by Richard McDermid. MC, Leiden, 14 December 2003
      V2.4.2: Use LONARR instead of INTARR to define the CLASS vector,
          to be able to deal with big images. Thanks to Tom Statler.
          MC, Leiden, 4 August 2004
      V2.4.3: Corrected bug introduced in version 2.3.1. It went undetected
          for a long time because it could only happen in special conditions.
          Now we recompute the index of the good bins after computing all
          centroids of the reassigned bins in reassign_bad_bins. Many thanks
          to Simona Ghizzardi for her clear analysis of the problem and
          the solution. MC, Leiden, 29 November 2004
      V2.4.4: Prevent division by zero for pixels with signal=0
          and noise=sqrt(signal)=0, as can happen from X-ray data.
          MC, Leiden, 30 November 2004
      V2.4.5: Added BIN2D prefix to internal routines to avoid possible
          naming conflicts. MC, Leiden, 3 December 2004
      V2.4.6: Added /NO_CVT keyword to optionally skip the CVT step of
          the algorithm. MC, Leiden, 27 August 2005
      V2.4.7: Verify that SIGNAL and NOISE are non negative vectors.
          MC, Leiden, 27 September 2005
      V2.4.8: Use geometric centroid of a bin during the bin-accretion stage,
          to allow the routine to deal with negative signal (e.g. in
          background-subtracted X-ray images). Thanks to Steven Diehl for
          pointing out the usefulness of dealing with negative signal.
          MC, Leiden, 23 December 2005
      V2.5.0: Added two new lines of code and the corresponding /WVT keyword
          to implement the nice modification to the algorithm proposed by
          Diehl & Statler (2006). MC, Leiden, 9 March 2006
      V2.5.1: Updated documentation. MC, Oxford, 3 November 2006
      V2.5.2: Print number of unbinned pixels. MC, Oxford, 28 March 2007
      V2.5.3: Fixed program stop, introduced in V2.5.0, with /NO_CVT keyword.
          MC, Oxford, 3 December 2007
      V2.5.4: Improved color shuffling for final plot.
          MC, Oxford, 30 November 2009
      V2.5.5: Added PIXSIZE keyword. MC, Oxford, 28 April 2010
      V2.5.6: Use IDL intrinsic function DISTANCE_MEASURE for
          automatic pixelSize, when PIXSIZE keyword is not given.
          MC, Oxford, 11 November 2011
      V2.5.7: Included safety termination criterion of Lloyd algorithm
          to prevent loops using /WVT. MC, Oxford, 24 March 2012
      V2.5.8: Update Voronoi tessellation at the exit of bin2d_cvt_equal_mass.
          This is only done when using /WVT, as DIFF may not be zero at the
          last iteration. MC, La Palma, 15 May 2012
      V2.6.0: Included new SN_FUNCTION to illustrate the fact that the user can
          define his own function to estimate the S/N of a bin if needed.
          MC, London, 19 March 2014
      V3.0.0: Translated from IDL into Python and tested against the original.
          MC, London, 19 March 2014
      V3.0.1: Support both Python 2.6/2.7 and Python 3. MC, Oxford, 25 May 2014
      V3.0.2: Avoid potential runtime warning while plotting.
          MC, Oxford, 2 October 2014
      V3.1: Allow generic signal-to-noise functions to be passed in from main
      level.  Also changed API for _sn_func.
          EWR, Bleakest Alberta, Jan 2016.
"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy import ndimage

#----------------------------------------------------------------------------

def _sn_func(signal, noise):
    """
    Generic function to calculate the S/N of a bin with spaxels.
    The Voronoi binning algorithm does not require this function to have a
    specific form and this generic one can be changed by the user if needed.

    The S/N returned by this function does not need to be an analytic function
    of S and N. There is no need for this function to return the actual S/N.
    Instead this function could return any quantity the user needs to equalize.

    For example _sn_func could be a procedure which uses ppxf to measure the
    velocity dispersion from the coadded spectrum of spaxels and
    returns the relative error in the dispersion.
    Of course an analytic approximation of S/N speeds up the calculation.

    """
    return  np.sum(signal)/np.sqrt(np.sum(noise**2))

#----------------------------------------------------------------------

def _weighted_centroid(x, y, density):
    """
    Computes weighted centroid of one bin.
    Equation (4) of Cappellari & Copin (2003)

    """
    mass = np.sum(density)
    xBar = np.sum(x*density)/mass
    yBar = np.sum(y*density)/mass

    return xBar, yBar

#----------------------------------------------------------------------

def _roundness(x, y, pixelSize):
    """
    Implements equation (5) of Cappellari & Copin (2003)

    """
    n = x.size
    equivalentRadius = np.sqrt(n/np.pi)*pixelSize
    xBar, yBar = np.mean(x), np.mean(y)  # Geometric centroid here!
    maxDistance = np.sqrt(np.max((x - xBar)**2 + (y - yBar)**2))
    roundness = maxDistance/equivalentRadius - 1.

    return roundness

#----------------------------------------------------------------------

def _accretion(x, y, signal, noise, targetSN, pixelSize, quiet, aggregator = _sn_func):
    """
    Implements steps (i)-(v) in section 5.1 of Cappellari & Copin (2003)

    """
    n = x.size
    classe = np.zeros(n, dtype=int)  # will contain the bin number of each given pixel
    good = np.zeros(n, dtype=bool)   # will contain 1 if the bin has been accepted as good

    # For each point, find the distance to all other points and select the minimum.
    # This is a robust but slow way of determining the pixel size of unbinned data.
    #
    if pixelSize is None:
        pixelSize = np.min(distance.pdist(np.column_stack([x, y])))

    currentBin = np.argmax(signal/noise)  # Start from the pixel with highest S/N
    SN = signal[currentBin]/noise[currentBin]

    # Rough estimate of the expected final bin number.
    # This value is only used to give an idea of the expected
    # remaining computation time when binning very big dataset.
    #
    w = signal/noise < targetSN
    maxnum = int(np.sum((signal[w]/noise[w])**2)/targetSN**2 + np.sum(~w))

    # The first bin will be assigned CLASS = 1
    # With N pixels there will be at most N bins
    #
    for ind in range(1, n+1):

        if not quiet:
            print(ind, ' / ', maxnum)

        classe[currentBin] = ind  # Here currentBin is still made of one pixel
        xBar, yBar = x[currentBin], y[currentBin]    # Centroid of one pixels

        while True:

            if np.all(classe):
                break  # Stops if all pixels are binned

            # Find the unbinned pixel closest to the centroid of the current bin
            #
            unBinned = np.where(classe == 0)[0]
            k = np.argmin((x[unBinned] - xBar)**2 + (y[unBinned] - yBar)**2)

            # (1) Find the distance from the closest pixel to the current bin
            #
            minDist = np.min((x[currentBin] - x[unBinned[k]])**2 + (y[currentBin] - y[unBinned[k]])**2)

            # (2) Estimate the `roundness' of the POSSIBLE new bin
            #
            nextBin = np.append(currentBin, unBinned[k])
            roundness = _roundness(x[nextBin], y[nextBin], pixelSize)

            # (3) Compute the S/N one would obtain by adding
            # the CANDIDATE pixel to the current bin
            #
            SNOld = SN
            SN = aggregator(signal[nextBin], noise[nextBin])

            # Test whether (1) the CANDIDATE pixel is connected to the
            # current bin, (2) whether the POSSIBLE new bin is round enough
            # and (3) whether the resulting S/N would get closer to targetSN
            #
            if (np.sqrt(minDist) > 1.2*pixelSize or roundness > 0.3
                or abs(SN - targetSN) > abs(SNOld - targetSN)):
                if SNOld > 0.8*targetSN:
                    good[currentBin] = 1
                break

            # If all the above 3 tests are negative then accept the CANDIDATE
            # pixel, add it to the current bin, and continue accreting pixels
            #
            classe[unBinned[k]] = ind
            currentBin = nextBin

            # Update the centroid of the current bin
            #
            xBar, yBar = np.mean(x[currentBin]), np.mean(y[currentBin])

        # Get the centroid of all the binned pixels
        #
        binned = classe > 0
        if np.all(binned):
            break  # Stop if all pixels are binned
        xBar, yBar = np.mean(x[binned]), np.mean(y[binned])

        # Find the closest unbinned pixel to the centroid of all
        # the binned pixels, and start a new bin from that pixel.
        #
        unBinned = np.where(classe == 0)[0]
        k = np.argmin((x[unBinned] - xBar)**2 + (y[unBinned] - yBar)**2)
        currentBin = unBinned[k]    # The bin is initially made of one pixel
        SN = signal[currentBin]/noise[currentBin]

    classe *= good  # Set to zero all bins that did not reach the target S/N

    return classe, pixelSize

#----------------------------------------------------------------------------

def _reassign_bad_bins(classe, x, y):
    """
    Implements steps (vi)-(vii) in section 5.1 of Cappellari & Copin (2003)

    """
    # Find the centroid of all successful bins.
    # CLASS = 0 are unbinned pixels which are excluded.
    #
    good = np.unique(classe[classe > 0])
    xnode = ndimage.mean(x, labels=classe, index=good)
    ynode = ndimage.mean(y, labels=classe, index=good)

    # Reassign pixels of bins with S/N < targetSN
    # to the closest centroid of a good bin
    #
    bad = classe == 0
    index = np.argmin((x[bad, None] - xnode)**2 + (y[bad, None] - ynode)**2, axis=1)
    classe[bad] = good[index]

    # Recompute all centroids of the reassigned bins.
    # These will be used as starting points for the CVT.
    #
    good = np.unique(classe)
    xnode = ndimage.mean(x, labels=classe, index=good)
    ynode = ndimage.mean(y, labels=classe, index=good)

    return xnode, ynode

#----------------------------------------------------------------------------

def _cvt_equal_mass(x, y, signal, noise, xnode, ynode, quiet, wvt, aggregator = _sn_func):
    """
    Implements the modified Lloyd algorithm
    in section 4.1 of Cappellari & Copin (2003).

    NB: When the keyword WVT is set this routine includes
    the modification proposed by Diehl & Statler (2006).

    """
    if wvt:
        dens = np.ones_like(signal)
    else:
        dens = (signal/noise)**2  # See beginning of section 4.1 of CC03
    scale = np.ones_like(xnode)   # Start with the same scale length for all bins

    for it in range(1, xnode.size):  # Do at most xnode.size iterations

        xnodeOld, ynodeOld = xnode.copy(), ynode.copy()

        # Computes (Weighted) Voronoi Tessellation of the pixels grid
        #
        classe = np.argmin(((x[:, None] - xnode)**2 + (y[:, None] - ynode)**2)/scale**2, axis=1)

        # Computes centroids of the bins, weighted by dens**2.
        # Exponent 2 on the density produces equal-mass Voronoi bins.
        # The geometric centroids are computed if WVT keyword is set.
        #
        good = np.unique(classe)
        for k in good:
            index = classe == k   # Find subscripts of pixels in bin k.
            xnode[k], ynode[k] = _weighted_centroid(x[index], y[index], dens[index]**2)
            if wvt:
                sn = aggregator(signal[index], noise[index])
                scale[k] = np.sqrt(index.sum()/sn)  # Eq. (4) of Diehl & Statler (2006)

        diff = np.sum((xnode - xnodeOld)**2 + (ynode - ynodeOld)**2)

        if not quiet:
            print('Iter: %4i  Diff: %.4g' % (it, diff))

        if diff == 0:
            break

    # If coordinates have changed, re-compute (Weighted) Voronoi Tessellation of the pixels grid
    #
    if diff > 0:
        classe = np.argmin(((x[:, None] - xnode)**2 + (y[:, None] - ynode)**2)/scale**2, axis=1)
        good = np.unique(classe)  # Check for zero-size Voronoi bins

    # Only return the generators and scales of the nonzero Voronoi bins

    return xnode[good], ynode[good], scale[good], it

#-----------------------------------------------------------------------

def _compute_useful_bin_quantities(x, y, signal, noise, xnode, ynode, scale, aggregator = _sn_func):
    """
    Recomputes (Weighted) Voronoi Tessellation of the pixels grid to make sure
    that the class number corresponds to the proper Voronoi generator.
    This is done to take into account possible zero-size Voronoi bins
    in output from the previous CVT (or WVT).

    """
    # classe will contain the bin number of each given pixel
    #
    classe = np.argmin(((x[:, None] - xnode)**2 + (y[:, None] - ynode)**2)/scale**2, axis=1)

    # At the end of the computation evaluate the bin luminosity-weighted
    # centroids (xbar, ybar) and the corresponding final S/N of each bin.
    #
    xbar = np.empty_like(xnode)
    ybar = np.empty_like(xnode)
    sn = np.empty_like(xnode)
    area = np.empty_like(xnode)
    good = np.unique(classe)
    for k in good:
        index = classe == k   # Find subscripts of pixels in bin k.
        xbar[k], ybar[k] = _weighted_centroid(x[index], y[index], signal[index])
        sn[k] = aggregator(signal[index], noise[index])
        area[k] = index.sum()

    return classe, xbar, ybar, sn, area

#-----------------------------------------------------------------------

def _display_pixels(x, y, counts, pixelSize):
    """
    Display pixels at coordinates (x, y) coloured with "counts".
    This routine is fast but not fully general as it assumes the spaxels
    are on a regular grid. This needs not be the case for Voronoi binning.

    """
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    nx = round((xmax - xmin)/pixelSize) + 1
    ny = round((ymax - ymin)/pixelSize) + 1
    img = np.full((nx, ny), np.nan)  # use nan for missing data
    j = np.round((x - xmin)/pixelSize).astype(int)
    k = np.round((y - ymin)/pixelSize).astype(int)
    img[j, k] = counts

    plt.imshow(np.rot90(img), interpolation='none', cmap='prism',
               extent=[xmin - pixelSize/2, xmax + pixelSize/2,
                       ymin - pixelSize/2, ymax + pixelSize/2])

#----------------------------------------------------------------------

def voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=True,
                         pixelsize=None, plot=True, quiet=True,
                         wvt=True, aggregator = _sn_func):
    """
    PURPOSE:
          Perform adaptive spatial binning of Integral-Field Spectroscopic
          (IFS) data to reach a chosen constant signal-to-noise ratio per bin.
          This method is required for the proper analysis of IFS
          observations, but can also be used for standard photometric
          imagery or any other two-dimensional data.
          This program precisely implements the algorithm described in
          section 5.1 of the reference below.

    EXPLANATION:
          Further information on VORONOI_2D_BINNING algorithm can be found in
          Cappellari M., Copin Y., 2003, MNRAS, 342, 345

    CALLING SEQUENCE:

        binNum, xBin, yBin, xBar, yBar, sn, nPixels, scale = \
            voronoi_2d_binning(x, y, signal, noise, targetSN,
                               plot=True, quiet=False, wvt=False,
                               cvt=True, pixelsize=None)

    """
    # This is the main program that has to be called from external programs.
    # It simply calls in sequence the different steps of the algorithms
    # and optionally plots the results at the end of the calculation.

    if not (x.size == y.size == signal.size == noise.size):
        raise ValueError('Input vectors (x, y, signal, noise) must have the same size')
    if not np.all((noise > 0) & np.isfinite(noise)):
        raise ValueError('NOISE must be a positive vector')

    # Perform basic tests to catch common input errors
    #
    if np.sum(signal)/np.sqrt(np.sum(noise**2)) < targetSN:
        raise ValueError("""Not enough S/N in the whole set of pixels.
            Many pixels may have noise but virtually no signal.
            They should not be included in the set to bin,
            or the pixels should be optimally weighted.
            See Cappellari & Copin (2003, Sec.2.1) and README file.""")
    if np.min(signal/noise) > targetSN:
        raise ValueError('All pixels have enough S/N and binning is not needed')

    # Prevent division by zero for pixels with signal=0 and
    # noise=sqrt(signal)=0 as can happen with X-ray data
    #
    noise = noise.clip(np.min(noise[noise > 0])*1e-9)

    print('Bin-accretion...')
    classe, pixelsize = _accretion(x, y, signal, noise, targetSN, pixelsize, quiet, aggregator = aggregator)
    print(np.max(classe), ' initial bins.')
    print('Reassign bad bins...')
    xNode, yNode = _reassign_bad_bins(classe, x, y)
    print(xNode.size, ' good bins.')
    if cvt:
        print('Modified Lloyd algorithm...')
        xNode, yNode, scale, it = _cvt_equal_mass(x, y, signal, noise, xNode, yNode,
                                                  quiet, wvt, aggregator = aggregator)
        print(it-1, ' iterations.')
    else:
        scale = 1.
    classe, xBar, yBar, sn, area = _compute_useful_bin_quantities(x, y, signal, noise, xNode, yNode, scale,
                                                                  aggregator = aggregator)
    w = area == 1
    print('Unbinned pixels: ', np.sum(w), ' / ', x.size)
    print('Fractional S/N scatter (%):', np.std(sn[~w] - targetSN, ddof=1)/targetSN*100)

    if plot:
        plt.clf()
        plt.subplot(211)
        rnd = np.argsort(np.random.random(xNode.size))  # Randomize bin colors
        _display_pixels(x, y, rnd[classe], pixelsize)
        plt.plot(xNode, yNode, '+w', scalex=False, scaley=False) # do not rescale after imshow()
        plt.xlabel('R (arcsec)')
        plt.ylabel('R (arcsec)')
        plt.title('Map of Voronoi bins')

        plt.subplot(212)
        rad = np.sqrt(xBar**2 + yBar**2)  # Use centroids, NOT generators
        plt.plot(rad[~w], sn[~w], 'or', label='Voronoi bins')
        plt.xlabel('R (arcsec)')
        plt.ylabel('Bin S/N')
        plt.axis([np.min(rad), np.max(rad), 0, np.max(sn)])  # x0, x1, y0, y1
        if np.sum(w) > 0:
            plt.plot(rad[w], sn[w], 'xb', label='single spaxels')
        plt.axhline(targetSN)
        plt.legend()
        plt.pause(0.01)  # allow plot to appear in certain cases

    return classe, xNode, yNode, xBar, yBar, sn, area, scale

#----------------------------------------------------------------------------
