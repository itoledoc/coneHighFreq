import distutils
import time as time_utilities
import socket
import re
import fnmatch, pickle, traceback, copy as python_copy
import scipy.odr
import distutils.spawn
import pylab as pb
import calDatabaseQuery
import os
from scipy import optimize
from types import NoneType
from xmlrpclib import ServerProxy
from pylab import *


def find_nearest_vlbi_sources(
        radec_string, radius=3, frequency=None, catalog='rfc_2015c_cat.txt',
        server='http://asa.alma.cl/sourcecat/xmlrpc', phasecal_flux_limit=0.150,
        checksource_flux_limit=0.05, ageLimit=30):
    """
    Reads a VLBI catalog from astrogeo.org and finds all sources within
    specified radius of the specified target. IT downloads the file first, if
    necessary.
    radecString: various acceptable formats, including '17:20:05.3 -03:20:15'
                 (see help au.radec2rad for a complete description)
                 Also accepts proper names such as J2219+1806 and 3C273
    radius: search radius (in degrees)
    frequency: if not None, then also query for the flux densities
    server: this is only needed to instantiate the calDatabaseQuery object
            that holds the catalog parsing logic
    phasecalFluxLimit: minimum flux required (in Janskys)
    checksourceFluxLimit: minimum flux required (in Janskys)
    ageLimit: maximum age (in days) to consider flux measurement current
    Args:
        radec_string:
        radius:
        frequency:
        catalog:
        server:
        phasecal_flux_limit:
        checksource_flux_limit:
        ageLimit:

    Returns:
    * if frequency is not specified: a list of source names
    * otherwise: a dictionary keyed by source name, with value equal to the
                 predicted flux density
    -Todd Hunter
    """

    sourcecat = ServerProxy(server)
    sourcelist, stypes = findnearestvlbi(
        radec_string, radius, catalog, sourcecat, return_dictionary=True)
    if frequency is None:
        sourceinfo = {}
        for src in sourcelist.keys():
            sourceinfo[src] = {'separation': sourcelist[src],
                               'type': stypes[src]}
        return sourceinfo
    else:
        print "Querying ALMA catalog for flux densities."
        mydict = {}
        maxflux = phasecal_flux_limit
        uncertainty = 0
        bestsrc = 'none'
        allflux = []
        allfluxsource = []
        allages = []
        for src in sourcelist.keys():
            mydict[src] = get_alma_flux(src, frequency, silent=True)
            if mydict[src] == -1:
                mydict[src] = None
            else:
                mydict[src]['separationDegrees'] = sourcelist[src]
                mydict[src]['type'] = stypes[src]
                allflux.append(mydict[src]['fluxDensity'])
                allfluxsource.append(src)
                allages.append(mydict[src]['meanAge'])
                if mydict[src]['fluxDensity'] > maxflux:
                    bestsrc = src
                    maxflux = mydict[src]['fluxDensity']
                    uncertainty = mydict[src]['fluxDensityUncertainty']
        print "Found %d objects with ALMA measurements." % (len(allflux))
        if bestsrc != 'none':
            print "The highest flux density object is %s with %.1f+-%.1f mJy at %.1f deg." % (
            bestsrc, maxflux * 1000, uncertainty * 1000,
            mydict[bestsrc]['separationDegrees'])
        totalSources = len(sourcelist.keys())
        allflux = np.array(allflux)
        allfluxsource = np.array(allfluxsource)
        allages = np.array(allages)
        idx = np.argsort(allflux)
        allflux = allflux[idx]
        allfluxsource = allfluxsource[idx]
        allages = allages[idx]
        phasecal = None
        checksource = None
        if len(allflux) < 1:
            print "No sources have millimeter flux density measurements"
            return mydict
        if allflux[-1] > phasecal_flux_limit:
            phasecal = allfluxsource[-1]
            phasecalFlux = allflux[-1] * 1000
            if (allflux[-1] > 2 * phasecal_flux_limit or allages[
                -1] > ageLimit):
                phasecalReobserve = False
            else:
                phasecalReobserve = True
            if len(allflux) > 1:
                if (allflux[
                        -2] > checksource_flux_limit):  # need to add angular separation check
                    checksource = allfluxsource[-2]
                    checksourceFlux = allflux[-2] * 1000
                    if (allflux[-2] > 2 * checksource_flux_limit or allages[
                        -2] > ageLimit):
                        checksourceReobserve = False
                    else:
                        checksourceReobserve = True
        elif allflux[-1] > checksource_flux_limit:
            checksource = allfluxsource[-1]
            checksourceFlux = allflux[-1] * 1000
            if (allflux[-1] > 2 * checksource_flux_limit or allages[
                -1] > ageLimit):
                checksourceReobserve = False
            else:
                checksourceReobserve = True

        if phasecal is None:
            print "No sources meet phasecal flux limit."
            if radius < 5:
                print "Try increasing search radius to 5 degrees to get more phasecal candidates."
            elif len(allflux) < totalSources:
                print "Schedule a cone search to measure flux densities of the other %d phasecal candidates." % (
                totalSources - len(allflux))
            else:
                print "Schedule a cone search to identify new sources for phase calibration."
        else:
            if phasecalReobserve:
                print "Recommended phase calibrator = %s, but need to reobserve it to confirm flux density (%.0f mJy)." % (
                phasecal, phasecalFlux)
            else:
                print "Recommended phase calibrator = %s (%.0f mJy)" % (
                phasecal, phasecalFlux)
        if checksource is None:
            if len(allflux) < 2:
                print "No potential check sources."
            else:
                print "None of the potential check sources meet check source flux limit."
            if radius < 5:
                print "Try increasing search radius to 5 degrees to get more check source candidates."
            elif len(allflux) < totalSources:
                print "Schedule a cone search to measure flux densities of the other %d check source candidates." % (
                totalSources - len(allflux))
            else:
                print "Schedule a cone search for new check sources."
        elif checksourceReobserve:
            print "Recommended check source = %s, but need to reobserve it to confirm flux density (%.0f mJy)." % (
            checksource, checksourceFlux)
        else:
            print "Recommended check source = %s (%.0f mJy)." % (
            checksource, checksourceFlux)

        return mydict


def findnearestvlbi(
        radec_string, radius, catalog, sourcecat, return_dictionary=True):
    if not os.path.exists(catalog):
        print "Downloading catalog %s from astrogeo.org" % catalog
        wget = distutils.spawn.find_executable(
            'wget', path=':'.join(sys.path) + ':' + os.environ['PATH'])
        cmd = wget + ' astrogeo.org/vlbi/solutions/%s/%s' % (
            catalog.replace('_cat.txt', ''), catalog)
        os.system(cmd)

    sname, stype, s_ra, s_dec, s_radeg, s_decdeg = read_vlbi_catalog(catalog)
    rarad, decrad = radec2rad(radec_string)
    sources = []
    mydict = {}
    typedict = {}
    for i in range(len(sname)):
        rav = np.radians(float(s_radeg[i]))
        decv = np.radians(float(s_decdeg[i]))
        separation = angularseparationradians(rarad, decrad, rav, decv)
        if separation < np.radians(radius):
            sources.append(sname[i])
            mydict[sname[i]] = np.degrees(separation)
            typedict[sname[i]] = stype[i]
    print "Found %d sources within %g degrees" % (len(sources), radius)
    if return_dictionary:
        return mydict, typedict
    else:
        return sources


def read_vlbi_catalog(catalog):
    file_obj = open(catalog, 'r')
    lines = file_obj.readlines()
    file_obj.close()

    sname = []
    stype = []
    s_ra = []
    s_dec = []
    s_radeg = []
    s_decdeg = []

    # loop over all sources in VLBI database
    for line in lines:
        if line[0] != '#':
            fields = line.split()
            if len(fields) != 25:
                print 'bad line ', line
                break
            stype.append(fields[0])
            sname.append(fields[2])
            ra = fields[3] + ':' + fields[4] + ':' + fields[5]
            s_ra.append(ra)
            s_radeg.append(dms2decimaldeg(ra, True))
            dec = fields[6] + ':' + fields[7] + ':' + fields[8]
            s_dec.append(dec)
            s_decdeg.append(dms2decimaldeg(dec, False))

    nvlbi = len(sname)
    print 'number of VLBI sources ', nvlbi
    return sname, stype, s_ra, s_dec, s_radeg, s_decdeg


def dms2decimaldeg(value, raflag):
    """
    Converts decimal degrees to Hours:Minutes:Seconds
    If ra is True, then we are converting a Ra measurment
    and we devide by 15 to go from 0-->360deg to 0--->24 hours.
    Else, we are converting a signed Dec measurement

    Args:
        value:
        raflag:

    """
    if raflag:
        ra = 15.
    else:
        ra = 1.

    vals = value.split(':')
    if vals[0][0] == '-':
        sgn = -1
        vals[0] = - float(vals[0])
    else:
        sgn = 1

    value = float(vals[0]) + float(vals[1]) / 60.0 + float(vals[2]) / 3600.

    return ra * sgn * value


def radec2rad(radecstring):
    """
    Convert a position from a single RA/Dec sexagesimal string to RA and
    Dec in radians.
    The RA and Dec portions can be separated by a comma or a space.
    The RA portion of the string must be colon-delimited, space-delimited,
        or 'h/m/s' delimited.
    The Dec portion of the string can be either ":", "." or space-delimited.
    If it is "." delimited, then it must have degrees, minutes, *and* seconds.
    See also rad2radec.
    -Todd Hunter

    Args:
        radecstring:
    """
    if radecstring.find('h') > 0 and radecstring.find('d') > 0:
        radecstring = radecstring.replace(
            'h', ':').replace('m', ':').replace('d', ':').replace('s', '')
    radec1 = radecstring.replace(',', ' ')
    tokens = radec1.split()
    if len(tokens) == 2:
        (ra, dec) = radec1.split()
    elif len(tokens) == 6:
        h, m, s, d, dm, ds = radec1.split()
        ra = '%s:%s:%s' % (h, m, s)
        dec = '%+f:%s:%s' % (float(d), dm, ds)
    else:
        print "Invalid format for RA/Dec string"
        return
    tokens = ra.strip().split(':')
    hours = 0
    for i, t in enumerate(tokens):
        hours += float(t) / (60. ** i)
    if dec.find(':') > 0:
        tokens = dec.lstrip().split(':')
    elif dec.find('.') > 0:
        try:
            (d, m, s) = dec.lstrip().split('.')
        except:
            (d, m, s, sfraction) = dec.lstrip().split('.')
            s = s + '.' + sfraction
        tokens = [d, m, s]
    else:  # just an integer
        tokens = [dec]
    dec1 = 0
    for i, t in enumerate(tokens):
        dec1 += abs(float(t) / (60. ** i))
    if dec.lstrip().find('-') == 0:
        dec1 = -dec1
    decrad = dec1 * np.pi / 180.
    ra1 = hours * 15
    rarad = ra1 * np.pi / 180.
    return rarad, decrad


def angularseparationradians(ra0, dec0, ra1, dec1, return_components=False):
    """
    Computes the great circle angle between two celestial coordinates.
    using the Vincenty formula (from wikipedia) which is correct for all
    angles, as long as you use atan2() to handle a zero denominator.
     See  http://en.wikipedia.org/wiki/Great_circle_distance
    Input and output are in radians.  It also works for the az,el coordinate
    system.

    returnComponents=True
    will return: [separation, raSeparation, decSeparation, raSeparationCosDec]
    See also angularSeparation()
    -- Todd Hunter

    Args:
        return_components:
        dec1:
        ra1:
        dec0:
        ra0:
    """
    result = angular_separation(
        ra0 * 180 / math.pi, dec0 * 180 / math.pi, ra1 * 180 / math.pi,
        dec1 * 180 / math.pi,
        return_components)
    if return_components:
        return np.array(result) * math.pi / 180.

    return result * math.pi / 180.


def angular_separation(ra0, dec0, ra1, dec1, return_components=False):
    """
    Computes the great circle angle_s between two celestial coordinates.
    using the Vincenty formula (from wikipedia) which is correct for all
    angles, as long as you use atan2() to handle a zero denominator.
     See  http://en.wikipedia.org/wiki/Great_circle_distance
    ra,dec must be given in degrees, as is the output.
    It also works for the az,el coordinate system.
    Component separations are field_0 minus field_1.
    See also angularSeparationRadians()
    -- Todd Hunter

    Args:
        return_components:
        dec1:
        ra1:
        dec0:
        ra0:
    """
    ra0 *= math.pi / 180.
    dec0 *= math.pi / 180.
    ra1 *= math.pi / 180.
    dec1 *= math.pi / 180.
    deltalong = ra0 - ra1
    argument1 = (((math.cos(dec1) * math.sin(deltalong)) ** 2) +
                 ((math.cos(dec0) * math.sin(dec1) - math.sin(dec0) * math.cos(
                     dec1) * math.cos(deltalong)) ** 2)) ** 0.5
    argument2 = math.sin(dec0) * math.sin(dec1) + math.cos(dec0) * math.cos(
        dec1) * math.cos(deltalong)
    angle_s = math.atan2(argument1, argument2) / (math.pi / 180.)
    if angle_s > 360:
        angle_s -= 360
    if return_components:
        cosdec = math.cos((dec1 + dec0) * 0.5)
        radegreescosdec = np.degrees(ra0 - ra1) * cosdec
        radegrees = np.degrees(ra0 - ra1)
        decdegrees = np.degrees(dec0 - dec1)
        if radegrees > 360:
            radegrees -= 360
        if decdegrees > 360:
            decdegrees -= 360
        retval = angle_s, radegrees, decdegrees, radegreescosdec
    else:
        retval = angle_s
    return retval


def get_alma_flux(
        sourcename, frequency, lowband=3, highband=7, ignorelowband=False,
        ignorehighband=False, simulatelowband=False, simulateHighBand=False,
        defaultSpectralIndex=-0.6, defaultSpectralIndexUncertainty=0.2,
        verbose=False, trials=10000, date='', searchAdjacentNames=False,
        server='', dayWindow=0, showplot=False, plotfile='', silent=False,
        separationThreshold=14, maximumSensibleSpectralIndex=0.0,
        lowbandFlux=None,
        highbandFlux=None, lowbandFrequency=None, highbandFrequency=None,
        lowbandUncertainty=None, highbandUncertainty=None,
        lowbandDate=None, highbandDate=None):
    """
    Queries the ALMA calibrator catalog for flux density measurements in two
    bands.
    It first computes the spectral index, then interpolates or extrapolates
    to the
    desired frequency.  It then performs a Monte-Carlo simulation to obtain
    the flux
    density uncertainty.  If the date is not specified, it simply gets the
    most recent
    measurements.  If the date is specified, then it finds the measurements
    in each
    band that are closest to that date, either in the past or future.

    frequency: the frequency for which to estimate the flux density.  It can
    be a
       floating point value in Hz or GHz, or string with units (Hz, kHz, MHz,
       GHz).
    lowband: the ALMA band to use for the low frequency measurement
    highband: the ALMA band to use for the high frequency measurement
    ignorelowband: ignore any measurements in the low frequency band
    ignorehighband: ignore any measurements in the high frequency band
    simulatelowband: set to 3.0 Jy (meant only for allow testing)
    simulateHighBand: set to 1.0 Jy (meant only for testing)
    defaultSpectralIndex: spectral index to use if only one band measurement
    is available
    defaultSpectralIndexUncertainty:  uncertainty to use
    trials: the number of trials to run in the Monte-Carlo simulation of the
    uncertainty
    date: the date in the past to begin looking for earlier measurements (
    default=today)
          formats accepted: '20120101' or '2012-01-01' or '2012/01/01'  where
          the
                             delimiter can be any non-integer character
    searchAdjacentNames: pass this flag to au.searchFlux
    server: pass this string to au.searchFlux (name of xmlrpc database URL)
            can be 'internal', 'external' or full URL
    dayWindow: if non-negative, then process all measurements within this
               many days of the first measurement found (per band)
    showplot: if True, then produce a plot with errorbars and model
    plotfile: write the plot to a file
    silent: if True, then don't print the normal status messages
    separationThreshold: in days, if Band3/7 measurements are further apart
    then this, or
          the mean age of all measurements used greater than this, then write
          a warning
    maximumSensibleSpectralIndex: if larger than this, then write a warning

    The following parameters allow you to use a different measurement from what
    is in the catalog:
    lowbandFlux: recently measured value to use instead of the database value
    (Jy)
    highbandFlux: recently measured value to use instead of the database
    value (Jy)
    lowbandFrequency: frequency of lowbandFlux
    highbandFrequency: frequency of highbandFlux
    lowbandUncertainty: uncertainty of lowbandFlux
    highbandUncertainty: uncertainty of highbandFlux
    lowbandDate: date for lowbandFlux
    highbandDate: date for highbandFlux

    Args:
        highbandDate:
        lowbandDate:
        highbandUncertainty:
        lowbandUncertainty:
        highbandFrequency:
        lowbandFrequency:
        highbandFlux:
        lowbandFlux:
        maximumSensibleSpectralIndex:
        separationThreshold:
        silent:
        plotfile:
        showplot:
        dayWindow:
        server:
        searchAdjacentNames:
        date:
        trials:
        verbose:
        defaultSpectralIndexUncertainty:
        defaultSpectralIndex:
        simulateHighBand:
        simulatelowband:
        ignorehighband:
        ignorelowband:
        highband:
        lowband:
        frequency:
        sourcename:
    """

    if date == '':
        date = getcurrentdate(delimiter='-')
    if (lowbandFlux is not None or
            lowbandFrequency is not None or
            lowbandUncertainty is not None or
            lowbandDate is not None):
        if (lowbandFlux is None or
                lowbandFrequency is None or
                lowbandUncertainty is None or
                lowbandDate is None):
            print("lowbandFlux, lowbandFrequency, lowbandUncertainty,"
                  "lowbandDate must all be specified.")
            return

        lowbandFrequency = parseFrequencyArgumentToHz(lowbandFrequency)

    if (highbandFlux is not None or
            highbandFrequency is not None or
            highbandUncertainty is not None or
            highbandDate is not None):
        if (highbandFlux is None or
                highbandFrequency is None or
                highbandUncertainty is None or
                highbandDate is None):
            print("highbandFlux, highbandFrequency, highbandUncertainty, "
                  "highbandDate must all be specified.")
            return
        highbandFrequency = parseFrequencyArgument(highbandFrequency)

    if type(frequency) == str:
        frequency = parseFrequencyArgumentToHz(frequency)

    if frequency < 2000:
        frequency *= 1e9

    allresults = {}
    noLowBandMeasurement = True
    noHighBandMeasurement = True

    if ignorelowband and ignorehighband:
        print("No measurements available.")
        return

    for band in [lowband, 6, highband, 9]:
        if band == lowband and ignorelowband:
            continue
        if band == highband and ignorehighband:
            continue
        if verbose:
            print(
                "Calling searchFlux('%s', band='%s', date='%s', "
                "returnMostRecent=True, verbose=%s, "
                "searchAdjacentNames=%s)" % (
                    sourcename, band, date, verbose, searchAdjacentNames))
        results = searchFlux(sourcename, band=band, date=date,
                             returnMostRecent=True,
                             verbose=verbose,
                             searchAdjacentNames=searchAdjacentNames,
                             server=server, dayWindow=dayWindow)
        if band == 9:
            band = highband

        if band == 6:
            band = lowband

        if verbose:
            print("Band %d results = %s" % (band, results))

        if results is None:
            return None

        if results == -1:
            if not silent:
                print("No band %s measurements found for %s in the catalog" % (
                    band, sourcename))
            continue  # source not found in catalog

        if dayWindow < 0:
            results = [results]

        allresults[band] = {}  # each key is an array of the same length
        allresults[band]['frequency'] = []
        allresults[band]['flux'] = []
        allresults[band]['uncertainty'] = []
        allresults[band]['age'] = []

        for result in results:
            if (type(result) == NoneType or
                    simulatelowband or
                    simulateHighBand):
                if type(result) == NoneType:
                    print "The calibrator database is not currently accessible."

                if simulatelowband or simulateHighBand:
                    if simulatelowband:
                        allresults[lowband]['frequency'].append(90e9)
                        allresults[lowband]['flux'].append(3)
                        allresults[lowband]['uncertainty'].append(0.3)
                        allresults[lowband]['age'].append(100)
                        noLowBandMeasurement = False

                    if simulateHighBand:
                        allresults[lowband]['frequency'].append(340e9)
                        allresults[lowband]['flux'].append(1)
                        allresults[lowband]['uncertainty'].append(0.1)
                        allresults[lowband]['age'].append(100)
                        noHighBandMeasurement = False

                else:
                    return None

            elif type(result) == int:
                print("No Band %d observation in the catalog for this source" %
                      band)
                if band == highband:
                    if noLowBandMeasurement:
                        print("No measurement in either band.  Cannot produce "
                              "a flux density.")
                        return None
            else:
                # Clean up the existing measurement, if necessary
                if result['uncertainty'] == 0.0:
                    print("No uncertainty in the catalog for this measurement, "
                          "assuming 10 percent.")
                    result['uncertainty'] = 0.1 * (result['flux'])
                if result['flux'] <= 0:
                    print("WARNING: Found flux density = %f, discarding" % (
                        result['flux']))
                    continue
                if band == lowband:
                    noLowBandMeasurement = False
                else:
                    noHighBandMeasurement = False
                if band == lowband and lowbandFlux is not None:
                    allresults[band]['frequency'].append(lowbandFrequency)
                    allresults[band]['flux'].append(lowbandFlux)
                    allresults[band]['age'].append(
                        computeIntervalBetweenTwoDays(date, lowbandDate))
                    allresults[band]['uncertainty'].append(lowbandUncertainty)
                elif band == highband and highbandFlux is not None:
                    allresults[band]['flux'].append(result['flux'])
                    allresults[band]['frequency'].append(highbandFlux)
                    allresults[band]['age'].append(
                        computeIntervalBetweenTwoDays(date, highbandDate))
                    allresults[band]['uncertainty'].append(highbandUncertainty)
                else:
                    allresults[band]['frequency'].append(result['frequency'])
                    allresults[band]['flux'].append(result['flux'])
                    allresults[band]['uncertainty'].append(
                        result['uncertainty'])
                    allresults[band]['age'].append(result['age'])
                if not silent:
                    print(
                        "Using Band %d measurement: %.3f +- %.3f "
                        "(age=%d days) %.1f GHz" % (
                            band, allresults[band]['flux'][-1],
                            allresults[band]['uncertainty'][-1],
                            allresults[band]['age'][-1],
                            allresults[band]['frequency'][-1] * 1e-9))

        if date != '' and dayWindow < 0:
            if type(result) != int:
                futureDate = computeFutureDate(date, result['age'])
            else:
                futureDate = computeFutureDate(date,
                                               365 * 10)  # look ahead 10 years
            futureResult = searchFlux(sourcename, band=band, date=futureDate,
                                      returnMostRecent=True, verbose=verbose,
                                      server=server)
            if type(futureResult) != int:
                if futureResult['date'] != result['date']:
                    # then the measurement after the observation date is closer
                    # in time than the one before, so use it
                    result = python_copy.deepcopy(futureResult)
                    if result['uncertainty'] == 0.0:
                        print("No uncertainty in the catalog for this "
                              "measurement, assuming 10 percent.")
                        result['uncertainty'] = 0.1 * (result['flux'])
                    allresults[band]['frequency'].append(result['frequency'])
                    allresults[band]['flux'].append(result['flux'])
                    allresults[band]['uncertainty'].append(
                        result['uncertainty'])
                    allresults[band]['age'].append(result['age'])
    # end 'for' loop over band
    if allresults == {}:
        return -1

    spectralIndex = defaultSpectralIndex
    spectralIndexUncertainty = defaultSpectralIndexUncertainty

    if noLowBandMeasurement or noHighBandMeasurement:
        if noLowBandMeasurement:
            band = highband
        else:
            band = lowband
        #        print "allresults.keys() = ", allresults.keys()
        freqs = [allresults[band]['frequency'][0]]
        fluxDensity = (
            allresults[band]['flux'][0] *
            (frequency / allresults[band]['frequency'][0]) ** spectralIndex)

        # The following does not account for uncertainty in the spectral index,
        # but prevents a crash when defining mydict later.

        fluxDensityUncertainty = (
            fluxDensity *
            allresults[band]['uncertainty'][0] / allresults[band]['flux'][0])

        intercept = np.log10(fluxDensity)

        # interceptUncertainty should simply be the flux density uncertainty on
        # the log scale
        interceptUncertainty = (0.434 * allresults[band]['uncertainty'][0] /
                                allresults[band]['flux'][0])
        meanAge = allresults[band]['age'][0]
        meanOfLogX = np.log10(allresults[band]['frequency'][0] * 1e-9)
        ageDifference = 0
    else:
        # compute the mean interval in days between the science observation and
        # the flux monitoring observation
        meanAge = np.mean(
            np.abs([allresults[lowband]['age'] + allresults[highband]['age']]))
        ageDifference = fabs(np.mean(allresults[lowband]['age']) - np.mean(
            allresults[highband]['age']))
        freqs = allresults[lowband]['frequency'] + allresults[highband][
            'frequency']
        fluxes = allresults[lowband]['flux'] + allresults[highband]['flux']
        errors = allresults[lowband]['uncertainty'] + allresults[highband][
            'uncertainty']
        if verbose:
            print(
                "Calling linfit().spectralindex(freqs=%s, fluxes=%s, errors=%s,"
                " showplot=%s, plotfile='%s', source='%s', silent=%s)" %
                (str(freqs), str(fluxes), str(errors), showplot, plotfile,
                 sourcename, silent))
        mydict = linfit().spectralindex(freqs=freqs, fluxes=fluxes,
                                        errors=errors, showplot=showplot,
                                        plotfile=plotfile, source=sourcename,
                                        silent=silent)
        meanOfLogX = mydict['meanOfLogX']
        spectralIndex = mydict['spectralIndex']
        spectralIndexUncertainty = mydict['spectralIndexUncertainty']
        intercept = mydict['intercept']
        interceptUncertainty = mydict['interceptUncertainty']
        fluxDensity = fluxes[0] * (frequency / freqs[0]) ** spectralIndex
    logfit = True
    if not noLowBandMeasurement and not noHighBandMeasurement:
        fluxDensityUncertainty, monteCarloFluxDensity = (
            linfit().computeStdDevMonteCarlo(
                spectralIndex,
                spectralIndexUncertainty,
                intercept,
                interceptUncertainty,
                frequency * 1e-9,
                trials, logfit, meanOfLogX,
                covar=mydict['covar'], silent=silent, returnMedian=True))
    else:
        # no monte carlo available, so just set it to the spectral index f.d.
        monteCarloFluxDensity = fluxDensity
    mydict = {'fluxDensity': fluxDensity, 'spectralIndex': spectralIndex,
              'spectralIndexUncertainty': spectralIndexUncertainty,
              'fluxDensityUncertainty': fluxDensityUncertainty,
              'meanAge': meanAge,
              'ageDifference': ageDifference,
              'monteCarloFluxDensity': monteCarloFluxDensity}
    if not silent:
        print "Result using spectral index of %.3f for %.6f GHz from %.6f GHz = %.6f +- %.6f Jy" % (
        spectralIndex, frequency * 1e-9, freqs[0] * 1e-9, mydict['fluxDensity'],
        mydict['fluxDensityUncertainty'])
        if meanAge > separationThreshold:
            print "WARNING: the mean time separation between the target date and the flux monitoring observations is %d days" % (
            meanAge)
        if ageDifference > separationThreshold:
            print "WARNING: the time separation between the Band %d and %d measurements is %d days" % (
            lowband, highband, ageDifference)
        if spectralIndex > maximumSensibleSpectralIndex:
            print 'WARNING: The spectral index of %+.1f is unusual for a quasar. There may be a problem with the flux monitoring data.' % (
            spectralIndex)
    return mydict


def getcurrentdate(delimiter='/'):
    """
    returns date in format: 'YYYY-MM-DD'
    -Todd Hunter

    Args:
        delimiter:
    """
    return time_utilities.strftime('%Y/%m/%d').replace('/', delimiter)


def parseFrequencyArgumentToHz(bandwidth):
    """
    Converts a frequency string into floating point in Hz, based on the units.
    If the units are not present, then the value is assumed to be GHz if less
    than 1000 (in contrast to parseFrequencyArgument).
    -Todd Hunter

    Args:
        bandwidth:
    """
    value = parseFrequencyArgumentToGHz(bandwidth) * 1e9
    return value


def parseFrequencyArgumentToGHz(bandwidth):
    """
    Converts a frequency string into floating point in GHz, based on the units.
    If the units are not present, then the value is assumed to be GHz if less
    than 1000.
    -Todd Hunter

    Args:
        bandwidth:
    """
    value = parseFrequencyArgument(bandwidth)
    if value > 1000:
        value *= 1e-9
    return value


def parseFrequencyArgument(bandwidth):
    """
    Converts a string frequency into floating point in Hz, based on the units.
    If the units are not present, then the value is simply converted to float.
    -Todd Hunter

    Args:
        bandwidth:
    """
    bandwidth = str(bandwidth)
    ghz = bandwidth.lower().find('ghz')
    mhz = bandwidth.lower().find('mhz')
    khz = bandwidth.lower().find('khz')
    hz = bandwidth.lower().find('hz')
    if ghz > 0:
        bandwidth = 1e9 * float(bandwidth[:ghz])
    elif mhz > 0:
        bandwidth = 1e6 * float(bandwidth[:mhz])
    elif khz > 0:
        bandwidth = 1e3 * float(bandwidth[:khz])
    elif hz > 0:
        bandwidth = float(bandwidth[:hz])
    else:
        bandwidth = float(bandwidth)
    return bandwidth


def searchFlux(sourcename=None, date='', band=None, fLower=1e9, fUpper=1e12,
               tunnel=False, maxrows=10, limit=1000, debug=False,
               server='', dateCriteria=0, verbose=True, measurements=None,
               returnMostRecent=False, searchAdjacentNames=False,
               showDateReduced=False, dayWindow=-1, showPolarization=False,
               types=[1, 4, 25], sourceBandLimit=100, showAllCoordinates=False,
               returnPosition=False):
    """
    For help on the main parameters, type:
        help au.calDatabaseQuery.CalibratorCatalogUpdate.searchFlux
    -Todd Hunter

    A common option is:
    date: string, YYYYMMDD, e.g. '20120101' or '2012-01-01'
          or '2012/01/01'  where delimiter can be any non-integer character
          or YYYYmonDD or YYYY/mon/DD or YYYY-mon-DD where mon can be
          Jan/JAN/jan etc.

    Additional options:
    searchAdjacentNames: if True, search nearby names (if given in format: [
    J]HHMM[+/-]DDM[M])
    server: '', 'external', 'internal', or full URL
    returnMostRecent: if True, return a dictionary describing the most
             recent measurement
    dayWindow: if non-negative, and returnMostRecent is True, then return a
          list of matches that are within this many days of the first find
    types: a list of integers: 25:grid_source, 4:line_source, 1:point_source,
              24: polarization_source, or a list of strings:
              ['grid','line','point','polarization']
    measurements: a dictionary of measurements (eg. as returned from wrapSearch)

    Args:
        returnPosition:
        showAllCoordinates:
        sourceBandLimit:
        types:
        showPolarization:
        dayWindow:
        showDateReduced:
        searchAdjacentNames:
        returnMostRecent:
        measurements:
        verbose:
        dateCriteria:
        server:
        debug:
        limit:
        maxrows:
        tunnel:
        fUpper:
        fLower:
        band:
        date:
        sourcename:
    """
    hostname = socket.gethostname()
    if server.find('internal') >= 0:
        server = 'http://sourcecat.osf.alma.cl/sourcecat/xmlrpc'
    if server.find('external') >= 0:
        server = 'http://asa.alma.cl/sourcecat/xmlrpc'
    if 'alma.cl' in hostname:
        if server != '':
            if debug: print "Using server = ", server
            ccu = calDatabaseQuery.CalibratorCatalogUpdate(server=server)
        else:
            server = 'http://sourcecat.osf.alma.cl/sourcecat/xmlrpc'
            # There was a time when internal server failed while external kept working, so always use it.

            # server = 'http://asa.alma.cl/sourcecat/xmlrpc'
            ccu = calDatabaseQuery.CalibratorCatalogUpdate(server=server)
    else:
        if server == '':
            server = 'http://asa.alma.cl/sourcecat/xmlrpc'
        if debug:
            print "Calling calDatabaseQuery.CalibratorCatalogUpdate(tunnel=%s,server='%s')" % (
            tunnel, server)
        ccu = calDatabaseQuery.CalibratorCatalogUpdate(tunnel=tunnel,
                                                       server=server)
    if ccu.connectionFailed:
        return None
    date = replaceMonth(date)
    date = fillZerosInDate(date)
    mytypes = convertSourceTypes(types)
    status = ccu.searchFlux(sourcename, date, band, fLower, fUpper, tunnel,
                            maxrows,
                            limit, debug, server, dateCriteria, verbose,
                            measurements,
                            returnMostRecent, showDateReduced,
                            sourceBandLimit=sourceBandLimit,
                            dayWindow=dayWindow,
                            showPolarization=showPolarization,
                            types=mytypes,
                            showAllCoordinates=showAllCoordinates,
                            returnPosition=returnPosition)
    if status == -1:
        if searchAdjacentNames:
            names = getAdjacentSourceNames(sourcename)
            if names == []:
                print "This name is not an allowed format to search for adjacent names.  Must be: [J/B]HHMM+/-DDM[M]"
            else:
                for sourcename in names:
                    status = ccu.searchFlux(sourcename, date, band, fLower,
                                            fUpper,
                                            tunnel, maxrows,
                                            limit, debug, server, dateCriteria,
                                            verbose, measurements,
                                            returnMostRecent, showDateReduced,
                                            sourceBandLimit=sourceBandLimit,
                                            dayWindow=dayWindow,
                                            showPolarization=showPolarization,
                                            types=mytypes,
                                            returnPosition=returnPosition)
                    if status != -1: break
        elif verbose:
            if maxrows > 0:
                print "You could try setting searchAdjacentNames=True"
            else:
                print " "
    return status


def replaceMonth(datestring):
    """
    Replaces a 3-character month string with its 2-character integer string.
    -Todd Hunter

    Args:
        datestring:
    """
    for i, month in enumerate(['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul',
                               'aug', 'sep', 'oct', 'nov', 'dec']):
        datestring = datestring.lower().replace(month, '%02d' % (i + 1))
    return datestring.upper()


def fillZerosInDate(datestring):
    """
    Converts a string like '2014-5-6' into '2014-05-06'.

    Args:
        datestring:
    """
    if datestring.find('-') > 0:
        key = '-'
    elif datestring.find('/') > 0:
        key = '/'
    else:
        return datestring
    datestring = datestring.split(key)
    datestring = ['%02d' % int(i) for i in datestring]
    datestring = key.join(datestring)
    return datestring


def convertSourceTypes(types):
    # Convert a list of calibrator catalog source types from string to
    # integer code.
    # Todd Hunter
    mytypes = []
    sourceTypes = {'grid': 25, 'line': 4, 'point': 1, 'polarization': 24}
    if type(types) == list or type(types) == np.ndarray:
        if type(types[0]) == str:
            for t in types:
                if t not in sourceTypes.keys():
                    print "Illegal source type = ", t
                    return
                mytypes.append(sourceTypes[t])
        else:
            mytypes = types
    elif type(types) == str:
        # Assume it is a comma-delimited string
        types = types.split(',')
        for t in types:
            if t not in sourceTypes.keys():
                print "Illegal source type = ", t
                return
            mytypes.append(sourceTypes[t])
    else:
        mytypes = [types]
    return mytypes


def getAdjacentSourceNames(sourcename):
    """
    Takes a source name of the format [J/B]1924[+/-]242[9]
    and returns the four adjacent names:  1925+2429, 1923+2429, 1924-2430,
    1924-2428
    -Todd Hunter

    Args:
        sourcename:
    """
    mymatch = re.match(r'.?[\d,%]{4}[+-]\d{3}', sourcename)
    sources = []
    if mymatch is not None:
        if re.match(r'[\d,%]{4}[+-]\d{3}', sourcename) is None:
            epoch = sourcename[0]
            radec = sourcename[1:]
        else:
            epoch = 'J'
            radec = sourcename
        if len(radec.split('-')) > 1:
            mysign = '-'
        else:
            mysign = '+'
        ra, dec = radec.split(mysign)
        decdigits = len(dec)
        if ra.find('%') < 0:
            ra2 = int(ra) + 1
            # convert 1960 to 2000
            if ra2 % 100 == 60: ra2 = (ra2 / 100 + 1) * 100
            if ra2 == 2400: ra2 = 0
            sources.append('%s%04d%c%s' % (epoch, ra2, mysign, dec))
            ra0 = int(ra) - 1
            # convert 1999 to 1959
            if ra0 % 100 == 99: ra0 = (ra0 / 100) * 100 + 59
            if ra0 < -1: ra0 = 2359
            sources.append('%s%04d%c%s' % (epoch, ra0, mysign, dec))
        if dec.find('%') < 0:
            dec2 = int(dec) + 1
            # convert 1960 to 2000
            if dec2 % 100 == 60:
                dec2 = (dec2 / 100 + 1) * 100
            dec0 = int(dec) - 1
            # convert 1999 to 1959
            if dec0 % 100 == 99: dec0 = (dec0 / 100) * 100 + 59
            sources.append('%s%s%c%0*d' % (epoch, ra, mysign, decdigits, dec2))
            sources.append('%s%s%c%0*d' % (epoch, ra, mysign, decdigits, dec0))
        if dec.find('%') < 0 and ra.find('%') < 0:
            sources.append(
                '%s%04d%c%0*d' % (epoch, ra2, mysign, decdigits, dec2))
            sources.append(
                '%s%04d%c%0*d' % (epoch, ra2, mysign, decdigits, dec0))
            sources.append(
                '%s%04d%c%0*d' % (epoch, ra0, mysign, decdigits, dec2))
            sources.append(
                '%s%04d%c%0*d' % (epoch, ra0, mysign, decdigits, dec0))
    return sources


def computeIntervalBetweenTwoDays(date1, date2):
    """
    Takes 2 strings of format 'YYYY-MM-DD' and returns the number of
    days between them.  Positive if date1 > date2.  See also subtractDays.
    -Todd Hunter

    Args:
        date2:
        date1:
    """
    date1 = date1.replace('-', '').replace('/', '')
    date2 = date2.replace('-', '').replace('/', '')
    delta = datetime.date(int(date1[0:4]), int(date1[4:6]), int(date1[6:])) - \
            datetime.date(int(date2[0:4]), int(date2[4:6]), int(date2[6:]))
    return delta.days


def computeFutureDate(dateString, age):
    """
    Takes a date string (returned by the ALMA calibrator database,
    i.e. YYYY-MM-DD)
    and age in days and computes a future date string to use that is the same
    number
    of days into the future from the original date.

    Args:
        age:
        dateString:
    """
    dateString = dateString.replace('/', '-')
    if dateString.find('-') < 0:
        dateString = dateString[:4] + '-' + dateString[4:6] + '-' + dateString[
                                                                    6:]
    mydate = datetime.datetime.strptime(dateString,
                                        "%Y-%m-%d") + datetime.timedelta(
        days=age)
    futureDateString = mydate.strftime('%Y-%m-%d')
    return futureDateString


class linfit:
    """
    Todd Hunter
    """

    def computeMultivariateNormalStats(self, slope, intercept, covar, x,
                                       trials=10000,
                                       logfit=False, meanOfLogX=0):
        y = []
        # If the slope and intercept are correlated, then this will pick correlated
        # random variables, which reduces the spread of the result.
        mymean = [intercept, slope]
        random_intercept, random_slope = np.random.multivariate_normal(mymean,
                                                                       covar,
                                                                       trials).T
        if logfit:
            for t in range(trials):
                y.append(10 ** (
                random_slope[t] * (np.log10(x) - meanOfLogX) + random_intercept[
                    t]))
        else:
            for t in range(trials):
                y.append(10 ** (random_slope[t] * x + random_intercept[t]))
        rms = np.std(y, axis=0)
        ymin = np.min(y, axis=0)
        ymax = np.max(y, axis=0)
        ymean = np.mean(y, axis=0)
        return rms, ymin, ymax, ymean

    def computeStdDevMonteCarlo(self, slope, slopeSigma, intercept,
                                interceptSigma, x, trials=10000,
                                logfit=False, meanOfLogX=0, covar=None,
                                silent=False, returnMedian=False):
        """
        Computes Y and its uncertainty given a fitted slope and intercept, uncertainies on them,
        and a value of X.  Returns the uncertainty.
        """
        y = []
        if covar is not None:
            # If the slope and intercept are correlated, then this will pick correlated
            # random variables, which reduces the spread of the result.
            mymean = [intercept, slope]
            random_intercept, random_slope = np.random.multivariate_normal(
                mymean, covar, trials).T
        if logfit:
            for t in range(trials):
                if covar is None:
                    y.append(10 ** ((slope + slopeSigma * pickRandomError()) * (
                    np.log10(
                        x) - meanOfLogX) + intercept + interceptSigma * pickRandomError()))
                else:
                    y.append(10 ** (
                    random_slope[t] * (np.log10(x) - meanOfLogX) +
                    random_intercept[t]))
        else:
            for t in range(trials):
                if covar is None:
                    y.append((
                             slope + slopeSigma * pickRandomError()) * x + intercept + interceptSigma * pickRandomError())
                else:
                    y.append(10 ** (random_slope[t] * x + random_intercept[t]))
        y = np.array(y)
        mad = MAD(y)
        median = np.median(y)
        # Clip out wild values beyond 6*sigma (using the MAD) when computing std
        rms = np.std(y[np.where(np.abs(y) < median + 6 * mad)])
        if not silent:
            print "Median Monte-Carlo result for %f = %f +- %f (scaled MAD = %f)" % (
            x, median, rms, mad)
        if returnMedian:
            return rms, median
        else:
            return rms

    def readFluxscaleResult(self, xfilename, yfilename, source, verbose=False,
                            maxpoints=0, spwlist=[], referenceFrame='TOPO',
                            debug=False):
        """
      Specific function to read CASA output files from listobs and fluxscale.
      It returns the log10() of the frequency and flux densities read. The
      flux density uncertainties returned are scaled by the flux density values
      (but the ratios are not logged).
      """
        fx = open(xfilename, 'r')
        fy = open(yfilename, 'r')
        lines = fy.readlines()
        fy.close()
        x = []
        y = []
        yerror = []
        skiplist = ''
        trueSource = ''
        sourcesFound = []
        ignoreSpw = []  # This will be a list of spws with "Insufficient data"
        for line in lines:
            if (line.find('Flux') >= 0 and (
                    source == '' or line.find(source) >= 0)):
                tokens = line.split()
                for t in range(len(tokens)):
                    if tokens[t].find('SpW') >= 0:
                        #  We have specified the spws to include
                        spw = int(tokens[t].split('=')[1])
                        if line.find('INSUFFICIENT DATA') >= 0:
                            ignoreSpw.append(spw)
                            print "Skipping spw %d due to INSUFFICIENT DATA result from fluxscale." % (
                            spw)
                            break
                    if tokens[t] == 'is:':
                        if spw not in spwlist and spwlist != []:
                            skiplist += str(spw) + ','
                        else:
                            # print "Using spw %d" % spw
                            if line.find('Hz) is') > 0:
                                moreTokens = 2  # casa >= 4.3
                            else:
                                moreTokens = 0  # casa <= 4.2
                            if tokens[t - 4 - moreTokens] != 'for':
                                # there is a blank in the name
                                trueSource = tokens[t - 4 - moreTokens] + ' ' + \
                                             tokens[t - 3 - moreTokens]
                            else:
                                trueSource = tokens[
                                    t - 3 - moreTokens]  # no blanks were found in the name
                            sourcesFound.append(trueSource)
                            if debug:
                                print "parsing flux from %s" % (tokens[t + 1])
                            y.append(float(tokens[t + 1]))
                            if debug:
                                print "parsing flux uncertainty from %s" % (
                                tokens[t + 3])
                            yerror.append(float(tokens[t + 3]))
                            break
                if len(y) == maxpoints and maxpoints > 0: break
        if len(y) == 0:
            if trueSource == '':
                print "Did not find any flux densities for source = %s" % (
                source)
            else:
                print "Did not find any flux densities for source = %s" % (
                trueSource)
            return [], [], [], [], []
        if len(skiplist) > 0:
            print "Skipping spw ", skiplist[0:-1]
        lines = fx.readlines()
        fx.close()
        bw = []
        freqUnits = ''
        spwsKept = []
        for line in lines:
            loc = line.find('Ch1(')
            if loc >= 0:
                freqUnits = line[loc + 4:loc + 7]
                if verbose:
                    print "Read frequency units = ", freqUnits
            if line.find(referenceFrame) >= 0:
                tokens = line.split()
                for t in range(len(tokens)):
                    if tokens[t] == referenceFrame:
                        try:
                            spw = int(tokens[t - 2])
                        except:
                            spw = int(tokens[t - 3])

                        if ((
                                    spw not in spwlist and spwlist != []) or spw in ignoreSpw):
                            #                            print "Skipping spw %d" % spw
                            continue
                        else:
                            spwsKept.append(spw)
                            bw.append(float(tokens[t + 3]))
                            x.append(
                                float(tokens[t + 1]) + bw[-1] * 0.5 * 0.001)
                        break
                    if len(x) == maxpoints and maxpoints > 0:
                        print "Stopping after readings %d points ---------------" % (
                        maxpoints)
                        break
        if len(y) != len(x):
            print "There is a mismatch between the number of spws in the %s frame (%d" % (
            referenceFrame, len(x))
            print "and the number of valid flux densities (%d)." % (len(y))
            if len(np.unique(sourcesFound)) > 1:
                print "Please limit the fit to one source using the 'source' parameter: ", np.unique(
                    sourcesFound)
            return [], [], [], [], [], []
        if verbose:
            print "Read %d x values for %s = " % (len(x), trueSource), x
            print "Read %d y values for %s = " % (len(y), trueSource), y
        else:
            print "Read %d x values for %s (avg=%.3f) = " % (
            len(x), trueSource, np.mean(x)), x
            print "Read %d y values for %s (avg=%.3f) = " % (
            len(y), trueSource, np.mean(y)), y
        if freqUnits.find('MHz') >= 0:
            x = np.array(x) * 0.001
            freqUnits = 'GHz'
        logx = np.log10(x)
        logyerror = list(np.array(yerror) / np.array(y))
        logy = np.log10(y)
        return logx, logy, logyerror, trueSource, freqUnits, spwsKept
        # end of readFluxscaleResult

    def parse_spw_argument(self, spw):
        """
      # returns an integer list of spws based on a string entered by the user
      #  e.g.  '1~7,9~15' yields [1,2,3,4,5,6,7,9,10,11,12,13,14,15]
      #        [1,2,3] yields [1,2,3]
      #        1 yields [1]
      """
        if type(spw) == list:
            return spw
        elif type(spw) == int:
            return [spw]
        sublists = spw.split(',')
        spwlist = []
        for s in sublists:
            spwrange = s.split('~')
            if len(spwrange) > 1:
                try:
                    firstSpw = int(spwrange[0])
                    try:
                        secondSpw = int(spwrange[1])
                        for w in range(int(spwrange[0]), int(spwrange[1]) + 1):
                            spwlist.append(w)
                    except:
                        print "Unrecognized element in spw string: %s" % (
                        spwrange[1])
                except:
                    print "Unrecognized element in spw string: %s" % (
                    spwrange[0])
            else:
                try:
                    spwlist.append(int(spwrange[0]))
                except:
                    print "Unrecognized element in spw string: %s" % (
                    spwrange[0])
        return spwlist

    def spectralIndexFitterMonteCarlo(self, filename, yfilename='', degree=1,
                                      source='',
                                      verbose=False, maxpoints=0, trials=2000,
                                      spw='',
                                      referenceFrame='TOPO', plaintext=False,
                                      lineNumbers=None, columns=None,
                                      yscale=1.0,
                                      freqs=[], fluxes=[], errors=[]):
        """
      If error bars are included with the data, then the return values are:
          p, covar contain the results of the weighted fit
          p2 contains the results of the unweighted fit
      If no error bars are included, then
          p contains the results of the unweighted fit
      lineNumbers: which lines to read, where 1 is the first line in the file
      columns: which columns to read for freq, flux, error (starting at 1)
      """
        spwlist = []
        nullReturn = ([], [], [], [], [], [], [], [], [], [], [], [])
        if spw != [] and spw != '':
            spwlist = self.parse_spw_argument(spw)
            if len(spwlist) < 2:
                print "I need more than 1 spws to fit a spectral index."
                return nullReturn
            else:
                print "Will use spw = ", spwlist
        if filename == '':
            x = np.array(freqs, dtype=np.float)
            y = np.array(fluxes, dtype=np.float)
            yerror = np.array(errors, dtype=np.float)
            if len(x) != len(y) or len(x) != len(yerror):
                print "Inconsistent array lengths of freqs, fluxes, errors"
                return nullReturn
            freqUnits = 'GHz'
            spwsKept = 0
            if type(x[0]) == str:
                newx = []
                for xf in x:
                    newx.append(parseFrequencyArgument(xf) * 1e-9)
                x = newx
            elif x[0] > 2000:
                x *= 1e-9
            yerror = list(yerror / y)
            x = np.log10(x)
            y = np.log10(y)
        elif plaintext == False:
            (x, y, yerror, source, freqUnits,
             spwsKept) = self.readFluxscaleResult(filename, yfilename, source,
                                                  verbose, maxpoints, spwlist,
                                                  referenceFrame=referenceFrame)
        else:
            result = self.readPlaintextFile(filename, lineNumbers, columns)
            if result is None: return nullReturn
            (x, y, yerror) = result
            source = 'source'
            freqUnits = 'GHz'
            spwsKept = 0
        y = np.array(y) * yscale
        if len(x) == 0 or len(y) == 0:
            return nullReturn

        # normalize the X values about 0 such that the uncertainties between slope
        # and intercept will be uncorrelated
        meanOfLogX = np.mean(x)
        x = x - meanOfLogX

        if len(yerror) > 0:
            fitfunc = lambda p, x: p[0] * x + p[1]
            errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
            pinit = [1.0, -1.0]
            slope = []
            yoffset = []
            errors = np.zeros(len(y))
            # first trial is with no errors
            for t in range(trials):
                ytrial = list(np.array(y) + np.array(yerror) * np.array(errors))
                out = optimize.leastsq(errfunc, pinit, args=(x, ytrial, yerror),
                                       full_output=1)
                p = out[0]
                if t == 0:
                    covar = out[1]
                slope.append(p[0])
                yoffset.append(p[1])
                errors = []
                for e in range(len(y)):
                    errors.append(pickRandomError())
            p = [slope[0], yoffset[0]]
            pmedian = [np.median(slope), np.median(yoffset)]
            perror = [np.std(slope), np.std(yoffset)]
            p2 = np.polyfit(x, y, degree)
        else:
            p = np.polyfit(x, y, degree)
            covar = []
            p2 = []
            perror = []
            pmedian = []
        # if you change the number of return values, be sure to change the definition of nullReturn above
        return (p, covar, x, y, yerror, p2, source, freqUnits, pmedian, perror,
                spwsKept, meanOfLogX)

    # end of spectralIndexFitterMonteCarlo

    def readPlaintextFile(self, filename, lineNumbers=None, columns=None):
        """
        Reads a file of the format:   Freq(GHz)   Flux   Error
        and returns the 3 columns as three arrays.
        lineNumbers: which lines to read, where 1 is the first line in the file
        columns: which columns to read for freq, flux, error (starting at 1)
        """
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        freq = []
        flux = []
        error = []
        linectr = 0
        for line in lines:
            linectr += 1
            if line[0] != '#':
                if lineNumbers is None or linectr in lineNumbers:
                    if columns is None:
                        a, b, c = line.split()
                        freq.append(float(a))
                        flux.append(float(b))
                        error.append(float(c))
                    else:
                        mycolumns = line.split()
                        n_mycolumns = len(mycolumns)
                        if len(columns) < 3:
                            print "invalid column list: need to specify at least 3"
                            return
                        if np.max(columns) >= n_mycolumns:
                            print "invalid column list: %d exceeds %d" % (
                            np.max(columns), n_mycolumns)
                            return
                            #                      print "parsing %s" % (mycolumns[columns[0]-1])
                        freq.append(float(mycolumns[columns[0] - 1]))
                        flux.append(float(mycolumns[columns[1] - 1]))
                        error.append(float(mycolumns[columns[2] - 1]))
                        #      print "Read %d lines" % (len(freq))
        logx = np.log10(freq)
        logyerror = list(np.array(error) / np.array(flux))
        logy = np.log10(flux)
        return logx, logy, logyerror

    def spectralIndexFitterCovarMatrix(self, filename='', yfilename='',
                                       degree=1, source='',
                                       verbose=False, maxpoints=0, spw='',
                                       referenceFrame='TOPO', plaintext=False,
                                       lineNumbers=None, columns=None,
                                       yscale=1.0,
                                       freqs=[], fluxes=[], errors=[]):
        """
      If error bars are included with the data, then the return values are:
          p, covar contain the results of the weighted fit
          p2 contains the results of the unweighted fit
      If no error bars are included, then
          p contains the results of the unweighted fit
      lineNumbers: which lines to read, where 1 is the first line in the file
      columns: which columns to read for freq, flux, error (starting at 1)
      """
        spwlist = []
        nullReturn = ([], [], [], [], [], [], [], [], [], [], [])
        if spw != [] and spw != '':
            spwlist = self.parse_spw_argument(spw)
            if len(spwlist) < 2:
                print "I need more than 1 spws to fit a spectral index."
                return nullReturn
            else:
                print "Will use spw = ", spwlist
        if filename == '':
            x = np.array(freqs, dtype=float)
            y = np.array(fluxes, dtype=float)
            yerror = np.array(errors, dtype=float)
            if len(x) != len(y) or len(x) != len(yerror):
                print "Inconsistent array lengths of freqs, fluxes, errors"
                return
            freqUnits = 'GHz'
            spwsKept = 0
            if type(x[0]) == str:
                newx = []
                for xf in x:
                    newx.append(parseFrequencyArgument(xf) * 1e-9)
                x = newx
            elif x[0] > 2000:
                x *= 1e-9
            yerror = list(yerror / y)
            x = np.log10(x)
            y = np.log10(y)
        elif plaintext == False:
            # yerror returned by readFluxscaleResult = y_uncertainty / y
            (x, y, yerror, source, freqUnits,
             spwsKept) = self.readFluxscaleResult(filename, yfilename, source,
                                                  verbose, maxpoints, spwlist,
                                                  referenceFrame=referenceFrame)
        else:
            result = self.readPlaintextFile(filename, lineNumbers, columns)
            if result is None: return
            (x, y, yerror) = result
            source = 'source'
            freqUnits = 'GHz'
            spwsKept = 0
        if len(x) == 0 or len(y) == 0:
            return nullReturn
        y = np.log10(np.array(10 ** y) * yscale)

        # normalize the X values about 0 such that the uncertainties between slope
        # and intercept will be uncorrelated
        meanOfLogX = np.mean(x)
        x = x - meanOfLogX

        # simple fit with no errors
        p2 = np.polyfit(x, y, degree)

        x = 10 ** x
        y = 10 ** y

        if len(yerror) > 0:
            fitfunc = lambda p, x: p[0] + p[1] * np.log10(x)
            errfunc = lambda p, x, y, err: (np.log10(y) - fitfunc(p, x)) / err
            pinit = [1.0, -1.0]
            out = optimize.leastsq(errfunc, pinit, args=(x, y, yerror),
                                   full_output=1)
            p = out[0]
            covar = out[1]
        else:
            print "No uncertainties were found for the flux density measurements."

        return (
        p, x, y, yerror, covar, source, freqUnits, spwsKept, p2, meanOfLogX)

    # end of spectralIndexFitterCovarMatrix

    def round_to_n(self, x, n):
        '''
        http://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python
        '''
        return np.round(x, -int(floor(log10(abs(x)))) + (n - 1))

    def calc_ticks(self, domain, tick_count, equidistant):
        if equidistant:
            ticks = np.logspace(np.log10(domain[0]), np.log10(domain[1]),
                                num=tick_count, base=10)
        else:
            ticks = np.linspace(domain[0], domain[1], num=tick_count)
        for n in range(1, 6):
            if len(set(
                    self.round_to_n(tick, n) for tick in ticks)) == tick_count:
                break
        return list(self.round_to_n(tick, n) for tick in ticks)

    def spectralIndexCovarMatrix(self, filename='', yfilename='', source='',
                                 verbose=False,
                                 maxpoints=0, spw='', plotdir='',
                                 labelspw=False, referenceFrame='TOPO',
                                 plaintext=False,
                                 lineNumbers=None, columns=None, yscale=1.0,
                                 plotunits='linear',
                                 freqs=[], fluxes=[], errors=[], plotfile='',
                                 showplot=True,
                                 xaxis=[], dashedCurvePoints=40, silent=False,
                                 axisEqual=True):
        """
      This function is designed to fit the spectral index to the results
      output from casa's fluxscale task, and make a summary plot. Currently,
      it requires the output from the casa listobs task to determine the center
      frequencies of each spectral window.  It calls scipy.optimize.leastsq()
      and uses the fractional covariance matrix that it returns to determine the
      uncertainty on the fitted slope on the basis of the error bars on each flux
      density, as long as there are more than 2 points.  This is the formula that
      Bryan Butler sent me in Feb 2013. However, for the case of only 2 points,
      the denominator becomes zero in his formula, so it instead uses the
      sqrt(covariance matrix) which seems to agree well with the Monte-Carlo
      method.
      -- Todd Hunter

      filename: contains a listobs output file
      yfilename: contains a fluxscale output file
      source: sourcename to choose from the (possibly) multi-source fluxscale file
      maxpoints: the maximum number of spws to select for the fit (0=no max)
      spw: the spws to use, e.g. ''=all, '1~3,5,6~8'=[1,2,3,5,6,7,8]
      plotdir: the directory in which to write the plotfile
      labelspw: draw spw numeric labels adjacent to each point
      referenceFrame: the frequency reference frame, if not 'TOPO'
      plaintext: if True, then read the freqs, flux densities and uncertainties
                 from a single plain text file, specified by filename
      lineNumbers: which lines to read, where 1 is the first line in the file
      columns: which columns to read for freq, flux, error (starting at 1)
      freqs,fluxes,errors: alternative to using filename
      showplot: if True, produce a plot showing error bars and model
      xaxis: a list of points for which to compute model y values
      dashedCurvePoints: number of points to define the min/max/rms model lines
      axisEqual: if True, then when plotting, set axis('equal')
      """
        # x & y are returned in linear units of frequency (normalized about 1.0 MHz)
        # and flux density (usually Jy)
        # yerror is the (uncertainty in y) divided by y
        result = self.spectralIndexFitterCovarMatrix(filename, yfilename,
                                                     degree=1, source=source,
                                                     verbose=verbose,
                                                     maxpoints=maxpoints,
                                                     spw=spw,
                                                     referenceFrame=referenceFrame,
                                                     plaintext=plaintext,
                                                     lineNumbers=lineNumbers,
                                                     columns=columns,
                                                     yscale=yscale, freqs=freqs,
                                                     fluxes=fluxes,
                                                     errors=errors)
        if result is None: return
        if result[0] == []: return
        (p, x, y, yerror, covar, source, freqUnits, spwsKept, p2,
         meanOfLogX) = result
        if plotdir == '':
            plotdir = './'
        if plotdir[-1] != '/':
            plotdir += '/'
        if p == []:
            # If the data import failed, then the fit will be blank, so return.
            return p

        # Get the result of the fit
        yoffset = p[0]
        slope = p[1]
        freq = x[0]
        freqUnNormalized = 10 ** (np.log10(freq) + meanOfLogX)
        amp = 10 ** (yoffset + np.log10(freq) * slope)

        # from Bryan Butler email on Feb 15/16, 2013
        summed_error = 0
        for ii in range(len(x)):
            model = yoffset + slope * x[ii]
            residual = (model - y[ii]) ** 2 / (yerror[ii] ** 2)
            summed_error += residual
        if len(x) > 2:
            residual_variance = summed_error / (len(x) - 2)
            slopeFitError = fabs(slope) / np.sqrt(
                covar[1][1] * residual_variance)
            yoffsetError = fabs(yoffset) / np.sqrt(
                covar[0][0] * residual_variance)
        else:
            # residual variance cannot be calculated due to divide by zero
            # but the following formula seem to agree with the Monte Carlo approach
            try:
                slopeFitError = fabs(slope) * np.sqrt(covar[1][1])
            except TypeError:
                slopeFitError = 10.
                mydict = {'spectralIndex': -0.6,
                          'spectralIndexUncertainty': 0.5,
                          'intercept': 0, 'interceptUncertainty': 0,
                          'meanOfLogX': 0, 'yaxis': [], 'yaxisUncertainty': [],
                          'covar': covar}

            yoffsetError = fabs(yoffset) * np.sqrt(covar[0][0])
        ampError = self.computeStdDevMonteCarlo(slope, slopeFitError, yoffset,
                                                yoffsetError, freqUnNormalized,
                                                logfit=True,
                                                meanOfLogX=meanOfLogX,
                                                covar=covar, silent=silent)
        mydict = {'spectralIndex': slope,
                  'spectralIndexUncertainty': slopeFitError,
                  'intercept': yoffset, 'interceptUncertainty': yoffsetError,
                  'meanOfLogX': meanOfLogX, 'yaxis': [], 'yaxisUncertainty': [],
                  'covar': covar}

        amp2 = (10.0 ** p2[1]) * (freq ** p2[0])
        if verbose:
            print "yoffset=%f, covar = " % yoffset, covar
        if not silent:
            if amp < 0.01:
                print "Error-weighted fit: Slope: %.3f+-%.3f  Flux D. @ %.3f%s: %.3f+-%.3f mJy" % (
                slope,
                slopeFitError, freqUnNormalized, freqUnits, amp * 1000,
                ampError * 1000)
                print "   Un-weighted fit: Slope: %.3f         Flux D. @ %.3f%s: %.3f mJy" % (
                p2[0], freqUnNormalized, freqUnits, amp2 * 1000)
            else:
                print "Error-weighted fit: Slope: %.3f+-%.3f  Flux D. @ %.3f%s: %.3f+-%.3f Jy" % (
                slope,
                slopeFitError, freqUnNormalized, freqUnits, amp, ampError)
                print "   Un-weighted fit: Slope: %.3f         Flux D. @ %.3f%s: %.3f Jy" % (
                p2[0], freqUnNormalized, freqUnits, amp2)

        if showplot:
            pb.clf()
            desc = pb.subplot(111)
        logx = np.log10(x)
        logy = np.log10(y)
        logxsorted = np.sort(logx)
        yfit = yoffset + logx * slope
        for myx in xaxis:
            yaxis = 10 ** (yoffset + (np.log10(myx) - meanOfLogX) * slope)
            mydict['yaxis'].append(yaxis)
            yaxisUncertainty = self.computeStdDevMonteCarlo(slope,
                                                            slopeFitError,
                                                            yoffset,
                                                            yoffsetError, myx,
                                                            logfit=True,
                                                            meanOfLogX=meanOfLogX,
                                                            covar=covar,
                                                            silent=silent)
            mydict['yaxisUncertainty'].append(yaxisUncertainty)
        dashedCurveX = 10 ** np.linspace(np.log10(np.min(x)),
                                         np.log10(np.max(x)), dashedCurvePoints)
        dashedCurveRms, dashedCurveMin, dashedCurveMax, dashedCurveMean = \
            self.computeMultivariateNormalStats(slope, yoffset, covar,
                                                dashedCurveX, trials=10000,
                                                logfit=True, meanOfLogX=0)
        minvalue = np.log10(dashedCurveMin)
        maxvalue = np.log10(dashedCurveMax)
        logxsorted = np.log10(dashedCurveX)

        if showplot:
            if plotunits == 'linear':
                if verbose:
                    print "Using linear plot units"
                lower = y * yerror
                upper = y * yerror
                x = 10 ** (
                logx + meanOfLogX)  # This x is in linear units of GHz
                pb.errorbar(x, y, yerr=[lower, upper], color='k', ls='None')
                pb.loglog(x, y, 'ko',
                          x, 10 ** yfit, 'r-',
                          10 ** (logxsorted + meanOfLogX),
                          dashedCurveMean - dashedCurveRms, 'r:',
                          10 ** (logxsorted + meanOfLogX),
                          dashedCurveMean + dashedCurveRms, 'r:'
                          )
                desc.set_xscale("log", subsx=range(10))
                desc.set_yscale("log", subsy=range(10))
                if freqUnits != '':
                    pb.xlabel('Frequency (%s)' % freqUnits)
                else:
                    pb.xlabel('Frequency')
                pb.ylabel('Flux Density (Jy)')
                x0 = np.min(x)
                x1 = np.max(x)
                xr = x1 - x0
                y0 = np.min(y - lower)
                y1 = np.max(y + upper)
                yr = y1 - y0
                pb.xlim([x0 - xr * 0.1, x1 + xr * 0.1])
                pb.ylim([y0 - yr * 0.1, y1 + yr * 0.1])
                if labelspw:
                    for i in range(len(x)):
                        pb.text(x[i], y[i], '%d' % (spwsKept[i]), size=10)
            else:
                logx += meanOfLogX
                logxsorted += meanOfLogX
                pb.plot(logx, logy, 'ko', logx, yfit, 'r-',
                        logxsorted, minvalue, 'r--', logxsorted, maxvalue,
                        'r--')
                lower = np.log10(y) - np.log10(y - y * yerror)
                upper = np.log10(y + y * yerror) - np.log10(y)
                pb.errorbar(logx, logy, yerr=[lower, upper], color='k',
                            ls='None')
                if freqUnits != '':
                    pb.xlabel('Log10(Frequency(%s))' % freqUnits)
                else:
                    pb.xlabel('Log10(Frequency)')
                pb.ylabel('Log10(FluxDensity(Jy))')
                if axisEqual:
                    pb.axis('equal')
                if labelspw:
                    for i in range(len(x)):
                        pb.text(logx[i], logy[i], '%d' % (spwsKept[i]), size=10)
            if amp < 0.01:
                pb.title(
                    source + '  spectral index=%+.3f+-%.3f, F(%.2f%s)=%.3f+-%.3f mJy' % (
                    slope, slopeFitError, freqUnNormalized, freqUnits,
                    amp * 1000, ampError * 1000), size=14)
            else:
                pb.title(
                    source + '  spectral index=%+.3f+-%.3f, F(%.2f%s)=%.3f+-%.3f Jy' % (
                    slope, slopeFitError, freqUnNormalized, freqUnits, amp,
                    ampError), size=14)
            desc.xaxis.grid(which='major')
            desc.yaxis.grid(which='major')
            desc.xaxis.set_major_formatter(
                matplotlib.ticker.ScalarFormatter(useOffset=False))
            desc.yaxis.set_major_formatter(
                matplotlib.ticker.ScalarFormatter(useOffset=False))
            desc.set_xticks(
                self.calc_ticks(pb.xlim(), tick_count=4, equidistant=True))
            desc.set_yticks(
                self.calc_ticks(pb.ylim(), tick_count=4, equidistant=True))
            if plaintext:
                yfilename = filename
            if yfilename != '':
                ytext = yfilename
                #              ytext = self.replaceUnderscores(ytext)
                pb.text(0.01, 0.01, 'data from %s' % ytext,
                        transform=desc.transAxes, size=10)
            if plotfile == True:
                if plaintext:
                    pngname = '%s%s.png' % (plotdir, yfilename)
                else:
                    pngname = '%s%s.%s.png' % (plotdir, yfilename, source)
            elif plotfile != '':
                pngname = plotfile
            if plotfile != '':
                if os.access('./', os.W_OK) == False:
                    fileWriteMessage = "Cannot write plot to current directory, therefore will write to /tmp."
                    pngname = '/tmp/' + os.path.basename(pngname)
                pb.savefig(pngname)
                print "Plot saved in %s" % pngname
            pb.draw()
        return mydict
        # end of spectralIndexCovarMatrix()

    def spectralindex(self, filename='', yfilename='', source='', verbose=False,
                      maxpoints=0, trials=0, spw='', plotdir='',
                      labelspw=False, referenceFrame='TOPO', plaintext=False,
                      lineNumbers=None, columns=None, yscale=1.0,
                      plotunits='linear',
                      freqs=[], fluxes=[], errors=[], plotfile='',
                      showplot=True, xaxis=[],
                      silent=False, axisEqual=True):
        """
      This function is designed to fit the spectral index to the results
      output from casa's fluxscale task. Currently, it requires the output
      from the listobs task to determine the center frequencies of each spw.
      Usage: spectralIndex(filename='',yfilename='',source='',verbose=False,
                           maxpoints=0,trials=0,spw='',plotdir='',plaintext=False)
      filename: contains a listobs output file
      yfilename: contains a fluxscale output file
      source: sourcename to choose from the (possibly) multi-source fluxscale file
      maxpoints: the maximum number of spws to select for the fit (0=no max)
      trials: if > 0, use a Monte-Carlo technique estimate the fit uncertainties,
         otherwise, use the sqrt(covar_matrix) from scipy.optimize.leastsq (default).
         There is a minimum number of 100 trials, and ~1000 is recommended.
      spw: the spws to use, e.g. ''=all, '1~3,5,6~8'=[1,2,3,5,6,7,8]
      plotdir: the directory in which to write the plotfile
      labelspw: draw spw numeric labels adjacent to each point
      referenceFrame: the frequency reference frame, if not 'TOPO'
      plaintext: if True, then read the freqs, flux densities and uncertainties
                 from a single plain text file, specified by filename
      lineNumbers: which lines to read, where 1 is the first line in the file
      columns: which columns to read for freq, flux, error (starting at 1)
      yscale: factor to scale the y values after reading from the file
      plotunits: 'linear' or 'log' (only used if trials==0)
      freqs,fluxes,errors: alternative to using filename
      xaxis: a list of points for which to compute model y values
      axisEqual: if True, then when plotting, set axis('equal')
      """
        if (plaintext == False and yfilename == '' and (
                    freqs == [] or fluxes == [] or errors == [])):
            print "When plaintext=False, you must also specify yfilename which is the output file from fluxscale."
            return
        if plotunits not in ['linear', 'log']:
            print "plotunits must be either 'linear' or 'log'"
            return
        if len(errors) == 2:
            # report the exact fit, and 1-sigma extrema fits
            si = np.log(fluxes[1] / float(fluxes[0])) / np.log(
                freqs[1] / float(freqs[0]))
            try:
                si1 = np.log(
                    (fluxes[1] + errors[1]) / (fluxes[0] - errors[0])) / np.log(
                    freqs[1] / float(freqs[0]))
            except ZeroDivisionError:
                print("ZeroDiv", freqs[1], float(freqs[0]), fluxes[0], errors[0],
                      np.log(freqs[1] / float(freqs[0])))
                si1 = np.log(
                    (fluxes[1] + errors[1]) / (fluxes[0])) / np.log(
                    freqs[1] / float(freqs[0]))
            si2 = np.log(
                (fluxes[1] - errors[1]) / (fluxes[0] + errors[0])) / np.log(
                freqs[1] / float(freqs[0]))
            if not silent:
                print "exact value: %f,  1-sigma extrema: %f, %f,  mean unc=%f" % (
                si, si1, si2, 0.5 * abs(si2 - si1))
        if trials > 0:
            if trials < 100:
                trials = 100
            return (
            self.spectralIndexMonteCarlo(filename, yfilename, source, verbose,
                                         maxpoints,
                                         trials, spw, plotdir, labelspw,
                                         referenceFrame, plaintext,
                                         lineNumbers, columns, yscale,
                                         plotunits, freqs, fluxes,
                                         errors, plotfile, showplot, xaxis,
                                         silent, axisEqual))
        else:
            return (
            self.spectralIndexCovarMatrix(filename, yfilename, source, verbose,
                                          maxpoints,
                                          spw, plotdir, labelspw,
                                          referenceFrame,
                                          plaintext, lineNumbers, columns,
                                          yscale, plotunits,
                                          freqs, fluxes, errors, plotfile,
                                          showplot, xaxis,
                                          silent, axisEqual))

    def spectralIndexMonteCarlo(self, filename='', yfilename='', source='',
                                verbose=False,
                                maxpoints=0, trials=2000, spw='', plotdir='',
                                labelspw=False, referenceFrame='TOPO',
                                plaintext=False, lineNumbers=None, columns=None,
                                yscale=1.0,
                                plotunits='linear', freqs=[], fluxes=[],
                                errors=[], plotfile='',
                                showplot=True, xaxis=[], silent=False,
                                axisEqual=True):
        """
        This function is designed to fit the spectral index to the results
        output from casa's fluxscale task. Currently, it requires the output
        from the listobs task to determine the center frequencies of each
        spectral window.  It runs a brief Monte-Carlo series of fits to
        determine the uncertainty on the fitted slope on the basis of the
        error bars on each flux density.    -- Todd Hunter

        filename: contains a listobs output file
        yfilename: contains a fluxscale output file
        source: sourcename to choose from the (possibly) multi-source fluxscale file
        maxpoints: the maximum number of spws to select for the fit (0=no max)
        trials: number of Monte-Carlo fits to run to estimate the fit uncertainties
        spw: the spws to use, e.g. ''=all, '1~3,5,6~8'=[1,2,3,5,6,7,8]
        plotdir: the directory in which to write the plotfile
        labelspw: draw spw numeric labels adjacent to each point
        referenceFrame: the frequency reference frame, if not 'TOPO'
        lineNumbers: which lines to read, where 1 is the first line in the file
        columns: which columns to read for freq, flux, error (starting at 1)
        xaxis: a list of points for which to compute model y values
        axisEqual: if True, then when plotting, set axis('equal')
        """
        (p, covar, x, y, yerror, p2, source, freqUnits, pmedian, perror,
         spwsKept, meanOfLogX) = \
            self.spectralIndexFitterMonteCarlo(filename, yfilename, degree=1,
                                               source=source,
                                               verbose=verbose,
                                               maxpoints=maxpoints,
                                               trials=trials, spw=spw,
                                               referenceFrame=referenceFrame,
                                               plaintext=plaintext,
                                               lineNumbers=lineNumbers,
                                               columns=columns, yscale=yscale,
                                               freqs=freqs, fluxes=fluxes,
                                               errors=errors)
        if plotdir == '':
            plotdir = './'
        if plotdir[-1] != '/':
            plotdir += '/'
        if p == []:
            return p
        if p2 != []:
            # Then we have two solutions, where the first one is the error-weighted fit.
            slope = p[0]
            yoffset = p[1]
            if verbose:
                print "Completed %d Monte-Carlo error trials." % trials
            slopeErr = perror[0]
            yoffsetErr = perror[1]
            freq = 10 ** x[0]
            amp = 10 ** (p[1] + x[0] * slope)
            ampError = amp * perror[1]
            mydict = {'spectralIndex': slope,
                      'spectralIndexUncertainty': slopeErr,
                      'intercept': yoffset, 'interceptUncertainty': yoffsetErr,
                      'meanOfLogX': meanOfLogX, 'yaxis': [],
                      'yaxisUncertainty': []}
            freqUnNormalized = 10 ** (x[0] + meanOfLogX)
            if not silent:
                print "Error-weighted fit: Slope: %.3f+-%.3f  Flux D. @ %.3f%s: %.3f+-%.3f" % (
                slope, slopeErr, freqUnNormalized, freqUnits, amp, ampError)
            if verbose:
                yfit = x * pmedian[0] + pmedian[1]
                print "Predicted values = ", yfit
        else:
            # We only have one solution, so put it into p2.
            p2 = p
        amp2 = 10.0 ** p2[1]
        amp2 *= freq ** p2[0]
        if not silent:
            print "   Un-weighted fit: Slope: %.3f         Flux D. @ %.3f%s: %.3f" % (
            p2[0], freqUnNormalized, freqUnits, amp2)
        yfit = x * p[0] + p[
            1]  # the fit result with zero Monte-Carlo errors added to the data
        for myx in xaxis:
            mydict['yaxis'].append(
                10 ** (p[1] + (np.log10(myx) - meanOfLogX) * p[0]))
            mydict['yaxisUncertainty'].append(
                self.computeStdDevMonteCarlo(slope, slopeErr, yoffset,
                                             yoffsetErr, myx,
                                             logfit=True, meanOfLogX=meanOfLogX,
                                             covar=covar, silent=silent))
        dashedCurveHigh = []
        dashedCurveLow = []
        dashedCurve = []
        xsorted = np.sort(x)
        for myx in xsorted:
            uncertainty = self.computeStdDevMonteCarlo(slope, slopeErr, yoffset,
                                                       yoffsetErr,
                                                       10 ** (myx + meanOfLogX),
                                                       logfit=True,
                                                       meanOfLogX=meanOfLogX,
                                                       covar=covar,
                                                       silent=silent)
            dashedCurveHigh.append(
                np.log10(10 ** (p[1] + myx * p[0]) + uncertainty))
            dashedCurveLow.append(
                np.log10(10 ** (p[1] + myx * p[0]) - uncertainty))
        if showplot:
            pb.clf()
            desc = pb.subplot(111)
            desc.xaxis.grid(which='major')
            desc.yaxis.grid(which='major')
            if pmedian != []:
                x += meanOfLogX
                xsorted += meanOfLogX
                pb.plot(x, y, 'ko', x, yfit, 'r-')
                pb.hold(True)
                pb.plot(xsorted, dashedCurveLow, 'r--',
                        xsorted, dashedCurveHigh, 'r--')
                pb.title(
                    source + '  spectral index=%+.3f+-%.3f, F(%.2f%s)=%.3f+-%.3f Jy' % (
                    slope, slopeErr,
                    freqUnNormalized, freqUnits, amp, ampError), size=14)
            else:
                pb.plot(x, y, 'ko', x, yfit, 'r-')
                pb.title(source + ' index=%.3f, F(%.2f%s)=%.3f' % (
                slope, freqUnNormalized, freqUnits, amp))
            if labelspw:
                for i in range(len(x)):
                    pb.text(x[i], y[i], '%d' % (spwsKept[i]), size=10)
            originalYerror = yerror * (10 ** y)
            lower = y - np.log10(10 ** y - originalYerror)
            upper = np.log10(10 ** y + originalYerror) - y
            pb.errorbar(x, y, yerr=[lower, upper], color='k', ls='None')
            if freqUnits != '':
                pb.xlabel('Log10(Frequency(%s))' % freqUnits)
            else:
                pb.xlabel('Log10(Frequency)')
            pb.ylabel('Log10(FluxDensity(Jy))')
            if axisEqual:
                pb.axis('equal')
            if plaintext:
                yfilename = filename
            if yfilename != '':
                ytext = yfilename
                #              ytext = self.replaceUnderscores(ytext)
                pb.text(0.01, 0.01, 'data from %s' % ytext,
                        transform=desc.transAxes, size=10)
            if plotfile == True:
                pngname = '%s%s.%s.png' % (plotdir, yfilename, source)
            elif plotfile != '':
                pngname = plotfile
            if plotfile != '':
                if os.access('./', os.W_OK) == False:
                    fileWriteMessage = "Cannot write plot to current directory, therefore will write to /tmp."
                    pngname = '/tmp/' + os.path.basename(pngname)
                pb.savefig(pngname)
                print "Plot saved in %s" % pngname
            pb.draw()
        return mydict

    # end of spectralIndexMonteCarlo

    def replaceUnderscores(self, y):
        """
      Replaced underscores with \_ to avoid strange effects when plotting text with pylab.text()
      """
        newy = ''
        for i in range(len(y)):
            if y[i] == '_':
                newy += '\_'
            else:
                newy += y[i]
        return newy

    def linfitFromFile(self, filename, xcol, ycol, yerrorcol=None, pinit=[0, 0],
                       plot=False, plotfile=None, xlabel=None, ylabel=None,
                       title=None, delimiter=None, residual=False):
        """
        Performs linear fit to data from the specified file, and data columns.
        xcol, ycol: columns to use from the file (starting at 0)
        yerrorcol: column to use for uncertainties (None -> 1% of ycol)
        pinit: contains the initial guess of [slope, intercept]
        plot: whether to generate a plot window
        plotfile: whether to generate a png file
        xlabel, ylabel, title: labels for the plot
        delimiter: column delimitier in your file
        residual: if True, show the residual of the fit in second panel
        Returns:
        the fitted coefficients (and residuals if requested)
        """
        x, y = getxyFromFile(filename, xcol, ycol, delimiter)
        if yerrorcol is None:
            yerror = y * 0.01
        elif type(yerrorcol) == np.float or type(yerrorcol) == float:
            yerror = y * yerrorcol
        else:
            yerror = y * 0.01
        p = self.linfit(x, y, yerror, pinit, plot, plotfile, xlabel, ylabel,
                        title, residual)
        return p

    def linfit(self, x, y, yerror, pinit=[0.57, 3.4], plot=False, plotfile=None,
               xlabel=None, ylabel=None, title=None, residual=False,
               excludeXrange=[], returnResiduals=False):
        """
        Basic linear function fitter with error bars in y-axis data points.
        Uses scipy.optimize.leastsq().  Accepts either lists or arrays.
        Example:
             lf = au.linfit()
             lf.linfit(x, y, yerror, pinit)
        Input:
             x, y: x and y axis values
             yerror: uncertainty in the y-axis values
             pinit contains the initial guess of [slope, intercept]
             residual: if True, show the residual of the fit in second panel
             excludeXrange: exclude a range of x-axis values when doing the fit
                 e.g. [109,116] to avoid 109-116 GHz points
             returnResiduals: if True, return [slope,intercept,data-fit]
        Output:
           The fit result as: [slope, y-intercept]
        - Todd Hunter
        """
        x = np.array(x)
        y = np.array(y)
        yerror = np.array(yerror)
        if excludeXrange != []:
            idx1 = np.where(x < excludeXrange[0])[0]
            idx2 = np.where(x > excludeXrange[1])[0]
            keep = np.union1d(idx1, idx2)
            x = x[keep]
            y = y[keep]
            yerror = yerror[keep]
        fitfunc = lambda p, x: p[1] + p[0] * x
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / (err ** 2)
        out = optimize.leastsq(errfunc, pinit, args=(x, y, yerror / y),
                               full_output=1)
        p = out[0]
        covar = out[1]
        residuals = y - (p[0] * x + p[1])
        if plot:
            pb.clf()
            if residual:
                desc = pb.subplot(211)
            else:
                desc = pb.subplot(111)
            pb.errorbar(x, y, yerr=yerror, fmt='o', color='b', ecolor='b')
            xline = np.array(pb.xlim())
            yline = p[0] * xline + p[1]
            pb.plot(xline, yline, 'k-')
            if title is None:
                pb.title('y = (%g)x %+g' % (p[0], p[1]))
            else:
                pb.title(title)
                pb.text(0.4, 0.92, 'y = (%g)x %+g' % (p[0], p[1]),
                        transform=desc.transAxes)
            if xlabel is not None:
                pb.xlabel(xlabel)
            if ylabel is not None:
                pb.ylabel(ylabel)
            if plotfile is None:
                plotfile = 'linfit.png'
            if residual:
                pb.subplot(212)
                pb.plot(x, residuals, 'b.')
                pb.ylabel('Residuals')
                if xlabel is not None:
                    pb.xlabel(xlabel)
            pb.savefig(plotfile)
            pb.draw()
        if returnResiduals:
            return [p[0], p[1], residuals]
        else:
            return p

    def odrFunction(self, B, x):
        return B[0] * x + B[1]

    def orthogonalDistanceRegression(self, x, y, xerror, yerror, loglog=False):
        """
        Performs a linear fit accounting for measurement uncertainties in both axes.
        loglog: if True, then take the log10 of the x and y vectors, and set
                xerror and yerror equal to xerror/x and yerror/y, respectively.
                Also, the intercept will be converted to the scaling coefficient.
        """
        if loglog:
            xerror = np.array(xerror) / np.array(x)
            yerror = np.array(yerror) / np.array(y)
            print "Taking log10 of x and y"
            x = np.log10(x)
            y = np.log10(y)
        mydata = scipy.odr.RealData(x, y, sx=xerror, sy=yerror)
        myodr = scipy.odr.ODR(mydata, scipy.odr.Model(self.odrFunction),
                              beta0=[1., 2.])
        myoutput = myodr.run()
        #        myoutput.pprint()
        beta = myoutput.beta
        sd_beta = myoutput.sd_beta
        if loglog:
            beta[1] = 10 ** (beta[1])  # convert intercept to scaling coefficent
            # treat sd_beta as a fractional standard deviation
            sd_beta[1] = 0.5 * (10 ** (beta[1] + sd_beta[1] * beta[1]) - 10 ** (
            beta[1] - sd_beta[1] * beta[1]))
        return beta, sd_beta


def pickRandomError():
    """
    Picks a random value from a Gaussian distribution with mean 0 and standard deviation = 1.
    """
    w = 1.0
    while w >= 1.0:
        x1 = 2.0 * random.random() - 1.0
        x2 = 2.0 * random.random() - 1.0
        w = x1 * x1 + x2 * x2

    w = np.sqrt((-2.0 * np.log(w)) / w)
    y1 = x1 * w
    y2 = x2 * w
    return y1


def MADoutliers(mylist, c=0.6745, sigma=3):
    """
    Find the outlier values, and their indices, in a list by computing
    its median and median absolute deviation (scaled to the expected rms).
    Returns: 2 lists
    * list of outlier values
    * list of their indices
    -Todd Hunter
    """
    mylist = np.array(mylist)
    mad = MAD(mylist, c)
    idx = np.where(abs(mylist - np.median(mylist)) > mad * sigma)[0]
    outliers = mylist[idx]
    return outliers, idx


def MAD(a, c=0.6745, axis=0):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """
    a = np.array(a)
    good = (a == a)
    a = np.asarray(a, np.float64)
    if a.ndim == 1:
        d = np.median(a[good])
        m = np.median(np.fabs(a[good] - d) / c)
    #        print  "mad = %f" % (m)
    else:
        d = np.median(a[good], axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = swapaxes(a[good], 0, axis)
        else:
            aswp = a[good]
        m = np.median(np.fabs(aswp - d) / c, axis=0)

    return m


def getxyFromFile(filename, xcol, ycol, delimiter=None, maxLines=None,
                  startAfter=None, stopAt=None, yhms=False, ydms=False,
                  verbose=False):
    """
    Gets two columns of numeric (or string) data from an ASCII file, as numpy arrays.
    xcol: the column number from which to read x-axis data
    ycol: the column number from which to read y-axis data
    If a colon is seen in the xcol & ycol data, and xcol!=ycol, then assume both
      columns are a sky position in sexagesimal format & convert them to radians.
    delimiter: the string that delimits columns within a row
               (default = None which means whitespace)
    startAfter: if defined, then skip all rows up to and including the first
        one containing this string
    stopAt: if defined, then skip all remaining rows when this string is seen
    yhms: if True, then read 3 columns starting at ycol as HH MM SS and convert to deg
    ydms: if True, then read 3 columns starting at ycol as DD MM SS and convert to deg
    Returns:
    Two lists of values.
    -Todd Hunter
    """
    f = open(filename, 'r')
    x = []
    y = []
    lines = f.readlines()
    for i, line in enumerate(lines):
        if startAfter is not None:
            if line.find(startAfter) >= 0:
                startAfter = None
                print "Starting at line ", i + 1
            #            print "Skipping line ", i
            continue
        if stopAt is not None and startAfter is None:
            if line.find(stopAt) >= 0:
                print "Stopping at line ", i
                break
        if len(line.strip()) == 0: continue
        if line.strip()[0] == '#' or line.strip()[0] == '/': continue
        if line.find('Source') >= 0: continue
        tokens = line.strip().split(delimiter)
        if len(tokens) < xcol or len(tokens) < ycol:
            print "Skipping row %d because there are only %d columns." % (
            i, len(tokens))
            continue
        # support either 0-based or 1-based column specification
        if len(tokens) == xcol:
            xcol -= 1
        if len(tokens) == ycol:
            ycol -= 1
        if yhms:
            if tokens[xcol].find(':') > 0:
                x.append(hmsToHours(tokens[xcol]))
            else:
                try:
                    x.append(float(tokens[xcol]))
                except:
                    x.append(tokens[xcol])  # allow strings
            y.append(15 * hmsToHours(
                tokens[ycol] + ':' + tokens[ycol + 1] + ':' + tokens[ycol + 2]))
        elif ydms:
            if tokens[xcol].find(':') > 0:
                x.append(hmsToHours(tokens[xcol]))
            else:
                try:
                    x.append(float(tokens[xcol]))
                except:
                    x.append(tokens[xcol])  # allow strings
            y.append(hmsToHours(
                tokens[ycol] + ':' + tokens[ycol + 1] + ':' + tokens[ycol + 2]))
        else:
            if (tokens[xcol].find(':') > 0 and tokens[ycol].find(
                    ':') > 0 and xcol != ycol):
                if verbose:
                    print "Running radec2rad('%s %s')" % (
                    tokens[xcol], tokens[ycol])
                rarad, decrad = radec2rad(tokens[xcol] + ' ' + tokens[ycol])
                x.append(rarad)
                y.append(decrad)
            else:
                if tokens[xcol].find(':') > 0:
                    x.append(hmsToHours(tokens[xcol]))
                else:
                    try:
                        x.append(float(tokens[xcol]))
                    except:
                        x.append(tokens[xcol])  # allow strings
                if tokens[ycol].find(':') > 0:
                    y.append(hmsToHours(tokens[ycol]))
                else:
                    try:
                        y.append(float(tokens[ycol]))
                    except:
                        y.append(tokens[ycol])  # allow strings
        if len(x) == maxLines:
            print "Stopping after reading %d data lines from %d file lines" % (
            maxLines, i + 1)
            break
    f.close()
    return np.array(x), np.array(y)


def hmsToHours(hms):
    """
    Converts string (or list of strings) of "HH", "HH:MM", "HH:MM:SS" or
    "HH:MM:SS.S" (with any number of trailing digits) into
    floating point hours.    See also hoursToHMS.
    Todd Hunter
    """
    if type(hms) == str:
        h = [hms]
    else:
        h = hms
    hours = []
    for x in h:
        hours.append(sum(float(t) * 60 ** (2 - i) for i, t in
                         enumerate(x.split(":"))) / 3600.)
    if type(hms) == str:
        hours = hours[0]
    return hours
