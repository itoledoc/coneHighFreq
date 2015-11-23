from pylab import *

bandDefinitions = {
    1: [31.3e9, 45e9],
    2: [67e9, 90e9],
    3: [84e9, 116e9],
    4: [125e9, 163e9],
    5: [163e9, 211e9],
    6: [211e9, 275e9],
    7: [275e9, 373e9],
    8: [385e9, 500e9],
    9: [602e9, 720e9],
    10: [787e9, 950e9]
}


def getBand(freq):
    """
    Converts a frequency into an ALMA band number.
    freq: can be given either as a floating point value in Hz, or a string
          with units at the end (GHz, MHz, kHz, or Hz).
    Todd Hunter

    Args:
        freq:
    """
    if type(freq) == str:
        freq = parseFrequencyArgument(freq)
    for band in bandDefinitions.keys():
        if ((freq <= bandDefinitions[band][1]) and (
                    freq >= bandDefinitions[band][0])):
            return band
    print "This frequency does not lie within any ALMA band."


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


def strDate2MJD(d):
    """
    Converts date in string format 20110809 or 2011x08x09 to MJD
    where 'x' can be any non-numeric character, like '-' or '/'

    Args:
        d:
    """
    #    print "strDate2MJD received: ", d
    hr = 0
    mn = 0
    sc = 0
    if ((d[4] == d[7]) and (d[4] < '0' or d[4] > '9') and (
                    d[7] < '0' or d[7] > '9')):
        # a delimiter is present
        year = d[0:4]
        month = d[5:7]
        day = d[8:10]
        if len(d) > 11:
            tokens = len(d[11:].split(':'))
            hr = d[11:].split(':')[0]
            if tokens > 1:
                mn = d[11:].split(':')[1]
            if tokens > 2:
                sc = d[11:].split(':')[2]
    else:
        year = d[0:4]
        month = d[4:6]
        day = d[6:8]
        if len(d) > 9:
            tokens = len(d[9:].split(':'))
            hr = d[9:].split(':')[0]
            if tokens > 1:
                mn = d[9:].split(':')[1]
            if tokens > 2:
                sc = d[9:].split(':')[2]
    year = int(year)
    month = int(month)
    day = int(day)
    hr = int(hr)
    mn = int(mn)
    sc = float(sc)
    mjd = ymdhmsToMJD(year, month, day, hr, mn, sc)
    return mjd


def ymdhmsToMJD(year, month, day, hour=0, minute=0, second=0.0):
    """
    converts 2010,12,31 to MJD on Dec 31, 2010 at UT=0
    converts 2010,12,31,9.50 to MJD on Dec 31, 2010 at UT=09:30
    converts 2010,12,31,9,0,5 to MJD on Dec 31, 2010 at UT=09:05
    required arguments: year, month, day
    optional arguments: hour, minute, second

    Args:
        second:
        minute:
        hour:
        day:
        month:
        year:
    """
    if month < 3:
        month += 12
        year -= 1
    a = floor(year / 100.)
    b = 2 - a + floor(a / 4.)
    UT = hour + minute / 60. + second / 3600.
    day += UT / 24.
    jd = floor(365.25 * (year + 4716)) + floor(
        30.6001 * (month + 1)) + day + b - 1524.5
    mjd = jdToMJD(jd)
    return mjd


def jdToMJD(JD):
    """
    Converts a JD value to MJD

    Args:
        JD:
    """
    MJD = JD - 2400000.5
    return MJD


def angularSeparation(ra0, dec0, ra1, dec1, returnComponents=False):
    """
  Computes the great circle angle_t between two celestial coordinates.
  using the Vincenty formula (from wikipedia) which is correct for all
  angles, as long as you use atan2() to handle a zero denominator.
     See  http://en.wikipedia.org/wiki/Great_circle_distance
  ra,dec must be given in degrees, as is the output.
  It also works for the az,el coordinate system.
  Component separations are field_0 minus field_1.
  See also angularSeparationRadians()
  -- Todd Hunter

    Args:
        returnComponents:
        dec1:
        ra1:
        dec0:
        ra0:
  """
    ra0 *= math.pi / 180.
    dec0 *= math.pi / 180.
    ra1 *= math.pi / 180.
    dec1 *= math.pi / 180.
    deltaLong = ra0 - ra1
    argument1 = (((math.cos(dec1) * math.sin(deltaLong)) ** 2) +
                 ((math.cos(dec0) * math.sin(dec1) - math.sin(dec0) * math.cos(
                     dec1) * math.cos(deltaLong)) ** 2)) ** 0.5
    argument2 = math.sin(dec0) * math.sin(dec1) + math.cos(dec0) * math.cos(
        dec1) * math.cos(deltaLong)
    angle_t = math.atan2(argument1, argument2) / (math.pi / 180.)
    if angle_t > 360:
        angle_t -= 360
    if returnComponents:
        cosdec = math.cos((dec1 + dec0) * 0.5)
        radegreesCosDec = np.degrees(ra0 - ra1) * cosdec
        radegrees = np.degrees(ra0 - ra1)
        decdegrees = np.degrees(dec0 - dec1)
        if radegrees > 360:
            radegrees -= 360
        if decdegrees > 360:
            decdegrees -= 360
        # positionAngle = -math.atan2(decdegrees*math.pi/180.,
        # radegreesCosDec*math.pi/180.)*180/math.pi
        retval = angle_t, radegrees, decdegrees, radegreesCosDec
    else:
        retval = angle_t
    return retval


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
        radecstring = radecstring.replace('h', ':').replace('m', ':').replace(
            'd', ':').replace('s', '')
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
        # noinspection PyBroadException
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


def angularSeparationRadians(ra0, dec0, ra1, dec1, returnComponents=False):
    """
  Computes the great circle angle between two celestial coordinates.
  using the Vincenty formula (from wikipedia) which is correct for all
  angles, as long as you use atan2() to handle a zero denominator.
     See  http://en.wikipedia.org/wiki/Great_circle_distance
  Input and output are in radians.  It also works for the az,el coordinate
  system.
  returnComponents=True will return: [separation, raSeparation,
  decSeparation, raSeparationCosDec]
  See also angularSeparation()
  -- Todd Hunter

    Args:
        returnComponents:
        dec1:
        ra1:
        dec0:
        ra0:
  """
    result = angularSeparation(ra0 * 180 / math.pi, dec0 * 180 / math.pi,
                               ra1 * 180 / math.pi, dec1 * 180 / math.pi,
                               returnComponents)
    if returnComponents:
        return np.array(result) * math.pi / 180.
    else:
        return result * math.pi / 180.
