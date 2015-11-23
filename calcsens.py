#!/usr/bin/env python
import sensitivity
import numpy as np
import sys
from optparse import OptionParser


def get_nominal_integration_time(calibration):
    if calibration == 'bandpass':
        return 300.
    elif calibration in ['phase']:
        return 60.
    elif calibration in ['pointing', 'focus']:
        return 10.
    elif calibration in ['amplitude']:
        return 120.
    elif calibration in ['check']:
        return 240.
    else:
        raise Exception('Unrecognized calibration type [%s]' % calibration)


def getNominalSNR(calibration):
    if calibration == 'bandpass':
        return 50.   # Antenna-based
    elif calibration in ['phase']:
        return 15.   # Antenna-based
    elif calibration in ['focus']:
        return 30.   # Antenna-based
    elif calibration in ['pointing']:
        return 20.   # Antenna-based
    elif calibration in ['amplitude']:
        return 100.  # Image
    elif calibration in ['check']:
        return 15.   # Image
    else:
        raise Exception('Unrecognized calibration type [%s]' % (calibration))


def getPWV(pwv, freq, el, s):
    if pwv == 'auto':
        pwv_ = s.selecticAutomaticPWV(freq, el)
    else:
        pwv_ = pwv
    print('[getPWV]    Assumed PWV=%5.3f [mm] %s' % (
        pwv_, "(auto)" if pwv is "auto" else ""))
    print('[getPWV]    f=%.1f [GHz] el=%.1f [deg]' % (
        freq / 1.0e9, el))
    return pwv_


def calcInterferometricSensitivity(
        freq, tint, BW, npol, el=60., pwv='auto', verbose=False,
        mode="image", nant=36, array='12m', returnFull=False):

    s = sensitivity.SensitivityCalculator(
        config_dir='/home/itoledo/Work/coneHighFreq/config/')
    pwv_ = getPWV(pwv, freq, el, s)

    if array in ['12m', 'TP']:
        D = 12.
        N = nant

    else:
        D = 7
        N = nant

    d = s.calcSensitivity(
        pwv_, freq, tint=tint, BW=BW, N_pol=npol, N=N, D=D,
        returnFull=True, el=el, mode=mode)
    if returnFull:
        return d
    return d['sensitivity']


def calcMinimumCalibratorFlux(
        calibration, spwList=None, SNR=None, tint=None,
        el=45., npol=1, returnFull=False, pwv='auto', verbose=False,
        nant=36, mode='image', array='12m'):

    """
    #   SNR_bcal : Minimum SNR for basndpass calibrator
    #   SNR_pcal : Minimum SNR for phase calibrator
    #   tint_bacl : Maximum integration time for bandpass calibrator
    #   tint_pcal : Maximum integration time for phase calibrator
    #   pmin_scan = 0.5
    """
    if spwList is None:
        raise Exception('Either spectralSpec or spwList should be specified')

    if SNR is None:
        SNR = getNominalSNR(calibration)

    if tint is None:
        raise Exception('Integration Time should be specified')

    dataList = []
    for spw in spwList:
        freq, BW = spw['freq'], spw['BW']
        if "nch" in spw:
            nch = spw["nch"]
        else:
            nch = 128
        # self.logInfo('spw=%s' % spw)
        # Determine parameters
        # if calibration == 'bandpass':
        #     # Channel bandwidth (Hz)
        #     # 128 MHz is recommended by CSV-2964.
        #     ch_BW = min(128.0e6, BW)
        #     # Smoothing window
        #     nsmooth = ch_BW / (BW / float(nch))
        #     BW = ch_BW
        #     mode = 'antenna'
        # elif calibration in ['pointing', 'phase', 'delay', 'focus']:
        #     mode = 'antenna'
        # elif calibration in ['amplitude', 'check']:
        #     mode = 'image'
        # else:
        #     raise Exception('Unsupported calibration type [%d]' % calibration)

        # Calculate required sensitivity
        d = calcInterferometricSensitivity(
            freq, tint, BW, npol, el=el, pwv=pwv, verbose=verbose, mode=mode,
            nant=nant, array=array, returnFull=True)

        d['SNR'] = SNR
        d['minFlux'] = d['sensitivity'] * SNR
        dataList.append(d)

        print(
            '[%s sensitivity] %8.3f [GHz] %8.3f [mJy] minSNR=%4.1f '
            'minFlux=%8.3f [mJy]' % (
                mode, freq / 1.0e9, d['sensitivity'] * 1.0e3, SNR,
                d['minFlux'] * 1.0e3))

    if returnFull:
        return dataList

    freqs = np.array([spw[0] for spw in spwList])
    minFluxList = np.array([d['minFlux'] for d in dataList])
    return freqs, minFluxList


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--bw", type=float, default=1.875e9,
        help="Bandwidth in Hz [default=1.875e9]")
    parser.add_option(
        "--npol", type=int, default=2,
        help="Number of polarizations [default=%default]")
    parser.add_option(
        "-c", "--calibration", default="phase",
        help="Calibratoin type [default=%default] (available types are phase, "
             "bandpass, check, and pointing.")
    parser.add_option(
        "-e", "--el", type=float, default=45.,
        help="Elevation in degrees [default=%default]")
    parser.add_option(
        "-C", "--arrayType", default="12m",
        help="Specify the type of array configuration ('TP', '7m', '12m'"
             "[default=%default]")
    parser.add_option(
        "-N", "--numAntennas", default=36, type=int,
        help="Specify the number of antennas. "
             "[default=%default]")
    parser.add_option(
        "-M", "--mode", default="image",
        help="Specify the sensitivy solution mode: "
             " 'image' - Imaging sensitivity; 'antenna' - Antenna based gain "
             "solution sensitivity. "
             "[default=%default]")
    parser.add_option(
        "-I", "--int", default=None,
        help="Specify integration time. Else is assumed to be 60s"
    )
    parser.usage = """
    %prog <frequency_in_Hz>

    # Example 1.
    PHASE freq=345 GHz, 12m array (12m x 36), BW=1.875 GHz, N(pol)=2
    getRequiredCalibratorFlux.py 345e9

    # Example 2.
    PHASE freq=345 GHz, 12m array (12m x 36), BW=1 GHz, N(pol)=1
    getRequiredCalibratorFlux.py 345e9 --bw=1e9 --npol=1

    # Example 3.
    PHASE freq=345 GHz, 7m array (7m x 10), BW=1 GHz, N(pol)=1
    getRequiredCalibratorFlux.py 345e9 --bw=1e9 --npol=1 -C 7m

    # Example 4.
    BANDPASS freq=800 GHz, 7m array (7m x 10)
    getRequiredCalibratorFlux.py 800e9 -C 7m -c bandpass

    # Example 5.
    CHECK_SOURCE freq=800 GHz, 7m array (7m x 10)
    getRequiredCalibratorFlux.py 800e9 -C 7m -c check"""

    opts, args = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    frequency = float(args[0])
    spwList = [dict(freq=frequency, BW=opts.bw)]
    tint = get_nominal_integration_time(opts.calibration)

    print "# Freq   = %.3f [GHz]" % (frequency / 1.0e9)
    print "# BW     = %.3f [GHz]" % (opts.bw / 1.0e9)
    print "# N(pol) = %d " % opts.npol
    print "# tint   = %.1f [seconds]" % tint
    print "# El     = %.1f [degrees]" % opts.el

    sensDataList = calcMinimumCalibratorFlux(
        opts.calibration,
        returnFull=True,
        spwList=spwList,
        tint=tint,
        el=opts.el,
        verbose=True,
        npol=opts.npol,
        nant=opts.numAntennas,
        mode=opts.mode,
        array=opts.arrayType)

    minFluxList = np.array([s["minFlux"] for s in sensDataList])
    for s in sensDataList:
        print "# Assumed PWV=%.3f [mm] Trx=%d [K] Tsys=%d [K] " % (
            s["pwv"], s["Trx"], s["Tsys"])
    print
    print "Required calibrator flux: %.2f [mJy]" % (minFluxList[0] * 1.e3)
    print
