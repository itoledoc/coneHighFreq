#!/usr/bin/env python
import calcsens as cs
import longbaseline_cal as lbl
import pandas as pd
import os

from optparse import OptionParser

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--bw", type=float, default=1.875e9,
        help="Bandwidth in Hz [default=1.875e9]")
    parser.add_option(
        "--npol", type=int, default=2,
        help="Number of polarizations [default=%default]")
    parser.add_option(
        "-s", "--snr", default=10, type=float,
        help="Signal to noise desired [defaul=%default]")
    parser.add_option(
        "-m", "--mode", default="antenna",
        help="Specify the sensitivy solution mode: "
             " 'image' - Imaging sensitivity; 'antenna' - Antenna based gain "
             "solution sensitivity. "
             "[default=%default]")
    parser.add_option(
        "-i", "--intTime", default=60., type=float,
        help="Integration Time on calibrator [default=%60 s] in seconds"
    )
    parser.add_option(
        "-e", "--el", type=float, default=45.,
        help="Elevation in degrees [default=%default]")
    parser.add_option(
        "-n", "--numAntennas", default=36, type=int,
        help="Specify the number of antennas. "
             "[default=%default]")
    parser.add_option(
        "-c", "--arrayType", default="12m",
        help="Specify the type of array configuration ('TP', '7m', '12m'"
             "[default=%default]")
    parser.add_option(
        "-r", "--radius", default=15., type=float,
        help="Search radius in degrees [default=%default]")

    parser.usage = """
    %prog <coord string [HH:mm:ss.ss DD:mm:ss.ss]> <frequency_in_Hz>
    """

    opts, args = parser.parse_args()
    coords = str(args[0])
    frequency = float(args[1])

    spwList = [dict(freq=frequency, BW=opts.bw)]
    sensDataList = cs.calcMinimumCalibratorFlux(
        'phase',
        returnFull=True,
        spwList=spwList,
        tint=opts.intTime,
        el=opts.el,
        verbose=True,
        npol=opts.npol,
        nant=opts.numAntennas,
        mode=opts.mode,
        array=opts.arrayType)

    minFluxList = cs.np.array([s["minFlux"] for s in sensDataList])
    for s in sensDataList:
        print "# Assumed PWV=%.3f [mm] Trx=%d [K] Tsys=%d [K] " % (
            s["pwv"], s["Trx"], s["Tsys"])
    print
    print "Required calibrator flux: %.2f [mJy]" % (minFluxList[0] * 1.e3)
    print

    minflux = minFluxList[0]

    a = lbl.find_nearest_vlbi_sources(
        coords, radius=opts.radius, frequency=frequency, phasecal_flux_limit=minflux)
    df = pd.DataFrame(a).transpose()

    df.query('fluxDensity > @minflux').sort(
        'separationDegrees').to_csv('calib_candidates_hf.csv')
    df.sort(
            'separationDegrees').to_csv('cone_hf.csv')
