import os,sys
import numpy as np
import ROOT

secsInDay = 3600.*24.

def nuCrsScnGhandi(energy):
    # take energy in eV
    # return xscn in m^2
    return (7.84e-40 * ROOT.TMath.Power(energy/1e9,0.363))

def nuCrsScnCTWCC(energy):
    # take energy in eV
    # return xscn in m^2
    #
    # A. Connolly, R. S. Thorne, and D. Waters, Phys. Rev.D 83, 113009 (2011).
    #Charge Current Only
    e = ROOT.TMath.Log10(energy/1e9)
    c = (-1.826, -17.31, -6.406, 1.431, -17.91) # nu, CC
    l = ROOT.TMath.Log(e - c[0])
    x = c[1] + (c[2] * l) + (c[3] * l*l) + (c[4]/l)
    s  = ROOT.TMath.Power(10.0, x)
    s *= 1e-4 # cm^2 -> m^2
    return s

def nuCrsScnCTW(energy):
    # take energy in eV
    # return xscn in m^2
    #
    # A. Connolly, R. S. Thorne, and D. Waters, Phys. Rev.D 83, 113009 (2011).
    e = np.log10(energy*1e-9)
    cCC = (-1.826, -17.31, -6.406, 1.431, -17.91) # nu, CC
    l = np.log(e - cCC[0])
    xCC = cCC[1] + (cCC[2] * l) + (cCC[3] * l*l) + (cCC[4]/l)
    sCC  = 10.0** xCC

    cNC = (-1.826, -17.31, -6.448, 1.431, -18.61) # nu, NC
    xNC = cNC[1] + (cNC[2] * l) + (cNC[3] * l*l) + (cNC[4]/l)
    sNC  = 10.0 ** xNC
    s = sCC+sNC
    s *= 1e-4 # cm^2 -> m^2
    return s

def intLen(energy,
           nuCrsScn=nuCrsScnCTW):
    # take energy in eV
    # return interaction length in m
    return (1.66e-27 / (917.0*nuCrsScn(energy)))

def getEventsPerFluxPerEnergy(energy, evol,
                              livetime, signalEff,
                              nuCrsScn=nuCrsScnCTW):
    # energy in eV.
    # energyBinsPerDecade=1 for decade bins, 2 for half-decade bins, etc
    # evol = effective volume (km^3) * solid angle (sr)
    # livetime in days
    #
    # calculates the term:  eff * Veff * Omega * T / L
    # can be called multiple times in a sum over stations (or flavor,
    # which is why the interaction length function can be specified)
    # prior to finding a limit or sensitivity
    #
    # value return has units cm^2 sr s

    global secInDays

    ret  = (evol * signalEff * 1e9) # km^3 -> m^3
    ret *= livetime * secsInDay
    ret /= intLen(energy, nuCrsScn)
    ret *= 1e4 # m^2 -> cm^2
    return ret

def calcULE2FluxFromEvtsPerFluxPerE(energy,
                                    evtsPerFluxPerEnergy,
                                    energyBinsPerDecade=1.000,
                                    upperLimOnEvents=2.300):
    # energyBinsPerDecade=1 for decade bins, 2 for half-decade bins, etc
    # upperLimOnEvents = 2.3 for Neyman UL w/ 0 background,
    #                    2.44 for F-C UL w/ 0 background, etc
    #
    # returns upper limit on E^2 * flux, in GeV cm^-2 s^-1 sr^-1

    ul  = upperLimOnEvents / evtsPerFluxPerEnergy
    ul *= energyBinsPerDecade / ROOT.TMath.Ln10()
    ul *= energy/1.0e9 # eV -> GeV
    return ul


def getUpperLimitE2Flux(energy, evol,
                        livetime, signalEff,
                        verbose=True,
                        energyBinsPerDecade=1.000,
                        upperLimOnEvents=2.300,
                        nuCrsScn=nuCrsScnCTW):
    # energy in eV.
    # evol = effective volume (km^3) * solid angle (sr)
    # livetime in days
    # energyBinsPerDecade=1 for decade bins, 2 for half-decade bins, etc
    # upperLimOnEvents = 2.3 for Neyman UL w/ 0 background,
    #                    2.44 for F-C UL w/ 0 background, etc
    #
    # returns upper limit on E^2 * flux, in GeV cm^-2 s^-1 sr^-1


    if (verbose):
        print "evol={0:e} at E={1:e}. eff={2}, livetime={3}".format(evol,energy,
                                                                    signalEff,
                                                                    livetime)

    return calcULE2FluxFromEvtsPerFluxPerE(
        energy,
        getEventsPerFluxPerEnergy(energy, evol, livetime, signalEff,
                                  nuCrsScn),
        energyBinsPerDecade,
        upperLimOnEvents)


def getMultiStnUpperLimitE2Flux(energy, evol,
                                liveTimeSignalEffPairs,
                                verbose=True,
                                energyBinsPerDecade=1.000,
                                upperLimOnEvents=2.300,
                                nuCrsScn=nuCrsScnCTW):
    # energy in eV.
    # evol = effective volume (km^3) * solid angle (sr)
    # livetimes in days
    # liveTimeSignalEffPairs = [ (T_0, eff_0), (T_1, eff_1), ... ]
    # energyBinsPerDecade=1 for decade bins, 2 for half-decade bins, etc
    # upperLimOnEvents = 2.3 for Neyman UL w/ 0 background,
    #                    2.44 for F-C UL w/ 0 background, etc
    #
    # returns upper limit on E^2 * flux, in GeV cm^-2 s^-1 sr^-1


    evtsPerFluxPerEnergy = 0.0
    for livetime, signalEff in liveTimeSignalEffPairs:
        if (verbose):
            print "evol={0:e} at E={1:e}. eff={2}, livetime={3}".format(
                evol,energy, signalEff, livetime)
        evtsPerFluxPerEnergy += getEventsPerFluxPerEnergy(
            energy, evol, livetime, signalEff, nuCrsScn)

    return calcULE2FluxFromEvtsPerFluxPerE(
        energy,
        evtsPerFluxPerEnergy,
        energyBinsPerDecade,
        upperLimOnEvents)


def getNumNusDetectedAt(energyExp, estepExp,
                        logVeffFcn, sigeffFcn,
                        logFluxFcn, livetime,
                        callWithEnergyExp,
                        eToEv, distToM,
                        doPrint=False,
                        nuCrsScn=nuCrsScnCTW):
    # see getNumNusDetected (below) for comments.
    # this function is probably not useful outside of the loop used in
    # getNumNusDetected
    #
    # returns the number of nus detected in the specific energy bin

    nnu = 0.0
    e = ROOT.TMath.Power(10.0, energyExp);
    evar = e
    if (callWithEnergyExp):
        evar = energyExp
    # extrapolate using the log
    flux = ROOT.TMath.Power(10.0, logFluxFcn(evar))
    veff = 0.0
    seff = 0.0
    if (flux>0.0):
        # try extrapolating using the log
        veff = ROOT.TMath.Power(10.0, logVeffFcn(evar))
        if (veff<0.0):
            print "Invalid Veff={0:e} at E={1:e}. Assuming Veff=0.".format(
                veff, evar)
            veff=0.0
        seff = sigeffFcn(evar)

        eev = e * eToEv
        nnu  = veff * seff * flux * livetime * distToM
        nnu /= intLen( eev, nuCrsScn )
        nnu *= eev * ROOT.TMath.Ln10() * estepExp

        if (doPrint):
            print "E={0:e}, flux={1:e}, veff={2:e}, sigeff={3:e}, "\
                "T={4:e}, iLen={5:e}, n={6:e}".format(
                    evar, flux, veff, seff, livetime, intLen(eev), nnu)



    else:
        print "Warning: flux={0:e} at E={1:e}. Ignoring extrapolation."\
            .format(flux, evar)

    return nnu, veff, seff, flux



def getNumNusDetected(eminExp, emaxExp, estepExp,
                      logVeffFcn, sigeffFcn,
                      logFluxFcn, livetime,
                      callWithEnergyExp, eToEv, distToM,
                      doPrint=False,
                      nuCrsScn=nuCrsScnCTW):
    #
    # NOTE: units are assumed to already be consistent across the board!!
    #
    # eminExp  = the minimum energy exponent (i.e. 15.0 for 10^15.0)
    # emaxExp  = the maximum energy exponent (i.e. 15.0 for 10^15.0)
    # estepExp = the size of the energy exponent step to use in the
    #            integration. finer steps in principle give better accuracy
    # logVeffFcn = a function that returns log10(Veff) for a given energy.
    #              It can be a python function or a function on an object
    #              (e.g. graphOfLogVeff.Eval)
    # sigeffFcn = a function that returns the signal efficiency for a given
    #             energy
    # fluxFcn    = a function that returns the flux for a given energy
    # logFluxFcn = a function that return log10(flux) for a given energy.
    # livetime  = the livetime over which to integrate
    # callWithEnergyExp = if true, call the functions with the exponent
    #                     (i.e. 15.0); if false, call the functions with
    #                     the energy value (i.e. 10^15.0)
    # eToEv     = how to scale the energy to units of eV (i.e. 1e9 for GeV->eV)
    # distToM   = how to scale the distance to units of m (i.e. 1e3 for km->m)
    # doPrint   = if true, print the energy, flux and veff being used at
    #             each step (default: false)
    # nuCrsScn  = the cross section to use when calculating the interaction
    #             length. see intLen above.
    #
    #
    # returns:
    #   a) the integrated number of neutrinos detected
    #   b) a dictionary containing graphs (x-axis = energy exponent)
    #      b["flux"]  : a graph of the flux used at each energy
    #      b["veff"]  : a graph of the veff * signal eff used at each energy
    #      b["num"]   : a graph of the number of neutrinos detected
    #                   in each energy step
    #      b["cumul"] : a graph of the cumulative number of neutrinos
    #                   detected vs energy

    graphs = { "flux" : None, "veff" : None, "num" : None, "cumul" : None }
    for k in graphs.iterkeys():
        g = ROOT.TGraph()
        g.SetName("g" + k.title() + "Used")
        g.SetMarkerStyle(33)
        graphs[k] = g

    def addToGr(g, x, y):
        g.SetPoint( g.GetN(), x, y )

    evals = np.arange(eminExp, emaxExp, estepExp)
    nnu = 0.0
    for e in evals:

        num, veff, seff, flux = getNumNusDetectedAt(e, estepExp,
                                                    logVeffFcn,
                                                    sigeffFcn,
                                                    logFluxFcn,
                                                    livetime,
                                                    callWithEnergyExp,
                                                    eToEv, distToM,
                                                    doPrint,
                                                    nuCrsScn)

        nnu += num

        addToGr(graphs["flux"],  e, flux)
        addToGr(graphs["veff"],  e, veff*seff)
        addToGr(graphs["num"],   e, num)
        addToGr(graphs["cumul"], e, nnu)

    return nnu, graphs
