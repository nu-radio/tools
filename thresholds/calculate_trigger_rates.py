import scipy.special
from scipy import stats
from NuRadioMC.utilities import units
from scipy import optimize as opt


def get_high_low_rate(sigma, dt=1 * units.ns, time_window_high_low=5 * units.ns):
    p1 = stats.norm.sf(sigma) / dt
    p2 = stats.norm.sf(3. / 4 * sigma) / dt
    return p1 * p2 * time_window_high_low


def get_global_trigger_rate(r_single, n_channels, n_coincidences, time_window):
    return n_coincidences * scipy.special.comb(n_channels, n_coincidences) * r_single**n_coincidences * time_window**(n_coincidences - 1)


def get_sigma(rate, dt=1 * units.ns, time_window_high_low=5 * units.ns):
    def obj(sigma):
        return get_high_low_rate(sigma, dt=dt, time_window_high_low=time_window_high_low) - rate
    res = opt.brentq(obj, 0, 10)
    return res


def get_threshold(r_global, n_channels, n_coincidences, time_window, dt=1 * units.ns, time_window_high_low=5 * units.ns):
    r_single = (r_global / n_coincidences / scipy.special.comb(n_channels, n_coincidences) /
                time_window**(n_coincidences - 1)) ** (1. / n_coincidences)
    sigma = get_sigma(r_single, dt, time_window_high_low=time_window_high_low)
    print("global rate {:.0f} mHz {}/{}, {:.0f} ns -> {:.1f}".format(r_global / units.Hz / units.milli, n_coincidences,
                                                                   n_channels, time_window, sigma))


if __name__ == "__main__":
    R = 0.01 * units.Hz
    s = get_threshold(R, 8, 3, time_window=200 * units.ns)
    s = get_threshold(R, 4, 4, time_window=40 * units.ns)
    s = get_threshold(R, 4, 2, time_window=40 * units.ns)
    s = get_threshold(R, 4, 4, time_window=5 * units.ns)
    s = get_threshold(R, 4, 2, time_window=5 * units.ns)
    
#     
#     R = .1 * units.Hz
#     s = get_threshold(R, 8, 3, time_window=200 * units.ns)
#     s = get_threshold(R, 4, 4, time_window=40 * units.ns)
#     s = get_threshold(R, 4, 2, time_window=40 * units.ns)
#     
#     R = 1 * units.Hz
#     s = get_threshold(R, 8, 3, time_window=200 * units.ns)
#     s = get_threshold(R, 4, 4, time_window=40 * units.ns)
#     s = get_threshold(R, 4, 2, time_window=40 * units.ns)
