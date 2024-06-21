import numpy as np
import matplotlib.pyplot as plt

path = "output/saxs_fitter/1ubq"

def s_to_q(s):
    return np.array(s)/(4*4*np.pi*np.pi)

def q_axis():
    return np.linspace(0.1, 1, 1000)

class five_gaussian_ff():
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def ff(self):
        q = q_axis()
        ff = np.zeros_like(q)
        for i in range(5):
            ff += self.a[i]*np.exp(-self.b[i]*q*q)
        return (ff + self.c)

    def ff_norm(self):
        v = self.ff()
        return v/v[0]

N       = five_gaussian_ff([11.893780,  3.277479,  1.858092, 0.858927, 0.912985],   s_to_q([0.000158, 10.232723, 30.344690, 0.656065, 0.217287]),   -11.804902  )
NH      = five_gaussian_ff([1.650531, 0.429639, 2.144736, 1.851894, 1.408921],      s_to_q([10.603730, 6.987283, 29.939901, 10.573859, 0.611678]),  0.510589    )
NH2     = five_gaussian_ff([1.904157, 1.942536,  2.435585,  0.730512, 1.379728],    s_to_q([10.803702, 10.792421, 29.610479, 6.847755, 0.709687]),  0.603738    )
CH_arom = five_gaussian_ff([2.168070, 1.275811,  1.561096,  0.742395, -6.151144],   s_to_q([12.642907, 18.420069, 41.768517, 1.535360, -0.045937]), 7.400917    )
CH_sp2  = five_gaussian_ff([2.909457, 0.484873,  1.515916,  0.207091, 1.541518],    s_to_q([13.934162, 23.229153, 41.991425, 4.983276, 0.679898]),  0.338296    )
CH_sp3  = five_gaussian_ff([2.909530, 0.485267,  1.516151,  0.206905, 1.541626],    s_to_q([13.933084, 23.221524, 41.990403, 4.974183, 0.679266]),  0.337670    )
CH2_sp3 = five_gaussian_ff([3.275723, 0.870037,  1.534606,  0.395078, 1.544562],    s_to_q([13.408502, 23.785175, 41.922444, 5.019072, 0.724439]),  0.377096    )
CH3_sp3 = five_gaussian_ff([3.681341, 1.228691,  1.549320,  0.574033, 1.554377],    s_to_q([13.026207, 24.131974, 41.869426, 4.984373, 0.765769]),  0.409294    )

def NH_plot():
    data_crysol = np.loadtxt(f"output/saxs_fitter/NH/crysol_aa.dat", skiprows=1)
    data_pepsi  = np.loadtxt(f"output/saxs_fitter/NH/pepsi_aa.dat",  skiprows=1)
    data_foxs   = np.loadtxt(f"output/saxs_fitter/NH/foxs_aa.dat",   skiprows=1)
    data_ausaxs = np.loadtxt(f"output/saxs_fitter/NH/ausaxs_aa.dat", skiprows=1)
    data_pepsi[:, 1] *= data_crysol[0, 1]/data_pepsi[0, 1]

    fig, ax = plt.subplots(figsize=(10, 6))
    plt.plot(data_crysol[:, 0], data_crysol[:, 1], label="crysol")
    plt.plot(data_pepsi[:, 0], data_pepsi[:, 1], label="pepsi")
    plt.plot(data_foxs[:, 0], data_foxs[:, 1], label="foxs")
    plt.plot(data_ausaxs[:, 0], data_ausaxs[:, 1], label="ausaxs")
    plt.plot(q_axis(), N.ff()*N.ff(), "--", label="$ff_N$")
    plt.plot(q_axis(), NH.ff()*NH.ff(), "--", label="$ff_{NH}$")
    plt.plot(q_axis(), NH2.ff()*NH2.ff()[0], "--", label="$ff_{NH2}$")
    plt.xlabel("q [1/A]")
    plt.ylabel("I(q)")
    plt.legend()
    plt.semilogx()
    plt.xlim(0.1, 1)
    # plt.ylim(0.8, 1.05)
    plt.show()

def CH_plot():
    data_crysol = np.loadtxt(f"output/saxs_fitter/CH/crysol_aa.dat", skiprows=1)
    data_pepsi  = np.loadtxt(f"output/saxs_fitter/CH/pepsi_aa.dat",  skiprows=1)
    data_foxs   = np.loadtxt(f"output/saxs_fitter/CH/foxs_aa.dat",   skiprows=1)
    data_ausaxs = np.loadtxt(f"output/saxs_fitter/CH/ausaxs_aa.dat", skiprows=1)

    data_pepsi[:, 1] /= data_pepsi[0, 1]
    data_ausaxs[:, 1] /= data_ausaxs[0, 1]

    fig, ax = plt.subplots(figsize=(10, 6))
    # plt.plot(data_crysol[:, 0], data_crysol[:, 1], label="crysol")
    plt.plot(data_pepsi[:, 0], data_pepsi[:, 1], label="pepsi")
    # plt.plot(data_foxs[:, 0], data_foxs[:, 1], label="foxs")
    plt.plot(data_ausaxs[:, 0], data_ausaxs[:, 1], label="ausaxs")
    plt.plot(q_axis(), CH_sp2.ff_norm()*CH_sp2.ff_norm(), "--", label="$ff_{CH SP2}$")
    plt.plot(q_axis(), CH_sp3.ff_norm()*CH_sp3.ff_norm(), "--", label="$ff_{CH SP3}$")
    plt.plot(q_axis(), CH2_sp3.ff_norm()*CH2_sp3.ff_norm(), "--", label="$ff_{CH2 SP3}$")
    plt.plot(q_axis(), CH3_sp3.ff_norm()*CH3_sp3.ff_norm(), "--", label="$ff_{CH3 SP3}$")
    plt.xlabel("q [1/A]")
    plt.ylabel("I(q)")
    plt.legend()
    plt.semilogx()
    plt.xlim(0.1, 1)
    # plt.ylim(0.8, 1.05)
    plt.show()

# NH_plot()
CH_plot()