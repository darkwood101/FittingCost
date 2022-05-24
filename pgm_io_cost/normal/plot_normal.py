import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

plt.rc('text', usetex = True)
params= {'text.latex.preamble' : [r'\usepackage{amsmath}',
                                  r'\usepackage[T1]{fontenc}',
                                  r'\usepackage[tt=false, type1=true]{libertine}',
                                  r'\usepackage[varqu]{zi4}',
                                  r'\usepackage[libertine]{newtxmath}']}
plt.rcParams.update(params)
plt.rcParams.update({'font.size': 9})

plt.rc('figure', figsize = (5, 4))

bar_width = 0.1

x_labels = [r"$4 \times 10^6$", r"$4 \times 10^7$", r"$4 \times 10^8$"]
exp = np.array([7812.5, 78125, 781250])
e32 = np.array([7821, 78189, 781434])
e64 = np.array([7817, 78141, 781331])
e128 = np.array([7816, 78132, 781281])

brexp = range(len(x_labels))
br32 = np.array([x + bar_width for x in brexp])
br64 = np.array([x + bar_width for x in br32])
br128 = np.array([x + bar_width for x in br64])


fig1, ax1 = plt.subplots()
#plt.title(r"Title", fontsize=9)
plt.xlabel(r"$n$")
plt.ylabel(r"IO cost")
plt.title(r"$x_i \sim \mathcal{{N}} \big( 0, 10^4 \big)$")
plt.bar(brexp, exp, color="blue", width=bar_width, edgecolor="grey", label="Predicted")
plt.bar(br32, e32, color="red", width=bar_width, edgecolor="grey", label=r"$\varepsilon = 32$")
plt.bar(br64, e64, color="yellow", width=bar_width, edgecolor="grey", label=r"$\varepsilon = 64$")
plt.bar(br128, e128, color="green", width=bar_width, edgecolor="grey", label=r"$\varepsilon = 128$")

plt.xticks([r + 1.5 * bar_width for r in brexp], x_labels)
plt.legend()

ax1.set_yscale("log")

#ax1.get_xaxis().set_major_formatter(ScalarFormatter())
plt.savefig("normal.pdf", bbox_inches="tight")

err32 = (e32 - exp) / exp
err64 = (e64 - exp) / exp
err128 = (e128 - exp) / exp

max_err = np.max(np.stack((err32, err64, err128)))
print(err32)
print(err64)
print(err128)
print(max_err)