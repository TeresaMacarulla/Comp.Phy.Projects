import numpy as np
import matplotlib.pyplot as plt

fname = "data/screen_prob_t2.dat"

data = np.loadtxt(fname, comments="#")
y = data[:, 0]
p = data[:, 1]

plt.figure(figsize=(6, 4))
plt.plot(y, p, lw=2)
plt.xlabel(r"$y$")
plt.ylabel(r"$p(y \mid x = 0.8, t = 0.002)$")
plt.title(r"Detection probability along screen at $x=0.8$")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plots/screen_probability_t2.pdf", dpi=200)
plt.show()
