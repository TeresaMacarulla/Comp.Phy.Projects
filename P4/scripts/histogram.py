# P4/scripts/eps_histograms.py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

TXT_DIR = Path(__file__).resolve().parents[1] / "txt"   # P4/txt

FILES = [
    ("eps_T_1.000000_burn-in_100.txt",     1.0),
    ("eps_T_2.400000_burn-in_100000.txt",  2.4),
]

def discrete_bins(values):
    """
    A partir de los valores discretos de epsilon, construye bordes de bins
    centrados en cada valor. Los bordes se ponen a mitad de camino
    entre valores consecutivos.
    """
    xu = np.array(sorted(np.unique(np.round(values, 12))), dtype=float)
    if xu.size == 1:
        w = 0.5
        return np.array([xu[0] - w, xu[0] + w]), xu
    mids = (xu[:-1] + xu[1:]) / 2.0
    first = xu[0] - (mids[0] - xu[0])
    last  = xu[-1] + (xu[-1] - mids[-1])
    edges = np.concatenate(([first], mids, [last]))
    return edges, xu

def make_histogram(txt_name, T):
    path = TXT_DIR / txt_name
    data = np.loadtxt(path)
    eps = data[:, 1]  # columna 2: epsilon = E/N

    edges, eps_levels = discrete_bins(eps)

    plt.figure(figsize=(7, 4.5))

    # Histograma normalizado manualmente (probabilidad de cada valor de epsilon)
    counts, edges = np.histogram(eps, bins=edges)
    probs = counts / counts.sum()  # normalización: suma total = 1
    centers = 0.5 * (edges[1:] + edges[:-1])

    plt.bar(centers, probs, width=np.diff(edges), align="center",
        alpha=0.7, edgecolor="black")


    # Marcas en los niveles discretos (opcional, útil para ver los centros)
    plt.scatter(eps_levels, np.zeros_like(eps_levels), s=12, c="k", zorder=3)

    plt.xlabel(r"$\epsilon = E/N$")
    plt.ylabel("Probability")
    # plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    PLOTS_DIR = Path(__file__).resolve().parents[1] / "plots"
    PLOTS_DIR.mkdir(exist_ok=True)  
    out_pdf = PLOTS_DIR / f"hist_{txt_name.replace('.txt','')}.pdf"
    plt.savefig(out_pdf)

    plt.bar(centers, probs, width=np.diff(edges), align="center",
        alpha=0.7, edgecolor="black",
        label=fr"$T={T:g}\, J/k_B$")
    print(f"{txt_name}: sum of probabilities = {probs.sum():.12f}")  # debería ser 1.0



if __name__ == "__main__":
    for fname, T in FILES:
        make_histogram(fname, T)
    plt.show()

