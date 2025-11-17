# P4/scripts/fit_Tc_vs_L.py

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    # --- Datos de la tabla ---
    L = np.array([40, 60, 80, 100], dtype=float)
    Tc = np.array([2.30, 2.295, 2.28, 2.27], dtype=float)  # [J/k_B]
    dTc = np.array([0.02, 0.005, 0.01, 0.01], dtype=float)  # errores

    # Variable de ajuste: x = 1/L
    x = 1.0 / L

    # --- Ajuste lineal Tc(L) = a * (1/L) + b ---
    # Usamos pesos inversos a la incertidumbre: w = 1/sigma
    coeffs, cov = np.polyfit(x, Tc, deg=1, w=1.0 / dTc, cov=True)
    a, b = coeffs          # a = pendiente, b = Tc(L->∞)
    da, db = np.sqrt(np.diag(cov))

    # Predicción en los puntos de datos (para calcular R^2)
    Tc_pred = a * x + b
    ss_res = np.sum((Tc - Tc_pred) ** 2)
    ss_tot = np.sum((Tc - Tc.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot

    # Curva suave para la recta ajustada
    x_fit = np.linspace(0.0, x.max() * 1.05, 200)
    Tc_fit = a * x_fit + b

    # --- Figura ---
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    # Datos con barras de error
    ax.errorbar(
        x, Tc, yerr=dTc,
        fmt="o", color="C0", ecolor="k", capsize=4,
        label=r"Data: $T_c(L)$"
    )

    # Recta ajustada
    ax.plot(
        x_fit, Tc_fit, "-",
        color="C1", label=r"Linear fit: $T_c(L) = a/L + T_c(\infty)$"
    )

    # Etiquetas y título
    ax.set_xlabel(r"$1/L$")
    ax.set_ylabel(r"$T_c(L)\;[J/k_B]$")
    ax.set_title(r"Finite-size scaling of $T_c(L)$ vs $1/L$")

    # Cuadro de texto con parámetros del ajuste
    textstr = "\n".join([
        rf"$a = {a:.3f} \pm {da:.3f}$",
        rf"$T_c(\infty) = {b:.4f} \pm {db:.4f}\, J/k_B$",
        rf"$R^2 = {r2:.4f}$"
    ])
    ax.text(
        0.05, 0.95, textstr,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()

    fig.tight_layout()

    # --- Guardar PDF en P4/plots ---
    script_dir = Path(__file__).resolve().parent
    plots_dir = script_dir.parent / "plots"
    plots_dir.mkdir(exist_ok=True)
    #out_pdf = plots_dir / "Tc_vs_invL_fit.pdf"
    #fig.savefig(out_pdf, bbox_inches="tight")
    #print(f"Saved figure to: {out_pdf}")
    plt.show()

if __name__ == "__main__":
    main()
