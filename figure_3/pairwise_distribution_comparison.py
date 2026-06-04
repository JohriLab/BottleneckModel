import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================================
# PARAMETERS - Edit these values to change plot settings
# ============================================================================
# Two panels side by side; each dict is one (k, d, m) triple. T_max and x-axis
# limit follow the same formulas as before unless you override per panel.
panels = [
    {'k': 500, 'd': 500, 'm': 0},
    {'k': 25, 'd': 975, 'm': 0},
]
N = 1e10
output_file = 'plots/distribution_comparison_analytical_presentation.png'  # output plot file
pdf = True
log_x_axis = False  # Set to True to use log scale for x-axis
# X-axis upper limit = x_axis_max_mult * (d + k); smaller values zoom in on early times
x_axis_max_mult = 4
# Figure width (in) ~ typical \\linewidth; use \\includegraphics[width=\\linewidth]{...}
# so matplotlib’s 12 pt text stays ~12 pt in the compiled PDF
scaling = 1.25
figure_width_in = 6.5*scaling
figure_height_in = 4.0*scaling


def panel_limits(d, k):
    """Default T_max and x-axis max from d, k."""
    T_max = 100 * (d + k)
    x_axis_max = x_axis_max_mult * (d + k)
    return T_max, x_axis_max


# ============================================================================

def full_model_distribution(m, N, d, k, T_max):
    """
    p_times[t - 1, i] is P(coalesce at time t | start in state i), t = 1, ..., T_max.
    """

    a = 1 - 1/k
    b = 1 - 1/N
    c = 1/d
    alpha = 2*m*(1-m)

    Q = np.array([[a*b*(1-alpha), a*b*alpha], [c*b, 1-c]])

    n = Q.shape[0]
    powers = np.empty((T_max, 2, 2), dtype=Q.dtype)
    powers[0] = np.eye(2)
    for t in range(1, T_max):
        powers[t] = powers[t - 1] @ Q

    I = np.eye(n)
    ones = np.ones(n, dtype=Q.dtype)
    iq_ones = (I-Q) @ ones
    # p_times = np.einsum("tij,j->ti", powers, iq_ones)
    # p_times = powers @ iq_ones

    mu = np.array([[c, 1-c]])
    p_times = mu @ powers @iq_ones

    q_T = powers[-1] @ Q
    tail_mass = mu @q_T @ ones

    return p_times, tail_mass[0]

def compute_expectation(m, N, d, k):
    """
    returns the expectation of T_T for fixed n, N, d, k
    """

    a = 1 - 1/k
    b = 1 - 1/N
    c = 1/d
    alpha = 2*m*(1-m)

    Q = np.array([[a*b*(1-alpha), a*b*alpha], [c*b, 1-c]])

    n = Q.shape[0]
    I = np.eye(n)
    F = np.linalg.inv(I-Q)

    mu = np.array([[c, 1-c]])
    ones = np.ones(n, dtype=Q.dtype)
    expectation = mu @ F @ ones

    return expectation[0]

def wf_distribution(Ne, T_max):
    """
    Computes a geometric distribution with p = Ne (i.e. the distribution of time to coalescence under the WF model).
    """
    t = np.arange(0, T_max)

    p = 1/Ne
    p_times = p * (1-p)**(t)
    return p_times

def plot_distributions(panel_series):
    """
    panel_series: list of dicts with keys model, WF, T_max, x_axis_max (optional).
    """
    os.makedirs("plots", exist_ok=True)
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 12,
        'axes.labelsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'lines.linewidth': 2.0,
        'axes.linewidth': 0.8,
        'grid.linewidth': 0.5,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    n = len(panel_series)
    if n != 2:
        raise ValueError(f'Expected exactly 2 panels, got {n}')

    fig, axes = plt.subplots(
        1,
        2,
        sharey=True,
        figsize=(figure_width_in, figure_height_in),
        layout='constrained',
    )

    handles = None
    panel_labels = ['a)', 'b)']
    for i, (ax, spec) in enumerate(zip(axes, panel_series)):
        label_x = -0.01#-0.02 if i == 0 else -0.02
        ax.text(
            label_x,
            1.05,
            panel_labels[i],
            transform=ax.transAxes,
            ha='left',
            va='bottom',
            fontsize=12,
            clip_on=False,
        )

        model = np.asarray(spec['model']).reshape(-1)
        wf = np.asarray(spec['WF']).reshape(-1)
        T_max = spec['T_max']
        x_max = spec.get('x_axis_max')

        t = np.arange(1, T_max + 1)
        h1, = ax.plot(t, model, linewidth=2.0, color='#009E73', marker='o',
                      markersize=3, label='Bottleneck Model')
        h2, = ax.plot(t, wf, linewidth=2.0, color='#D55E00', marker='s',
                      markersize=3, label='Wright-Fisher')
        if handles is None:
            handles = (h1, h2)

        ax.set_xlabel('Time to coalescence')
        if log_x_axis:
            ax.set_xscale('log')
        if x_max is not None:
            ax.set_xlim(left=0, right=x_max)
        else:
            ax.set_xlim(left=0)
        ax.grid(True, alpha=0.3)
        ax.set_box_aspect(1)

    axes[0].set_ylabel('Probability')

    fig.legend(
        handles,
        [h.get_label() for h in handles],
        loc='outside lower center',
        ncol=2,
        frameon=True,
    )

    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    if pdf:
        try:
            pdf_file = output_file.rsplit('.', 1)[0] + '.pdf'
            plt.savefig(pdf_file, bbox_inches='tight')
        except Exception:
            pass
    plt.close()

def main():
    panel_series = []
    for i, p in enumerate(panels):
        k, d, m = p['k'], p['d'], p['m']
        T_max, x_axis_max = panel_limits(d, k)
        if 'T_max' in p:
            T_max = p['T_max']
        if 'x_axis_max' in p:
            x_axis_max = p['x_axis_max']

        print(f"Panel {i + 1}: computing model with m={m}, N={N}, d={d}, k={k}, T_max={T_max}.")
        model, tail = full_model_distribution(m, N, d, k, T_max)
        print(f"  Tail mass: {tail:.2e}")
        Ne = compute_expectation(m, N, d, k)
        print(f"  Computing WF distribution with Ne={Ne:.2e}")
        wf = wf_distribution(Ne, T_max)
        panel_series.append({
            'model': model,
            'WF': wf,
            'T_max': T_max,
            'x_axis_max': x_axis_max,
        })

    print("Plotting")
    plot_distributions(panel_series)
    print(f"Plot saved at {output_file}")

if __name__ == '__main__':
    main()
