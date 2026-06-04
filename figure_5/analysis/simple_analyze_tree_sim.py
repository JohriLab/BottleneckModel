import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import argparse

from replicate_utils import load_jsonl_data

# Single size for axis labels, ticks, panel tags, and placeholder text.
FIG_FONT_SIZE = 14


def _savefig(path):
    """Vector formats avoid dpi; raster outputs use 300 dpi."""
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.pdf', '.svg', '.eps', '.ps'):
        plt.savefig(path, format=ext[1:], bbox_inches='tight')
    else:
        plt.savefig(path, dpi=300, bbox_inches='tight')


def compute_kingman_sfs_folded(x, num_samples, theta):
    """
    Compute folded Kingman SFS.
    
    For folded SFS, the formula is:
    SFS[i] = θ × (1/i + 1/(n-i)) / (1 + δ(i, n-i))
    where δ(i, n-i) = 1 when i = n-i, 0 otherwise
    
    Args:
        x: Array of allele frequencies [1, 2, ..., L] where L is the SFS length
        num_samples: Number of samples (n)
        theta: Population-scaled mutation rate (θ = 2*N*u)
    
    Returns:
        numpy array of SFS values matching the length of x
    """
    kingman_sfs = np.zeros_like(x, dtype=float)
    
    for idx, freq in enumerate(x):
        if freq == num_samples - freq:  # i = n-i case
            kingman_sfs[idx] = theta / freq
        else:  # i ≠ n-i case
            kingman_sfs[idx] = theta * (1/freq + 1/(num_samples - freq))
    
    return kingman_sfs


def compute_mean_sfs(df_filtered):
    """
    Compute mean SFS from filtered dataframe.
    
    Args:
        df_filtered: DataFrame filtered to specific (d, k, m) values
    
    Returns:
        mean_sfs: numpy array of mean SFS values
        sd_sfs: numpy array of standard error of SFS values (SD / sqrt(n_reps))
        u: mutation rate
    """
    if df_filtered.empty:
        return np.array([]), np.array([]), None
    
    u = df_filtered['u'].iloc[0]
    sfs_list = df_filtered['sfs'].tolist()
    
    # Stack all SFS arrays and scale by u
    sfs_matrix = np.stack(sfs_list) * u
    mean_sfs = sfs_matrix.mean(axis=0)
    n_reps = sfs_matrix.shape[0]
    sd_sfs = sfs_matrix.std(axis=0, ddof=1) / np.sqrt(n_reps)
    
    return mean_sfs, sd_sfs, u


def calc_pi_from_sfs(sfs, num_samples):
    """
    Calculate π (nucleotide diversity) from folded SFS.
    
    For folded SFS, π = Σ(i=1 to L) SFS[i] * i * (n-i) / (n choose 2)
    where L is the length of the folded SFS (floor(n/2))
    
    Args:
        sfs: Folded SFS array
        num_samples: Number of samples (n)
    
    Returns:
        pi: Nucleotide diversity
    """
    if len(sfs) == 0:
        return 0.0
    
    x = np.arange(1, len(sfs) + 1)
    pi = np.sum(sfs * x * (num_samples - x)) / (num_samples * (num_samples - 1) / 2)
    return pi


def _panel_label_outside_top_left(ax, label, fontsize):
    """Panel tag above the axes, outside the data area (top-left)."""
    ax.text(
        0.0,
        1.02,
        label,
        transform=ax.transAxes,
        fontsize=fontsize,
        ha='left',
        va='bottom',
        clip_on=False,
    )


def _format_float_no_sci(x, max_decimals=8):
    """
    Format a float in fixed-point notation (no scientific notation),
    trimming trailing zeros.
    """
    s = format(float(x), f".{int(max_decimals)}f")
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    # Avoid "-0"
    return "0" if s in ("-0", "-0.0") else s


def _param_combo_label(d_val, k_val, m_val, tv_distance=None):
    """Per-panel parameter label (m, d, k), with optional TV distance."""
    m_str = _format_float_no_sci(m_val, max_decimals=8)
    base = f"m={m_str}\nd={d_val}\nk={k_val}"
    if tv_distance is None:
        return base
    return f"Pairwise TV={float(tv_distance):.4f}\n{base}"


def _annotate_param_combo(ax, label, fontsize):
    ax.text(
        0.98,
        0.97,
        label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=fontsize,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.8, edgecolor="none"),
    )


def _plot_sfs_on_ax(
    ax,
    mean_sfs_sim,
    sd_sfs_sim,
    kingman_sfs,
    x,
    fontsize=FIG_FONT_SIZE,
    ylabel=True,
    xlabel=True,
    yscale="linear",
):
    """Draw simulation vs Kingman SFS on a matplotlib Axes."""
    cap = max(2, int(round(fontsize / 6)))
    ax.bar(
        x, mean_sfs_sim,
        yerr=sd_sfs_sim,
        color='#2E8B57',
        width=0.8,
        error_kw={'capsize': cap, 'elinewidth': 1, 'alpha': 0.9, 'ecolor': 'grey'},
        alpha=0.7,
    )
    ax.errorbar(
        x, kingman_sfs,
        color='#FF8C00',
        linewidth=0,
        alpha=1.0,
        marker='o',
        markerfacecolor='none',
        markeredgecolor='#FF8C00',
        markeredgewidth=1.8,
        markersize=7.5,
        capsize=cap,
        elinewidth=1,
    )
    if xlabel:
        ax.set_xlabel('Derived Allele Count', fontsize=fontsize)
    else:
        ax.set_xlabel('')
    if ylabel:
        ax.set_ylabel('Frequency', fontsize=fontsize)
    else:
        ax.set_ylabel('')
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.ticklabel_format(axis='y', style='plain', useOffset=False)
    if yscale == "log":
        ax.set_yscale("log", nonpositive="clip")
        y_all = np.concatenate([np.asarray(mean_sfs_sim).ravel(), np.asarray(kingman_sfs).ravel()])
        y_pos = y_all[np.isfinite(y_all) & (y_all > 0)]
        if y_pos.size:
            ax.set_ylim(bottom=max(0.8 * float(np.min(y_pos)), 1e-12))


def plot_sfs_comparison(
    mean_sfs_sim,
    sd_sfs_sim,
    kingman_sfs,
    x,
    d_val,
    k_val,
    m_val,
    num_samples,
    u,
    output_file,
    yscale="linear",
):
    """
    Create a comparison plot of simulation SFS vs Kingman SFS.
    
    Args:
        mean_sfs_sim: Mean SFS from simulations
        sd_sfs_sim: Standard error of simulation SFS
        kingman_sfs: Expected Kingman SFS
        x: Array of allele counts [1, 2, ..., L]
        d_val, k_val, m_val: Parameter values
        num_samples: Number of samples
        u: Mutation rate
        output_file: Output path (e.g. .pdf)
    """
    fontsize = FIG_FONT_SIZE
    fig, ax = plt.subplots(figsize=(6, 4.5))
    _plot_sfs_on_ax(ax, mean_sfs_sim, sd_sfs_sim, kingman_sfs, x, fontsize=fontsize, yscale=yscale)
    _annotate_param_combo(ax, _param_combo_label(d_val, k_val, m_val), fontsize)
    plt.tight_layout()
    _savefig(output_file)
    plt.close()
    print(f"Plot saved: {output_file}")


def plot_sfs_grid(df, panels, output_file, input_label='', yscale="linear"):
    """
    2x2 grid of SFS comparisons; each panel is one (d, k, m) triple.

    panels: list of four (d, k, m) tuples, in row-major order (top-left,
            top-right, bottom-left, bottom-right).
    """
    fontsize = FIG_FONT_SIZE
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    axes_flat = axes.ravel()
    panel_labels = [f'{chr(ord("a") + i)})' for i in range(4)]

    # Row-major: a=0, b=1, c=2, d=3 — omit y-axis label on b and d; x-axis label only on c and d.
    for idx, (ax, panel, panel_letter) in enumerate(zip(axes_flat, panels, panel_labels)):
        if len(panel) == 4:
            d_val, k_val, m_val, tv_distance = panel
        else:
            d_val, k_val, m_val = panel
            tv_distance = None
        df_filtered = df[
            (df['k'] == k_val) & (df['d'] == d_val) & (df['m'] == m_val)
        ]
        if df_filtered.empty:
            ax.text(
                0.5, 0.5,
                'No data',
                ha='center', va='center', transform=ax.transAxes, fontsize=fontsize
            )
            ax.set_axis_off()
            _panel_label_outside_top_left(ax, panel_letter, fontsize)
            continue

        mean_sfs_sim, sd_sfs_sim, u = compute_mean_sfs(df_filtered)
        if mean_sfs_sim.size == 0:
            ax.text(
                0.5, 0.5, 'No SFS data',
                ha='center', va='center', transform=ax.transAxes, fontsize=fontsize
            )
            ax.set_axis_off()
            _panel_label_outside_top_left(ax, panel_letter, fontsize)
            continue

        num_samples = df_filtered['num_samples'].iloc[0]
        L = len(mean_sfs_sim)
        x = np.arange(1, L + 1)
        pi_sim = calc_pi_from_sfs(mean_sfs_sim, num_samples)
        theta = pi_sim
        kingman_sfs = compute_kingman_sfs_folded(x, num_samples, theta)

        ylabel = idx not in (1, 3)
        xlabel = idx not in (0, 1)
        _plot_sfs_on_ax(
            ax, mean_sfs_sim, sd_sfs_sim, kingman_sfs, x,
            fontsize=fontsize, ylabel=ylabel, xlabel=xlabel, yscale=yscale,
        )
        _panel_label_outside_top_left(ax, panel_letter, fontsize)
        _annotate_param_combo(ax, _param_combo_label(d_val, k_val, m_val, tv_distance=tv_distance), fontsize)

    # Shared legend centered below the grid (proxy artists for bar + points).
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor='#2E8B57', edgecolor='none', alpha=0.7, label='Bottleneck Model Simulation'),
        Line2D(
            [0],
            [0],
            marker='o',
            linestyle='None',
            color='#FF8C00',
            markerfacecolor='none',
            markeredgecolor='#FF8C00',
            markeredgewidth=1.8,
            markersize=7.5,
            label='Kingman',
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.075),
        ncol=2,
        frameon=True,
        fontsize=fontsize,
    )

    _savefig(output_file)
    plt.close()
    print(f"Grid plot saved: {output_file}{' (' + input_label + ')' if input_label else ''}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare simulation SFS to expected Kingman SFS for one (d, k, m) triple, "
            "or a 2x2 grid of four triples (--grid with four --panel)."
        )
    )
    parser.add_argument('--k', type=int, default=None, help='k value (single-plot mode)')
    parser.add_argument('--d', type=int, default=None, help='d value (single-plot mode)')
    parser.add_argument('--m', type=float, default=0.0, help='m value (default: 0.0)')
    parser.add_argument(
        '--panel',
        nargs='+',
        action='append',
        metavar='VAL',
        help=(
            'One panel as either: d k m (no TV label) OR d k m tv_distance. '
            'Use with --grid; repeat exactly four times.'
        ),
    )
    parser.add_argument(
        '--grid',
        action='store_true',
        help='Write one 2x2 figure from four --panel d k m lines (matplotlib subplots).',
    )
    parser.add_argument('--input-file', type=str, default='simulation_results.jsonl',
                       help='Input JSONL file with simulation results')
    parser.add_argument('--output', type=str, default=None,
                       help='Output path, e.g. .pdf (required for --grid; default for single plot)')
    parser.add_argument(
        '--yscale',
        choices=('linear', 'log'),
        default='linear',
        help="Y-axis scale for SFS plots (default: linear).",
    )
    
    args = parser.parse_args()

    if args.grid:
        if not args.panel or len(args.panel) != 4:
            parser.error('--grid requires exactly four --panel arguments, each: d k m')
        if not args.output:
            parser.error('--grid requires --output path for the combined figure')
        panels = []
        for p in args.panel:
            if len(p) == 3:
                panels.append((int(p[0]), int(p[1]), float(p[2])))
            elif len(p) == 4:
                panels.append((int(p[0]), int(p[1]), float(p[2]), float(p[3])))
            else:
                parser.error('--panel must have 3 or 4 values: d k m [tv_distance]')
        df = load_jsonl_data(args.input_file)
        out_dir = os.path.dirname(args.output)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        plot_sfs_grid(df, panels, args.output, input_label=args.input_file, yscale=args.yscale)
        return

    if args.k is None or args.d is None:
        parser.error('single-plot mode requires --k and --d (or use --grid with four --panel)')

    # Load data
    df = load_jsonl_data(args.input_file)
    
    # Filter to specified (k, d, m) values
    df_filtered = df[
        (df['k'] == args.k) & 
        (df['d'] == args.d) & 
        (df['m'] == args.m)
    ]
    
    if df_filtered.empty:
        print(f"No data found for k={args.k}, d={args.d}, m={args.m}")
        return
    
    # Compute mean SFS from simulation
    mean_sfs_sim, sd_sfs_sim, u = compute_mean_sfs(df_filtered)
    
    if mean_sfs_sim.size == 0:
        print("No SFS data found")
        return
    
    # Get num_samples
    num_samples = df_filtered['num_samples'].iloc[0]
    
    # Determine x from the length of the SFS array
    L = len(mean_sfs_sim)
    x = np.arange(1, L + 1)
    
    # Calculate π from simulation SFS
    pi_sim = calc_pi_from_sfs(mean_sfs_sim, num_samples)
    
    # Compute expected Kingman SFS with matching π
    # For Kingman: π = θ, so theta = pi_sim
    theta = pi_sim
    kingman_sfs = compute_kingman_sfs_folded(x, num_samples, theta)
    
    # Ensure the plots directory exists before saving plots
    if not os.path.exists("plots"):
        os.makedirs("plots")

    # Determine output filename
    if args.output:
        output_file = args.output
    else:
        output_file = f"plots/sfs_d{args.d}_k{args.k}_m{args.m}.pdf"
    
    # Create plot
    plot_sfs_comparison(
        mean_sfs_sim, sd_sfs_sim, kingman_sfs,
        x, args.d, args.k, args.m, num_samples, u, output_file, yscale=args.yscale
    )


if __name__ == "__main__":
    main()
