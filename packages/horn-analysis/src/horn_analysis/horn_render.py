"""3D horn geometry rendering with surface-of-revolution visualization.

Produces a side-by-side figure: 3D surface-of-revolution plot (left)
and 2D wall profile cross-section (right). Uses the shared plot theme
for consistent styling across the analysis package.
"""

import argparse
import base64
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 – registers 3D projection

from horn_analysis import plot_theme


# -- Radius profile functions ------------------------------------------------

def _radius_profile(
    z: np.ndarray,
    throat_radius: float,
    mouth_radius: float,
    length: float,
    profile: str,
) -> np.ndarray:
    """Compute radius along the horn axis for a given flare profile.

    Parameters
    ----------
    z : array-like
        Axial positions (m), expected in [0, length].
    throat_radius, mouth_radius, length : float
        Horn dimensions (m).
    profile : str
        One of ``"conical"``, ``"exponential"``, ``"hyperbolic"``.

    Returns
    -------
    numpy array of radii at each z position.
    """
    z = np.asarray(z, dtype=float)
    if profile == "conical":
        return throat_radius + (mouth_radius - throat_radius) * z / length
    elif profile == "exponential":
        m = np.log(mouth_radius / throat_radius)
        return throat_radius * np.exp(m * z / length)
    elif profile == "hyperbolic":
        m = np.arccosh(mouth_radius / throat_radius)
        return throat_radius * np.cosh(m * z / length)
    else:
        raise ValueError(
            f"Unknown profile '{profile}'. "
            f"Choose from: conical, exponential, hyperbolic"
        )


# -- 3D rendering ------------------------------------------------------------

def render_horn_3d(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    profile: str,
    output_file: str,
    *,
    n_z: int = 80,
    n_theta: int = 60,
    show_profile: bool = True,
    figsize: tuple = (14, 6),
):
    """Render a 3D surface-of-revolution horn with optional 2D profile panel.

    Parameters
    ----------
    throat_radius, mouth_radius, length : float
        Horn dimensions (m).
    profile : str
        Flare profile name.
    output_file : str
        Path to save the rendered image.
    n_z : int
        Number of axial samples for the surface mesh.
    n_theta : int
        Number of angular samples around the axis.
    show_profile : bool
        If True, add a 2D cross-section panel on the right.
    figsize : tuple
        Figure size (width, height) in inches.
    """
    plot_theme.apply_theme()

    z_vals = np.linspace(0, length, n_z)
    r_vals = _radius_profile(z_vals, throat_radius, mouth_radius, length, profile)
    theta = np.linspace(0, 2 * np.pi, n_theta)

    Z, Theta = np.meshgrid(z_vals, theta)
    R = np.meshgrid(r_vals, theta)[0]
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)

    ncols = 2 if show_profile else 1
    fig = plt.figure(figsize=figsize)

    # -- 3D surface panel --
    ax3d = fig.add_subplot(1, ncols, 1, projection="3d")
    ax3d.plot_surface(
        Z, X, Y,
        cmap="coolwarm",
        alpha=0.85,
        edgecolor="none",
        rstride=1,
        cstride=1,
    )
    ax3d.set_xlabel("Axial position (m)")
    ax3d.set_ylabel("X (m)")
    ax3d.set_zlabel("Y (m)")
    ax3d.set_title(f"{profile.capitalize()} Horn — 3D View")
    ax3d.view_init(elev=20, azim=-60)

    # Equal aspect ratio for all three axes
    max_range = max(length, 2 * mouth_radius) / 2
    mid_z = length / 2
    ax3d.set_xlim(mid_z - max_range, mid_z + max_range)
    ax3d.set_ylim(-max_range, max_range)
    ax3d.set_zlim(-max_range, max_range)

    # -- 2D profile panel --
    if show_profile:
        ax2d = fig.add_subplot(1, 2, 2)
        r_mm = r_vals * 1000
        z_mm = z_vals * 1000

        ax2d.fill_between(z_mm, -r_mm, r_mm, alpha=0.15, color=plot_theme.COLORS["primary"])
        ax2d.plot(z_mm, r_mm, color=plot_theme.COLORS["primary"], linewidth=1.4, label="Wall")
        ax2d.plot(z_mm, -r_mm, color=plot_theme.COLORS["primary"], linewidth=1.4)

        # Annotate throat and mouth radii
        ax2d.annotate(
            f"Throat: {throat_radius * 1000:.1f} mm",
            xy=(0, throat_radius * 1000),
            xytext=(length * 1000 * 0.15, mouth_radius * 1000 * 0.9),
            fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8),
            color="#555555",
        )
        ax2d.annotate(
            f"Mouth: {mouth_radius * 1000:.1f} mm",
            xy=(length * 1000, mouth_radius * 1000),
            xytext=(length * 1000 * 0.65, mouth_radius * 1000 * 0.9),
            fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8),
            color="#555555",
        )

        ax2d.set_xlabel("Axial position (mm)")
        ax2d.set_ylabel("Radius (mm)")
        ax2d.set_title(f"{profile.capitalize()} Horn — Wall Profile")
        ax2d.set_aspect("equal", adjustable="datalim")
        plot_theme.setup_grid(ax2d)

    plot_theme.save_figure(fig, output_file)


def fig_to_b64_3d(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    profile: str,
    **kwargs,
) -> str:
    """Render a 3D horn to a base64 data-URI string for HTML embedding.

    Accepts the same keyword arguments as :func:`render_horn_3d`
    (except ``output_file``).

    Returns
    -------
    str
        ``data:image/png;base64,...`` encoded image.
    """
    plot_theme.apply_theme()

    buf = io.BytesIO()
    # Render to a temporary buffer by writing to BytesIO
    n_z = kwargs.pop("n_z", 80)
    n_theta = kwargs.pop("n_theta", 60)
    show_profile = kwargs.pop("show_profile", True)
    figsize = kwargs.pop("figsize", (14, 6))

    z_vals = np.linspace(0, length, n_z)
    r_vals = _radius_profile(z_vals, throat_radius, mouth_radius, length, profile)
    theta = np.linspace(0, 2 * np.pi, n_theta)

    Z, Theta = np.meshgrid(z_vals, theta)
    R = np.meshgrid(r_vals, theta)[0]
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)

    ncols = 2 if show_profile else 1
    fig = plt.figure(figsize=figsize)

    ax3d = fig.add_subplot(1, ncols, 1, projection="3d")
    ax3d.plot_surface(
        Z, X, Y,
        cmap="coolwarm",
        alpha=0.85,
        edgecolor="none",
        rstride=1,
        cstride=1,
    )
    ax3d.set_xlabel("Axial position (m)")
    ax3d.set_ylabel("X (m)")
    ax3d.set_zlabel("Y (m)")
    ax3d.set_title(f"{profile.capitalize()} Horn — 3D View")
    ax3d.view_init(elev=20, azim=-60)

    max_range = max(length, 2 * mouth_radius) / 2
    mid_z = length / 2
    ax3d.set_xlim(mid_z - max_range, mid_z + max_range)
    ax3d.set_ylim(-max_range, max_range)
    ax3d.set_zlim(-max_range, max_range)

    if show_profile:
        ax2d = fig.add_subplot(1, 2, 2)
        r_mm = r_vals * 1000
        z_mm = z_vals * 1000
        ax2d.fill_between(z_mm, -r_mm, r_mm, alpha=0.15, color=plot_theme.COLORS["primary"])
        ax2d.plot(z_mm, r_mm, color=plot_theme.COLORS["primary"], linewidth=1.4)
        ax2d.plot(z_mm, -r_mm, color=plot_theme.COLORS["primary"], linewidth=1.4)
        ax2d.set_xlabel("Axial position (mm)")
        ax2d.set_ylabel("Radius (mm)")
        ax2d.set_title(f"{profile.capitalize()} Horn — Wall Profile")
        ax2d.set_aspect("equal", adjustable="datalim")
        plot_theme.setup_grid(ax2d)

    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


# -- CLI ----------------------------------------------------------------------

def main():
    """Command-line interface: ``horn-render``."""
    parser = argparse.ArgumentParser(
        description="Render a 3D horn geometry visualization."
    )
    parser.add_argument("output", help="Output image file path (e.g. horn.png)")
    parser.add_argument("--throat-radius", type=float, required=True,
                        help="Throat radius in metres")
    parser.add_argument("--mouth-radius", type=float, required=True,
                        help="Mouth radius in metres")
    parser.add_argument("--length", type=float, required=True,
                        help="Horn length in metres")
    parser.add_argument("--profile", type=str, default="conical",
                        choices=["conical", "exponential", "hyperbolic"],
                        help="Horn flare profile (default: conical)")
    parser.add_argument("--no-profile-panel", action="store_true",
                        help="Omit the 2D cross-section panel")
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    render_horn_3d(
        throat_radius=args.throat_radius,
        mouth_radius=args.mouth_radius,
        length=args.length,
        profile=args.profile,
        output_file=args.output,
        show_profile=not args.no_profile_panel,
    )


if __name__ == "__main__":
    main()
