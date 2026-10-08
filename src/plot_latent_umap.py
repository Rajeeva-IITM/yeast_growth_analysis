"""UMAP and PCA of the per-condition chemistry embedding.
"""

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import umap  # noqa: E402
from dotenv import load_dotenv  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.metrics import silhouette_score  # noqa: E402

load_dotenv()

matplotlib.rcParams["svg.fonttype"] = "none"  # keep labels as editable text, not paths
matplotlib.rcParams["font.size"] = 10

N_LATENT = 256
MARKER_CLEARANCE = 5.0  # half-size of a scatter marker, in points
CHAR_WIDTH = 0.62  # mean glyph width as a fraction of font size, for label box estimates
LATENT_COLUMNS = [f"latent_{i}" for i in range(N_LATENT)]

# Okabe-Ito colours; markers differ too so the figure survives greyscale printing.
STYLE = {
    "salt": ("#0072B2", "o", "Salts"),
    "carbon_source": ("#009E73", "s", "Carbon sources"),
    "other": ("#999999", "^", "Others"),
}


def load(path: Path) -> tuple[pl.DataFrame, np.ndarray]:
    """Read the compound table and its latent matrix.

    Parameters
    ----------
    path : Path
        ``latent_conditions.csv``.

    Returns
    -------
    tuple[pl.DataFrame, np.ndarray]
        Metadata frame and the ``(n_compounds, 256)`` matrix.
    """
    frame = pl.read_csv(path)
    matrix = frame.select(LATENT_COLUMNS).to_numpy().astype(np.float64)
    return frame.drop(LATENT_COLUMNS), matrix


def collapse(frame: pl.DataFrame, matrix: np.ndarray) -> tuple[np.ndarray, list[str], list[str]]:
    """Reduce compounds to distinct latent vectors, merging labels that coincide.

    Parameters
    ----------
    frame : pl.DataFrame
        Metadata, carrying ``VectorID``, ``Key`` and ``Class``.
    matrix : np.ndarray
        Latent matrix aligned to ``frame``.

    Returns
    -------
    tuple[np.ndarray, list[str], list[str]]
        Distinct vectors, one display label per vector, and one class per vector.
    """
    vectors, labels, classes = [], [], []
    for (_,), group in frame.with_row_index("row").group_by(["VectorID"], maintain_order=True):
        members = group["Key"].to_list()
        present = set(group["Class"].to_list())
        if len(present) > 1:
            raise ValueError(
                f"compounds {members} share a latent vector but span classes {sorted(present)}; "
                "the figure cannot colour them"
            )

        vectors.append(matrix[group["row"].to_list()[0]])
        labels.append(" / ".join(key.replace("_", " ") for key in sorted(members)))
        classes.append(present.pop())

    return np.vstack(vectors), labels, classes


def boxes(
    labels: list[str], xlim: tuple, ylim: tuple, size: tuple, fontsize: float
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Estimated half-width and half-height of each text label, in data units.

    Parameters
    ----------
    labels : list[str]
        Label text; only its length is used.
    xlim, ylim : tuple
        Axis limits.
    size : tuple
        Axes width and height in inches.
    fontsize : float
        Label font size in points.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, float, float]
        Half-widths, half-heights, and the data units per point on each axis.
    """
    unit_x = (xlim[1] - xlim[0]) / (size[0] * 72)
    unit_y = (ylim[1] - ylim[0]) / (size[1] * 72)
    half_w = np.array([len(text) * CHAR_WIDTH * fontsize * unit_x / 2 for text in labels])
    half_h = np.full(len(labels), fontsize * 1.25 * unit_y / 2)
    return half_w, half_h, unit_x, unit_y


def nudge(
    points: np.ndarray, labels: list[str], xlim: tuple, ylim: tuple, size: tuple, fontsize: float
) -> np.ndarray:
    """Place each text label near its marker without colliding, deterministically.

    Greedy candidate selection rather than a repulsion loop: for every label, 24 candidate
    offsets around its point are scored against the markers and the labels already placed,
    and the cheapest is taken. This converges by construction -- the earlier physics loop
    oscillated, because pushing a label off a marker pushed it into a neighbouring label and
    back again. Crowded points are placed first, so they get first pick of the free space.

    Parameters
    ----------
    points : np.ndarray
        ``(n, 2)`` marker coordinates in data units.
    labels : list[str]
        Label text.
    xlim, ylim : tuple
        Axis limits.
    size : tuple
        Axes width and height in inches.
    fontsize : float
        Label font size in points.

    Returns
    -------
    np.ndarray
        ``(n, 2)`` label anchor coordinates in data units.
    """
    half_w, half_h, unit_x, unit_y = boxes(labels, xlim, ylim, size, fontsize)
    marker_w, marker_h = MARKER_CLEARANCE * unit_x, MARKER_CLEARANCE * unit_y

    directions = [
        (0.0, 1.0), (0.0, -1.0), (1.0, 0.0), (-1.0, 0.0),
        (0.75, 0.75), (-0.75, 0.75), (0.75, -0.75), (-0.75, -0.75),
    ]
    radii = (1.0, 1.7, 2.6)

    def overlap(ax, ay, aw, ah, bx, by, bw, bh) -> float:
        """Intersection area of two axis-aligned boxes, in squared points."""
        dx = (aw + bw) - abs(ax - bx)
        dy = (ah + bh) - abs(ay - by)
        if dx <= 0 or dy <= 0:
            return 0.0
        return (dx / unit_x) * (dy / unit_y)

    # Most crowded first: a label with many markers nearby has the fewest escape routes.
    span = max(np.ptp(points[:, 0]) / unit_x, np.ptp(points[:, 1]) / unit_y)
    separation = np.hypot(
        (points[:, None, 0] - points[None, :, 0]) / unit_x,
        (points[:, None, 1] - points[None, :, 1]) / unit_y,
    )
    crowding = (separation < 0.12 * span).sum(axis=1)
    order = sorted(range(len(labels)), key=lambda i: (-crowding[i], -half_w[i], i))

    anchors = np.zeros_like(points)
    placed: list[tuple] = []
    for i in order:
        best, best_score = None, None
        for radius in radii:
            for ux, uy in directions:
                x = points[i, 0] + ux * (half_w[i] + marker_w) * radius
                y = points[i, 1] + uy * (half_h[i] + marker_h) * radius * 1.7

                score = sum(
                    overlap(x, y, half_w[i], half_h[i], px, py, marker_w, marker_h)
                    for px, py in points
                ) + sum(
                    overlap(x, y, half_w[i], half_h[i], qx, qy, qw, qh)
                    for qx, qy, qw, qh in placed
                )
                # Tie-break toward the marker, so labels stay attributable.
                score += 0.6 * radius
                # A label running off the axes is worse than one sitting near a neighbour.
                outside = (
                    max(0.0, (xlim[0] + half_w[i]) - x) + max(0.0, x - (xlim[1] - half_w[i]))
                ) / unit_x + (
                    max(0.0, (ylim[0] + half_h[i]) - y) + max(0.0, y - (ylim[1] - half_h[i]))
                ) / unit_y
                score += 40.0 * outside

                if best_score is None or score < best_score:
                    best, best_score = (x, y), score

        # Final guard: no label may leave the axes even if every candidate did.
        anchors[i] = (
            min(max(best[0], xlim[0] + half_w[i]), xlim[1] - half_w[i]),
            min(max(best[1], ylim[0] + half_h[i]), ylim[1] - half_h[i]),
        )
        placed.append((anchors[i][0], anchors[i][1], half_w[i], half_h[i]))

    return anchors


def count_overlaps(
    anchors: np.ndarray, labels: list[str], xlim: tuple, ylim: tuple, size: tuple, fontsize: float
) -> int:
    """Number of label pairs whose estimated boxes still intersect.

    Parameters
    ----------
    anchors : np.ndarray
        Label anchors from :func:`nudge`.
    labels : list[str]
        Label text.
    xlim, ylim, size, fontsize
        As in :func:`nudge`.

    Returns
    -------
    int
        Remaining overlapping pairs.
    """
    half_w, half_h, _, _ = boxes(labels, xlim, ylim, size, fontsize)

    total = 0
    for i in range(len(anchors)):
        for j in range(i + 1, len(anchors)):
            if (
                abs(anchors[j, 0] - anchors[i, 0]) < half_w[i] + half_w[j]
                and abs(anchors[j, 1] - anchors[i, 1]) < half_h[i] + half_h[j]
            ):
                total += 1
    return total


def draw(
    coords: np.ndarray,
    labels: list[str],
    classes: list[str],
    counts: dict[str, int],
    axis_labels: tuple[str, str],
    title: str,
    out_path: Path,
    fontsize: float = 9.0,
) -> int:
    """Draw one 2-D embedding to SVG.

    Parameters
    ----------
    coords : np.ndarray
        ``(n, 2)`` embedding.
    labels : list[str]
        One label per point.
    classes : list[str]
        One class per point.
    counts : dict[str, int]
        Compounds per class, for the legend. Not derivable from ``classes``: two pairs of
        compounds share a marker, so the marker count is 2 below the compound count.
    axis_labels : tuple[str, str]
        x and y axis text.
    title : str
        Figure title.
    out_path : Path
        Destination ``.svg``.
    fontsize : float
        Compound label size.

    Returns
    -------
    int
        Label pairs still overlapping after nudging.
    """
    figure, axis = plt.subplots(figsize=(9.5, 7.6))

    span = max(np.ptp(coords[:, 0]), np.ptp(coords[:, 1]))
    pad = (0.12 + 0.012 * fontsize) * span
    xlim = (coords[:, 0].min() - pad, coords[:, 0].max() + pad)
    ylim = (coords[:, 1].min() - pad, coords[:, 1].max() + pad * 1.6)
    axis.set_xlim(xlim)
    axis.set_ylim(ylim)

    classes_array = np.array(classes)
    for name, (colour, marker, legend) in STYLE.items():
        mask = classes_array == name
        if not mask.any():
            continue
        axis.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=colour,
            marker=marker,
            s=46,
            linewidths=0.6,
            edgecolors="white",
            label=f"{legend} (n={counts[name]})",
            zorder=3,
        )

    box = axis.get_position()
    size = (figure.get_size_inches()[0] * box.width, figure.get_size_inches()[1] * box.height)
    anchors = nudge(coords, labels, xlim, ylim, size, fontsize)

    for point, anchor, text in zip(coords, anchors, labels):
        if np.hypot(*(anchor - point)) > 0.02 * (xlim[1] - xlim[0]):
            axis.plot(
                [point[0], anchor[0]],
                [point[1], anchor[1]],
                color="#bbbbbb",
                linewidth=0.4,
                zorder=1,
            )
        axis.annotate(
            text,
            xy=anchor,
            fontsize=fontsize,
            ha="center",
            va="center",
            color="#222222",
            zorder=4,
        )

    axis.set_xlabel(axis_labels[0])
    axis.set_ylabel(axis_labels[1])
    axis.set_title(title, fontsize=10)
    axis.grid(alpha=0.2, linewidth=0.5)
    axis.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=9)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)

    figure.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(figure)

    return count_overlaps(anchors, labels, xlim, ylim, size, fontsize)


def main() -> None:
    """Entry point."""
    project_dir = Path(os.environ["PROJECT_DIR"])

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-dir", type=Path, default=project_dir / "results/latent_space")
    parser.add_argument("--out-dir", type=Path, default=project_dir / "results/latent_space")
    parser.add_argument("--n-neighbors", type=int, default=8)
    parser.add_argument("--min-dist", type=float, default=0.15)
    parser.add_argument("--metric", default="cosine")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--fontsize", type=float, default=9.0)
    args = parser.parse_args()

    frame, matrix = load(args.in_dir / "latent_conditions.csv")
    vectors, labels, classes = collapse(frame, matrix)
    counts = dict(frame.group_by("Class").len().iter_rows())  # pyrefly: ignore
    print(f"{frame.height} compounds -> {len(labels)} distinct latent vectors")

    embedding = np.asarray(
        umap.UMAP(
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist,
            metric=args.metric,
            n_components=2,
            random_state=args.seed,
        ).fit_transform(vectors)
    )

    pca = PCA(n_components=2)
    principal = pca.fit_transform(vectors)
    variance = pca.explained_variance_ratio_

    args.out_dir.mkdir(parents=True, exist_ok=True)
    umap_overlaps = draw(
        embedding,
        labels,
        classes,
        counts,
        ("UMAP 1", "UMAP 2"),
        f"Chemistry embedding, UMAP (n_neighbors={args.n_neighbors}, "
        f"min_dist={args.min_dist}, metric={args.metric})",
        args.out_dir / "latent_umap.svg",
        args.fontsize,
    )
    pca_overlaps = draw(
        principal,
        labels,
        classes,
        counts,
        (f"PC1 ({variance[0]:.1%})", f"PC2 ({variance[1]:.1%})"),
        "Chemistry embedding, PCA",
        args.out_dir / "latent_pca.svg",
        args.fontsize,
    )

    # One row per compound, so the collided pairs share coordinates rather than vanishing.
    position = {}
    for index, label in enumerate(labels):
        for member in label.split(" / "):
            position[member] = index

    rows = []
    for record in frame.iter_rows(named=True):
        index = position[record["Key"].replace("_", " ")]
        rows.append(
            {
                "Key": record["Key"],
                "Name": record["Name"],
                "Class": record["Class"],
                "VectorID": record["VectorID"],
                "UMAP1": float(embedding[index, 0]),
                "UMAP2": float(embedding[index, 1]),
                "PC1": float(principal[index, 0]),
                "PC2": float(principal[index, 1]),
            }
        )
    pl.DataFrame(rows).write_csv(args.out_dir / "latent_embedding_coords.csv")

    silhouette = silhouette_score(matrix, frame["Class"].to_numpy(), metric="cosine")
    print(
        f"PCA explained variance: PC1 {variance[0]:.1%}, PC2 {variance[1]:.1%}, "
        f"sum {variance.sum():.1%}"
    )
    print(f"silhouette of the 3 classes in raw 256-D space (cosine): {silhouette:.3f}")
    print(f"residual label overlaps: UMAP {umap_overlaps}, PCA {pca_overlaps}")
    print(f"-> {args.out_dir}")


if __name__ == "__main__":
    main()
