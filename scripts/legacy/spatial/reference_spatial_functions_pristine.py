def rectangular_limits(xy: np.ndarray, pad_fraction: float = 0.03) -> tuple[float, float, float, float]:
    x0, x1 = float(xy[:, 0].min()), float(xy[:, 0].max())
    y0, y1 = float(xy[:, 1].min()), float(xy[:, 1].max())
    xpad = max(x1 - x0, 1.0) * pad_fraction
    ypad = max(y1 - y0, 1.0) * pad_fraction
    return x0 - xpad, x1 + xpad, y0 - ypad, y1 + ypad

def format_square(ax, limits: tuple[float, float, float, float]) -> None:
    xmin, xmax, ymin, ymax = limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_box_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.5)

def format_zoom_out(ax, limits: tuple[float, float, float, float]) -> None:
    xmin, xmax, ymin, ymax = limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

def add_zoom_roi_box(ax, roi: tuple[float, float, float, float]) -> None:
    xmin, xmax, ymin, ymax = roi
    ax.add_patch(Rectangle(
        (xmin, ymin), xmax - xmin, ymax - ymin, fill=False,
        edgecolor="black", linewidth=0.8, linestyle="--", zorder=9,
    ))

def plot_distribution(ax, xy: np.ndarray, target: np.ndarray, limits, target_color: str, zoom: bool = False) -> None:
    # Distribution 保持清爽的灰/红点，无 density 线
    size_other = 0.12 if not zoom else 0.7
    size_target = 0.45 if not zoom else 2.0
    
    ax.scatter(xy[~target, 0], xy[~target, 1], s=size_other, c="#E0E0E0", alpha=0.6, edgecolors="none", rasterized=True, zorder=1)
    ax.scatter(xy[target, 0], xy[target, 1], s=size_target, c=target_color, alpha=1.0, edgecolors="none", rasterized=True, zorder=4)

    (format_square if zoom else format_zoom_out)(ax, limits)

def nice_scale_bar(roi_side_coordinate: float, pixel_to_um: float) -> float:
    roi_um = roi_side_coordinate * pixel_to_um
    candidates = [20, 50, 100, 200, 500, 1000]
    valid = [x for x in candidates if x <= roi_um * 0.28]
    return float(valid[-1] if valid else candidates[0])

def add_scale_bar(ax, roi, pixel_to_um: float) -> float:
    xmin, xmax, ymin, ymax = roi
    width = xmax - xmin
    height = ymax - ymin
    bar_um = nice_scale_bar(width, pixel_to_um)
    bar_coord = bar_um / pixel_to_um
    x0 = xmin + width * 0.06
    y0 = ymin + height * 0.10
    ax.plot([x0, x0 + bar_coord], [y0, y0], color="black", linewidth=1.5, solid_capstyle="butt", zorder=8)
    ax.text(x0 + bar_coord / 2, y0 + height * 0.055, f"{int(bar_um)} µm", ha="center", va="bottom", fontsize=8, color="black", zorder=8)
    return bar_um