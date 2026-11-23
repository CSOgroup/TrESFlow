"""Observed paired-fragment distributions for the bulk report."""
from collections import Counter
from html import escape
from math import ceil

from .plots import DNA, GRID, svg_figure


def fragment_figure(lengths, scope, full_range=False):
    """Plot a frequency polygon with all eligible fragments as the denominator."""
    lengths = Counter({int(length): int(count) for length, count in lengths.items() if int(count)})
    total = sum(lengths.values())
    maximum = max(lengths, default=0)
    limit = max(1000, ceil(maximum / 100) * 100) if full_range else 1000
    # The main view keeps 5 bp resolution, including short-fragment structure.
    # Bound full-range plotting work for unusually long inserts.
    bin_bp = max(5, ceil(limit / 1000 / 5) * 5) if full_range else 5
    bins = Counter()
    for length, count in lengths.items():
        if length <= limit:
            bins[min(length // bin_bp, (limit-1)//bin_bp)] += count
    outside = sum(n for length, n in lengths.items() if length > limit)
    median = lower_middle = None
    cumulative = 0
    for length, count in sorted(lengths.items()):
        cumulative += count
        if lower_middle is None and cumulative >= total / 2:
            lower_middle = length
        if cumulative > total / 2:
            median = (lower_middle + length) / 2
            break

    width, height = 1000, 370
    left, right, top, bottom = 80, 965, 76, 280
    peak = max(bins.values(), default=0) / total * 100 if total else 0
    ymax = ceil(peak / 5) * 5 if peak >= 5 else (ceil(peak * 2) / 2 or 1)
    parts = [
        '<text class="panel-title" x="24" y="28">Fragment length distribution</text>',
        f'<text class="axis muted" x="24" y="51">{escape(scope)}</text>',
        f'<text class="axis" x="{left}" y="69">Fragments (%)</text>',
    ]
    for tick in range(5):
        value = ymax * tick / 4
        y = bottom - (bottom-top)*tick/4
        parts.append(f'<path d="M{left} {y}H{right}" stroke="{GRID}" fill="none"/>')
        parts.append(f'<text class="axis muted" x="{left-10}" y="{y+4}" text-anchor="end">{value:g}</text>')
    for tick in range(6):
        value = limit*tick/5
        x = left+(right-left)*tick/5
        parts.append(f'<text class="axis muted" x="{x}" y="{bottom+23}" text-anchor="middle">{value:,.0f}</text>')
    parts.append(f'<path d="M{left} {top}V{bottom}H{right}" stroke="#aeb4b4" fill="none"/>')
    coordinates = []
    for index in range(ceil(limit/bin_bp)):
        x = left+(index+.5)*bin_bp/limit*(right-left)
        pct = bins[index]/total*100 if total else 0
        y = bottom-pct/ymax*(bottom-top)
        coordinates.append(f"{x:.2f},{y:.2f}")
        if bins[index]:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="{DNA}">'
                         f'<title>{index*bin_bp}–{min((index+1)*bin_bp, limit)} bp: '
                         f'{bins[index]:,} fragments ({pct:.2f}%)</title></circle>')
    points = " ".join(coordinates)
    parts.append(f'<polygon points="{left},{bottom} {points} {right},{bottom}" fill="{DNA}" opacity=".12"/>')
    parts.append(f'<polyline points="{points}" fill="none" stroke="{DNA}" stroke-width="2" stroke-linejoin="round"/>')
    median_text = f"{median:,.1f}".removesuffix(".0") + " bp" if median is not None else "unavailable"
    if median is not None and median <= limit:
        x = left+median/limit*(right-left)
        parts.append(f'<path d="M{x:.2f} {top}V{bottom}" stroke="#59666b" stroke-dasharray="5 5" fill="none">'
                     f'<title>Median {median_text}</title></path>')
    parts.append(f'<text class="axis" x="{(left+right)/2}" y="{bottom+49}" text-anchor="middle">Fragment length (bp)</text>')
    tail = f"{100*outside/total:.2f}%" if total else "unavailable"
    parts.append(f'<text class="axis muted" x="24" y="351">n = {total:,} final aligned fragments · Median {median_text} · {bin_bp} bp bins · {tail} beyond {limit:,} bp</text>')
    title = f"Fragment lengths · {scope}" + (" · full range" if full_range else "")
    return dict(kind="fragment_lengths_full" if full_range else "fragment_lengths", run=scope,
                svg=svg_figure(width, height, "".join(parts), title))
