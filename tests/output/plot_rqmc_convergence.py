import csv
import math
from pathlib import Path

csv_path = Path(r'tests/output/rqmc_convergence.csv')

out_svg = csv_path.with_suffix('.svg')
decomp_svg = csv_path.with_name('rqmc_error_decomposition.svg')
effdim_svg = csv_path.with_name('rqmc_effdim_vs_slope.svg')

rows = []
with csv_path.open(newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if not row or row.get('include_in_plot') in (None, ''):
            continue
        required = ['dimension', 'bridge', 'N', 'error', 'slope', 'sampling_error', 'discretization_bias', 'effective_dim_90']
        if any(row.get(key) in (None, '') for key in required):
            continue
        if int(row['include_in_plot']) != 1:
            continue
        row['dimension'] = int(row['dimension'])
        row['bridge'] = int(row['bridge'])
        row['N'] = int(row['N'])
        row['error'] = float(row['error'])
        row['slope'] = float(row['slope'])
        row['sampling_error'] = float(row['sampling_error'])
        row['discretization_bias'] = float(row['discretization_bias'])
        row['effective_dim_90'] = int(row['effective_dim_90'])
        rows.append(row)

dimensions = sorted({row['dimension'] for row in rows})
plot_order = ['MC', 'QMC', 'RQMC (Digital Shift)', 'RQMC (Owen Scramble)']
colors = {
    'MC': '#222222',
    'QMC': '#1f77b4',
    'RQMC (Digital Shift)': '#2ca02c',
    'RQMC (Owen Scramble)': '#d62728',
}

panel_w = 620
panel_h = 280
left_pad = 80
right_pad = 20
top_pad = 45
bottom_pad = 55
plot_w = panel_w - left_pad - right_pad
plot_h = panel_h - top_pad - bottom_pad
outer_pad = 20
svg_w = outer_pad * 3 + panel_w * 2
svg_h = outer_pad * (len(dimensions) + 1) + panel_h * len(dimensions)

def esc(text):
    return (str(text)
        .replace('&', '&amp;')
        .replace('<', '&lt;')
        .replace('>', '&gt;'))

def log_map(value, lo, hi, start, span):
    return start + (math.log10(value) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * span

def grouped(panel_rows):
    by_label = {}
    for label in plot_order:
        series = [row for row in panel_rows if row['plot_label'] == label]
        if series:
            by_label[label] = sorted(series, key=lambda row: row['N'])
    return by_label

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 20px; font-weight: bold; }',
    '.subtitle { font-size: 12px; fill: #444; }',
    '.panel-title { font-size: 14px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '.guide { stroke-width: 1.2; fill: none; opacity: 0.85; }',
    '.series { stroke-width: 2.2; fill: none; }',
    '</style>',
]

parts.append(f'<text x="{svg_w / 2}" y="26" text-anchor="middle" class="title">MC vs QMC vs RQMC convergence</text>')
parts.append(f'<text x="{svg_w / 2}" y="44" text-anchor="middle" class="subtitle">Geometric-Asian call under GBM, log-log error vs total paths N</text>')

for row_idx, dim in enumerate(dimensions):
    for col_idx, bridge in enumerate([0, 1]):
        x0 = outer_pad + col_idx * (panel_w + outer_pad)
        y0 = outer_pad + row_idx * (panel_h + outer_pad)
        plot_x0 = x0 + left_pad
        plot_y0 = y0 + top_pad

        panel_rows = [row for row in rows if row['dimension'] == dim and row['bridge'] == bridge]
        series_map = grouped(panel_rows)
        all_x = sorted({row['N'] for row in panel_rows})
        all_y = [row['error'] for row in panel_rows if row['error'] > 0.0]
        if not all_x or not all_y:
            continue

        y_min = min(all_y)
        y_max = max(all_y)
        y_min = 10 ** math.floor(math.log10(y_min))
        y_max = 10 ** math.ceil(math.log10(y_max))
        if y_min == y_max:
            y_min /= 10.0
            y_max *= 10.0

        parts.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="#ffffff" stroke="#dddddd"/>')
        bridge_label = 'ON' if bridge else 'OFF'
        parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + 20}" text-anchor="middle" class="panel-title">d = {dim}, Brownian bridge {bridge_label}</text>')

        for decade in range(int(math.log10(y_min)), int(math.log10(y_max)) + 1):
            tick_val = 10 ** decade
            yy = plot_y0 + plot_h - log_map(tick_val, y_min, y_max, 0, plot_h)
            parts.append(f'<line x1="{plot_x0}" y1="{yy:.2f}" x2="{plot_x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
            parts.append(f'<text x="{plot_x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">1e{decade}</text>')

        for xv in all_x:
            xx = log_map(xv, min(all_x), max(all_x), plot_x0, plot_w)
            parts.append(f'<line x1="{xx:.2f}" y1="{plot_y0}" x2="{xx:.2f}" y2="{plot_y0 + plot_h}" class="grid"/>')
            parts.append(f'<text x="{xx:.2f}" y="{plot_y0 + plot_h + 18}" text-anchor="middle" class="axis">{xv}</text>')

        parts.append(f'<rect x="{plot_x0}" y="{plot_y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
        parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + panel_h - 12}" text-anchor="middle" class="axis">N (total paths)</text>')
        parts.append(f'<text x="{x0 + 18}" y="{y0 + panel_h / 2}" text-anchor="middle" transform="rotate(-90 {x0 + 18} {y0 + panel_h / 2})" class="axis">Absolute error</text>')

        if 'MC' in series_map:
            mc = series_map['MC']
            x_ref = mc[0]['N']
            y_ref = mc[0]['error']
            points = []
            for xv in all_x:
                yv = y_ref * (xv / x_ref) ** (-0.5)
                xx = log_map(xv, min(all_x), max(all_x), plot_x0, plot_w)
                yy = plot_y0 + plot_h - log_map(yv, y_min, y_max, 0, plot_h)
                points.append(f'{xx:.2f},{yy:.2f}')
            parts.append(f'<polyline points="{" ".join(points)}" class="guide" stroke="#666666" stroke-dasharray="6,4"/>')

        if 'QMC' in series_map:
            qmc = series_map['QMC']
            x_ref = qmc[0]['N']
            y_ref = qmc[0]['error']
            points = []
            for xv in all_x:
                yv = y_ref * (xv / x_ref) ** (-1.0)
                xx = log_map(xv, min(all_x), max(all_x), plot_x0, plot_w)
                yy = plot_y0 + plot_h - log_map(yv, y_min, y_max, 0, plot_h)
                points.append(f'{xx:.2f},{yy:.2f}')
            parts.append(f'<polyline points="{" ".join(points)}" class="guide" stroke="#8888cc" stroke-dasharray="2,4"/>')

        for label in plot_order:
            if label not in series_map:
                continue
            series = series_map[label]
            pts = []
            for item in series:
                xx = log_map(item['N'], min(all_x), max(all_x), plot_x0, plot_w)
                yy = plot_y0 + plot_h - log_map(item['error'], y_min, y_max, 0, plot_h)
                pts.append((xx, yy))
            point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
            parts.append(f'<polyline points="{point_str}" class="series" stroke="{colors[label]}"/>')
            for xx, yy in pts:
                parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3.2" fill="{colors[label]}"/>')

        legend_x = plot_x0 + 8
        legend_y = plot_y0 + 16
        legend_items = []
        for label in plot_order:
            if label in series_map:
                slope = series_map[label][0]['slope']
                legend_items.append((label, f'{label} (slope={slope:.2f})'))
        legend_items.append(('guide_mc', 'MC guide N^-1/2'))
        legend_items.append(('guide_qmc', 'Ideal QMC guide N^-1'))

        for idx, item in enumerate(legend_items):
            key, text = item
            yy = legend_y + idx * 16
            if key == 'guide_mc':
                parts.append(f'<line x1="{legend_x}" y1="{yy}" x2="{legend_x + 18}" y2="{yy}" stroke="#666666" stroke-width="1.2" stroke-dasharray="6,4"/>')
            elif key == 'guide_qmc':
                parts.append(f'<line x1="{legend_x}" y1="{yy}" x2="{legend_x + 18}" y2="{yy}" stroke="#8888cc" stroke-width="1.2" stroke-dasharray="2,4"/>')
            else:
                parts.append(f'<line x1="{legend_x}" y1="{yy}" x2="{legend_x + 18}" y2="{yy}" stroke="{colors[key]}" stroke-width="2.2"/>')
                parts.append(f'<circle cx="{legend_x + 9}" cy="{yy}" r="3.0" fill="{colors[key]}"/>')
            parts.append(f'<text x="{legend_x + 24}" y="{yy + 4}" class="legend">{esc(text)}</text>')

parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')

decomp_rows = [row for row in rows if row['N'] == max(r['N'] for r in rows)]
bar_w = 180
bar_h = 240
svg_w2 = outer_pad * 3 + bar_w * 2
svg_h2 = outer_pad * (len(dimensions) + 1) + bar_h * len(dimensions)
parts2 = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w2}" height="{svg_h2}" viewBox="0 0 {svg_w2} {svg_h2}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 20px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.label { font-size: 11px; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
]
parts2.append(f'<text x="{svg_w2 / 2}" y="26" text-anchor="middle" class="title">QMC Error Decomposition</text>')
max_err = max(max(r['sampling_error'], r['discretization_bias']) for r in decomp_rows)
max_err = 1.1 * max_err if max_err > 0 else 1.0
for row_idx, dim in enumerate(dimensions):
    for col_idx, bridge in enumerate([0, 1]):
        panel = [r for r in decomp_rows if r['dimension'] == dim and r['bridge'] == bridge and r['plot_label'] == 'QMC']
        if not panel:
            continue
        item = panel[0]
        x0 = outer_pad + col_idx * (bar_w + outer_pad)
        y0 = outer_pad + row_idx * (bar_h + outer_pad)
        parts2.append(f'<rect x="{x0}" y="{y0}" width="{bar_w}" height="{bar_h}" fill="#ffffff" stroke="#dddddd"/>')
        parts2.append(f'<text x="{x0 + bar_w / 2}" y="{y0 + 18}" text-anchor="middle" class="label">d={dim}, bridge {"ON" if bridge else "OFF"}</text>')
        plot_y0 = y0 + 30
        plot_h2 = 140
        base_y = plot_y0 + plot_h2
        sampling_h = plot_h2 * item['sampling_error'] / max_err
        bias_h = plot_h2 * item['discretization_bias'] / max_err
        parts2.append(f'<rect x="{x0 + 30}" y="{base_y - sampling_h:.2f}" width="36" height="{sampling_h:.2f}" fill="#1f77b4"/>')
        parts2.append(f'<rect x="{x0 + 88}" y="{base_y - bias_h:.2f}" width="36" height="{bias_h:.2f}" fill="#d62728"/>')
        parts2.append(f'<line x1="{x0 + 18}" y1="{base_y}" x2="{x0 + bar_w - 18}" y2="{base_y}" class="frame"/>')
        parts2.append(f'<text x="{x0 + 48}" y="{base_y + 16}" text-anchor="middle" class="axis">sampling</text>')
        parts2.append(f'<text x="{x0 + 106}" y="{base_y + 16}" text-anchor="middle" class="axis">bias</text>')
        parts2.append(f'<text x="{x0 + bar_w / 2}" y="{y0 + 195}" text-anchor="middle" class="label">effective dim 90% = {item["effective_dim_90"]}</text>')
        parts2.append(f'<text x="{x0 + bar_w / 2}" y="{y0 + 212}" text-anchor="middle" class="label">ordering impact shown by early-dim concentration</text>')
parts2.append('</svg>')
decomp_svg.write_text('\n'.join(parts2), encoding='utf-8')
print(f'Wrote {decomp_svg}')

slope_rows = []
seen = set()
for row in rows:
    key = (row['dimension'], row['bridge'], row['plot_label'])
    if key in seen:
        continue
    seen.add(key)
    slope_rows.append(row)

qmc_like = [row for row in slope_rows if row['plot_label'] != 'MC']
if qmc_like:
    min_x = min(row['effective_dim_90'] for row in qmc_like)
    max_x = max(row['effective_dim_90'] for row in qmc_like)
    if min_x == max_x:
        min_x -= 1
        max_x += 1

    min_y = min(min(row['slope'] for row in qmc_like), -1.0)
    max_y = max(max(row['slope'] for row in qmc_like), -0.5)
    y_pad = 0.08 * max(0.2, max_y - min_y)
    min_y -= y_pad
    max_y += y_pad

    svg_w3 = 760
    svg_h3 = 520
    left3 = 90
    right3 = 30
    top3 = 55
    bottom3 = 70
    plot_w3 = svg_w3 - left3 - right3
    plot_h3 = svg_h3 - top3 - bottom3

    def x_map(x):
        return left3 + (x - min_x) / (max_x - min_x) * plot_w3

    def y_map(y):
        return top3 + plot_h3 - (y - min_y) / (max_y - min_y) * plot_h3

    marker_map = {
        (0, 'QMC'): 'circle',
        (1, 'QMC'): 'square',
        (0, 'RQMC (Digital Shift)'): 'triangle',
        (1, 'RQMC (Digital Shift)'): 'diamond',
        (0, 'RQMC (Owen Scramble)'): 'triangle_down',
        (1, 'RQMC (Owen Scramble)'): 'hex',
    }

    parts3 = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w3}" height="{svg_h3}" viewBox="0 0 {svg_w3} {svg_h3}">',
        '<style>',
        'text { font-family: Arial, sans-serif; fill: #111; }',
        '.title { font-size: 20px; font-weight: bold; }',
        '.axis { font-size: 12px; fill: #444; }',
        '.legend { font-size: 11px; }',
        '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
        '.frame { stroke: #444; stroke-width: 1; fill: none; }',
        '.guide-qmc { stroke: #1f77b4; stroke-width: 1.5; stroke-dasharray: 4,4; opacity: 0.8; }',
        '.guide-mc { stroke: #666666; stroke-width: 1.5; stroke-dasharray: 6,4; opacity: 0.8; }',
        '</style>',
        f'<text x="{svg_w3 / 2}" y="28" text-anchor="middle" class="title">Effective Dimension vs QMC Convergence Slope</text>',
        f'<text x="{svg_w3 / 2}" y="46" text-anchor="middle" class="axis">As EffDim90 rises, QMC slope typically degrades from ideal -1 toward MC-like -0.5</text>',
    ]

    x_ticks = sorted(set(row['effective_dim_90'] for row in qmc_like))
    for x_tick in x_ticks:
        xx = x_map(x_tick)
        parts3.append(f'<line x1="{xx:.2f}" y1="{top3}" x2="{xx:.2f}" y2="{top3 + plot_h3}" class="grid"/>')
        parts3.append(f'<text x="{xx:.2f}" y="{top3 + plot_h3 + 20}" text-anchor="middle" class="axis">{x_tick}</text>')

    for y_tick in [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5]:
        if y_tick < min_y or y_tick > max_y:
            continue
        yy = y_map(y_tick)
        parts3.append(f'<line x1="{left3}" y1="{yy:.2f}" x2="{left3 + plot_w3}" y2="{yy:.2f}" class="grid"/>')
        parts3.append(f'<text x="{left3 - 10}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{y_tick:.1f}</text>')

    parts3.append(f'<rect x="{left3}" y="{top3}" width="{plot_w3}" height="{plot_h3}" class="frame"/>')
    parts3.append(f'<text x="{svg_w3 / 2}" y="{svg_h3 - 18}" text-anchor="middle" class="axis">Effective Dimension (EffDim90)</text>')
    parts3.append(f'<text x="24" y="{svg_h3 / 2}" text-anchor="middle" transform="rotate(-90 24 {svg_h3 / 2})" class="axis">Log-log convergence slope</text>')

    y_qmc = y_map(-1.0)
    y_mc = y_map(-0.5)
    parts3.append(f'<line x1="{left3}" y1="{y_qmc:.2f}" x2="{left3 + plot_w3}" y2="{y_qmc:.2f}" class="guide-qmc"/>')
    parts3.append(f'<line x1="{left3}" y1="{y_mc:.2f}" x2="{left3 + plot_w3}" y2="{y_mc:.2f}" class="guide-mc"/>')
    parts3.append(f'<text x="{left3 + plot_w3 - 8}" y="{y_qmc - 6:.2f}" text-anchor="end" class="legend">Ideal QMC slope -1</text>')
    parts3.append(f'<text x="{left3 + plot_w3 - 8}" y="{y_mc - 6:.2f}" text-anchor="end" class="legend">MC slope -0.5</text>')

    for row in qmc_like:
        xx = x_map(row['effective_dim_90'])
        yy = y_map(row['slope'])
        label = row['plot_label']
        color = colors.get(label, '#000000')
        marker = marker_map.get((row['bridge'], label), 'circle')
        if marker == 'circle':
            parts3.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="5" fill="{color}"/>')
        elif marker == 'square':
            parts3.append(f'<rect x="{xx - 5:.2f}" y="{yy - 5:.2f}" width="10" height="10" fill="{color}"/>')
        elif marker == 'triangle':
            parts3.append(f'<polygon points="{xx:.2f},{yy - 6:.2f} {xx - 6:.2f},{yy + 5:.2f} {xx + 6:.2f},{yy + 5:.2f}" fill="{color}"/>')
        elif marker == 'diamond':
            parts3.append(f'<polygon points="{xx:.2f},{yy - 6:.2f} {xx - 6:.2f},{yy:.2f} {xx:.2f},{yy + 6:.2f} {xx + 6:.2f},{yy:.2f}" fill="{color}"/>')
        elif marker == 'triangle_down':
            parts3.append(f'<polygon points="{xx - 6:.2f},{yy - 5:.2f} {xx + 6:.2f},{yy - 5:.2f} {xx:.2f},{yy + 6:.2f}" fill="{color}"/>')
        else:
            parts3.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3.5" fill="{color}"/>')
            parts3.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="6" fill="none" stroke="{color}" stroke-width="1.5"/>')
        bridge_label = 'bridge ON' if row['bridge'] == 1 else 'bridge OFF'
        parts3.append(f'<text x="{xx + 8:.2f}" y="{yy - 6:.2f}" class="legend">{esc(label)}, {bridge_label}, d={row["dimension"]}</text>')

    parts3.append('</svg>')
    effdim_svg.write_text('\n'.join(parts3), encoding='utf-8')
    print(f'Wrote {effdim_svg}')
