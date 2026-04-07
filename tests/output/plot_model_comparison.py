import csv
from collections import defaultdict
from pathlib import Path

base = Path(r'tests/output')

slice_rows = list(csv.DictReader((base / 'model_comparison_slices.csv').open(newline='')))
rmse_rows = list(csv.DictReader((base / 'model_comparison_rmse.csv').open(newline='')))
colors = {'Rough Bergomi': '#111111', 'Local Vol (Dupire)': '#1f77b4', 'Heston': '#d62728', 'Black-Scholes': '#2ca02c'}

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def render_line_panel(title, rows, kind, out_path):
    width, height = 720, 420
    left, right, top, bottom = 70, 20, 45, 55
    plot_w = width - left - right
    plot_h = height - top - bottom
    x0, y0 = left, top
    groups = defaultdict(lambda: {'x': [], 'y': []})
    for row in rows:
        if row['kind'] != kind:
            continue
        groups[row['model']]['x'].append(float(row['y']))
        groups[row['model']]['y'].append(float(row['implied_vol']))
    all_x = [x for g in groups.values() for x in g['x']]
    all_y = [y for g in groups.values() for y in g['y']]
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    if xmax == xmin: xmax = xmin + 1.0
    if ymax == ymin: ymax = ymin + 1e-6
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text { font-family: Arial, sans-serif; fill: #111; } .title { font-size: 18px; font-weight: bold; } .axis { font-size: 11px; fill: #444; } .legend { font-size: 11px; } .grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; } .frame { stroke: #444; stroke-width: 1; fill: none; }</style>',
        f'<text x="{width/2}" y="24" text-anchor="middle" class="title">{esc(title)}</text>',
    ]
    for i in range(5):
        xv = xmin + (xmax - xmin) * i / 4.0
        xx = x0 + (xv - xmin) / (xmax - xmin) * plot_w
        parts.append(f'<line x1="{xx:.2f}" y1="{y0}" x2="{xx:.2f}" y2="{y0 + plot_h}" class="grid"/>')
        parts.append(f'<text x="{xx:.2f}" y="{y0 + plot_h + 18}" text-anchor="middle" class="axis">{xv:.3g}</text>')
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4.0
        yy = y0 + plot_h - (yv - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<line x1="{x0}" y1="{yy:.2f}" x2="{x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{yv:.3g}</text>')
    parts.append(f'<rect x="{x0}" y="{y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
    for model, vals in groups.items():
        pts = []
        for x, y in zip(vals['x'], vals['y']):
            xx = x0 + (x - xmin) / (xmax - xmin) * plot_w
            yy = y0 + plot_h - (y - ymin) / (ymax - ymin) * plot_h
            pts.append((xx, yy))
        point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{colors[model]}" stroke-width="2.2"/>')
        for xx, yy in pts:
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{colors[model]}"/>')
    lx, ly = x0 + 8, y0 + 14
    for i, model in enumerate(groups.keys()):
        yy = ly + i * 16
        parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{colors[model]}" stroke-width="2.2"/>')
        parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{colors[model]}"/>')
        parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(model)}</text>')
    parts.append('</svg>')
    out_path.write_text('\n'.join(parts), encoding='utf-8')

def render_bar_panel(title, rows, out_path):
    width, height = 520, 360
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text { font-family: Arial, sans-serif; fill: #111; } .title { font-size: 18px; font-weight: bold; } .axis { font-size: 11px; fill: #444; } .frame { stroke: #444; stroke-width: 1; fill: none; }</style>',
        f'<text x="{width/2}" y="24" text-anchor="middle" class="title">{esc(title)}</text>',
    ]
    vals = [(row['model'], float(row['rmse'])) for row in rows]
    ymax = max(v for _, v in vals) * 1.1
    x0, y0, w, h = 60, 50, 420, 250
    parts.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" class="frame"/>')
    bw = w / max(1, len(vals)) * 0.55
    for i, (model, value) in enumerate(vals):
        cx = x0 + (i + 0.5) * w / len(vals)
        bh = h * value / ymax if ymax > 0 else 0
        parts.append(f'<rect x="{cx - bw/2:.2f}" y="{y0 + h - bh:.2f}" width="{bw:.2f}" height="{bh:.2f}" fill="{colors[model]}"/>')
        parts.append(f'<text x="{cx:.2f}" y="{y0 + h + 18}" text-anchor="middle" class="axis">{esc(model)}</text>')
        parts.append(f'<text x="{cx:.2f}" y="{y0 + h - bh - 6:.2f}" text-anchor="middle" class="axis">{value:.4f}</text>')
    parts.append('</svg>')
    out_path.write_text('\n'.join(parts), encoding='utf-8')

render_line_panel('Smile Comparison', slice_rows, 'smile', base / 'model_comparison_smile.svg')
render_line_panel('ATM Term Structure Comparison', slice_rows, 'term', base / 'model_comparison_term.svg')
render_bar_panel('Surface Calibration Error', rmse_rows, base / 'model_comparison_rmse.svg')
print('Wrote model comparison SVGs to', base)
