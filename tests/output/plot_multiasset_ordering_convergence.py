import csv
from collections import defaultdict
from pathlib import Path

csv_path = Path(r'tests/output/multiasset_ordering_convergence.csv')

rows = list(csv.DictReader(csv_path.open(newline='')))
out_svg = csv_path.with_suffix('.svg')
colors = {
    'MC': '#222222',
    'QMC asset-first': '#1f77b4',
    'QMC time-first': '#ff7f0e',
    'QMC PCA+bridge': '#2ca02c',
    'QMC adaptive': '#d62728',
}
series = defaultdict(lambda: {'x': [], 'y': [], 'slope': None})
for row in rows:
    series[row['label']]['x'].append(float(row['paths']))
    series[row['label']]['y'].append(float(row['error']))
    series[row['label']]['slope'] = float(row['slope']) if row['slope'] not in ('', 'nan') else None

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

labels = ['MC', 'QMC asset-first', 'QMC time-first', 'QMC PCA+bridge', 'QMC adaptive']
all_x = [x for s in labels for x in series[s]['x']]
all_y = [y for s in labels for y in series[s]['y']]
lx = [__import__('math').log10(v) for v in all_x]
ly = [__import__('math').log10(v) for v in all_y]
xmin, xmax = min(lx), max(lx)
ymin, ymax = min(ly), max(ly)
if xmax == xmin: xmax = xmin + 1
if ymax == ymin: ymax = ymin + 1

svg_w, svg_h = 960, 560
left, right, top, bottom = 86, 24, 60, 64
pw, ph = svg_w - left - right, svg_h - top - bottom
px, py = left, top

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 18px; font-weight: bold; }',
    '.subtitle { font-size: 12px; fill: #555; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
    '<text x="480" y="24" text-anchor="middle" class="title">Multi-Asset QMC Dimension Ordering Convergence</text>',
    '<text x="480" y="42" text-anchor="middle" class="subtitle">2-asset spread option under correlated GBM, log-log error vs paths</text>',
]
for i in range(5):
    xv = xmin + (xmax - xmin) * i / 4.0
    xx = px + (xv - xmin) / (xmax - xmin) * pw
    parts.append(f'<line x1="{xx:.2f}" y1="{py}" x2="{xx:.2f}" y2="{py + ph}" class="grid"/>')
    parts.append(f'<text x="{xx:.2f}" y="{py + ph + 18}" text-anchor="middle" class="axis">10^{xv:.2g}</text>')
for i in range(5):
    yv = ymin + (ymax - ymin) * i / 4.0
    yy = py + ph - (yv - ymin) / (ymax - ymin) * ph
    parts.append(f'<line x1="{px}" y1="{yy:.2f}" x2="{px + pw}" y2="{yy:.2f}" class="grid"/>')
    parts.append(f'<text x="{px - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">10^{yv:.2g}</text>')
parts.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" class="frame"/>')

for label in labels:
    pts = []
    for x, y in zip(series[label]['x'], series[label]['y']):
        xx = px + (__import__('math').log10(x) - xmin) / (xmax - xmin) * pw
        yy = py + ph - (__import__('math').log10(y) - ymin) / (ymax - ymin) * ph
        pts.append((xx, yy))
    point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
    parts.append(f'<polyline points="{point_str}" fill="none" stroke="{colors[label]}" stroke-width="2.2"/>')
    for xx, yy in pts:
        parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{colors[label]}"/>')

guide_specs = [(-0.5, 'MC guide N^-1/2', '#888888'), (-1.0, 'QMC guide N^-1', '#666666')]
anchor_x = max(all_x)
anchor_y = min(all_y)
for slope, guide_label, guide_color in guide_specs:
    guide_pts = []
    c = anchor_y / (anchor_x ** slope)
    for x in sorted(set(all_x)):
        y = c * (x ** slope)
        xx = px + (__import__('math').log10(x) - xmin) / (xmax - xmin) * pw
        yy = py + ph - (__import__('math').log10(y) - ymin) / (ymax - ymin) * ph
        guide_pts.append((xx, yy))
    point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in guide_pts)
    parts.append(f'<polyline points="{point_str}" fill="none" stroke="{guide_color}" stroke-width="1.4" stroke-dasharray="6,5"/>')

parts.append(f'<text x="{svg_w/2}" y="{svg_h - 16}" text-anchor="middle" class="axis">Paths N</text>')
parts.append(f'<text x="20" y="{svg_h/2}" text-anchor="middle" transform="rotate(-90 20 {svg_h/2})" class="axis">Absolute / RMS error</text>')
lx, ly0 = 110, 90
for i, label in enumerate(labels):
    yy = ly0 + i * 16
    slope = series[label]['slope']
    suffix = f' (slope={slope:.2f})' if slope is not None else ''
    parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{colors[label]}" stroke-width="2.2"/>')
    parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{colors[label]}"/>')
    parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(label + suffix)}</text>')
parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')
