import csv
from pathlib import Path

csv_path = Path(r'tests/output/parallel_scaling.csv')

rows = list(csv.DictReader(csv_path.open(newline='')))
out_svg = csv_path.with_suffix('.svg')

studies = ['strong', 'weak']
colors = {'strong': '#1f77b4', 'weak': '#d62728'}
series = {s: {'threads': [], 'speedup': [], 'eff': [], 'time': [], 'paths': []} for s in studies}
for row in rows:
    s = row['study']
    if s not in series:
        continue
    series[s]['threads'].append(int(row['threads']))
    series[s]['speedup'].append(float(row['speedup']))
    series[s]['eff'].append(float(row['efficiency']))
    series[s]['time'].append(float(row['mean_time_sec']))
    series[s]['paths'].append(int(row['paths']))

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def panel(parts, title, ylabel, key, x0, y0, w, h, guides=None):
    left, right, top, bottom = 64, 20, 38, 55
    pw = w - left - right
    ph = h - top - bottom
    px = x0 + left
    py = y0 + top
    xs = [x for s in studies for x in series[s]['threads']]
    ys = [y for s in studies for y in series[s][key]]
    if not xs or not ys:
        return
    if guides:
        ys = ys + guides
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax = xmin + 1
    if ymax == ymin:
        ymax = ymin + 1
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    parts.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="#fff" stroke="#dddddd"/>')
    parts.append(f'<text x="{x0 + w/2}" y="{y0 + 20}" text-anchor="middle" class="panel">{esc(title)}</text>')
    for i in range(5):
        xv = xmin + (xmax - xmin) * i / 4.0
        xx = px + (xv - xmin) / (xmax - xmin) * pw
        parts.append(f'<line x1="{xx:.2f}" y1="{py}" x2="{xx:.2f}" y2="{py + ph}" class="grid"/>')
        parts.append(f'<text x="{xx:.2f}" y="{py + ph + 18}" text-anchor="middle" class="axis">{xv:.0f}</text>')
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4.0
        yy = py + ph - (yv - ymin) / (ymax - ymin) * ph
        parts.append(f'<line x1="{px}" y1="{yy:.2f}" x2="{px + pw}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{px - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{yv:.3g}</text>')
    if guides:
        for gv in guides:
            yy = py + ph - (gv - ymin) / (ymax - ymin) * ph
            parts.append(f'<line x1="{px}" y1="{yy:.2f}" x2="{px + pw}" y2="{yy:.2f}" stroke="#888" stroke-dasharray="5,4" stroke-width="1.2"/>')
    parts.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" class="frame"/>')
    for s in studies:
        pts = []
        for x, y in zip(series[s]['threads'], series[s][key]):
            xx = px + (x - xmin) / (xmax - xmin) * pw
            yy = py + ph - (y - ymin) / (ymax - ymin) * ph
            pts.append((xx, yy))
        if not pts:
            continue
        point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{colors[s]}" stroke-width="2.2"/>')
        for xx, yy in pts:
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{colors[s]}"/>')
    parts.append(f'<text x="{x0 + w/2}" y="{y0 + h - 14}" text-anchor="middle" class="axis">Threads / cores</text>')
    parts.append(f'<text x="{x0 + 18}" y="{y0 + h/2}" text-anchor="middle" transform="rotate(-90 {x0 + 18} {y0 + h/2})" class="axis">{esc(ylabel)}</text>')

parts = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="780" viewBox="0 0 1120 780">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 18px; font-weight: bold; }',
    '.panel { font-size: 13px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
    '<text x="560" y="24" text-anchor="middle" class="title">Parallel Scaling Study</text>',
]
panel(parts, 'Speedup vs Cores', 'Speedup', 'speedup', 20, 50, 520, 320, guides=[1.0])
panel(parts, 'Efficiency vs Cores', 'Efficiency', 'eff', 580, 50, 520, 320, guides=[1.0])
panel(parts, 'Runtime vs Cores', 'Mean wall time (s)', 'time', 300, 400, 520, 320)
lx = 90
ly = 86
labels = {'strong': 'Strong scaling', 'weak': 'Weak scaling'}
for i, s in enumerate(studies):
    yy = ly + i * 16
    parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{colors[s]}" stroke-width="2.2"/>')
    parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{colors[s]}"/>')
    parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(labels[s])}</text>')
parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')
