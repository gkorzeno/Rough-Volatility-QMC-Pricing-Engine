import csv
from pathlib import Path

base = Path(r'tests/output')


def read_rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def line_points(xs, ys, x0, y0, w, h):
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0e-6
    pts = []
    for x, y in zip(xs, ys):
        xx = x0 + (x - xmin) / (xmax - xmin) * w
        yy = y0 + h - (y - ymin) / (ymax - ymin) * h
        pts.append((xx, yy))
    return pts, xmin, xmax, ymin, ymax

def render_panel(title, xlabel, ylabel, series, out_path, width=700, height=420):
    margin_l = 70
    margin_r = 20
    margin_t = 45
    margin_b = 55
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b
    x0 = margin_l
    y0 = margin_t
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>',
        'text { font-family: Arial, sans-serif; fill: #111; }',
        '.title { font-size: 18px; font-weight: bold; }',
        '.axis { font-size: 11px; fill: #444; }',
        '.legend { font-size: 11px; }',
        '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
        '.frame { stroke: #444; stroke-width: 1; fill: none; }',
        '</style>',
        f'<text x="{width/2}" y="24" text-anchor="middle" class="title">{esc(title)}</text>',
    ]
    all_x = []
    all_y = []
    computed = []
    for label, color, xs, ys in series:
        pts, xmin, xmax, ymin, ymax = line_points(xs, ys, x0, y0, plot_w, plot_h)
        computed.append((label, color, pts))
        all_x.extend(xs)
        all_y.extend(ys)
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0e-6
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
    for label, color, pts in computed:
        point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        for xx, yy in pts:
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{color}"/>')
    parts.append(f'<text x="{width/2}" y="{height - 12}" text-anchor="middle" class="axis">{esc(xlabel)}</text>')
    parts.append(f'<text x="18" y="{height/2}" text-anchor="middle" transform="rotate(-90 18 {height/2})" class="axis">{esc(ylabel)}</text>')
    lx = x0 + 8
    ly = y0 + 14
    for idx, item in enumerate(computed):
        label, color, _ = item
        yy = ly + idx * 16
        parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{color}" stroke-width="2.2"/>')
        parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{color}"/>')
        parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(label)}</text>')
    parts.append('</svg>')
    out_path.write_text('\n'.join(parts), encoding='utf-8')

smile = read_rows(base / 'rough_smile.csv')
term = read_rows(base / 'rough_term_structure.csv')
svi = read_rows(base / 'rough_vs_svi.csv')
var = read_rows(base / 'rough_mc_variance.csv')
calib = read_rows(base / 'rough_calibration_noise.csv')

render_panel(
    'Rough Bergomi Smile Diagnostics',
    'Strike',
    'Implied Vol',
    [('Rough Smile', '#1f77b4',
      [float(r['strike']) for r in smile],
      [float(r['implied_vol']) for r in smile])],
    base / 'rough_smile.svg'
)

render_panel(
    'ATM Term Structure',
    'Maturity',
    'Implied Vol',
    [('ATM IV', '#d62728',
      [float(r['expiry']) for r in term],
      [float(r['implied_vol']) for r in term])],
    base / 'rough_term_structure.svg'
)

svi_groups = {}
for row in svi:
    key = row['expiry']
    svi_groups.setdefault(key, {'x': [], 'rough': [], 'svi': []})
    svi_groups[key]['x'].append(float(row['strike']))
    svi_groups[key]['rough'].append(float(row['rough_iv']))
    svi_groups[key]['svi'].append(float(row['svi_iv']))
first_key = sorted(svi_groups.keys())[0]
render_panel(
    f'Rough Surface vs SVI Fit (T={first_key})',
    'Strike',
    'Implied Vol',
    [
        ('Rough', '#1f77b4', svi_groups[first_key]['x'], svi_groups[first_key]['rough']),
        ('SVI', '#2ca02c', svi_groups[first_key]['x'], svi_groups[first_key]['svi']),
    ],
    base / 'rough_vs_svi.svg'
)

render_panel(
    'Monte Carlo Variance Analysis',
    'Paths',
    'Surface RMSE',
    [('RMSE', '#9467bd',
      [float(r['paths']) for r in var],
      [float(r['rmse']) for r in var])],
    base / 'rough_mc_variance.svg'
)

if calib:
    by_sampler = {}
    for row in calib:
        by_sampler.setdefault(row['sampler'], {'x': [], 'prmse': [], 'ormse': []})
        by_sampler[row['sampler']]['x'].append(float(row['paths']))
        by_sampler[row['sampler']]['prmse'].append(float(row['parameter_rmse_mean']))
        by_sampler[row['sampler']]['ormse'].append(float(row['objective_rmse_mean']))

    render_panel(
        'Calibration Under Noise: Parameter RMSE',
        'Paths per evaluation',
        'Parameter RMSE',
        [(label, color, vals['x'], vals['prmse']) for label, color, vals in [
            ('MC', '#222222', by_sampler.get('MC', {'x': [], 'prmse': []})),
            ('QMC', '#1f77b4', by_sampler.get('QMC', {'x': [], 'prmse': []})),
            ('RQMC', '#d62728', by_sampler.get('RQMC', {'x': [], 'prmse': []})),
        ] if vals['x']],
        base / 'rough_calibration_noise.svg'
    )

print('Wrote rough-vol SVG diagnostics to', base)
