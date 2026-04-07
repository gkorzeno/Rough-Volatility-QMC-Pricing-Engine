import csv
from pathlib import Path

csv_path = Path(r'tests/output/rough_calibration_noise.csv')

out_svg = csv_path.with_suffix('.svg')
rows = []
with csv_path.open(newline='') as f:
    rows = list(csv.DictReader(f))

samplers = ['MC', 'QMC', 'RQMC']
colors = {'MC': '#222222', 'QMC': '#1f77b4', 'RQMC': '#d62728'}
series = {s: {'x': [], 'prmse': [], 'ormse': [], 'pstd': [], 'prmse_slope': None, 'ormse_slope': None, 'pstd_slope': None} for s in samplers}
for row in rows:
    s = row['sampler']
    if s not in series:
        continue
    series[s]['x'].append(float(row['paths']))
    series[s]['prmse'].append(float(row['parameter_rmse_mean']))
    series[s]['ormse'].append(float(row['objective_rmse_mean']))
    series[s]['pstd'].append(float(row['parameter_std_norm']))
    series[s]['prmse_slope'] = float(row['parameter_rmse_slope']) if row['parameter_rmse_slope'] not in ('', 'nan') else None
    series[s]['ormse_slope'] = float(row['objective_rmse_slope']) if row['objective_rmse_slope'] not in ('', 'nan') else None
    series[s]['pstd_slope'] = float(row['parameter_std_norm_slope']) if row['parameter_std_norm_slope'] not in ('', 'nan') else None

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def render_panel(title, ylabel, key, slope_key, x0, y0, panel_w, panel_h, parts):
    left = 70
    right = 18
    top = 38
    bottom = 55
    plot_w = panel_w - left - right
    plot_h = panel_h - top - bottom
    px = x0 + left
    py0 = y0 + top
    all_x = [x for s in samplers for x in series[s]['x']]
    all_y = [y for s in samplers for y in series[s][key]]
    if not all_x or not all_y:
        return
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    if xmax == xmin:
        xmax = xmin + 1
    if ymax == ymin:
        ymax = ymin + 1e-6
    parts.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="#ffffff" stroke="#dddddd"/>')
    parts.append(f'<text x="{x0 + panel_w/2}" y="{y0 + 20}" text-anchor="middle" class="panel-title">{esc(title)}</text>')
    for i in range(5):
        xv = xmin + (xmax - xmin) * i / 4.0
        xx = px + (xv - xmin) / (xmax - xmin) * plot_w
        parts.append(f'<line x1="{xx:.2f}" y1="{py0}" x2="{xx:.2f}" y2="{py0 + plot_h}" class="grid"/>')
        parts.append(f'<text x="{xx:.2f}" y="{py0 + plot_h + 18}" text-anchor="middle" class="axis">{xv:.0f}</text>')
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4.0
        yy = py0 + plot_h - (yv - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<line x1="{px}" y1="{yy:.2f}" x2="{px + plot_w}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{px - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{yv:.3g}</text>')
    parts.append(f'<rect x="{px}" y="{py0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
    for s in samplers:
        if not series[s]['x']:
            continue
        pts = []
        for x, y in zip(series[s]['x'], series[s][key]):
            xx = px + (x - xmin) / (xmax - xmin) * plot_w
            yy = py0 + plot_h - (y - ymin) / (ymax - ymin) * plot_h
            pts.append((xx, yy))
        point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{colors[s]}" stroke-width="2.2"/>')
        for xx, yy in pts:
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{colors[s]}"/>')
    parts.append(f'<text x="{x0 + panel_w/2}" y="{y0 + panel_h - 14}" text-anchor="middle" class="axis">Paths per evaluation</text>')
    parts.append(f'<text x="{x0 + 18}" y="{y0 + panel_h/2}" text-anchor="middle" transform="rotate(-90 {x0 + 18} {y0 + panel_h/2})" class="axis">{esc(ylabel)}</text>')
    lx = px + 8
    ly = py0 + 14
    for i, s in enumerate(samplers):
        yy = ly + i * 16
        parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{colors[s]}" stroke-width="2.2"/>')
        parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{colors[s]}"/>')
        slope = series[s][slope_key]
        slope_txt = f' (slope={slope:.2f})' if slope is not None else ''
        parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(s + slope_txt)}</text>')

parts = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="780" viewBox="0 0 1120 780">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 18px; font-weight: bold; }',
    '.panel-title { font-size: 13px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
    '<text x="560" y="24" text-anchor="middle" class="title">Calibration Under Noise: MC vs QMC vs RQMC</text>',
]
render_panel('Recovered Parameter RMSE', 'Parameter RMSE', 'prmse', 'prmse_slope', 20, 40, 520, 320, parts)
render_panel('Objective RMSE', 'Surface RMSE', 'ormse', 'ormse_slope', 580, 40, 520, 320, parts)
render_panel('Recovered Parameter StdDev Norm', 'StdDev RMS(xi0, eta, rho, H)', 'pstd', 'pstd_slope', 300, 390, 520, 320, parts)
parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')
