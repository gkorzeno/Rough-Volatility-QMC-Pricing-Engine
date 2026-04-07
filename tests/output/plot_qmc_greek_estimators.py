import csv
import math
from pathlib import Path

csv_path = Path(r'qmc_greek_estimators.csv')

out_svg = csv_path.with_suffix('.svg')
rows = []
with csv_path.open(newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if not row or row.get('N') in (None, ''):
            continue
        row['N'] = int(row['N'])
        for key in ['delta_rmse', 'vega_rmse', 'delta_stability', 'vega_stability']:
            row[key] = float(row[key]) if row.get(key) not in (None, '', 'nan') else float('nan')
        for key in ['delta_slope', 'vega_slope']:
            row[key] = float(row[key]) if row.get(key) not in (None, '', 'nan') else float('nan')
        rows.append(row)

metrics = [
    ('delta_rmse', 'Delta RMSE vs Black-Scholes'),
    ('vega_rmse', 'Vega RMSE vs Black-Scholes'),
    ('delta_stability', 'Delta Stability (cross-run stddev)'),
    ('vega_stability', 'Vega Stability (cross-run stddev)'),
]
estimators = ['pathwise', 'likelihood_ratio', 'aad']
sampler_colors = {
    'MC': '#222222',
    'QMC': '#1f77b4',
    'RQMC (Owen)': '#d62728',
}
dash_map = {
    'pathwise': '',
    'likelihood_ratio': '6,4',
    'aad': '2,4',
}
label_map = {
    'pathwise': 'Pathwise',
    'likelihood_ratio': 'Likelihood Ratio',
    'aad': 'AAD',
}

panel_w = 540
panel_h = 280
outer = 22
left = 72
right = 18
top = 38
bottom = 56
plot_w = panel_w - left - right
plot_h = panel_h - top - bottom
svg_w = outer * 3 + panel_w * 2
svg_h = outer * 3 + panel_h * 2

def esc(text):
    return (str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 20px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.panel-title { font-size: 13px; font-weight: bold; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '.series { stroke-width: 2.0; fill: none; }',
    '</style>',
]
parts.append(f'<text x="{svg_w / 2}" y="26" text-anchor="middle" class="title">Pathwise vs Likelihood Ratio vs AAD under MC / QMC / RQMC</text>')

all_n = sorted({row['N'] for row in rows})
for idx, (metric, title) in enumerate(metrics):
    col = idx % 2
    row_idx = idx // 2
    x0 = outer + col * (panel_w + outer)
    y0 = outer + row_idx * (panel_h + outer)
    plot_x0 = x0 + left
    plot_y0 = y0 + top
    panel_rows = [r for r in rows if math.isfinite(r[metric])]
    if not panel_rows:
        continue
    vals = [r[metric] for r in panel_rows if r[metric] > 0]
    y_min = 10 ** math.floor(math.log10(min(vals)))
    y_max = 10 ** math.ceil(math.log10(max(vals)))
    if y_min == y_max:
        y_min /= 10.0
        y_max *= 10.0

    def x_map(x):
        return plot_x0 + (math.log10(x) - math.log10(min(all_n))) / (math.log10(max(all_n)) - math.log10(min(all_n))) * plot_w
    def y_map(y):
        return plot_y0 + plot_h - (math.log10(y) - math.log10(y_min)) / (math.log10(y_max) - math.log10(y_min)) * plot_h

    parts.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="#ffffff" stroke="#dddddd"/>')
    parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + 20}" text-anchor="middle" class="panel-title">{esc(title)}</text>')

    for decade in range(int(math.log10(y_min)), int(math.log10(y_max)) + 1):
        yy = y_map(10 ** decade)
        parts.append(f'<line x1="{plot_x0}" y1="{yy:.2f}" x2="{plot_x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{plot_x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">1e{decade}</text>')

    for n in all_n:
        xx = x_map(n)
        parts.append(f'<line x1="{xx:.2f}" y1="{plot_y0}" x2="{xx:.2f}" y2="{plot_y0 + plot_h}" class="grid"/>')
        parts.append(f'<text x="{xx:.2f}" y="{plot_y0 + plot_h + 18}" text-anchor="middle" class="axis">{n}</text>')

    parts.append(f'<rect x="{plot_x0}" y="{plot_y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
    parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + panel_h - 14}" text-anchor="middle" class="axis">N (paths)</text>')

    legend_y = plot_y0 + 16
    legend_x = plot_x0 + 8
    legend_i = 0
    for sampler in ['MC', 'QMC', 'RQMC (Owen)']:
        for estimator in estimators:
            series = [r for r in rows if r['sampler'] == sampler and r['estimator'] == estimator and math.isfinite(r[metric]) and r[metric] > 0]
            if not series:
                continue
            series = sorted(series, key=lambda r: r['N'])
            pts = []
            for item in series:
                pts.append((x_map(item['N']), y_map(item[metric])))
            point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
            dash = dash_map[estimator]
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ''
            color = sampler_colors[sampler]
            parts.append(f'<polyline points="{point_str}" class="series" stroke="{color}"{dash_attr}/>')
            for xx, yy in pts:
                parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{color}"/>')

            ly = legend_y + legend_i * 15
            parts.append(f'<line x1="{legend_x}" y1="{ly}" x2="{legend_x + 18}" y2="{ly}" stroke="{color}" stroke-width="2"{dash_attr}/>')
            slope_vals = [r[metric.replace("rmse", "slope")] for r in rows if r['sampler'] == sampler and r['estimator'] == estimator and 'rmse' in metric and math.isfinite(r.get(metric.replace("rmse", "slope"), float('nan')))]
            slope_txt = f' slope={slope_vals[0]:.2f}' if slope_vals else ''
            parts.append(f'<text x="{legend_x + 24}" y="{ly + 4}" class="legend">{esc(sampler + " / " + label_map[estimator] + slope_txt)}</text>')
            legend_i += 1

parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')
