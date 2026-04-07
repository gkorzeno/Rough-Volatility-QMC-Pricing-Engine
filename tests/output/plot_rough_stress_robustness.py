import csv
from pathlib import Path

base = Path(r'tests/output')

pricing_rows = list(csv.DictReader((base / 'rough_stress_pricing.csv').open(newline='')))
summary_rows = list(csv.DictReader((base / 'rough_stress_summary.csv').open(newline='')))

scenario_order = [row['scenario'] for row in summary_rows]
family_by_scenario = {row['scenario']: row['family'] for row in summary_rows}
sampler_colors = {'MC': '#222222', 'QMC': '#1f77b4', 'RQMC': '#d62728'}
family_colors = {'baseline': '#9ecae1', 'correlation': '#fdd0a2', 'vol-of-vol': '#fdae6b', 'roughness': '#c7e9c0'}

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def grouped_bars(parts, panel_x, panel_y, panel_w, panel_h, title, values_by_key, y_label, guide_lines=None):
    left, right, top, bottom = 64, 20, 40, 78
    plot_w = panel_w - left - right
    plot_h = panel_h - top - bottom
    x0 = panel_x + left
    y0 = panel_y + top
    values = [v for v in values_by_key.values() if v is not None]
    if not values:
        return
    if guide_lines:
        values.extend(guide_lines)
    ymin = min(values)
    ymax = max(values)
    if ymax == ymin:
        ymax = ymin + 1.0
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    parts.append(f'<rect x="{panel_x}" y="{panel_y}" width="{panel_w}" height="{panel_h}" fill="#fff" stroke="#dddddd"/>')
    parts.append(f'<text x="{panel_x + panel_w/2}" y="{panel_y + 20}" text-anchor="middle" class="legend">{esc(title)}</text>')
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4.0
        yy = y0 + plot_h - (yv - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<line x1="{x0}" y1="{yy:.2f}" x2="{x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{yv:.3g}</text>')
    if guide_lines:
        for gv in guide_lines:
            yy = y0 + plot_h - (gv - ymin) / (ymax - ymin) * plot_h
            parts.append(f'<line x1="{x0}" y1="{yy:.2f}" x2="{x0 + plot_w}" y2="{yy:.2f}" stroke="#888" stroke-dasharray="5,4" stroke-width="1.2"/>')
    parts.append(f'<rect x="{x0}" y="{y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
    group_w = plot_w / max(1, len(scenario_order))
    bar_w = group_w * 0.18
    offsets = {'MC': -bar_w * 1.2, 'QMC': 0.0, 'RQMC': bar_w * 1.2}
    for i, scenario in enumerate(scenario_order):
        gx = x0 + (i + 0.5) * group_w
        bg = family_colors.get(family_by_scenario.get(scenario, ''), '#f4f4f4')
        parts.append(f'<rect x="{x0 + i * group_w + 2:.2f}" y="{y0}" width="{group_w - 4:.2f}" height="{plot_h}" fill="{bg}" opacity="0.25"/>')
        for sampler in ['MC', 'QMC', 'RQMC']:
            value = values_by_key.get((scenario, sampler))
            if value is None:
                continue
            yy = y0 + plot_h - (value - ymin) / (ymax - ymin) * plot_h
            base = y0 + plot_h - (0.0 - ymin) / (ymax - ymin) * plot_h if ymin <= 0.0 <= ymax else (y0 + plot_h)
            top_y = min(yy, base)
            height = abs(base - yy)
            parts.append(f'<rect x="{gx + offsets[sampler] - bar_w/2:.2f}" y="{top_y:.2f}" width="{bar_w:.2f}" height="{height:.2f}" fill="{sampler_colors[sampler]}"/>')
        parts.append(f'<text x="{gx:.2f}" y="{y0 + plot_h + 18}" text-anchor="middle" class="axis">{esc(scenario)}</text>')
    parts.append(f'<text x="{panel_x + panel_w/2}" y="{panel_y + panel_h - 14}" text-anchor="middle" class="axis">Stress scenario</text>')
    parts.append(f'<text x="{panel_x + 18}" y="{panel_y + panel_h/2}" text-anchor="middle" transform="rotate(-90 {panel_x + 18} {panel_y + panel_h/2})" class="axis">{esc(y_label)}</text>')

pricing_slopes = {}
final_rmse = {}
for row in pricing_rows:
    pricing_slopes[(row['scenario'], row['sampler'])] = float(row['slope']) if row['slope'] not in ('', 'nan') else None
    key = (row['scenario'], row['sampler'])
    current = final_rmse.get(key)
    if current is None or int(row['paths']) > current[0]:
        final_rmse[key] = (int(row['paths']), float(row['rmse']))
final_rmse = {k: v[1] for k, v in final_rmse.items()}

improvement = {}
severity = {}
for row in summary_rows:
    scenario = row['scenario']
    improvement[(scenario, 'MC')] = 0.0
    improvement[(scenario, 'QMC')] = float(row['qmc_minus_mc_slope'])
    improvement[(scenario, 'RQMC')] = float(row['rqmc_minus_mc_slope'])
    severity[(scenario, 'MC')] = float(row['atm_iv'])
    severity[(scenario, 'QMC')] = float(row['left_skew'])
    severity[(scenario, 'RQMC')] = float(row['rho_abs_error'])

parts = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1260" height="940" viewBox="0 0 1260 940">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 19px; font-weight: bold; }',
    '.subtitle { font-size: 12px; fill: #555; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #d0d0d0; stroke-width: 1; opacity: 0.5; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
    '<text x="630" y="24" text-anchor="middle" class="title">Rough Bergomi Stress Robustness</text>',
    '<text x="630" y="42" text-anchor="middle" class="subtitle">Extreme correlation, high vol-of-vol, and very small H show where QMC drifts away from ideal behavior</text>',
]

grouped_bars(parts, 20, 60, 1220, 270, 'Pricing convergence slope by scenario', pricing_slopes, 'Log-log slope of pricing RMSE vs paths', [-1.0, -0.5])
grouped_bars(parts, 20, 350, 1220, 250, 'Final-N pricing RMSE by scenario', final_rmse, 'RMSE at largest path count')
grouped_bars(parts, 20, 620, 600, 260, 'QMC / RQMC slope advantage over MC', improvement, 'Slope minus MC slope')
grouped_bars(parts, 640, 620, 600, 260, 'Regime severity proxies', severity, 'ATM IV / left skew / |rho error|')

lx, ly = 90, 78
for idx, sampler in enumerate(['MC', 'QMC', 'RQMC']):
    yy = ly + idx * 16
    parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{sampler_colors[sampler]}" stroke-width="2.2"/>')
    parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{sampler}</text>')

fx, fy = 980, 78
families = [('baseline', 'Baseline'), ('correlation', 'Extreme corr'), ('vol-of-vol', 'High vol-of-vol'), ('roughness', 'Very small H')]
for idx, (key, label) in enumerate(families):
    yy = fy + idx * 16
    parts.append(f'<rect x="{fx}" y="{yy - 9}" width="14" height="10" fill="{family_colors[key]}" opacity="0.55"/>')
    parts.append(f'<text x="{fx + 22}" y="{yy}" class="legend">{esc(label)}</text>')

parts.append('</svg>')
(base / 'rough_stress_robustness.svg').write_text('\n'.join(parts), encoding='utf-8')
print('Wrote', base / 'rough_stress_robustness.svg')
