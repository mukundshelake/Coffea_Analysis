from jinja2 import Environment, FileSystemLoader
from coffea.util import load, save
# Configure Jinja with custom block start and end strings
env = Environment(loader=FileSystemLoader('.'), variable_start_string='((*', variable_end_string='*))', block_start_string='((%', block_end_string='%))', comment_start_string='<!--', comment_end_string='-->')

template = env.get_template('unfilled.tex')

# Data to fill in the template
data = {}
for era in ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']:
    for region in ['A', 'B', 'C', 'D']:
        for param in ['Fu', 'Fd', 'Fu_ud', 'Fd_ud', 'Du', 'Dd', 'FuDu', 'FdDd', 'Fu_vs_Fd', 'Du_vs_Dd', 'Fu_ud_vs_Fd_ud', 'FuDu_vs_FdDd']:
            keyName = f'{param}_{era}_{region}'
            valueName = f'../plots/paramDistributions/{era}/region_{region}/{param}.png'
            data[keyName] = valueName

out = load('../outputs/Fq_results.coffea')
for era in ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']:
    for param in ['Fu', 'Fd', 'Fu_ud', 'Fd_ud', 'Du', 'Dd', 'FuDu', 'FdDd']:
        key_withFlow = f"{param}_{era}_withFlow"
        value_withFlow = round(out[era]['Full_y0_range']['with_flow'][param],4)
        data[key_withFlow] = value_withFlow
        key_withoutFlow = f"{param}_{era}_withoutFlow"
        value_withoutFlow = round(out[era]['Full_y0_range']['without_flow'][param],4)
        data[key_withoutFlow] = value_withoutFlow


# Render the template with data
output = template.render(data)

# Save the rendered LaTeX code
with open('filled.tex', 'w') as f:
    f.write(output)
