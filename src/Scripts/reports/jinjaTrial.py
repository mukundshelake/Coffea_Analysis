from jinja2 import Environment, FileSystemLoader

# Set up the environment
env = Environment(loader=FileSystemLoader('.'))
template = env.get_template('trial.tex')

# Data to fill in the template
data = {
    'name': 'John Doe',
    'age': 30
}

# Render the template with data
output = template.render(data)

# Save the rendered LaTeX code
with open('output.tex', 'w') as f:
    f.write(output)
