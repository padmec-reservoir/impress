import yaml

with open('variable_settings.yml', 'r') as f:
    config = yaml.load(f)

nodes = config['nodes']
edges = config['edges']
faces = config['faces']
volumes = config['volumes']
pressure = config['pressure']
