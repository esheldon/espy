import os

def get_config_dir():
    dir=os.environ['ESPY_DIR']
    return os.path.join(dir, 'gmix_sdss', 'config')

def get_config_path(run):
    dir=get_config_dir()
    return os.path.join(dir, '%s.yaml' % run)

def read_config(run):
    import yaml

    path=get_config_path(run)
    conf=yaml.load(open(path))

    if conf['run'] != run:
        mess="run does not match itself in config: %s instead of  %s"
        mess=mess % (conf['run'],run)
        raise ValueError(mess)

    return conf


