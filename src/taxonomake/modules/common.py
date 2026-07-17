import os.path

def _make_absolute(dir, *paths):
    return os.path.normpath(os.path.join(dir, *paths))

def _make_absolute_if_path(dir, path):
    return _make_absolute(dir, path) if os.path.dirname(path) != "" else path

def config_dir(config):
    return os.path.dirname(config["configfilepath"])

def config_sample_reads1(config):
    return [_make_absolute(config_dir(config), s) for s in config["samples"]["reads1"]]

def config_sample_reads2(config):
    return [_make_absolute(config_dir(config), s) for s in config["samples"]["reads2"]]

def config_sample_names(config):
    return config["samples"]["names"]

def config_truth(config):
    return _make_absolute(config_dir(config), config["truth"])

def config_genomes_list(config):
    return _make_absolute(config_dir(config), config["genomes_list"])

def config_taxonomy(config):
    return _make_absolute(config_dir(config), config["taxonomy"])

def config_readsim_bin(config):
    return _make_absolute_if_path(config_dir(config), config["readsim"]["binary"])

def config_gtdbtk_dir(config, *paths):
    return _make_absolute(config_dir(config), config["gtdbtk"]["directory"], *paths)

def config_gtdbtk_data(config):
    return _make_absolute(config_dir(config), config["gtdbtk"]["data"])