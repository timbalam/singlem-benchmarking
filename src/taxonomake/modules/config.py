
import tomllib

def load_config(file):
    config = tomllib.loads(file)

    return Config(config)

class Config:
    def __init__(self, config):
        self.samples = {s:ReadsFile(p) for s, p in config["samples"].items()}
        self.truths = {s:ProfileFile(p) for s, p in config["truth"].items()}

        # optional
        self.annotated_genomes = {s:SequenceFile(p) for s, p in config["genomes"].items()}
        self.readsim_tool = ReadSimTool(config["readsim"])


class ReadsFile:
    def __init__(self, path):
        self.path = path

class ProfileFile:
    def __init__(self, path):
        self.path = path

class SequenceFile:
    def __init__(self, path):
        self.path = path

class ReadSimTool:
    def __init__(self, read_sim_options):
        self._read_sim_options = read_sim_options