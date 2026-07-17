
import os
import logging
import subprocess
import shutil
import sys
from ruamel.yaml import YAML
import shlex
from snakemake.utils import validate

SNAKEFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Snakefile")
SCHEMA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.schema.yaml")

def process_community_description(file, *, prefix, cores = 8, snakemake_args):
    yaml = YAML()
    with open(file) as f:
        conf = yaml.load(f)
    validate(conf, SCHEMA)
    cmd = [
        shutil.which("snakemake"),
        "--snakefile",
        SNAKEFILE,
        "--directory",
        f"{prefix}",
        "--rerun-incomplete",
        "--keep-going",
        "--configfile",
        f"{file}", 
        "--nolock",
        "--cores",
        f"{cores}",
        "--config",
        f"configfilepath={os.path.abspath(file)}"
    ]
    if snakemake_args != "":
        cmd += snakemake_args

    logging.debug(f"Command: {shlex.join(cmd)}")
    logging.info("Executing: %s" % shlex.join(cmd))
    proc = subprocess.Popen(
        cmd,
        #shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1
    )
    
    proc.wait()

    for line in proc.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()
    
    for line in proc.stderr:
        sys.stderr.write(line)
        sys.stderr.flush()
    
    if proc.returncode == 0:
        logging.info("Finished")
    else:
        sys.exit(1)

def get_samples(toml, dir):
    try:
        conf_samples = toml["samples"]
    except KeyError:
        raise Exception("'samples' missing")

    try:
        conf_sample_names = conf_samples["names"]
    except KeyError:
        raise InvalidCommunityDescription("'samples.names' missing")
    
    try:
        conf_sample_reads1 = [make_absolute(dir, p) for p in conf_samples["reads1"]]
        conf_sample_reads2 = [make_absolute(dir, p) for p in conf_samples["reads2"]]
        
        return PairedSamples(
            samples = conf_sample_names,
            reads1 = conf_sample_reads1,
            reads2 = conf_sample_reads2
        )
    except KeyError:
        raise InvalidCommunityDescription("'samples.reads1' or 'samples.reads2' missing")

class PairedSamples:
    def __init__(self, *, samples, reads1, reads2):
        self.samples = samples
        self.reads1 = reads1
        self.reads2 = reads2
    
    def config(self, *, readsim_tool, **args):
        return readsim_tool.paired_config(samples = self.samples,
                                          reads1 = self.reads1,
                                          reads2 = self.reads2,
                                          **args)

def get_truth(toml, dir):
    try:
        truth = toml["truth"]
        return make_absolute(dir, truth)
    except KeyError:
        raise InvalidCommunityDescription("'truth' missing")

def get_genomes_list(toml, dir):
    try:
        genomes_list = toml["genomes_list"]
        return make_absolute(dir, genomes_list)
    except KeyError:
        return None

def get_taxonomy(toml, dir):
    try:
        taxonomy = toml["taxonomy"]
        return make_absolute(dir, taxonomy)
    except KeyError:
        raise InvalidCommunityDescription("'taxonomy' missing.")
    
def get_gtdbtk(toml, dir):
    try:
        conf_gtdbtk = toml["gtdbtk"]
    except KeyError:
        return None
    
    try:
        conf_gtdbtk_dir = conf_gtdbtk["dir"]
        conf_gtdbtk_dir = make_absolute(conf_gtdbtk_dir, dir)
    except KeyError:
        raise InvalidCommunityDescription("'gtdbtk.dir' missing.")
    
    conf_gtdbtk_release = conf_gtdbtk.get("release")
    try:
        conf_gtdbtk_data = conf_gtdbtk["data"]
        conf_gtdbtk_data = make_absolute(config_gtdbtk_data, dir)
    except KeyError:
        conf_gtdbtk_data = None
    
    return GtdbTkAssignTaxonomy(
        dir = conf_gtdbtk_dir,
        release = conf_gtdbtk_release,
        data = conf_gtdbtk_data
    )

class GtdbTkAssignTaxonomy:
    def __init__(self, *, dir, release, data):
        self.dir = dir
        self.release = release
        self.data = data
   
def get_readsim_tool(toml, dir):
    try:
        conf_readsim = toml["readsim"]
    except KeyError:
        return None
    
    try:
        conf_readsim_tool = conf_readsim["tool"]
    except KeyError:
        raise InvalidCommunityDescription("'readsim.tool' missing")
    
    if conf_readsim_tool == "art":
        try:
            conf_readsim_art_bin = conf_readsim["bin"]
        except KeyError:
            raise InvalidCommunityDescription("'readsim.bin' missing")

        return ArtSimTool(
            read_length = conf_readsim.get("read_length"),
            bin = make_absolute(dir, conf_readsim_art_bin)
                if is_path(conf_readsim_art_bin)
                else conf_readsim_art_bin
        )
    else:
        raise InvalidCommunityDescription(f"Unknown 'readsim.tool' option: {conf_readsim_tool}")

def is_path(name):
    return os.path.dirname(name) != ""

def get_config(*, samples, **args):
    return samples.config(**args)

class ArtSimTool:
    def __init__(self, *, read_length, bin):
        self.read_length = read_length
        self.bin = bin

    def paired_config(self, *, samples, reads1, reads2,
                      coverage_file, genomes_list,
                      taxonomy, gtdbtk_data, gtdbtk_dir,
                      gtdbtk_release, threads):
        config = {
            "samples": samples,
            "reads1": reads1,
            "reads2": reads2,
            "coverage_file": coverage_file,
            "genomes_list": genomes_list,
            "taxonomy": taxonomy,
            "threads": threads,
            "read_length": self.read_length,
            "art_bin": self.bin,
            "gtdbtk_data": gtdbtk_data,
            "gtdbtk_dir": gtdbtk_dir,
            "gtdbtk_release": gtdbtk_release
        }
        return config

class InvalidCommunityDescription(Exception):
    pass