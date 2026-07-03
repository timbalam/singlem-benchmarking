
import tomllib
import polars as pl
import os
import logging
import subprocess
import shutil
import sys
from ruamel.yaml import YAML

def load_community_description(file):
    with open(file, "rb") as f:
        toml = tomllib.load(f)

    dir = os.path.dirname(os.path.abspath(file))
    return CommunityDescription(
        samples = get_samples(toml, dir),
        truth = get_truth(toml, dir),
        taxonomy = get_taxonomy(toml, dir),
        genomes_list = get_genomes_list(toml, dir),
        readsim_tool = get_readsim_tool(toml, dir),
        gtdbtk = get_gtdbtk(toml, dir)
    )

def make_absolute(dir, *paths):
    return os.path.normpath(os.path.join(dir, *paths))

class CommunityDescription:
    def __init__(self, *, samples, truth, taxonomy = None, genomes_list = None,
                 readsim_tool = None, gtdbtk = None):
        self.samples = samples
        self.truth = truth
        self.taxonomy = taxonomy
        self.genomes_list = genomes_list
        self.readsim_tool = readsim_tool
        self.gtdbtk = gtdbtk

    def process(self, *, prefix, cores = 8, snakemake_args):  
        output_config = os.path.join(prefix, 'config.yaml')

        conf, workflow = get_config_and_workflow(
            samples = self.samples,
            readsim_tool = self.readsim_tool,
            coverage_file = self.truth,
            taxonomy = self.taxonomy,
            genomes_list = self.genomes_list,
            gtdbtk_data = self.gtdbtk_data if self.gtdbk is not None,
            gtdbtk_dir = self.gtdbtk_dir if self.gtdbk is not None,
            gtdbtk_release = self.gtdbtk_release if self.gtdbk is not None,
            threads = cores
        )
        
        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False    

        with open(output_config, "w") as f:
            yaml.dump(conf, f)
        logging.info(f"Configuration file written to {output_config}")

        cmd = (
            "{snakemake} --snakefile {snakefile} --directory {prefix} "
            "--rerun-incomplete --keep-going "
            "--configfile {config_file} --nolock "
            "--cores {cores} "
            "{snakemake_args} "
            "{workflow}"
        ).format(
            snakemake = shutil.which("snakemake"),
            snakefile = get_snakefile(),
            prefix = prefix,
            cores = cores,
            config_file = output_config,
            snakemake_args = snakemake_args,
            workflow = workflow
        )

        logging.debug(f"Command: {cmd}")
        logging.info("Executing: %s" % cmd)
        proc = subprocess.Popen(
            cmd,
            shell=True,
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
            logging.info("Finished: %s" % workflow)
        else:
            sys.exit(1)


def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

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
    
    def config_and_workflow(self, *, readsim_tool, **args):
        return readsim_tool.paired_config_and_workflow(samples = self.samples,
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
    
    return GtdbTKAssignTaxonomy(
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

def get_config_and_workflow(*, samples, **args):
    return samples.config_and_workflow(**args)

class ArtSimTool:
    def __init__(self, *, read_length, bin):
        self.read_length = read_length
        self.bin = bin

    def paired_config_and_workflow(self, *, samples, reads1, reads2,
                                   coverage_file, genomes_list,
                                   taxonomy, gtdbtk_data, gtdbtk_dir,
                                   gtdbtk_release,
                                   threads):
        workflow = "simulate_art_paired_reads"
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
        return config, workflow

class InvalidCommunityDescription(Exception):
    pass