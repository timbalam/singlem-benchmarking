
import tomllib
import polars as pl
import os
import logging
import subprocess
from ruamel.yaml import YAML

def load_community_description(file):
    toml = tomllib.loads(file)

    return CommunityDescription(toml)

class CommunityDescription:
    def __init__(self, *, toml):
        self.samples = get_samples(toml)
        
        self.truth = get_truth(toml)

        #optional
        self.annotated_genomes = get_genomes(toml)
        
        self.readsim_tool = get_readsim_tool(toml)

    def process(self, *, prefix, snakemake_args):  
        output_config = os.path.join(prefix, 'config.yaml')

        conf = self.readsim_tool.config(
            reads1 = self.reads1,
            reads2 = self.reads2,
            coverage_file = self.truth,
            genomes_file = self.annotated_genomes,
            threads = 8
        )
        workflow = self.readsim_tool.workflow
        
        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False    

        with open(output_config, "w") as f:
            yaml.dump(conf, f)
        logging.info(f"Configuration file written to {output_config}")

        cmd = (
            "snakemake --snakefile {snakefile} --directory {working_dir} "
            "--rerun-incomplete --keep-going "
            "--configfile {config_file} --nolock "
            "{snakemake_args} "
            "{target_rule}"
        ).format(
            snakefile=get_snakefile(),
            working_dir=prefix,
            config_file=output_config,
            snakemake_args=snakemake_args,
            target_rule=workflow
        )

        logging.debug(f"Command: {cmd}")
        logging.info("Executing: %s" % cmd)
        subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1
        )

def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def get_samples(toml):

    try:
        conf_samples = toml["samples"]
        return pl.DataFrame({
            sample: conf_samples.keys(),
            path: conf_samples.values()
        })
    except KeyError:
        raise Exception("'samples' missing")

def get_truth(toml):
    try:
        truth = toml["truth"]
        return truth
    except KeyError:
        raise Exception("'truth' missing")

def get_genomes(toml):
    toml.get("genomes_file")
   
def get_readsim_tool(toml):
    try:
        conf_readsim = toml["readsim"]
    except KeyError:
        return None
    
    try:
        conf_readsim_tool = conf_readsim["tool"]
    except KeyError:
        raise Exception("'readsim.tool' missing")
    
    if conf_readsim_tool == "art":
        try:
            conf_readsim_art_bin = conf_readsim["bin"]
        except KeyError:
            raise Exception("'readsim.bin' missing")

        return ArtSimTool(
            read_length = toml.get("read_length"),
            bin = conf_readsim_art_bin
        )
    else:
        raise Exception(f"Unknown 'readsim.tool' option: {conf_readsim_tool}")

class ArtSimTool:
    def __init__(self, *, read_length, binary):
        self.read_length = read_length
        self.binary = binary
        self.workflow = "simulate_art_reads"

    def config(self, *, reads1, reads2, coverage_file, genomes_file, threads):
        return {
            reads1: reads1,
            reads2: reads2,
            coverage_file: truth,
            genomes_file: annotated_genomes,
            threads: threads,
            read_length: self.read_length,
            art_bin: self.binary
        }