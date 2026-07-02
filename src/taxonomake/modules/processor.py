
import os
import logging
import subprocess
from ruamel.yaml import YAML


def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def process(config, prefix):
    
    output = os.path.join(prefix, 'config.yaml')
    
    yaml = YAML()
    yaml.version = (1, 1)
    yaml.default_flow_style = False    

    with open(output, "w") as f:
        yaml.dump(config, f)
    logging.info(f"Configuration file written to {output}")

    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir} "
        "{jobs} --rerun-incomplete --keep-going "
        "--configfile {config_file} --nolock "
        "{profile} "
        "{target_rule}"
    ).format(
        snakefile=get_snakefile(),
        working_dir=output,
        jobs="--cores {}".format(cores) if cores is not None else "--jobs 1",
        config_file=config,
        profile="" if not profile else "--profile {}".format(profile),
        target_rule=workflow if workflow != "None" else ""
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
        bufsize=1,
    )

    


