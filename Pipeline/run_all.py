# Usage:
# python3 run_all.py
# make sure put some reads in Reads folder and define correspondent sample IDs in config.yaml
# also specify memory and threads

import os.path
import subprocess

# running from external directory workaround
file_path = os.path.realpath(__file__)
snakefile = os.path.dirname(file_path) + "/full_pipeline.smk"

def main():
    # runs snakemake pipeline
    subprocess.call(f"snakemake -s {snakefile} --nolock", shell=True)

if __name__ == '__main__':
    main()
