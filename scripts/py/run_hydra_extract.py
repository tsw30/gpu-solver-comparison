"""Run Hydra extract script"""

import os
import subprocess

RUN_DIR = "/rds/project/rds-Z81G6YisvY4/tsw30/hydra/LSF5/CFD/Airmeter/" \
          "urans_ar60p5_3"
PROBES = ['L00S', 'BPIN', 'ESIN']
OUT_DIR = RUN_DIR + "/avgd"
SCRIPT_PATH = "/home/tsw30/python/gpu-solver-comparison/scripts/py/" \
              "hydra_extract.py"

def main():
    """Use subprocesses for each domain processed."""
    for probe in PROBES:
        output_file = os.path.join(OUT_DIR, probe + '_data.json')
        command = ['python ' + SCRIPT_PATH + ' --probe_name ' + probe +
                   ' --output_file ' + output_file]
        subprocess.call(command)


if __name__ == "__main__":
    main()
