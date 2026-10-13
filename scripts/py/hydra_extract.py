"""Extract radial profile and averaged data from a hydra simulation output
file.
"""

import numpy as np
from typing import Union, cast
import argparse
import os
import json
from vtk_utils.hydra.io import HydraVTKReader
from vtk_utils.proc.general import calc_streamwise_cfr
from vtk_utils.utils.general import extract_volumes
from vtk_utils.proc.surface import SurfaceWrapper, isoSurface, average
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
SAMPLED_VARIABLES = ["p0", "T0", "rho", "vx"]
AVG_VARIABLES = ["p0", "T0"]
PROBE__ZONE = {"L00S": "L01R-C", "ESIN": "L01S-C", "BPIN": "L01B-C"}


def load_measurement_coords(file_path: Union[str, None] = None) -> dict:
    """Load measurement coordinates from json file."""
    if file_path is None:
        file_path = os.path.join(SCRIPT_PATH, "..", "measurement_coords.json")
    with open(file_path, "r") as f:
        coords_dict = json.load(f)
    return coords_dict


def load_domain(
        inpt_file: str, zone: Union[str, None] = None) -> vtkMultiBlockDataSet:
    """Load domain flow solution."""
    if not os.path.isfile(inpt_file):
        raise FileNotFoundError(f"Input file {inpt_file} not found.")
    reader = HydraVTKReader(inpt_file)
    domain = reader.read(zones=[zone] if zone is not None else None,
                         read_mean_flow=True,
                         read_probes=False
                         )
    return domain


def extract_radial_profile(
        domain: vtkMultiBlockDataSet,
        coords_dict: dict,
        probe_name: str,
        num_points: int = 101) -> SurfaceWrapper:
    """Extract radial profile from domain at given probe location."""
    zone = cast(vtkMultiBlockDataSet, domain.GetBlock(0))
    volume = extract_volumes(zone)
    if probe_name == "L01R":
        iso_var = "x"
        iso_value = np.mean(coords_dict[probe_name], axis=0)[0]
    else:
        iso_var = "Streamwise CFR"
        iso_value = -0.025
        calc_streamwise_cfr(zone=zone, store=True)
    surface = isoSurface(volume, var=iso_var, value=iso_value)()
    averaged_surface = average(surface)
    mo_profile = averaged_surface.mixed_out_average(
        dim="pitchwise", nDivs=num_points-1, distribution="cosine")
    return mo_profile


def area_weights(coords: np.ndarray) -> np.ndarray:
    """Calculate area weights for given coordinates."""
    r = coords[:, 1]
    weights = np.zeros_like(r)
    weights[0] = 0.25 * (r[0] + r[1]) ** 2 - r[0] ** 2
    weights[-1] = r[-1] ** 2 - 0.25 * (r[-1] + r[-2]) ** 2
    weights[1:-1] = (0.25 * (r[2:] + r[1:-1]) ** 2 -
                     0.25 * (r[1:-1] + r[:-2]) ** 2)
    return weights


def rig_sampling(
        coords_dict: dict, probe_name: str,
        profile: SurfaceWrapper) -> dict:
    """Sample rig coordinates from radial profile, based on span."""
    profile.calculator(variables=SAMPLED_VARIABLES + ["r"], store=True)
    r_profile = profile.get_cell_centred_values("r")
    span_profile = ((r_profile - r_profile.min()) /
                    (r_profile.max() - r_profile.min()))
    r_rig = np.array(coords_dict[probe_name])[:, 1]
    span_rig = (r_rig - r_rig.min()) / (r_rig.max() - r_rig.min())

    sampled_data = {"r": r_rig}
    for var in SAMPLED_VARIABLES:
        var_profile = profile.get_cell_centred_values(var)
        sampled_data[var] = np.interp(span_rig, span_profile, var_profile)
        # Extrapolate endwall points from real measurments
        if probe_name == "L00S":
            sampled_data[var][0] = sampled_data[var][1]
            sampled_data[var][-1] = sampled_data[var][-2]
        else:
            sampled_data[var][0] = sampled_data[var][1] + \
                (r_rig[0] - r_rig[1]) * (sampled_data[var][2] -
                                         sampled_data[var][1]) / \
                (r_rig[2] - r_rig[1])
            sampled_data[var][-1] = sampled_data[var][-2] + \
                (r_rig[-1] - r_rig[-2]) * (sampled_data[var][-2] -
                                           sampled_data[var][-3]) / \
                (r_rig[-2] - r_rig[-3])

    return sampled_data


def rig_average(
        sampled_data: dict, weights: np.ndarray) -> tuple[dict, dict]:
    """Calculate area-weighted average of sampled data."""

    averaged_data = {}
    mass_averaged_data = {}
    rhovx = sampled_data["rho"] * sampled_data["vx"]
    mass_weights = weights * rhovx
    for var in AVG_VARIABLES:
        data = sampled_data[var]
        averaged_data[var] = np.average(data, weights=weights)
        mass_averaged_data[var] = np.average(data, weights=mass_weights)

    return averaged_data, mass_averaged_data


def process_probe(
        inpt_file: str, coords_dict: dict,
        probe_name: str) -> tuple[dict, dict, dict, dict]:
    """Process a single probe: extract radial profile, sample rig coordinates,
    and calculate area-weighted averages."""

    zone = PROBE__ZONE[probe_name]
    domain = load_domain(inpt_file, zone=zone)
    profile = extract_radial_profile(domain, coords_dict, probe_name)
    sampled_data = rig_sampling(coords_dict, probe_name, profile)
    weights = area_weights(np.array(coords_dict[probe_name]))
    averaged_data, mass_averaged_data = rig_average(sampled_data, weights)

    # Extract full CFD MMMO values for AVG_VARIABLES
    cfd_averaged_data = {}
    cfd_avg = average(profile)
    for var in AVG_VARIABLES:
        v = cfd_avg.mass_average(var, dim="full")
        print(type(v[0]), type(v[1]))
        cfd_averaged_data[var] = v[0]
        cfd_averaged_data["mass_flow"] = np.absolute(np.sum(v[1]))

    return averaged_data, mass_averaged_data, sampled_data, cfd_averaged_data


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract radial profile and averaged data from a hydra"
                    " simulation output file.")
    parser.add_argument("inpt_file", type=str, help="Path to the input VTK"
                        " file.")
    parser.add_argument("--num_points", type=int, default=101,
                        help="Number of points in the radial profile "
                             "(default: 101).")
    parser.add_argument("--coords_file", type=str, default=None,
                        help="Path to the measurement coordinates JSON file "
                             "(default: None).")
    parser.add_argument("--output_file", type=str,
                        default="./extracted_data.json",
                        help="Path to the output JSON file "
                        "(default: ./extracted_data.json).")
    parser.add_argument("--probe_name", type=str, default="BPIN",
                        help="Name of the probe to process (default: BPIN).")
    return parser.parse_args()


def main():
    """Main function to extract radial profile and averaged data from a hydra
    simulation output file."""
    args = parse_args()
    coords_dict = load_measurement_coords(args.coords_file)

    results = {}
    sampled = {}
    if args.probe_name:
        probe_names = [args.probe_name]
    else:
        raise ValueError("Valid probe name must be specified.")

    for probe_name in probe_names:
        (
            averaged_data,
            mass_averaged_data,
            sampled_data,
            cfd_averaged_data,
        ) = process_probe(args.inpt_file, coords_dict, probe_name)
        results[probe_name] = {
            "averaged_data": averaged_data,
            "mass_averaged_data": mass_averaged_data,
            "cfd_MMMO_data": cfd_averaged_data
        }
        for k, v in sampled_data.items():
            if isinstance(v, np.ndarray):
                sampled_data[k] = v.tolist()
        sampled[probe_name] = sampled_data

    # Save results to JSON file
    output_dir = os.path.dirname(args.output_file)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    output_base = os.path.splitext(args.output_file)
    avgd_file = output_base[0] + ".json"
    radial_file = output_base[0] + "_radial.json"
    with open(avgd_file, "w") as f:
        json.dump(results, f, indent=4)
    with open(radial_file, "w") as f:
        json.dump(sampled, f, indent=4)

    print(f"Results saved to {output_base[0]}")


if __name__ == "__main__":
    main()
