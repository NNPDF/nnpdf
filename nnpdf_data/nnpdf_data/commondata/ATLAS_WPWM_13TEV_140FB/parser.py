# author: JK
import numpy as np
import yaml
import glob
import os
from nnpdf_data.filter_utils.utils import prettify_float


yaml.add_representer(float, prettify_float)

with open("sources.yaml") as f:
    definitions_yaml = yaml.safe_load(f)
definitions = definitions_yaml["definitions"]

MTW_BIN_MAP = {
    0: (200, 300),
    1: (300, 425),
    2: (425, 600),
    3: (600, 900),
    4: (900, 2000),
}
SQRTS = 13000

def symerr(delta_plus, delta_minus):
    semi_diff = (delta_plus - delta_minus) * 0.5
    average = (delta_plus + delta_minus) * 0.5
    sigma = np.sqrt(average ** 2 + 2 * semi_diff ** 2)
    return semi_diff, sigma

def parse_files_singlediff(pattern, tag):
    data, bins, uncertainties, unc_values = [], [], {}, []
    for fname in sorted(glob.glob(pattern)):
        with open(fname) as f:
            content = yaml.safe_load(f)
        mtw_bin = None
        indep = content["independent_variables"][0]["values"]
        dep = content["dependent_variables"][0]["values"]

        for mtw_bin, entry in zip(indep, dep):
            val = entry["value"]
            bins.append({
                "m_t": {
                    "min": mtw_bin["low"],
                    "mid": 0.5*(mtw_bin["low"]+mtw_bin["high"]),
                    "max": mtw_bin["high"]
                },
                "sqrts": {
                    "min": None,
                    "mid": SQRTS,
                    "max": None
                }
            })
            unc_bin = {}
            
            lumi_label = "Luminosity_error"
            if lumi_label not in uncertainties:
                uncertainties[lumi_label] = {
                    "description": "Luminosity error",
                    "treatment": "MULT",
                    "type": "ATLASLUMI15"
                }
            unc_bin[lumi_label] = 0.0083 * val
            sum_delta = 0 

            for err in entry["errors"]:
                label= err["label"]
                if label in definitions:
                    uncertainties[label] = definitions[label]
                if 'asymerror' in err:
                    minus = err['asymerror']['minus']
                    plus = err['asymerror']['plus']
                elif 'symerror' in err:
                    plus = minus = err['symerror']
                else: 
                    raise ValueError(f"No error found for {pattern}, {tag}")
            
                delta, sym_err = symerr(plus, minus) 
                unc_bin[label] = sym_err
                sum_delta += delta
                unc_values.append(unc_bin)
                
            shifted_data = val + sum_delta 
            data.append(shifted_data) 
    
    return data, bins, uncertainties, unc_values

def parse_files_doublediff(pattern, tag):
    data, bins, uncertainties, unc_values, deltas = [], [], {}, [], []
    for fname in sorted(glob.glob(pattern)):
        with open(fname) as f:
            content = yaml.safe_load(f)

        mtw_bin = None
        for k in MTW_BIN_MAP:
            if f"mtw{k}" in fname:
                mtw_bin = k
                break
        if mtw_bin is None:
            raise RuntimeError(f"Cannot determine mtw bin from {fname}")

        indep = content["independent_variables"][0]["values"]
        dep = content["dependent_variables"][0]["values"]

        mtw_low, mtw_high = MTW_BIN_MAP[mtw_bin]
        mtw_mid = 0.5 * (mtw_low + mtw_high)

        for eta_bin, entry in zip(indep, dep):
            val = entry["value"]
            

            bins.append({
                "eta": {
                    "min": eta_bin["low"],
                    "mid": 0.5 * (eta_bin["low"] + eta_bin["high"]),
                    "max": eta_bin["high"]
                },
                "m_t": {
                    "min": mtw_low,
                    "mid": mtw_mid,
                    "max": mtw_high
                },
                "sqrts": {
                    "min": None,
                    "mid": SQRTS,
                    "max": None
                }
            })
            unc_bin, delta_bin = {}, {}
            lumi_label = "Luminosity uncertainty"
            if lumi_label not in uncertainties:
                uncertainties[lumi_label] = {
                    "description": "Luminosity error",
                    "treatment": "MULT",
                    "type": "ATLASLUMI15"
                }
            unc_bin[lumi_label] = 0.0083 * val
            delta_bin[lumi_label] = 0.0

            sum_delta = 0

            for err in entry["errors"]:
                label = err["label"]
                
                if label in definitions:
                    uncertainties[label] = definitions[label]
                if 'asymerror' in err:
                    minus = err['asymerror']['minus']
                    plus = err['asymerror']['plus']
                elif 'symerror' in err:
                    plus = minus = err['symerror']
                else: 
                    raise ValueError(f"No error found for {pattern}, {tag}")
                
                delta, sym_err = symerr(plus, minus) 
                unc_bin[label] = sym_err
                sum_delta += delta
                

            unc_values.append(unc_bin)
            shifted_data = val + sum_delta
            data.append(shifted_data)

    return data, bins, uncertainties, unc_values


def mergedata(dataset):
    data, bins, uncertainties, unc_values = [], [], {}, []
    for d, b, u, uv in dataset:
        data += d
        bins += d
        unc_values += uv
        uncertainties.update(u)
    return data, bins, uncertainties, unc_values


def write_yaml(filename, content):
    with open(filename, "w") as f:
        yaml.dump(content, f, sort_keys=False)

def write_all(base, data, bins, uncertainties, unc_values):
    write_yaml(f"data_{base}.yaml", {"data_central": data})
    write_yaml(f"kinematics_{base}.yaml", {"bins": bins})
    write_yaml(f"uncertainties_{base}.yaml", {"definitions": uncertainties, "bins": unc_values})


def main():
#    ## doublediff
#    # WP
#    data, bins, uncertainties, unc_values = parse_files_doublediff("./rawdata/lep_physical_plus_absetamtw_mtw*.yaml", "WP")
#    write_yaml("data_WP_DIF_MTW-ABSETA.yaml", {"data_central": data})
#    write_yaml("kinematics_WP_DIF_MTW-ABSETA.yaml", {"bins": bins})
#    write_yaml("uncertainties_WP_DIF_MTW-ABSETA.yaml", {"definitions": uncertainties, "bins": unc_values})
#
#    ## WM
#    data, bins, uncertainties, unc_values = parse_files_doublediff("./rawdata/lep_physical_minus_absetamtw_mtw*.yaml", "WM")
#    write_yaml("data_WM_DIF_MTW-ABSETA.yaml", {"data_central": data})
#    write_yaml("kinematics_WM_DIF_MTW-ABSETA.yaml", {"bins": bins})
#    write_yaml("uncertainties_WM_DIF_MTW-ABSETA.yaml", {"definitions": uncertainties, "bins": unc_values})
#    
#    ## singlediff
#    # WP
#    data, bins, uncertainties, unc_values = parse_files_singlediff("./rawdata/lep_physical_plus_mtw.yaml", "WP")
#    write_yaml("data_WP_DIF_MTW.yaml", {"data_central": data})
#    write_yaml("kinematics_WP_DIF_MTW.yaml", {"bins": bins})
#    write_yaml("uncertainties_WP_DIF_MTW.yaml", {"definitions": uncertainties, "bins": unc_values})
#    # WM
#    data, bins, uncertainties, unc_values = parse_files_singlediff("./rawdata/lep_physical_minus_mtw.yaml", "WM")
#    write_yaml("data_WM_DIF_MTW.yaml", {"data_central": data})
#    write_yaml("kinematics_WM_DIF_MTW.yaml", {"bins": bins})
#    write_yaml("uncertainties_WM_DIF_MTW.yaml", {"definitions": uncertainties, "bins": unc_values})


    # 1D, muon + electron summed, W+ and W- in sequence
    combined_WP = parse_files_singlediff("./rawdata/lep_physical_plus_mtw.yaml", "WP")
    combined_WM = parse_files_singlediff("./rawdata/lep_physical_minus_mtw.yaml", "WM")
    combined_W = mergedata([combined_WP, combined_WM])
    write_all("WPWM_LEP_SDIF_MTW", *combined_W)    

    # 1D, muon only, W+ and W- in sequence
    muo_WP = parse_files_singlediff("./rawdata/muo_plus_mtw.yaml", "WP")
    muo_WM = parse_files_singlediff("./rawdata/muo_minus_mtw.yaml", "WM")
    muo_W = mergedata([muo_WP, muo_WM])
    write_all("WPWM_MUON_SDIF_MTW", *muo_W)

    # 2D, muon + electron summed, W+ and W- in sequence
    combined_WP_2D = parse_files_doublediff("./rawdata/lep_physical_plus_absetamtw_mtw*.yaml", "WP")
    combined_WM_2D = parse_files_doublediff("./rawdata/lep_physical_minus_absetamtw_mtw*.yaml", "WM")
    combined_W_2D = mergedata([combined_WP_2D, combined_WM_2D])
    write_all("WPWM_LEP_DDIF_MTW-ABSETA", *combined_W_2D)

    # 2D muon only, W+ and W- in sequence
    muo_WP_2D = parse_files_doublediff("./rawdata/muo_plus_absetamtw_mtw*.yaml", "WP")
    muo_WM_2D = parse_files_doublediff("./rawdata/muo_minus_absetamtw_mtw*.yaml", "WM")
    muo_W_2D = mergedata([muo_WP_2D, muo_WM_2D])
    write_all("WPWM_MUON_DDIF_MTW-ABSETA", *muo_W_2D) 

if __name__ == "__main__":
    main()



