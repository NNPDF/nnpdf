# author: JK
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


def parse_files_doublediff(pattern, tag):
    data, bins, uncertainties, unc_values = [], [], {}, []
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
            data.append(val)
            

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
            unc_bin = {}
            lumi_label = "Luminosity uncertainty"
            if lumi_label not in uncertainties:
                uncertainties[lumi_label] = {
                    "description": "Luminosity error",
                    "treatment": "MULT",
                    "type": "ATLASLUMI15"
                }
            unc_bin[lumi_label] = 0.0083 * val


            for err in entry["errors"]:
                label = err["label"]
                
                if label in definitions:
                    uncertainties[label] = definitions[label]
                unc_bin[label] = err["symerror"]

            unc_values.append(unc_bin)

    return data, bins, uncertainties, unc_values

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
            data.append(val)
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
            
            for err in entry["errors"]:
                label= err["label"]
                if label in definitions:
                    uncertainties[label] = definitions[label]
                unc_bin[label] = err["symerror"]
            unc_values.append(unc_bin)
    
    
    return data, bins, uncertainties, unc_values


def write_yaml(filename, content):
    with open(filename, "w") as f:
        yaml.dump(content, f, sort_keys=False)


def main():
    ## doublediff
    # WP
    data, bins, uncertainties, unc_values = parse_files_doublediff("./rawdata/lep_physical_plus_absetamtw_mtw*.yaml", "WP")
    write_yaml("data_WP_DIF_MTW-ABSETA.yaml", {"data_central": data})
    write_yaml("kinematics_WP_DIF_MTW-ABSETA.yaml", {"bins": bins})
    write_yaml("uncertainties_WP_DIF_MTW-ABSETA.yaml", {"definitions": uncertainties, "bins": unc_values})

    ## WM
    data, bins, uncertainties, unc_values = parse_files_doublediff("./rawdata/lep_physical_minus_absetamtw_mtw*.yaml", "WM")
    write_yaml("data_WM_DIF_MTW-ABSETA.yaml", {"data_central": data})
    write_yaml("kinematics_WM_DIF_MTW-ABSETA.yaml", {"bins": bins})
    write_yaml("uncertainties_WM_DIF_MTW-ABSETA.yaml", {"definitions": uncertainties, "bins": unc_values})
    
    ## singlediff
    # WP
    data, bins, uncertainties, unc_values = parse_files_singlediff("./rawdata/lep_physical_plus_mtw.yaml", "WP")
    write_yaml("data_WP_DIF_MTW.yaml", {"data_central": data})
    write_yaml("kinematics_WP_DIF_MTW.yaml", {"bins": bins})
    write_yaml("uncertainties_WP_DIF_MTW.yaml", {"definitions": uncertainties, "bins": unc_values})
    # WM
    data, bins, uncertainties, unc_values = parse_files_singlediff("./rawdata/lep_physical_minus_mtw.yaml", "WM")
    write_yaml("data_WM_DIF_MTW.yaml", {"data_central": data})
    write_yaml("kinematics_WM_DIF_MTW.yaml", {"bins": bins})
    write_yaml("uncertainties_WM_DIF_MTW.yaml", {"definitions": uncertainties, "bins": unc_values})

if __name__ == "__main__":
    main()



