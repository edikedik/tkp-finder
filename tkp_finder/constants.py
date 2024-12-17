from pathlib import Path

PSKD_NAME = "PSKD"
HIT_SEP = "~"
MSA_PATH = Path(__file__).parent / "resources" / "sbmsa_latest.fasta"
PLANT_HMM_URL = (
    "https://raw.githubusercontent.com/edikedik/tkp-finder/master"
    "/Appendix_4/Plant_Pkinase_fam.hmm"
)
PLANT_HMMS_PATH = Path(__file__).parent / "resources" / "plant_hmms"
RESOURCES_PATH = Path(__file__).parent / "resources"


if __name__ == "__main__":
    raise RuntimeError
