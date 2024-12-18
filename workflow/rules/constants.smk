INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
BMK_DIR = config.get("benchmarks_dir", "benchmarks")
LOG_DIR = config.get("logs_dir", "logs")
HMM_PROFILE = config["hmm_profile"]
RGX_CHR = re.compile(r"(chr[0-9XY]+)")
IDENT_THR = config.get("identity_thr", 90.0)
