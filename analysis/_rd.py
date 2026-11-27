import pyreadr, sys
for p in ["/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/data/topolska/ext/pre_processed_data/scaled_variants_aPCA.rds",
          "/mnt/c/Users/armin/Desktop/AXZ/insertion-lab/data/topolska/ext/additional_files/tsuboyama_sec_struc.rds"]:
    print("=====", p.split("/")[-1])
    try:
        r = pyreadr.read_r(p)
        for k, v in r.items():
            print(" key:", k, "shape:", getattr(v, "shape", None))
            print(" cols:", list(v.columns)[:40])
            print(v.head(3).to_string()[:1500])
    except Exception as e:
        print(" ERROR", type(e).__name__, e)
