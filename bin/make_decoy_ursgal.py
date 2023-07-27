#!/usr/bin/env python
import sys
import ursgal

target = sys.argv[1]
output = sys.argv[2]

params = {
    "enzyme": "trypsin",
    "decoy_generation_mode": "reverse_protein"
}
uc = ursgal.UController(params=params)
uc.execute_misc_engine(
    input_file=target,
    engine="generate_target_decoy_1_0_0",
    output_file_name=output
)
