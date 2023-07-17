#!/bin/bash

input="$1"
model="$2"
outdir="$3"
export CUDA_VISIBLE_DEVICES=0

python $SMSNET_M --model_dir $model --inference_input_file $input
