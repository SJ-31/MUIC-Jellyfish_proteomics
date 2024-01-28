#!/usr/bin/env python

import anc2vec.train as builder
import numpy as np

go_obo = "/mnt/data/shannc/nf/data/reference/go.obo"
save_to = "/mnt/data/shannc/nf/data/reference/go_embedded.npz"
es = builder.fit(
    go_obo,
    embedding_sz=200,
    batch_sz=64,
    num_epochs=100,
)
np.save(save_to, embds=es)

# to load do
# np.load(save_to, allow_pickle=True)["embds"].item() will give you the dictionary
