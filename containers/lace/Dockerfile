# Custom Lace 1.14.1 with NetworkX 2.x fix
FROM docker.io/mambaorg/micromamba:1.5.8

# Install Lace with NetworkX pinned to 2.x (not 3.x which breaks .node attribute)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    lace=1.14.1 \
    "networkx>=2.8,<3.0" \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
