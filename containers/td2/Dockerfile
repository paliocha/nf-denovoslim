FROM mambaorg/micromamba:2.5.0

# Install Python + Perl via conda
RUN micromamba install -y -n base -c conda-forge \
       python=3.11 \
       perl \
    && micromamba clean --all --yes

# Enable conda env for subsequent RUN steps
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install TD2 from PyPI (includes PSAURON)
RUN pip install --no-cache-dir TD2

# Verify
RUN TD2.LongOrfs --help && TD2.Predict --help

LABEL maintainer="martin.paliocha@nmbu.no"
LABEL description="TD2 (TransDecoder 2) with PSAURON for ORF prediction"
