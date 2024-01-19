FROM papaemmelab/docker-hmftools:v1.0.0

# Install dependencies
RUN apt-get update && apt-get install -y tar curl python3

# Install Google Cloud SDK
RUN curl -sSL https://sdk.cloud.google.com | bash
ENV PATH $PATH:/root/google-cloud-sdk/bin

# Set version as an environment variable
ENV GENOME_VERSION=37
ENV HMFTOOLS_VERSION=5.28
ENV HMFTOOLS_VERSION_UNDERSCORE 5_28
ENV REF_DIR=hmf_pipeline_resources.${GENOME_VERSION}_v${HMFTOOLS_VERSION}

# Download the file using gsutil
RUN \
    mkdir -p /data && \
    gsutil cp -r gs://hmf-public/HMFtools-Resources/dna_pipeline/v${HMFTOOLS_VERSION_UNDERSCORE}/${REF_DIR}.gz /data/ && \
    tar -xzf /data/${REF_DIR}.gz -C /data && \
    rm /data/${REF_DIR}.gz

COPY . /app
WORKDIR /app

ENTRYPOINT ["bash"]