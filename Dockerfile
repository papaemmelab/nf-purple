FROM nextflow/nextflow

# Install dependencies
RUN yum install -y tar

# Install https://github.com/askimed/nf-test
RUN mkdir -p /opt/bin && \
    cd /opt/bin && \
    curl -fsSL https://code.askimed.com/install/nf-test | bash && \
    export PATH=/opt/bin:$PATH && \
    nf-test version

COPY . /app
WORKDIR /app

ENTRYPOINT ["nextflow", "-version"]
