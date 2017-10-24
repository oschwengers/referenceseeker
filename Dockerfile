
FROM ubuntu:xenial

LABEL version "1.0"
LABEL maintainer "oliver.schwengers@computational.bio.uni-giessen.de"

RUN apt-get -y update && apt-get -y install \
    libidn11\
    python3\
    python3-pip
RUN pip3 install biopython

COPY refseekr /
COPY share/ /share

ENV REFERENCE_SEEKER_HOME=/

WORKDIR /data

ENTRYPOINT [ "python3", "/refseekr.py", "--db", "/db" ]
CMD ["--help"]
