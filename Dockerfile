FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update


# -----------------------------------------

# dependencies for vcf file import, plink, gemma
RUN mkdir /kb/deps
COPY ./deps /kb/deps
ENV PATH=$PATH:/kb/module/deps/bin

RUN curl -J -L https://easygwas.ethz.ch/down/dataset/download/1/ -o /kb/deps/testdata/AtPolyDB.zip

# GEMMA tools
#ARG gemmaurl=https://github.com/genetics-statistics/GEMMA/releases/download/0.98.1/gemma-0.98.1-linux-static.gz
#ARG plinkurl=http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip

#RUN curl $plinkurl -o /kb/module/deps/plink.zip
#RUN curl $gemmaurl -o /kb/module/deps/gemma.zip

# python libraries
# RUN apt-get update -y \
#     && apt-get install -y python3-minimal python3-pip

RUN pip install --upgrade pip \
	&& pip install --upgrade requests \
    && pip install -q pyvcf

RUN mkdir /data

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
