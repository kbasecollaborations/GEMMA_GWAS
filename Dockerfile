FROM kbase/sdkbase2:python
MAINTAINER UT Plant Sciences
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

RUN apt-get update -y \
    && apt-get install -y vim \
	&& apt-get install -y bcftools \
	&& apt-get install -y time

RUN pip install --upgrade pip \
	&& pip install --upgrade requests \
	&& pip install pandas \
    && pip install -q pyvcf

RUN mkdir /data

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
