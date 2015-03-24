FROM ubuntu

RUN apt-get update
RUN apt-get install -y gfortran build-essential make gcc build-essential python python-dev wget libgsl0ldbl gsl-bin libgsl0-dev python-pip

# blas
ADD blas.sh /tmp/blas.sh
RUN chmod 755 /tmp/blas.sh
RUN /tmp/blas.sh
ENV BLAS /usr/local/lib/libfblas.a

# lapack
ADD lapack.sh /tmp/lapack.sh
RUN chmod 755 /tmp/lapack.sh
RUN /tmp/lapack.sh
ENV LAPACK /usr/local/lib/liblapack.a

RUN pip install numpy scipy
RUN easy_install -U ete2

WORKDIR /opt

ADD * /opt/

RUN g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`
RUN chmod a+rwx ../opt

#
# Download and unpack Mutect
#
#RUN wget https://github.com/broadinstitute/mutect/releases/download/1.1.5/muTect-1.1.5-bin.zip
#RUN wget https://github.com/broadinstitute/picard/releases/download/1.122/picard-tools-1.122.zip

#RUN unzip muTect-1.1.5-bin.zip
#RUN unzip picard-tools-1.122.zip

# Link the picard tools to /opt/picard
#RUN ln -s picard-tools-1.122 picard