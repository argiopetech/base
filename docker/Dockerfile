FROM alpine as build

WORKDIR /base

RUN buildDeps="gcc \
               g++ \
               git \
               cmake \
               gsl-dev \
               boost-dev \
               openssl \
               make \
              " \
 && apk add --no-progress --no-cache $buildDeps \
 && wget https://github.com/argiopetech/base/archive/master.tar.gz \
 && tar xf master.tar.gz \
 && rm master.tar.gz \
 && cd base-master \
 && ./build_local.sh


FROM alpine

WORKDIR /base

RUN apk add --no-cache --no-progress gsl libstdc++

COPY --from=build /base/base-master/BUILD/bin/* /usr/local/bin/
