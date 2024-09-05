FROM rockylinux:9 AS BUILDER

RUN yum install -y make gcc-c++ zlib-devel boost-{devel,iostreams,program-options}
WORKDIR /suprdupr
COPY * .
RUN make all

FROM rockylinux:9 AS RUNNER
RUN yum install -y boost-{iostreams,program-options} procps && yum clean all
COPY --from=BUILDER /suprdupr/suprDUPr /usr/bin/suprDUPr
COPY --from=BUILDER /suprdupr/suprDUPr.read_id /usr/bin/suprDUPr.read_id
COPY --from=BUILDER /suprdupr/filterfq /usr/bin/filterfq
