FROM ubuntu:latest
WORKDIR /app
COPY . /app
RUN apt-get update && apt-get install -y --no-install-recommends wget
RUN wget --no-check-certificate https://github.com/Augroup/miniQuant/releases/latest/download/miniQuant_linux_latest.tar.gz && tar -zxvf miniQuant_linux_latest.tar.gz && cd miniQuant_linux && chmod +x miniQuant
CMD ["miniQuant"]