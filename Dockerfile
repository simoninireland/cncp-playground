ARG BASE_IMAGE=jupyter/minimal-notebook
FROM $BASE_IMAGE

MAINTAINER simon.dobson@st-andrews.ac.uk

EXPOSE 8050

WORKDIR /tmp
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY newman-ziff-extended.ipynb .
