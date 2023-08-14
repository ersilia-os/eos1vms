FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install rdkit-pypi
RUN pip install onnxruntime

RUN pip install pandas==1.1.2

WORKDIR /repo
COPY . /repo
