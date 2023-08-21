FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN pip install rdkit==2023.3.3
RUN pip install onnxruntime==1.15.1
RUN pip install pandas==2.0.3

WORKDIR /repo
COPY . /repo
