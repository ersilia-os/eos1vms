FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge rdkit
RUN pip install scikit-learn==0.22

WORKDIR /repo
COPY . /repo
