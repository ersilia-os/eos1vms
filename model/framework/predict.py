from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import joblib
import sys

infile=sys.argv[1]
outfile=sys.argv[2]
checkpoints=sys.argv[3]

mdl = joblib.load("{0}/10uM/mNB_10uM_all.pkl".format(checkpoints))
chemid = list(mdl.targets)

TOP = 50

def predict(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    res = np.zeros(len(fp), np.int32)
    DataStructs.ConvertToNumpyArray(fp, res)
    probs=mdl.predict_proba(res.reshape(1,-1))[0]
    idxs=np.argsort(-probs)
    preds = []
    for i in range(TOP):
        idx = idxs[i]
        preds += [(chemid[idx], probs[idx])]
    return preds

with open(infile, "r") as f:
    smiles_list = []
    for line in f:
        smi = line.rstrip()
        smiles_list += [smi]

with open(outfile, "w") as f:
    f.write("{0}\t{1}\t{2}\n".format("smiles", "chembl_id", "probability"))
    for smi in smiles_list:
        preds = predict(smi)
        for pred in preds:
            f.write("{0}\t{1}\t{2}\n".format(smi, pred[0], pred[1]))
