"""Validate the delivered assets and all step 7 states, including atom identities."""
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
steps = json.loads((ROOT / "steps.json").read_text(encoding="utf-8"))
summary = json.loads((ROOT / "model-validation.json").read_text())
expected_states = [
    "COc1ccc(N)c(N)c1.S=C=S",
    "COc1ccc(NC(=S)S)c(N)c1",
    "COc1ccc2c(c1)[NH2+]C([S-])(S)N2",
    "COc1ccc2c(c1)NC([S-])([SH2+])N2",
    "COc1ccc2[nH]c(=S)[nH]c2c1.S",
    "COc1ccc2[nH]c(S)nc2c1.S",
]


def canon(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=False)


def parent(mol, idx):
    return tuple(a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors())


checked = []
for step in steps:
    entry = next(r for r in summary["results"] if r["step"] == step["num"])
    assert entry["frames"] == step["frames"] and entry["atom_count"] == step["atoms"]
    blocks = [b.lstrip("\n") for b in (ROOT / f"models/step-{step['num']:02d}.sdf")
              .read_text().split("$$$$") if b.strip()]
    assert len(blocks) == step["frames"], step["num"]
    assert (ROOT / step["diagram"]).is_file()
    if step.get("stages"):
        covered = []
        for stage in step["stages"]:
            covered += list(range(stage["start"], stage["end"] + 1))
            assert (ROOT / stage["diagram"]).is_file()
        assert covered == list(range(1, step["frames"] + 1))
    if step["num"] == 7:
        assert step["frames"] == 226 and step["atoms"] == 23
        assert [s["start"] for s in step["stages"]] == [1, 46, 91, 136, 181]
        assert all("7a" in step["stages"][i]["name"] for i in (1, 2))
        assert all("7b" in step["stages"][i]["name"] for i in (3, 4))
        models = [Chem.MolFromMolBlock(b, removeHs=False) for b in blocks]
        assert all(m is not None and m.GetNumAtoms() == 23 for m in models)
        atoms = [a.GetAtomicNum() for a in models[0].GetAtoms()]
        for mol in models:
            assert [a.GetAtomicNum() for a in mol.GetAtoms()] == atoms
            assert rdMolDescriptors.CalcMolFormula(mol) == "C8H10N2OS2"
            assert Chem.GetFormalCharge(mol) == 0
            assert np.isfinite(mol.GetConformer().GetPositions()).all()
        endpoints = [models[i] for i in [0, 45, 90, 135, 180, 225]]
        for mol, expected in zip(endpoints, expected_states):
            assert canon(mol) == canon(Chem.MolFromSmiles(expected))
        mapping = json.loads((ROOT / "step07-validation.json").read_text())["atom_mapping"]
        retained, leaving = mapping["retained_sulfur"], mapping["departing_sulfur"]
        cc = mapping["CS2_carbon"]
        for mol in endpoints[4:]:
            assert mol.GetBondBetweenAtoms(cc, retained) is not None
            assert mol.GetBondBetweenAtoms(cc, leaving) is None
            neighbours = mol.GetAtomWithIdx(leaving).GetNeighbors()
            assert len(neighbours) == 2 and all(a.GetAtomicNum() == 1 for a in neighbours)
        for left, right, donor, acceptor in [
            (0, 1, mapping["first_nitrogen"], leaving),
            (2, 3, mapping["ring_closing_nitrogen"], leaving),
            (4, 5, mapping["ring_closing_nitrogen"], retained),
        ]:
            moved = [(parent(endpoints[left], i), parent(endpoints[right], i))
                     for i, atomic in enumerate(atoms)
                     if atomic == 1 and parent(endpoints[left], i) != parent(endpoints[right], i)]
            assert moved == [((donor,), (acceptor,))], (left, right, moved)
        spt = (ROOT / step["script"]).read_text()
        assert spt.count("$$$$") == 226
        assert "PARTIAL" in spt
        checked.append("226 valid frames; six native states; retained/departing sulfur and three proton transfers")
    if step["num"] == 4:
        report = json.loads((ROOT / "activation-validation.json").read_text())
        for n, state in zip([0, 45, 90, 135, 180], report["stages"]):
            assert canon(Chem.MolFromMolBlock(blocks[n])) == canon(Chem.MolFromSmiles(state["smiles"]))
        checked.append("Step 4: five existing states")
    if step["num"] == 9:
        mol = Chem.MolFromMolBlock(blocks[-1], removeHs=False)
        mol.GetConformer().Set3D(True)
        Chem.AssignStereochemistryFrom3D(mol)
        rdCIPLabeler.AssignCIPLabels(mol)
        assert next(a for a in mol.GetAtoms() if a.GetSymbol() == "S").GetProp("_CIPCode") == "S"
        checked.append("Step 9: S configuration retained")

manifest = json.loads((ROOT / "diagrams/manifest.json").read_text())
for item in manifest["images"]:
    image = ROOT / "diagrams" / (item["name"] + ".png")
    assert hashlib.sha256(image.read_bytes()).hexdigest() == item["sha256"], image
print(json.dumps({"status": "PASS", "steps": len(steps), "checks": checked, "diagram_hashes": len(manifest["images"])}))
