"""Build the step 7 teaching sequence from an explicit native ChemDraw export.

Usage: python tools/build_step07.py --source-cdxml <native-export.cdxml>
Requires RDKit and NumPy. The input document is read only.
"""
from argparse import ArgumentParser
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
FRAMES_PER_STAGE = 45
NAMES = [
    "7 · CS₂ addition",
    "7a · Ring closure",
    "7a · Proton relay",
    "7b · H₂S loss",
    "7b · Thione ⇌ thiol",
]
FILES = ["addition", "cyclisation", "proton-relay", "sulfur-loss", "tautomerisation"]
DESCRIPTIONS = [
    "The first amine adds to CS₂; a condensed proton transfer gives the neutral dithiocarbamic acid shown.",
    "The neighbouring amine closes the five-membered ring, giving the cyclic tetrahedral N⁺/S⁻ intermediate.",
    "A solvent-assisted proton relay transfers an N–H proton to the departing sulfur, giving S–H₂⁺.",
    "Collapse restores C=S and releases the second sulfur as H₂S, giving benzimidazole-2-thione.",
    "Proton transfer and bond rearrangement connect the thione with the thiol tautomer, intermediate I.",
]
SCOPE = (
    "Neutral-form teaching sequence matching ChemDraw; protonation depends on the basic "
    "medium and work-up. Solvent and base are omitted from the proton relays. "
    "Motion is illustrative, not a calculated reaction pathway."
)


def canonical(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=False)


def finish(edit):
    mol = edit.GetMol()
    mol.UpdatePropertyCache(strict=True)
    Chem.SanitizeMol(mol)
    return mol


def bond(edit, a, b, order):
    edit.GetBondBetweenAtoms(a, b).SetBondType(
        {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE}[order]
    )


def hydrogens(edit, atom, number):
    edit.GetAtomWithIdx(atom).SetNumExplicitHs(number)
    edit.GetAtomWithIdx(atom).SetNoImplicit(True)


def charge(edit, atom, value):
    edit.GetAtomWithIdx(atom).SetFormalCharge(value)


def keyed_hydrogens(mol, previous=None):
    """Keep each heavy atom and each transferred proton's identity across states."""
    expanded = Chem.AddHs(mol)
    keys = {i: f"A{i}" for i in range(mol.GetNumAtoms())}
    parents = {}
    for atom in expanded.GetAtoms():
        if atom.GetAtomicNum() == 1:
            parents.setdefault(atom.GetNeighbors()[0].GetIdx(), []).append(atom.GetIdx())
    current, spare, needed = {}, [], []
    if previous is None:
        serial = 0
        for parent, indices in parents.items():
            current[parent] = []
            for idx in indices:
                keys[idx] = f"H{serial}"
                current[parent].append(keys[idx])
                serial += 1
    else:
        for parent in sorted(set(previous) | set(parents)):
            old, indices = previous.get(parent, []), parents.get(parent, [])
            current[parent] = []
            for key, idx in zip(old, indices):
                keys[idx] = key
                current[parent].append(key)
            spare += old[len(indices):]
            needed += [(parent, idx) for idx in indices[len(old):]]
        assert len(spare) == len(needed), "A proton was created or lost"
        for key, (parent, idx) in zip(spare, needed):
            keys[idx] = key
            current[parent].append(key)
    order = sorted(keys, key=lambda idx: (keys[idx][0] == "H", int(keys[idx][1:])))
    return Chem.RenumberAtoms(expanded, order), current


def align(moving, fixed):
    mc, fc = moving.mean(0), fixed.mean(0)
    u, _, vt = np.linalg.svd((moving - mc).T @ (fixed - fc))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return rotation, fc - mc @ rotation


def conform(mol, seed):
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) < 0:
        params.useRandomCoords = True
        assert AllChem.EmbedMolecule(mol, params) >= 0
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return np.array(mol.GetConformer().GetPositions())


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--source-cdxml", required=True, type=Path)
    args = parser.parse_args()
    xml = ET.parse(args.source_cdxml).getroot()

    def fragment(fid):
        original = xml.find(f".//fragment[@id='{fid}']")
        assert original is not None, f"Missing native fragment {fid}"
        root = ET.Element("CDXML")
        ET.SubElement(root, "page").append(deepcopy(original))
        mol = Chem.MolsFromCDXML(ET.tostring(root, encoding="unicode"))[0]
        mol.RemoveAllConformers()
        return mol

    acid = fragment(10227)
    carbon = next(a.GetIdx() for a in acid.GetAtoms()
                  if sum(n.GetAtomicNum() == 16 for n in a.GetNeighbors()) == 2)
    sulfurs = [a.GetIdx() for a in acid.GetAtomWithIdx(carbon).GetNeighbors()
               if a.GetAtomicNum() == 16]
    retained = next(a for a in sulfurs if acid.GetBondBetweenAtoms(carbon, a).GetBondTypeAsDouble() == 2)
    leaving = next(a for a in sulfurs if a != retained)
    first_n = next(a.GetIdx() for a in acid.GetAtomWithIdx(carbon).GetNeighbors() if a.GetAtomicNum() == 7)
    second_n = next(a.GetIdx() for a in acid.GetAtoms() if a.GetAtomicNum() == 7 and a.GetIdx() != first_n)
    ring = [a.GetIdx() for a in acid.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic()]
    assert len(ring) == 6

    edit = Chem.RWMol(acid)
    edit.RemoveBond(first_n, carbon)
    bond(edit, carbon, leaving, 2)
    hydrogens(edit, leaving, 0)
    hydrogens(edit, first_n, 2)
    reactants = finish(edit)

    edit = Chem.RWMol(acid)
    edit.AddBond(second_n, carbon, Chem.BondType.SINGLE)
    bond(edit, carbon, retained, 1)
    charge(edit, retained, -1)
    charge(edit, second_n, 1)
    cyclic = finish(edit)

    edit = Chem.RWMol(cyclic)
    hydrogens(edit, second_n, 1)
    charge(edit, second_n, 0)
    hydrogens(edit, leaving, 2)
    charge(edit, leaving, 1)
    protonated = finish(edit)

    edit = Chem.RWMol(protonated)
    edit.RemoveBond(carbon, leaving)
    bond(edit, carbon, retained, 2)
    charge(edit, retained, 0)
    charge(edit, leaving, 0)
    thione = finish(edit)

    edit = Chem.RWMol(thione)
    Chem.Kekulize(edit, clearAromaticFlags=True)
    bond(edit, carbon, retained, 1)
    bond(edit, carbon, second_n, 2)
    hydrogens(edit, second_n, 0)
    hydrogens(edit, retained, 1)
    thiol = finish(edit)

    states = [reactants, acid, cyclic, protonated, thione, thiol]
    native_ids = [[1113, 1135], [10227], [10261], [10301], [10375], [10407]]
    expected = [
        [fragment(1113), fragment(1135)], [acid], [fragment(10261)],
        [fragment(10301)], [fragment(10375), Chem.MolFromSmiles("S")],
        [fragment(10407), Chem.MolFromSmiles("S")],
    ]
    reports = []
    elements = Counter(a.GetSymbol() for a in Chem.AddHs(states[0]).GetAtoms())
    for index, (state, wanted) in enumerate(zip(states, expected)):
        actual_graphs = sorted(canonical(m) for m in Chem.GetMolFrags(state, asMols=True))
        assert actual_graphs == sorted(canonical(m) for m in wanted), (index, actual_graphs)
        assert Counter(a.GetSymbol() for a in Chem.AddHs(state).GetAtoms()) == elements
        assert Chem.GetFormalCharge(state) == 0
        reports.append({
            "state": index, "frame": index * FRAMES_PER_STAGE + 1,
            "smiles": canonical(state), "formula": rdMolDescriptors.CalcMolFormula(state),
            "net_charge": 0, "native_fragment_ids": native_ids[index], "native_graph_match": True,
        })

    models, previous = [], None
    for state in states:
        mol, previous = keyed_hydrogens(state, previous)
        models.append(mol)
    count = models[0].GetNumAtoms()
    coords = [conform(mol, 7200 + i) for i, mol in enumerate(models)]
    base = coords[1]
    center = base[ring].mean(0)
    _, _, axes = np.linalg.svd(base[ring] - center)
    if np.linalg.det(axes) < 0:
        axes[-1] *= -1
    coords = [(p - center) @ axes.T for p in coords]
    base = coords[1]
    # Give the functional groups a clear right-facing view.
    if base[carbon, 0] < base[ring, 0].mean():
        coords = [p * np.array([-1, -1, 1]) for p in coords]
        base = coords[1]
    for i, mol in enumerate(models):
        rotation, translation = align(coords[i][ring], base[ring])
        coords[i] = coords[i] @ rotation + translation
        for component in Chem.GetMolFrags(mol):
            if first_n in component:
                continue
            indices = list(component)
            main_indices = list(next(c for c in Chem.GetMolFrags(mol) if first_n in c))
            points = coords[i][indices]
            x = coords[i][main_indices, 0].max() + 2.7
            coords[i][indices] += np.array([x - points[:, 0].min(), -points[:, 1].mean(), -points[:, 2].mean()])

    def bonds(mol):
        return {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))): b.GetBondTypeAsDouble()
                for b in mol.GetBonds()}

    bond_sets = [bonds(mol) for mol in models]
    blocks, connections = [], []
    frames = 5 * FRAMES_PER_STAGE + 1
    for frame in range(frames):
        segment = min(frame // FRAMES_PER_STAGE, 4)
        local = (frame - segment * FRAMES_PER_STAGE) / FRAMES_PER_STAGE
        progress = max(0, min(1, (local - .12) / .76))
        fraction = progress * progress * (3 - 2 * progress)
        left, right = segment, segment + 1
        mol = Chem.Mol(models[left if fraction < .5 else right])
        mol.RemoveAllConformers()
        conf = Chem.Conformer(count)
        conf.Set3D(True)
        points = coords[left] + (coords[right] - coords[left]) * fraction
        for index, xyz in enumerate(points):
            conf.SetAtomPosition(index, tuple(xyz))
        mol.AddConformer(conf)
        mol.SetProp("_Name", f"Step 7 stage {segment + 1}; frame {frame + 1}/{frames}")
        blocks.append(Chem.MolToMolBlock(mol, kekulize=False) + "\n$$$$\n")
        if .02 < fraction < .98:
            for pair in sorted(set(bond_sets[left]) | set(bond_sets[right])):
                before, after = bond_sets[left].get(pair, 0), bond_sets[right].get(pair, 0)
                if before == after:
                    continue
                style = "PARTIALDOUBLE" if min(before, after) >= 1 and max(before, after) >= 2 else "PARTIAL"
                a, b = (frame * count + value for value in pair)
                connections.append(f"connect (atomIndex={a}) (atomIndex={b}) {style};")
    sdf = "".join(blocks)
    setup = "\n".join(connections) + """
select all; spacefill 22%; wireframe 0.10; color atoms CPK;
background white; set antialiasDisplay true; set frank true;
set perspectiveDepth false; set zoomLarge false; center {*}; zoom 105;
rotate x 12; rotate y -10;
select formalCharge != 0; label %[formalCharge]; color labels [x16405a]; font label 15; select all;
animation fps 10; animation mode ONCE; frame 1;
"""
    (ROOT / "models/step-07.sdf").write_text(sdf, encoding="utf-8", newline="\r\n")
    (ROOT / "models/step-07.spt").write_text(
        'set autobond false;\nload DATA "model"\n' + sdf + 'END "model";\n' + setup, encoding="utf-8", newline="\r\n")
    steps = json.loads((ROOT / "steps.json").read_text(encoding="utf-8"))
    step = next(s for s in steps if s["num"] == 7)
    step.update(frames=frames, atoms=count, description=SCOPE)
    step["phases"] = ["Diamine + CS₂", "Addition, ring closure, proton relay and sulfur loss", "Thiol I + H₂S"]
    step["stages"] = [
        {"name": name, "start": i * FRAMES_PER_STAGE + 1,
         "end": (i + 1) * FRAMES_PER_STAGE if i < 4 else frames,
         "diagram": f"diagrams/step-07-{FILES[i]}.png", "description": DESCRIPTIONS[i]}
        for i, name in enumerate(NAMES)
    ]
    (ROOT / "steps.json").write_text(json.dumps(steps, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\r\n")
    report = {
        "status": "PASS", "rdkit": rdBase.rdkitVersion,
        "source_cdxml_sha256": hashlib.sha256(args.source_cdxml.read_bytes()).hexdigest(),
        "frames": frames, "atoms_per_frame": count, "states": reports,
        "atom_mapping": {"CS2_carbon": carbon, "retained_sulfur": retained,
                         "departing_sulfur": leaving, "first_nitrogen": first_n,
                         "ring_closing_nitrogen": second_n},
        "scope": SCOPE,
        "labels": "Player 7a and 7b subdivide the upper and lower rows of native ChemDraw panel 07A.",
    }
    (ROOT / "step07-validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    summary_path = ROOT / "model-validation.json"
    summary = json.loads(summary_path.read_text())
    entry = next(r for r in summary["results"] if r["step"] == 7)
    entry.update(frames=frames, intermediate_graphs="All six states match native ChemDraw",
                 atom_mapping="Retained/departing sulfurs and three proton transfers verified")
    summary["step7_revision"] = "2026-10-20"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\r\n")
    print(json.dumps({"status": "PASS", "frames": frames, "atoms": count, "native_states": len(states)}))


if __name__ == "__main__":
    main()
