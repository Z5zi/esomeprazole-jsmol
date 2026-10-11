# Esomeprazole synthesis in JSmol

Nine genuine interactive JSmol animations paired side by side with native ChemDraw mechanism diagrams, embedded in Microsoft Sway through a GeoGebra Activity Web element.

Use the numbered buttons to select a reaction. Drag to rotate, scroll to zoom, play or pause, and scrub through the frames. Most steps contain 61 frames; step 4 now contains 181 frames across addition, collapse, alpha deprotonation and the formal rearrangement. Stage buttons jump to each event and the matching 2D mechanism changes with the timeline. The `?step=1` through `?step=9` query selects the initial step.

The diagrams were exported as SVG directly from the native ChemDraw document, then rendered to PNG at four times the vector viewport dimensions. Full panels are 5,820 pixels wide. SVG copies are retained, and `diagrams/manifest.json` records dimensions and hashes. Unicode subscripts use system-font fallback. The original on-screen footer was removed at the user's request; scientific scope remains documented here and in the native mechanism captions.

## Scientific scope

These are illustrative interpolations between prepared 3D structures, not calculated intrinsic reaction coordinates or transition-state trajectories. Dashed bonds mark schematic changes. Multistage steps and proton transfers are condensed as described on screen. The detailed titanium catalyst cycle is not represented.

The endpoint molecular graphs and net charges were checked against the source models with RDKit. Step 4 additionally matches all five native structural states in panel 04A and panel 04. Its atom mapping follows the formal six-electron drawing: the former acyl carbonyl oxygen forms the bond to the exo methylene, and the former N-oxide oxygen becomes the ester carbonyl oxygen. See `activation-validation.json`. The final esomeprazole sulfoxide was independently reassigned from 3D coordinates with the accurate CIP labeler and confirmed as **S**.

JSmol is loaded from its official St Olaf distribution. Molecular models and exported mechanism diagrams are served from this repository. The native editable CDX, account information and credentials are not included.

## Sources

- https://patents.google.com/patent/EP0484265A1/en
- https://patents.google.com/patent/EP0130729A2/en
- https://patents.google.com/patent/WO2010134099A1/en
- Cotton et al., Tetrahedron: Asymmetry 11 (2000), 3819–3825. https://doi.org/10.1016/S0957-4166(00)00352-9
- Kurzawa et al. (2023), mechanistic discussion of O-acylation, deprotonation and alternative Boekelheide pathways: https://doi.org/10.1002/chem.202204015
- JSmol documentation: https://wiki.jmol.org/index.php/Jmol_JavaScript_Object

## Hosting

GitHub Pages serves the root of the `main` branch. `.nojekyll` prevents Jekyll processing. The page can be embedded with an HTTPS iframe; Sway uses a supported GeoGebra Activity as the outer container.
