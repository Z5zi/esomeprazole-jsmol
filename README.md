# Esomeprazole synthesis in JSmol

Nine genuine interactive JSmol animations for teaching, intended to be embedded in Microsoft Sway through a GeoGebra Activity Web element.

Use the numbered buttons to select a reaction. Drag to rotate, scroll to zoom, play or pause, and scrub through 61 frames per reaction. The `?step=1` through `?step=9` query selects the initial step.

## Scientific scope

These are illustrative interpolations between prepared 3D structures, not calculated intrinsic reaction coordinates or transition-state trajectories. Dashed bonds mark schematic changes. Multistage steps and proton transfers are condensed as described on screen. The detailed titanium catalyst cycle is not represented.

The 18 endpoint molecular graphs and their net charges were checked against the source models with RDKit. The final sulfoxide was independently reassigned from 3D coordinates using the accurate CIP labeler and confirmed as **S**. See `model-validation.json`.

JSmol is loaded from its official St Olaf distribution. Molecular model data are served from this repository. No original ChemDraw files, OneNote diagrams, student account information, or credentials are included.

## Sources

- https://patents.google.com/patent/EP0484265A1/en
- https://patents.google.com/patent/EP0130729A2/en
- https://patents.google.com/patent/WO2010134099A1/en
- Cotton et al., Tetrahedron: Asymmetry 11 (2000), 3819–3825. https://doi.org/10.1016/S0957-4166(00)00352-9
- JSmol documentation: https://wiki.jmol.org/index.php/Jmol_JavaScript_Object

## Hosting

GitHub Pages serves the root of the `main` branch. `.nojekyll` prevents Jekyll processing. The page can be embedded with an HTTPS iframe; Sway uses a supported GeoGebra Activity as the outer container.
