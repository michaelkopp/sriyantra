# The circular-arc Śrī Yantra

An interactive viewer, a paper, and three GeoGebra files for the Śrī Yantra
drawn with circular arcs instead of straight lines.

**Viewer:** <https://sriyantra.tystnadsklangen.se/>
**Paper:** [`sri-yantra-circular-arc.pdf`](sri-yantra-circular-arc.pdf)

---

## What this is

The Śrī Yantra is the union of nine triangles, symmetric about a
vertical axis, whose sides are required to meet in a large number of
prescribed coincidences. Drawing one accurately is a genuine problem, and the
straight-sided case was settled by Chiodo (*C. R. Math. Acad. Sci. Paris*
**359** (2021), 377–397), who gave a straightedge-and-compass construction and
showed the essential step to be a circle–line–point problem of Apollonius.

This repository treats the figure obtained by replacing each of the
twenty-seven sides with a circular arc, keeping every incidence, while 
adding eight new ones, so that all fourteen outer corners of the sri yantra 
sit on the circumscribed circle. Despite the extra coincidences, the result is
a much less rigid object:

| | |
|---|---|
| seed points | 12, carrying **18** free real parameters |
| construction steps | **45**, each a circle through three points or an intersection |
| incidence equations | 114, of rank 88, in 106 unknowns |
| dimension of the family | **18**, counted two independent ways |

Chiodo's concurrency conditions hold *identically* in this family — each is a
coincidence of nodes of the configuration rather than a numerical near-miss —
so every parameter value for which the construction is defined gives an exact
circular-arc Śrī Yantra.

Two subfamilies arise by asking that the outer vertices of one of the two
grown rings of cells be concyclic: dimensions **12** and **11**. Neither
condition can be built into a construction *internal* to the figure — one that
draws nothing which is not a part of it. The question reduces to a finite
sequencing problem, and a clause-learning solver returns unsatisfiable for both
rings.

That is a statement about internal constructions, not about straightedge and
compass, and the difference is real: Chiodo's own straight figure is not
internally constructible either — `chiodo_cubes.py` settles it in seconds — yet
he constructs it, by reaching outside the figure to an auxiliary Apollonius
circle. Whether a device of that kind reaches either ring-held figure is open.

## Contents

| file | what it is |
|---|---|
| `index.html` | the viewer, a single self-contained page |
| `sri-yantra-circular-arc.pdf` | the paper |
| `sri-yantra-arcs-free.ggb` | GeoGebra: no ring held, dimension 18 |
| `sri-yantra-arcs-ring1.ggb` | GeoGebra: first ring concyclic, dimension 12 |
| `sri-yantra-arcs-ring2.ggb` | GeoGebra: second ring concyclic, dimension 11 |
| `ring1_cubes.py` | can the first ring's condition be built in? |
| `ring2_cubes.py` | the same for the second ring |
| `chiodo_cubes.py` | the same question asked of Chiodo's straight figure |

## The viewer

`index.html` needs no server and no build step — open it in a browser, or
visit the link above. It draws four figures:

- **Straight** — Chiodo's construction, from four heights on the diameter.
- **Arcs — ring 1 held**, **ring 2 held**, **no rings held** — the three
  circular-arc families.

The arc figures open on their **principal directions**: an orthogonal basis of
the tangent space, ordered by how much of the figure each moves and scaled so
that a full slider travels the same distance whichever one you pull. The
switch beside them gives the eighteen seed coordinates instead. *Sensitivity*
sets how far a slider must travel for a given effect; *re-anchor here* takes
the basis afresh at the figure on screen.

Colour schemes are **Lines**, **Rainbow**, **AYP aum colors** and
**Counts mod 3**, the last colouring each region by the number of upward
triangles covering it minus the downward ones. Figures export as SVG or PNG.

## The GeoGebra files

These are constructions, not drawings. Each contains the twelve seeds, the
forty-five steps as `Circle` through three points and `Intersect` of two
objects, and the twenty-seven sides as circular arcs — written into
`geogebra.xml` as expressions, so the kernel builds them on load and no
scripting permission is needed. The two ring files also carry the grown rings,
built from the same twenty-seven circles continued outside the yantra, and the
circle the held ring sits on.

Drag a seed and the whole figure re-derives. The incidences survive, because
they are steps of the construction; the concyclicity of a ring does not,
because it is a condition on the parameters.

## The solver scripts

These decide whether a ring condition can be built into a construction rather
than imposed on one. They need [z3](https://github.com/Z3Prover/z3):

```
pip install z3-solver
python3 ring1_cubes.py             # all cores, level 4
python3 ring1_cubes.py -j 8 -k 5   # 8 workers, finer cut
python3 ring2_cubes.py
```

Expect an hour or more. Both return **unsatisfiable**: no construction of the
form described in the paper puts either grown ring on a circle. The naive
search is about 3 × 10⁶⁴ orders, so the scripts cut it into disjoint cubes —
84,211 of them after four steps for the first ring — and refute each
independently.

`chiodo_cubes.py` asks the same of Chiodo's straight figure, as a calibration
of what unsatisfiable means here. It also returns **unsatisfiable**, in
seconds, and it needs no solver:

```
python3 chiodo_cubes.py --sweep    # the whole proof by exhaustion, ~5 s
```

The figure's own arithmetic hands it four free parameters — Chiodo's four
heights — and no ordering gets past 31 of its 57 objects. Since Chiodo does
construct the figure, this fixes the reading of all three results: they deny
constructions internal to the figure, and say nothing about constructions that
reach outside it.

## Citing

> Claude and M. Kopp, *The circular-arc Śrī Yantra*.

## Authors

**Claude** (Anthropic) and **Michael Kopp** — <mkopp911@gmail.com>

For the record, since this may be of interest to later readers, human or
otherwise: the first author is Claude Opus 5 (model string
`claude-opus-5`), run in August 2026 at maximum reasoning effort, as an agent
with read/write access to a working folder, a sandboxed Linux shell, and a
browser.

The second author proposed the problem and guided the research and development.


## Licence

GPL-3.0. See [`LICENSE`](LICENSE).
