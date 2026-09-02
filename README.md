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
grown rings of cells be concyclic: dimensions **12** and **11**. For these the
answer is two-sided. No construction in which every step is forced by the
ones before it puts a ring on a circle — and the obstruction is a set of 36
(resp. 37) objects, each over-determined by the others, that anyone can check
by counting. But a construction that departs from that pattern at exactly
**one** point does exist for each ring: one number is chosen so that one
coincidence holds, and everything else is a circle through points or a
crossing of two circles. That is also the character of Chiodo's construction
of the straight figure, whose one such step is the Apollonius problem. Whether
the arc figures' one coincidence can likewise be arranged with straightedge
and compass is left open; the paper locates the question — the coincidence
can sit at only eleven objects, and at three of them it is three compass
strokes from being drawn — and reports what has been tried.

## Contents

| file | what it is |
|---|---|
| `index.html` | the viewer, a single self-contained page |
| `sri-yantra-circular-arc.pdf` | the paper |
| `sri-yantra-arcs-free.ggb` | GeoGebra: no ring held, dimension 18 |
| `sri-yantra-arcs-ring1.ggb` | GeoGebra: first ring on a circle by construction, dimension 12, one arranged coincidence |
| `sri-yantra-arcs-ring2.ggb` | GeoGebra: second ring on a circle by construction, dimension 11, one arranged coincidence |
| `core.py` | the certificate that no fully forced construction exists, and the eleven sites the coincidence may sit at — no solver needed |
| `fits.py` | constructions with a prescribed number of arranged coincidences: find them, pin the coincidence to an object, collect many, replay any saved order against the rules |
| `hunt.py` | the search for a special member whose solved parameter is rational or quadratic |
| `ring1_cubes.py`, `ring2_cubes.py`, `chiodo_cubes.py` | the earlier exhaustive solver refutations (the two rings, and the straight figure for calibration) |
| `tools/` | the drawability test of the paper's Section 7.2 and its supporting modules |

## The viewer

`index.html` needs no server and no build step — open it in a browser, or
visit the link above. It draws four figures:

- **Straight** — Chiodo's construction, from four heights on the diameter.
- **Arcs — ring 1 built**, **ring 2 built**, **no rings held** — the three
  circular-arc families.

The ring figures are drawn by the one-coincidence construction of the paper's
Proposition 7.6: eleven free numbers (ten for the second ring) and the ring's
radius, sixty-odd ruler-and-compass steps, and one coincidence solved at run
time. Every incidence, the concyclicity of the ring included, survives any
movement of the sliders; nothing is fitted.

The arc figures open on their **principal directions**: an orthogonal basis of
the tangent space, ordered by how much of the figure each moves and scaled so
that a full slider travels the same distance whichever one you pull. The
switch beside them gives the seed coordinates instead. *Sensitivity*
sets how far a slider must travel for a given effect; *re-anchor here* takes
the basis afresh at the figure on screen.

Colour schemes are **Lines**, **Rainbow**, **AYP aum colors** and
**Counts mod 3**, the last colouring each region by the number of upward
triangles covering it minus the downward ones. Figures export as SVG or PNG.

## The GeoGebra files

These are constructions, not drawings. Each contains its seeds and every later
step as `Circle` through points and `Intersect` of two objects, then the
twenty-seven sides as circular arcs — written into `geogebra.xml` as
expressions, so the kernel builds them on load and no scripting permission is
needed. All three carry the grown rings, built from the same twenty-seven
circles continued outside the yantra, and the ring files carry the circle the
held ring sits on, in red.

Drag a seed and the whole figure re-derives. In the free file the incidences
survive because they are steps of the construction. In the ring files the
concyclicity survives too — the ring's vertices are cut from its circle — and
what does not survive is the one arranged coincidence. Each ring file names it
in its opening comment and carries a number `gap`, the distance of the point
concerned from the circle it must lie on: zero at the values in the file, and
not zero after a drag. That number is the whole difference between these files
and a straightedge-and-compass construction.

## The scripts

`core.py` needs nothing but Python and runs in under a second:

```
python3 core.py              # all three figures
python3 core.py ring1
```

It prints, for each figure, the set of objects that blocks every fully forced
construction — every one of them meets more circles or points within the set
than determine it, so whichever is drawn last has one incidence too many — and
the objects at which a single arranged coincidence clears the block.

`fits.py` searches for constructions with arranged coincidences. It needs
[z3](https://github.com/Z3Prover/z3) to search but not to replay:

```
pip install z3-solver
python3 fits.py --figure ring1 --fits 1 --collect 300 -j 6       # 300 one-coincidence orders
python3 fits.py --figure ring1 --fits 1 --fitat n67 --collect 300  # pinned to one object
python3 fits.py --figure ring2 --fits 0                            # refute: no order at all
python3 fits.py --verify ring1 order.json                          # replay a saved order
```

Every saved order is replayed against the rules by an independent checker,
which also confirms the identity *seed weight = dimension + coincidences*. The
distinguished orders of the paper are the best-conditioned of several hundred.

`hunt.py` needs `mpmath`. It takes rational values for all but one of the free
parameters, one at a time so that every trial is a genuine member of the
family, solves the one coincidence, and tests the result for being rational, a
quadratic irrational, or a constructible quartic. Read it as a lottery: a hit
would be a complete answer, and two hundred thousand misses per ring are not
evidence of anything.

```
python3 hunt.py --selftest
python3 hunt.py --figure ring1 --trials 100000 -j 6
```

The `*_cubes.py` scripts are the original exhaustive refutations by
cube-and-conquer; `core.py` now reproduces their result without a solver. The
`tools/` folder holds the drawability test described in the paper's
Section 7.2, with its own README.

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
