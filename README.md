# The circular-arc Śrī Yantra

An interactive viewer, two papers and five GeoGebra files for the Śrī Yantra
drawn with circular arcs instead of straight lines, and for its extension by
one more ring of triangles.

**Viewer:** <https://sriyantra.tystnadsklangen.se/>
**Papers:** [`sri-yantra-circular-arc.pdf`](sri-yantra-circular-arc.pdf), *The circular-arc Śrī Yantra*;
[`sri-yantra-extended.pdf`](sri-yantra-extended.pdf), *Extended circular-arc Śrī Yantras*

---

## What this is

The Śrī Yantra is the union of nine triangles, symmetric about a
vertical axis, whose sides are required to meet in a large number of
prescribed coincidences. Drawing one accurately is a genuine problem, and the
straight-sided case was settled by Chiodo (*C. R. Math. Acad. Sci. Paris*
**359** (2021), 377–397), who gave a straightedge-and-compass construction and
showed the essential step to be a circle–line–point problem of Apollonius.

The first paper treats the figure obtained by replacing each of the
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

## The extended figure

The second paper asks whether the triangular structure can be continued by one
more ring. Read from the centre outward, the classical figure is built ring by
ring inside nested *kites*, the quadrilaterals bounded by the legs of one
upward and one downward triangle: each ring has its top and bottom vertices at
the apexes of one kite and its other outer vertices on the legs of the next,
and these incidences are exactly Chiodo's twelve concurrences. A further ring
needs a fifth kite, whose legs carry the six outer corners of the ring of
fourteen: an upward triangle U0 through the corners of D1, D2, D3 and a
downward triangle D0 through those of U1, U2, U3. The figure then has eleven
triangles, and the extension is asked to be *laminar*: no leg of a triangle
meets a leg of another triangle of the same family.

The paper proves two results.

- **No straight extension.** No laminar straight Śrī Yantra has both triples
  of outer corners collinear with a laminar extension. The question reduces to
  thirteen points on the projective line, and a computer-assisted certificate
  decides it over the 58,786 cones given by the triangulations of a 13-gon,
  with interval arithmetic and linear-programming certificates checked in exact
  rational arithmetic, and an asymptotic analysis at the one vertex of the
  associahedron where the equations and the inequalities approach each other.
- **Two constructed arc extensions.** With circular arcs the extension exists.
  Two such figures, whose new rings have ten and fourteen cells with every
  outer node on the circumscribed circle, are given by straightedge-and-compass
  constructions with 22 and 20 free parameters in which every incidence holds
  by construction and no coincidence is arranged. Their laminarity is certified
  in interval arithmetic.

## Contents

| file | what it is |
|---|---|
| `index.html` | the viewer, a single self-contained page |
| `sri-yantra-circular-arc.pdf` | the first paper, *The circular-arc Śrī Yantra* |
| `sri-yantra-extended.pdf` | the second paper, *Extended circular-arc Śrī Yantras* |
| `sri-yantra-arcs-free.ggb` | GeoGebra: no ring held, dimension 18 |
| `sri-yantra-arcs-ring1.ggb` | GeoGebra: first ring on a circle by construction, dimension 12, one arranged coincidence |
| `sri-yantra-arcs-ring2.ggb` | GeoGebra: second ring on a circle by construction, dimension 11, one arranged coincidence |
| `sri-yantra-arcs-extended-2.ggb` | GeoGebra: extended figure with a ring of ten cells, 22 free parameters, no arranged coincidence |
| `sri-yantra-arcs-extended-3.ggb` | GeoGebra: extended figure with a ring of fourteen cells, 20 free parameters, no arranged coincidence |
| `core.py` | the certificate that no fully forced construction exists, and the eleven sites the coincidence may sit at — no solver needed |
| `fits.py` | constructions with a prescribed number of arranged coincidences: find them, pin the coincidence to an object, collect many, replay any saved order against the rules |
| `hunt.py` | the search for a special member whose solved parameter is rational or quadratic |
| `ring1_cubes.py`, `ring2_cubes.py`, `chiodo_cubes.py` | the earlier exhaustive solver refutations (the two rings, and the straight figure for calibration) |
| `tools/` | the drawability test of the first paper's Section 7.2 and its supporting modules |
| `extended/` | the programs and data of the second paper: the certifiers of the straight case, the checks, and the paper's source |

## The viewer

`index.html` needs no server and no build step — open it in a browser, or
visit the link above. It draws six figures:

- **Straight** — Chiodo's construction, from four heights on the diameter.
- **Arcs — ring 1 built**, **ring 2 built**, **no rings held** — the three
  circular-arc families.
- **Arcs extended 2**, **Arcs extended 3** — the two constructed extended
  figures of the second paper, with eleven triangles and a ring of ten or
  fourteen cells on the circle; 22 and 20 free numbers, and nothing solved
  while the sliders move.

The ring figures are drawn by the one-coincidence construction of the first
paper's Proposition 7.6: eleven free numbers (ten for the second ring) and the ring's
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
needed. The three files of the first paper carry the grown rings, built from
the same twenty-seven circles continued outside the yantra, and the ring files
carry the circle the held ring sits on, in red. The two extended files carry
the eleven triangles and the ring they cut, every outer node of the ring on the
circle by construction.

Drag a seed and the whole figure re-derives. In the free file the incidences
survive because they are steps of the construction. In the ring files the
concyclicity survives too — the ring's vertices are cut from its circle — and
what does not survive is the one arranged coincidence. Each ring file names it
in its opening comment and carries a number `gap`, the distance of the point
concerned from the circle it must lie on: zero at the values in the file, and
not zero after a drag. That number is the whole difference between these files
and a straightedge-and-compass construction. The two extended files have no
arranged coincidence, so every incidence survives a drag.

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
`tools/` folder holds the drawability test described in the first paper's
Section 7.2, with its own README.

### The second paper

Everything the second paper uses is in `extended/`, described in its
Section 10 and in the READMEs of the subfolders:

| folder | what it is |
|---|---|
| `extended/certifier/cones/` | the certificate of the straight case: the 58,786 cones, the interval test, the linear-programming certificates verified in exact rationals, a resumable parallel run, and the records of the run (Python with NumPy and SciPy) |
| `extended/certifier/fan/` | the certified regions near the vertex where the equations and inequalities meet (Node.js) |
| `extended/certifier/cells/` | the second, independent certification by a cell decomposition, with its records |
| `extended/check/` | exact checks of the algebra, the tables of the paper, the figures, and `arcs/`, the interval-arithmetic check of the two constructions |
| `extended/paper/` | the LaTeX source and figures of the paper |

The check of the two arc constructions needs nothing but Python and runs in
about a second:

```
python3 extended/check/arcs/certify.py
```

The cone certifier runs its own soundness tests before it starts; its records
took 71 hours of processor time on ten cores.

## Citing

> Claude and M. Kopp, *The circular-arc Śrī Yantra*, 2026.
>
> Claude and M. Kopp, *Extended circular-arc Śrī Yantras*, 2026.

## Authors

**Claude** (Anthropic) and **Michael Kopp** — <mkopp911@gmail.com>

For the record, since this may be of interest to later readers, human or
otherwise: two Claude models contributed, and the work divides between them as
follows.

- *The circular-arc Śrī Yantra*, the viewer, the GeoGebra files and the scripts
  at the top level are the work of Claude Opus 5 (model string
  `claude-opus-5`), run in August 2026 at maximum reasoning effort, as an agent
  with read/write access to a working folder, a sandboxed Linux shell, and a
  browser. So are the constructions of the two extended figures, in the viewer
  and in their GeoGebra files.
- *Extended circular-arc Śrī Yantras*, its proofs and the programs in
  `extended/` are largely the work of Claude Opus 5.5 (model string
  `claude-opus-5-5`), run in September and October 2026 at extra reasoning effort, run as an agent with file,
  shell and browser access. It also revised the explanatory text of the viewer
  slightly.

The second author proposed the problems and guided the research and
development.


## Licence

GPL-3.0. See [`LICENSE`](LICENSE).
