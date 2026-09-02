# tools

The scripts behind Section 7.2 of the paper — whether the one arranged
coincidence can be drawn. They need `mpmath` and `numpy`, and the two
`*_cubes.py` files from the repository root beside them (for the ring-vertex
tables); copy or symlink those in.

    python3 scanall.py ring1 orders_ring1_f1_n67.json ...   # the drawability test
    python3 concurrent.py ring2 orders_ring2_f1_n31.json    # the same, one file, verbose
    python3 dscan.py ring1 orders_ring1_f1.json             # depth of every parameter

The order files come from `fits.py --collect` in the repository root, for
example `python3 fits.py --figure ring1 --fits 1 --fitat n67 --collect 300`.

| file | role |
|---|---|
| `concurrent.py` | the test itself: stands-still measured by displacement, both shapes, every hit re-verified |
| `scanall.py` | `concurrent.py` over many collections, one line each |
| `dcheck.py`, `dscan.py`, `degree.py` | the depth of a parameter — number of square roots between it and the coincidence — measured, and the graph version it replaced |
| `score_orders.py`, `emit_ringbuild.py`, `ring_construct.py`, `recipe.json` | compile an order into an executable recipe; `ringbuild.json` in the viewer is their output |
| `flipexec.py` | the executor in arbitrary precision, with a branch flip |
