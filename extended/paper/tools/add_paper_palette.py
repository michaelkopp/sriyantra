"""Write a copy of the viewer with the paper palette added (the ring tones of check/fig1/palette.py), for exporting
Figure 3. The published viewer is left unchanged.  Usage: python3 add_paper_palette.py index.html index_paper.html"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'check', 'fig1'))
from palette import RING, PALE, hexc
s = open(sys.argv[1]).read()
css = "  body.pal-paper{--bg:#ffffff;--ink:#2a2a2a;--gold:#2a2a2a;--dim:#777777;\n    --x-bg:#ffffff;--x-ink:#2a2a2a;\n"
css += "    " + "".join("--x-f-c%d:%s;" % (i, hexc(RING[i])) for i in range(6)) + "\n"
css += "    " + "".join("--x-f-lc%d:%s;" % (i, hexc(PALE[i])) for i in range(6)) + "\n"
css += "    --x-f-grow:%s;--x-f-grow2:%s;--x-f-grow-lt:%s;--x-f-grow2-lt:%s;--x-ring:#2a2a2a}\n" % (hexc(RING[5]), hexc(RING[5]), hexc(PALE[5]), hexc(PALE[5]))
css += "  body.pal-paper .tri{stroke:var(--x-ink);stroke-width:.8}\n"
a = "  body.pal-rainbow2{"; assert a in s; s = s.replace(a, css + a, 1)
b = "  rainbow2:{ scheme: 'rainbow',"; assert b in s
s = s.replace(b, "  paper:{ scheme: 'rainbow', dot: 0.013, dotOff: true, frame: false, field: true, fieldIn: '--x-f-lc4', outerIn: '--x-f-grow-lt', fixedField: true },\n" + b, 1)
open(sys.argv[2], 'w').write(s)
print('written', sys.argv[2])
