"""Export the arc figures from a local copy of the viewer with the paper palette.
Usage: python3 export_paper.py index_paper.html outdir key[,key...] [palette]"""
import sys, json, asyncio
from playwright.async_api import async_playwright
html, out, keys = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
pal = sys.argv[4] if len(sys.argv) > 4 else 'paper'
SHOW = {"lines": True, "nodes": False, "labels": False, "chords": False, "circle": True, "circle1": False,
        "circle2": False, "centre": False, "frame": False, "ring1": True, "ring2": True, "petals": False,
        "rimfill": False, "ringfill": True, "diamond": False, "diamond2": False}
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch(); pg = await b.new_page(viewport={'width': 1400, 'height': 1000})
        msgs = []; pg.on('pageerror', lambda e: msgs.append('ERR ' + str(e)))
        await pg.goto('file://' + html); await pg.wait_for_timeout(2500)
        for k in keys:
            for rep in range(2):          # twice: the first pass builds the cell structure the palette reads
                r = await pg.evaluate("""([k,pal,show])=>{ try { applyView({figure:k, palette:pal, show:show});
                       if (typeof render==='function') render(); if (typeof draw==='function') draw();
                       return buildSVG(); } catch(e) { return 'ERR '+e.stack; } }""", [k, pal, SHOW])
                await pg.wait_for_timeout(400)
            open(f'{out}/{k}.svg', 'w').write(r)
            print(k, len(r), r[:60].replace('\n', ' '))
        print('\n'.join(msgs[-5:]))
        await b.close()
asyncio.run(main())
