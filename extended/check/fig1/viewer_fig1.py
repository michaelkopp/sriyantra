"""Step 8a. The nine triangles of the viewer's straight construction at the parameters of preset Fig1paper2 -> T9.json
Usage: python3 viewer_fig1.py path/to/index.html   (the viewer of github.com/michaelkopp/sriyantra; needs Playwright)"""
import sys, os, json, asyncio
from playwright.async_api import async_playwright
HERE=os.path.dirname(os.path.abspath(__file__))
pre=json.load(open(os.path.join(HERE,'fig1paper2_preset.json')))['presets']['Fig1paper2']
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await b.new_page()
        msgs=[]; pg.on('pageerror', lambda e: msgs.append(str(e)))
        await pg.goto('file://'+os.path.abspath(sys.argv[1])); await pg.wait_for_timeout(2500)
        r=await pg.evaluate("(q)=>exactFigure(q[0],q[1],q[2],q[3])",pre['exactp'])
        json.dump(r,open(os.path.join(HERE,'T9.json'),'w'),indent=1)
        print('T9.json written', msgs[:5])
        await b.close()
asyncio.run(main())
