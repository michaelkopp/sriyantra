import sys, asyncio, re
from playwright.async_api import async_playwright
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch()
        for f in sys.argv[1:]:
            s=open(f).read()
            m=re.search(r'width="([\d.]+)" height="([\d.]+)"',s); w,h=float(m.group(1)),float(m.group(2))
            pg=await b.new_page(viewport={'width':int(w)+2,'height':int(h)+2})
            await pg.set_content(f'<html><body style="margin:0">{s[s.index("<svg"):]}</body></html>')
            await pg.screenshot(path=f[:-4]+'.png', clip={'x':0,'y':0,'width':w,'height':h})
            # pdf at the svg's aspect: width in inches
            await pg.pdf(path=f[:-4]+'.pdf', width=f'{w}px', height=f'{h}px', print_background=True, page_ranges='1')
            print(f, w, h)
        await b.close()
asyncio.run(main())
