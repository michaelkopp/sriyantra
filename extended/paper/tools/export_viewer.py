import sys, json, asyncio
from playwright.async_api import async_playwright
html=sys.argv[1]; keys=sys.argv[2].split(','); out=sys.argv[3]; pal=sys.argv[4] if len(sys.argv)>4 else 'none'
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch()
        pg=await b.new_page(viewport={'width':1400,'height':1000})
        msgs=[]
        pg.on('console', lambda m: msgs.append(m.text))
        pg.on('pageerror', lambda e: msgs.append('ERR '+str(e)))
        await pg.goto('file://'+html)
        await pg.wait_for_timeout(2500)
        clicks=sys.argv[5].split(',') if len(sys.argv)>5 and sys.argv[5] else []
        for k in keys:
            await pg.evaluate("k=>setDataset(k)",k)
            await pg.wait_for_timeout(300)
            for c in clicks:
                await pg.click('#'+c)
                await pg.wait_for_timeout(200)
            clicks=[]
            r=await pg.evaluate("""([k,pal])=>{ try{ if(typeof setPalette==='function') setPalette(pal); else {palette=pal;}
                 if (typeof render==='function') render(); return buildSVG(); } catch(e){ return 'ERR '+e.stack; } }""",[k,pal])
            open(f'{out}/{k}.svg','w').write(r)
            print(k, len(r), r[:80].replace('\n',' '))
        print('\n'.join(msgs[-10:]))
        await b.close()
asyncio.run(main())
