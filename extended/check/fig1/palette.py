"""Subdued ring palette: rainbow order from the centre outward, lightness rising monotonically outward so that the
rings stay distinct in greyscale and the large outer cells stay light. OKLCH -> sRGB."""
import math
def oklch(L,C,h):
    a=C*math.cos(math.radians(h)); b=C*math.sin(math.radians(h))
    l_=L+0.3963377774*a+0.2158037573*b; m_=L-0.1055613458*a-0.0638541728*b; s_=L-0.0894841775*a-1.2914855480*b
    l,m,s=l_**3,m_**3,s_**3
    r= 4.0767416621*l-3.3077115913*m+0.2309699292*s
    g=-1.2684380046*l+2.6097574011*m-0.3413193965*s
    bb=-0.0041960863*l-0.7034186147*m+1.7076147010*s
    f=lambda x: 12.92*x if x<=0.0031308 else 1.055*x**(1/2.4)-0.055
    return tuple(min(1,max(0,f(v))) for v in (r,g,bb))
HUES=[305,250,165,90,58,28]          # centre, 8, 10, 10, 14, new ring
LS  =[0.560,0.630,0.700,0.770,0.835,0.890]
CS  =[0.090,0.085,0.080,0.085,0.075,0.060]
RING=[oklch(L,C,h) for L,C,h in zip(LS,CS,HUES)]
PALE=[oklch(0.955,0.022,h) for h in HUES]
def hexc(c): return '#%02x%02x%02x'%tuple(round(255*v) for v in c)
def lum(c):
    lin=lambda v: v/12.92 if v<=0.04045 else ((v+0.055)/1.055)**2.4
    r,g,b=map(lin,c); return 0.2126*r+0.7152*g+0.0722*b
if __name__=='__main__':
    for i,(c,p) in enumerate(zip(RING,PALE)):
        print(i,hexc(c),'Y %.3f'%lum(c),'pale',hexc(p),'Y %.3f'%lum(p))
