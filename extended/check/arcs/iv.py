"""Interval arithmetic in decimal floating point at 60 significant digits. Every operation is computed with
correct rounding to the nearest and the result is then widened by one unit in the last digit on each side, which
contains the exact result. The high precision keeps the enclosures narrow through the sixty-odd steps of a
construction."""
import decimal
from decimal import Decimal as D, localcontext
from fractions import Fraction
PREC=60
CTX=decimal.Context(prec=PREC,rounding=decimal.ROUND_HALF_EVEN,Emin=-999999,Emax=999999)
def _dn(x): return CTX.next_minus(x)
def _up(x): return CTX.next_plus(x)
def _frac_lo(q): return _dn(CTX.divide(D(q.numerator),D(q.denominator)))
def _frac_hi(q): return _up(CTX.divide(D(q.numerator),D(q.denominator)))
class I:
    __slots__=('lo','hi')
    def __init__(s,lo,hi=None):
        if hi is None: hi=lo
        if isinstance(lo,Fraction): lo=_frac_lo(lo) if lo.denominator!=1 else D(lo.numerator)
        if isinstance(hi,Fraction): hi=_frac_hi(hi) if hi.denominator!=1 else D(hi.numerator)
        if isinstance(lo,(int,float)): lo=D(lo)
        if isinstance(hi,(int,float)): hi=D(hi)
        s.lo=lo; s.hi=hi
        assert s.lo<=s.hi
    @staticmethod
    def c(x): return x if isinstance(x,I) else I(Fraction(x) if not isinstance(x,D) else x)
    def __add__(s,o): o=I.c(o); return I(_dn(CTX.add(s.lo,o.lo)),_up(CTX.add(s.hi,o.hi)))
    __radd__=__add__
    def __neg__(s): return I(-s.hi,-s.lo)
    def __sub__(s,o): o=I.c(o); return I(_dn(CTX.subtract(s.lo,o.hi)),_up(CTX.subtract(s.hi,o.lo)))
    def __rsub__(s,o): return I.c(o)-s
    def __mul__(s,o):
        o=I.c(o); p=[CTX.multiply(a,b) for a in (s.lo,s.hi) for b in (o.lo,o.hi)]
        return I(_dn(min(p)),_up(max(p)))
    __rmul__=__mul__
    def __truediv__(s,o):
        o=I.c(o)
        if o.lo<=0<=o.hi: raise ZeroDivisionError('interval division by an interval containing 0')
        p=[CTX.divide(a,b) for a in (s.lo,s.hi) for b in (o.lo,o.hi)]
        return I(_dn(min(p)),_up(max(p)))
    def __rtruediv__(s,o): return I.c(o)/s
    def sqrt(s):
        if s.lo<0: raise ValueError('sqrt of an interval reaching below 0')
        return I(max(D(0),_dn(CTX.sqrt(s.lo))),_up(CTX.sqrt(s.hi)))
    def pos(s): return s.lo>0
    def neg(s): return s.hi<0
    def mid(s): return float((s.lo+s.hi)/2)
    def w(s): return float(s.hi-s.lo)
    def __repr__(s): return '[%s, %s]'%(s.lo,s.hi)
def cos_sin(t):
    """enclosures of cos t and sin t for a rational t with |t| < 4, from the Taylor series in exact rational
    arithmetic, the remainder being bounded by the first omitted term"""
    t=Fraction(t); c=Fraction(0); s=Fraction(0); term=Fraction(1); n=0
    while True:
        if n%4==0: c+=term
        elif n%4==1: s+=term
        elif n%4==2: c-=term
        else: s-=term
        n+=1; term=term*t/n
        if n>8 and abs(term)<Fraction(1,10**70): break
    r=abs(term)
    return I(c-r,c+r),I(s-r,s+r)
