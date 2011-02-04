from sage.all import multinomial as mn, ordered_partitions as op, factorial
import re

def s2111(n):
    return sum( 3*mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if a > b ) + \
           sum( 3*mn(a,b,c,d,e)/2 for a,b,c,d,e in op(n-1,5) if a == b)

def s2100(n):
    return sum( 2*mn(a,b,c) for a,b,c in op(n-1,3) if a > b > 1 ) + \
           sum( mn(a,b,c) for a,b,c in op(n-1,3) if a > b == 1 ) + \
           sum( mn(a,b,c) for a,b,c in op(n-1,3) if a == b > 1 )

def s2110(n):
    return sum( 2*mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a > b ) + \
           sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a == b )

def s2101(n):
    return sum( 3*mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a > b > 1 ) + \
           sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a > b == 1 ) + \
           sum( 3*mn(a,b,c,d)/2 for a,b,c,d in op(n-1,4) if a == b > 1 )

def s1111(n):
    return sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a>1 )

def s1200(n):
    return sum( mn(a,b,c) for a,b,c in op(n-1,3) if min(a,c) > 1 )

def s1210(n):
    return sum( 3*mn(a,b,c,d)/2 for a,b,c,d in op(n-1,4) if a==b and c>1) + \
           sum( mn(a,b,c,d) for a,b,c,d in op(n-1, 4) if c == 1) + \
           sum( 3*mn(a,b,c,d) for a,b,c,d in op(n-1, 4) if c>1 and a>b )

def s1201(n):
    return sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if min(b,c) >= 2)

def s1211(n):
    return sum( mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if c == 1) + \
           sum( 3*mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if c>1 and a > b) + \
           sum( 3*mn(a,b,c,d,e)/2 for a,b,c,d,e in op(n-1,5) if a==b and c>1 )

def s0300(n):
    return sum( mn(a,b,c)/6 for a,b,c in op(n-1,3) if min(a,b,c) > 1 )

def s0301(n):
    return sum( mn(a,b,c,d)/6 for a,b,c,d in op(n-1,4) if min(a,b,c) > 1 )

def s0210(n):
    return sum( mn(a,b,c)/2 for a,b,c in op(n-1,3) if a == b and c > 1 ) + \
           sum( mn(a,b,c) for a,b,c in op(n-1,3) if a > b and c > 1 )

def s0211(n):
    return sum( mn(a,b,c,d)/2 for a,b,c,d in op(n-1,4) if a == b and c > 1 ) + \
           sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a > b and c > 1 )

def s1300(n):
    return sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a>1 and b>c ) + \
           sum( mn(a,b,c,d)/2 for a,b,c,d in op(n-1,4) if a>1 and b==c )

def s1301(n):
    return sum( mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if a>1 and b>c ) + \
           sum( mn(a,b,c,d,e)/2 for a,b,c,d,e in op(n-1,5) if a>1 and b==c )

def s2001(n):
    return sum( mn(a,b,c) for a,b,c in op(n-1,3) if a > b > 1 ) + \
           sum( mn(a,b,c)/2 for a,b,c in op(n-1,3) if a == b > 1 )

def s2011(n):
    return sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a > b ) + \
           sum( mn(a,b,c,d)/2 for a,b,c,d in op(n-1,4) if a == b )

def s2200(n):
    return sum( 3*mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a > b > 1 ) + \
           sum( 2*mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a > b == 1 ) + \
           sum( 3*mn(a,b,c,d)/2 for a,b,c,d in op(n-1,4) if a == b > 1 ) + \
           sum( mn(a,b,c,d) for a,b,c,d in op(n-1,4) if a == b == 1 and c > d) + \
           sum( mn(a,b,c,d)/2 for a,b,c,d in op(n-1,4) if a == b == 1 and c == d)

def s2210(n):
    return s2111(n)

def s2201(n):
    return sum( 3*mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if a > b > 1 ) + \
           sum( 2*mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if a > b == 1 ) + \
           sum( 3*mn(a,b,c,d,e)/2 for a,b,c,d,e in op(n-1,5) if a == b > 1 ) + \
           sum( mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if a == b == 1 and c > d) + \
           sum( mn(a,b,c,d,e)/2 for a,b,c,d,e in op(n-1,5) if a == b == 1 and c == d)

def s2211(n):
    return sum( 3*mn(a,b,c,d,e,f) for a,b,c,d,e,f in op(n-1,6) if a > b ) + \
           sum( 3*mn(a,b,c,d,e,f)/2 for a,b,c,d,e,f in op(n-1,6) if a == b)

def s2300(n):
    return sum( mn(a,b,c,d,e) for a,b,c,d,e in op(n-1,5) if a>b ) + \
           sum( mn(a,b,c,d,e)/2 for a,b,c,d,e in op(n-1,5) if a==b )

def s2301(n):
    return sum( mn(a,b,c,d,e,f) for a,b,c,d,e,f in op(n-1,6) if a>b ) + \
           sum( mn(a,b,c,d,e,f)/2 for a,b,c,d,e,f in op(n-1,6) if a==b )

def triangles(n):
    return sum(f(n) for name,f in globals().iteritems() if re.match("s\d{4}", name))
