# algorithms of primes, prime powers, etc.

from cwc.utils.logging import is_debug, logger, set_debug
from itertools import permutations, combinations
from random import randint
from math import gcd

import sys
import galois

from core import binom

def moebius(n):
    """Compute the Möbius function mu(n) for a positive integer n."""
    if n <= 0:
        raise ValueError("n must be a positive integer")
    if n == 1:
        return 1
    p = 2
    count = 0
    while p * p <= n:
        if n % p == 0:
            n //= p
            count += 1
            if n % p == 0:
                return 0
        p += 1
    if n > 1:
        count += 1
    return -1 if count % 2 == 1 else 1

def LengthLimitedLyndonWords(s, n):
    """Duval's algorithm to generate nonempty Lyndon words of length <= n over an s-symbol alphabet."""
    w = [-1]
    while w:
        w[-1] += 1
        yield list(w)
        m = len(w)
        while len(w) < n:
            w.append(w[-m])
        while w and w[-1] == s - 1:
            w.pop()

def LyndonWordsWithLength(s, n):
    """Generate Lyndon words of length exactly n over an s-symbol alphabet."""
    if n == 0:
        yield []
    else:
        for w in LengthLimitedLyndonWords(s, n):
            if len(w) == n:
                yield w

def CountLyndonWords(s, n):
    """Count the number of length-n Lyndon words over s symbols using Möbius inversion."""
    if n == 0:
        return 1
    total = 0
    for i in range(1, n + 1):
        if n % i == 0:
            total += moebius(n // i) * (s ** i)
    return total // n

# computes factor n*(n-1)*(n-2)*...*(n-k+1) 
def factor(n, k=-1):
    """
    A fast way to calculate factor
    """
    if k==-1:
        k=n
    logger.debug('factor =',n,'!/',(n-k),'!',)
    if 0 <= k <= n:
        ntok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            n -= 1
            logger.debug(ntok)
        logger.debug(ntok)
        return ntok 
    else:
        return 0

def primePowersSieve(limit):
    a = [True] * limit                          # Initialize the primality list
    a[0] = a[1] = False
    for (i, isprime) in enumerate(a):
        if isprime:
            yield i
            for n in range(i, limit//i):     # Mark factors non-prime
                #print('prime ',i,' n=',n,' where ',n%i) 
                if n%i!=0:
                    a[i*n] = False

def nextPrimePower(n):
    for q in primePowersSieve(2*n):
        if q>=n:
            return q

def prevPrimePower(n):
    ret=2
    for q in primePowersSieve(n+2):
        logger.debug('prime',q)
        if q>n:
            return ret
        else:
            ret=q

#convert an integer to a polinom of  base
def int2poly(x, base, max_deg=-1):
    digits = []
    xx=x
    while x:
        digits.append(x % base)
        x = int(x / base)
    
    if max_deg!=-1:
        while len(digits)<max_deg:
            digits.append(0)

#    digits.reverse()
    logger.debug(xx, 'in base',base,'is',digits)
    return digits

# prime power predicate
def findWitness(n, k=5): # miller-rabin
    s, d = 0, n-1
    while d % 2 == 0:
        s, d = s+1, d/2
    for i in range(k):
        a = randint(2, n-1)
        x = pow(a, int(d), n)
        if x == 1 or x == n-1: continue
        for r in range(1, s):
            x = (x * x) % n
            if x == 1: return a
            if x == n-1: break
        else: return a
    return 0

# returns p,k such that n=p**k, or 0,0
# assumes n is an integer greater than 1
def primePower(n):
    def checkP(n, p):
        k = 0
        while n > 1 and n % p == 0:
            n, k = n / p, k + 1
        if n == 1: return p, k
        else: return 0, 0
    if n % 2 == 0: return checkP(n, 2)
    q = n
    while True:
        a = findWitness(q)
        if a == 0: return checkP(n, q)
        d = gcd(pow(a,q,n)-a, q)
        if d == 1 or d == q: return 0, 0
        q = d

def irrednum(q,deg):
    return CountLyndonWords(q,deg)

"""irrednum= {
 2: [1, 2, 3, 6, 9, 18, 30, 56, 99, 186, 335, 630, 1161, 2182, 4080, 7710, 14532, 27594, 52377],
 3: [3, 8, 18, 48, 116, 312, 810, 2184, 5880],
 4: [6, 20, 60, 204, 670, 2340, 8160, 29120, 104754],
 5: [10, 40, 150, 624, 2580, 11160, 48750, 217000, 976248],
 7: [21, 112, 588, 3360, 19544, 117648, 720300, 4483696, 28245840],
 8: [28, 168, 1008, 6552, 43596, 299592, 2096640, 14913024, 107370900],
 9: [36, 240, 1620, 11808, 88440, 683280, 5380020, 43046640, 348672528],
 11: [55, 440, 3630, 32208, 295020, 2783880, 26793030, 261994040, 261994040],
 13: [78, 728, 7098, 74256, 804076, 8964072, 101962770, 1178277464],
 16: [120, 1360, 16320, 209712, 2795480, 38347920, 536862720, 7635496960],
 17: [136, 1632, 20808, 283968, 4022064, 58619808, 871959240, 13176430176],
 19: [361, 7201, 137161, 2613241, 49651921, 943523641, 17926956361, 340614647281, 6471678428641],
 23: [253, 4048, 69828, 1287264, 24670536, 486403632, 9788838180, 200128072144],
 25: [300, 5200, 97500, 1953120, 40687400, 871930800, 19073437500, 423855250000],
 27: [351, 6552, 132678, 2869776, 64566684, 1494336168, 35303625630, 847288607256]}

irredpol={
	2: [['0', '1'], 
        ['01'], 
        ['001', '011'],
        ['0001', '0011', '0111'], 
        ['00001', '00011', '00101', '00111', '01011', '01111'],
        ['000001', '000011', '000101', '000111', '001011', '001101', '001111', '010111', '011111'],
        ['0000001', '0000011', '0000101', '0000111', '0001001', '0001011', '0001101', '0001111', '0010011']]
}"""


class FieldElementWrapper:
    def __init__(self, val, field_wrapper):
        self.val = val
        self.field_wrapper = field_wrapper

    @property
    def exp_coefs(self):
        if self.field_wrapper.n == 1:
            return [int(self.val)]
        else:
            return list(self.val.vector())[::-1]

    def symbol(self):
        return int(self.val)

    def strPoly(self):
        coefs = self.exp_coefs
        ret = ''
        for i, c in enumerate(coefs):
            if c != 0:
                rr = ''
                if c != 1 or i == 0:
                    rr += str(c)
                if i > 0:
                    rr += "x"
                if i > 1:
                    rr += "^" + str(i)
                if ret != '' and rr != '':
                    ret = rr + " + " + ret
                else:
                    ret += rr
        if ret == '':
            ret = '0'
        return ret

    def __hash__(self):
        return hash(int(self.val))

    def __eq__(self, other):
        if isinstance(other, FieldElementWrapper):
            return self.val == other.val
        return self.val == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        other_val = other.val if isinstance(other, FieldElementWrapper) else other
        return FieldElementWrapper(self.val + other_val, self.field_wrapper)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        other_val = other.val if isinstance(other, FieldElementWrapper) else other
        return FieldElementWrapper(self.val - other_val, self.field_wrapper)

    def __rsub__(self, other):
        other_val = other.val if isinstance(other, FieldElementWrapper) else other
        return FieldElementWrapper(other_val - self.val, self.field_wrapper)

    def __mul__(self, other):
        other_val = other.val if isinstance(other, FieldElementWrapper) else other
        return FieldElementWrapper(self.val * other_val, self.field_wrapper)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other_val = other.val if isinstance(other, FieldElementWrapper) else other
        return FieldElementWrapper(self.val / other_val, self.field_wrapper)

    def __rtruediv__(self, other):
        other_val = other.val if isinstance(other, FieldElementWrapper) else other
        return FieldElementWrapper(other_val / self.val, self.field_wrapper)

    def __pow__(self, power):
        return FieldElementWrapper(self.val ** power, self.field_wrapper)

    def __str__(self):
        return str(self.val)

    def __repr__(self):
        return f"GFElement({int(self.val)})"


class GaloisField:
    def __init__(self, p, n=1, coefs=None):
        self.p = p
        self.n = n
        self.dim = p ** n
        if n > 1 and coefs:
            base_field = galois.GF(p)
            poly = galois.Poly(coefs[::-1], field=base_field)
            self.GF = galois.GF(self.dim, irreducible_poly=poly)
        else:
            self.GF = galois.GF(self.dim)
            
        self._elements = []
        for idx in range(self.dim):
            self._elements.append(FieldElementWrapper(self.GF(idx), self))


    def __getitem__(self, idx):
        if idx < self.dim and idx >= -self.dim:
            return self._elements[idx]
        else:
            raise IndexError("Error, element out of bounds.")

    def __len__(self):
        return self.dim

    def __iter__(self):
        return iter(self._elements)

    def print(self):
        print(f"Galois Field of order {self.dim}")
        print(f"Base prime: {self.p}")
        print(f"Extension degree: {self.n}")
        if self.n > 1:
            print(f"Irreducible polynomial: {self.GF.irreducible_poly}")
 

# q is a prime
# FR(q^i)
# from http://fchabaud.free.fr/English/Poly/triform.php?FDeg=2&LDeg=1000&Output=HGF%28p%29%2Fmalcolm#D11
def getGF(i, p):
    if i==1:
        return GaloisField(p)
    elif p==2: 
        if i==2: # 1 + x + x^2 
            return GaloisField(2, 2, [1, 1, 1])
        elif i==3: # 1 + x + x^3
            return GaloisField(2, 3, [1, 1, 0, 1])
        elif i==4: # 1 + x + x^4 
            return GaloisField(2, 4, [1, 1, 0, 0, 1])
        elif i==5: # 1 + x^2 + x^5 
            return GaloisField(2, 5, [1, 0, 1, 0, 0, 1])
        elif i==6: # 1 + x + x^6 
            return GaloisField(2, 6, [1, 1, 0, 0, 0, 0, 1])
        elif i==7: # 1 + x^6 + x^7 
            return GaloisField(2, 7, [1, 0, 0, 0, 0, 0, 1, 1])
        elif i==8: # X8+X7+X6+X+1
            return GaloisField(2, 8, [1, 1, 0, 0, 0, 0, 1, 1, 1])
        elif i==9: # 1 + x^5 + x^9 
            return GaloisField(2, 9, [1, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        elif i==10: # 1 + x^7 + x^10 
            return GaloisField(2, 10, [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 ])
    elif p==3:
        if i==2:  # 2 + x + x^2 
            return GaloisField(3, 2, [2, 1, 1])
        elif i==3: # 1 + 2x + x^3
            return GaloisField(3, 3, [1, 2, 0, 1])
        elif i==4: # 2 + x + x^4
            return GaloisField(3, 4, [2, 1, 0, 0, 1])
    elif p==5:
        if i==2:
            return GaloisField(5, 2, [3, 2, 1])
        if i==3: 
            return GaloisField(5, 3, [2, 3, 0, 1])
        if i==4: # 2 + 3 x + 3 x^2 + 2 x^3 + x^4
            return GaloisField(5, 4, [2, 3, 3, 2, 1])
    elif p==7:
        if i==2:# 3 + x + x^2
            return GaloisField(7, 2, [3, 1, 1])
    elif p==11:
        if i==2:# 7 + x + x^2
            return GaloisField(11, 2, [7, 1, 1])
    elif p==13:
        if i==2:# 2 + x + x^2
            return GaloisField(13, 2, [2, 1, 1])
    elif p==17:
        if i==2:# 3 + x + x^2
            return GaloisField(17, 2, [3, 1, 1])


def poly2list(poly,k,GF):
    block_ids=[]
    for x in range(k):
        digit=GF[0]
        for j,c in enumerate(poly):
            aa=pow(GF[x],j)*GF[c]
            logger.debug('digit:',digit.symbol(),'+(',pow(GF[x],j).symbol(),'=',GF[x].symbol(),'^',j,')*',GF[c].symbol(),'=',aa.symbol(),')') 
            digit+=aa
        #digit= self.GF.evaluate(poly, self.GF[x]) 
        logger.debug('poly=',poly,'item=',GF[x].symbol(),'->',digit.symbol())
        block_ids.append(digit.symbol())
    return block_ids

def is_trivial_params(m, D, w):
    """
    Checks if a code with parameters (m, D, w) is a trivial case
    that can be generated on demand without saving a physical file.
    """
    if w > m:
        return True
    if w + D // 2 > m and m >= w:
        return True
    if w == D // 2:
        return True
    if D == 2:
        return True
    if m < 2 * w:
        return True
    return False

