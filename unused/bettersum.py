def msum(iterable):
    "Full precision summation using multiple floats for intermediate values"

    partials = []               # sorted, non-overlapping partial sums
    for x in iterable:
        newpartials = []
        for y in partials:
            if abs(x) < abs(y):
                x, y = y, x
            hi = x + y
            lo = y - (hi - x)
            if lo:
                newpartials.append(lo)
            x = hi
        partials = newpartials + [x]
    return sum(partials, 0.0)


from math import frexp

def lsum(iterable):
    "Full precision summation using long integers for intermediate values"
    # Transform (exactly) a float to m * 2 ** e where m and e are integers.
    # Adjust (tmant,texp) and (mant,exp) to make texp the common exponent.
    # Given a common exponent, the mantissas can be summed directly.

    tmant, texp = 0, 0
    for x in iterable:
        mant, exp = frexp(x)
        mant, exp = int(mant * 2.0 ** 53), exp-53
        if texp < exp:
            mant <<= exp - texp
        elif texp > exp:
            tmant <<= texp - exp
            texp = exp
        tmant += mant
    return float(str(tmant)) * 2.0 ** texp


from decimal import getcontext, Decimal, Inexact
getcontext().traps[Inexact] = True

def dsum(iterable):
    "Full precision summation using Decimal objects for intermediate values"
    # Transform (exactly) a float to m * 2 ** e where m and e are integers.
    # Convert (mant, exp) to a Decimal and add to the cumulative sum.
    # If the precision is too small for exact conversion and addition,
    # then retry with a larger precision.

    total = Decimal(0)
    for x in iterable:
        mant, exp = frexp(x)
        mant, exp = int(mant * 2.0 ** 53), exp-53
        while True:
            try:
                newtotal = total + mant * Decimal(2) ** exp
            except Inexact:
                getcontext().prec += 1
            else:
                total = newtotal
                break
    return float(total)


from random import random, normalvariate, shuffle

def test(nvals):
    for j in range(1000):
        vals = [7, 1e100, -7, -1e100, -9e-20, 8e-20] * 10
        s = 0
        for i in range(nvals):
            v = normalvariate(0, random())**7 - s
            s += v
            vals.append(v)
        shuffle(vals)
        assert msum(vals) == lsum(vals) == dsum(vals)
        print('.', end=' ')
    print('Tests Passed')

if __name__ == '__main__':
    test(200)
    
