from math import log10 , floor

def round_it(x, sig):
    return round(x, sig-int(floor(log10(abs(x))))-1)
