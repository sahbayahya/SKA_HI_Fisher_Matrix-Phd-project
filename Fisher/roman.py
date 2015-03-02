# Roman numerals
# I = 1
# V = 5
# X = 10
# L = 50
# C = 100
# D = 500
# M = 1000

def roman(x):
    s = ''
    n = x / 1000
    s = 'M' * n

    x = x % 1000
    n = x / 100
    if n == 9:
        s += 'CM'
    elif n == 4:
        s += 'CD'
    else:
        if n >= 5:
            s += 'D'
            n -= 5
        s += 'C' * n

    x = x % 100
    n = x / 10
    if n == 9:
        s += 'XC'
    elif n == 4:
        s += 'XL'
    else:
        if n >= 5:
            s += 'L'
            n -= 5
        s += 'X' * n

    x = x % 10
    n = x / 1
    if n == 9:
        s += 'IX'
    elif n == 4:
        s += 'IV'
    else:
        if n >= 5:
            s += 'V'
            n -= 5
        s += 'I' * n

    return s
