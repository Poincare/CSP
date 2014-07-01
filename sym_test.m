x = sym('x', 'real')
y = sym('y', 'real')
expr = (x+y)^100
pretty(expand(expr))