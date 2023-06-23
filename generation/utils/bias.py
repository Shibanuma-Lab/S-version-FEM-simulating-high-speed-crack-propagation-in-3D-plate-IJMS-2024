import numpy as np

def bias(num_elements, bias, length):
    rb = 0
    xsol = 0
    
    if bias == 1 or bias == 1.:
        xsol = [length/num_elements]
        rb = 1.
    else:
        rb = bias**(1/(num_elements - 1))
        xsol = np.roots([(1-rb**num_elements)/(1-rb), -length])
        xsol = [x.real for x in xsol if np.isclose(x.imag, 0, atol=1e-8)]
        
    return xsol, rb