import math
import numpy as np

def my_equal(x,y):
    return abs(x-y) < 1e-8

def decode_glazer(s):
    """
    :Input: Glazer notation, e.g., :math:`a^-a^-a^-`
    :Output: Essential information about how to set the tiltings,
             i.e., :math:`\mathbf{k}_{\omega_{x,y,z}}`.
             And the relationship between tilting angles.

    """

    # s is a string of Glazer notation.
    if(len(s) != 6):
        print("Glazer notation has the form 'a-a-a-', your input seems incorect!")
        exit()

    k_ox = [1,1,1]
    if ( s[1] == '-' ):
       k_ox = [1,1,1]
    elif ( s[1] == '+' ):
       k_ox = [0,1,1]
    elif ( s[1] == '0' ):
       k_ox = [0,0,0]
    else: 
       print("Glazer notation has the form 'a-a-a-', your input seems incorect!")
       exit()

    k_oy = [1,1,1]
    if ( s[3] == '-' ):
       k_oy = [1,1,1]
    elif ( s[3] == '+' ):
       k_oy = [1,0,1]
    elif ( s[3] == '0' ):
       k_oy = [0,0,0]
    else:
       print("Glazer notation has the form 'a-a-a-', your input seems incorect!")
       exit()

    k_oz = [1,1,1]
    if ( s[5] == '-' ):
       k_oz = [1,1,1]
    elif ( s[5] == '+' ):
       k_oz = [1,1,0]
    elif ( s[5] == '0' ):
       k_oz = [0,0,0]
    else:
       print("Glazer notation has the form 'a-a-a-', your input seems incorect!")
       exit()

    k = np.array([k_ox,k_oy,k_oz]) * math.pi

    a = s[0]
    b = s[2]
    c = s[4]

    # the angles of rotation should be consistent with the letters of glazer notation
    relat1 = [a==b, b==c, c==a]

    return (k,relat1)
