##############
# Exercise 2.7
##############


def isCharged(x):
    if(x== 'D' or x== 'E' or x== 'R' or x== 'H' or x== 'K'):
        return True
    else:
        return False

#LIFWVMCYA

def isHydrophobic(x):
    if(x== 'L' or x== 'I' or x== 'F' or x== 'W' or x== 'V' or x== 'M' or x== 'Y' or x== 'A'):
        return True
    else:
        return False

def isNegativelyCharged(x):
    if(x== 'D' or x== 'E'):
        return True
    else:
        return False

def isAcid(x):
    if(x== 'D' or x== 'E'):
        return True
    else:
        return False

def isBasic(x):
    if(x== 'R' or x== 'H' or x== 'K'):
        return True
    else:
        return False

def isPositivelyCharged(x):
    if(x== 'R' or x== 'K' or x== 'H'):
        return True
    else:
        return False

def isAromatic(x):
    if(x== 'F' or x== 'W' or x== 'Y' or x== 'H'):
        return True
    else:
        return False

def isProline(x):
    if(x== 'P'):
        return True
    else:
        return False

def containsSulfur(x):
    if(x== 'C' or x== 'M'):
        return True
    else:
        return False

#DERKHNQSTY

def isPolar(x):
    if(x== 'D' or x== 'E' or x== 'R' or x== 'K' or x== 'H' or x== 'N' or x== 'Q' or x== 'S' or x== 'T' or x== 'Y'):
        return True
    else:
        return False
    