import math

class Constants:
    COPPER_RESISTIVITY = 1.72E-8  # Ohm * m
    #E0 = 8.85418782E-12  # C^2 / (N * m^2)
    E0 = 8.85418782E-12  # C^2 / (N * m^2)
    M0 = (4 * math.pi) * 1E-7  # (T * m) / A

    K = 1 / (4 * math.pi * E0)  # (N * m^2) / C^2; koulomb's konstant kekeke

    E = 1.60217662E-19  # C; this is the elementary charge. Not to be confused with E0, which is epsilon naught.

    MASS_PROTON = 1.6726219E-27  # Kg
    MASS_ELECTRON = 9.10938356E-31  # Kg
    MASS_NEUTRON = 1.674927471E-27  # Kg

    #COPPER_CUBE_SIZE = 1  # m (fake)
    #COPPER_ION_RADIUS = 0.2  # m (fake)
    COPPER_MASS = 1.055E-25 # Kg

    COPPER_CUBE_SIZE = 2.6E-8 # m
    COPPER_ION_RADIUS = 9.1E-11  # m