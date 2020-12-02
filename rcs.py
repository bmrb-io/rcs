from scipy.special import ellipk, ellipe

def get_ring_current_shift(ring_center, ring_normal, proton_position):

    r_ring = 1.39 #1.18 for 5 atoms
    zs_proton, rho_proton = convert_to_polar(
        ring_center, ring_normal, proton_position, r_ring
    )
    total_shift = 0
    for z_proton in zs_proton:
        if check_proximity(rho_proton, z_proton):
            if rho <= 1.0:
                total_shift += 10
            else:
                total_shift -= 10
        else:
            dist_a = (1.0 + rho_proton) ** 2 + z_proton ** 2
            dist_b = (1.0 - rho_proton) ** 2 + z_proton ** 2
            ksq = 4.0 * rho_proton / dist_a
            term = ellipk(ksq) + ellipe(ksq) * (
                -1 + 2 * (1 - rho_proton) / dist_b
            )
            total_shift += term / (dist_a ** 0.5)
    
    return total_shift


    


def convert_to_polar(ring_center, ring_normal, proton_position, r_ring):

    proton_vector = [proton_position[i] - ring_center[i] for i in range(3)]
    z_proton = sum([proton_vector[i] * ring_normal[i] for i in range(3)])
    rho_proton = (
        sum([proton_position[i] ** 2 for i in range(3)]) - z_proton ** 2
    )
    rho_proton = rho_proton ** 0.5

    z_top_proton = (z_proton - 0.64) / r_ring
    z_bot_proton = (z_proton + 0.64) / r_ring
    zs_proton = [z_top_proton, z_bottom_proton]
    rho_proton = rho_proton / r_ring

    return zs_proton, rho_proton

def get_rcfac(ring_type):
    rc_fac_dict = {
        'PHE': 1.00
        'TRP5': 0.56
        'TRP6': 1.04
        'TYR': 0.94
    }

    return rc_fac_dict[ring_type]

def check_proximity(rho_proton, z_proton):

    dist_b = (1.0 - rho_proton) ** 2 + z_proton ** 2
    return (dist_b <= 0.1)