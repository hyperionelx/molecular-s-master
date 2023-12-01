import numpy as np


def SI2au(SI, PhysicalQuantity):
    conversion_factors = {
        'm': 9.1093837015e-31,  # 质量 mass
        'e': 1.602176634e-19,   # 电荷 charge
        'a0': 5.29177210903e-11,  # 长度 length
        'E': 4.3597447222071e-18,  # 能量 energy
        'epsilon': 1.11265005545e-10,  # 介电常数 permittivity
        't': 2.4188843265857e-17,  # 时间 time
        'f': 1 / 2.4188843265857e-17,  # 频率 frequency
        'omega': (1 / 2.4188843265857e-17) * 2 * np.pi,  # 角频率
        'hbar': 1.054571817e-34,  # 约化普朗克常数/角动量
        'P': 1.99285191410e-24,  # 动量
        'ED': 8.4783536255e-30,  # 电偶极矩：e*a0
        'EQ': 4.4865515246e-40,  # 电四极矩：e*a0^2
        'MD': 1.85480201566e-23,  # 磁偶极矩：e*a0
        'EF': 5.14220674763e11,  # 电场强度
        'EP': 27.211386245988,  # 电位
        'I': 6.4364099007e19,  # 光强
        'v': 2.18769126364e6   # 速度
    }

    if PhysicalQuantity in conversion_factors:
        return SI/conversion_factors[PhysicalQuantity]
    else:
        raise ValueError('PhysicalQuantity 请输入: m, e, a0, E, epsilon, t, f, omega, hbar, P, ED, EQ, MD, EF, EP, I, v')
