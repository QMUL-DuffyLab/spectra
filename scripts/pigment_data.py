#!/usr/bin/env python3

# note : usage is pigment_data[ligand_code][state][variable]
pigment_data = {
        'CLA': { # note that chlorophyll Qy := S1, Qx := S2
            'D' : 4.0,
            'S1': {
                'energy': 14900.0,
                'lifetime': 4.0,
                'reorg': 0.0, # actually will need to read this in
                },
            'S2': {
                # 'energy': 16300.0,
                'energy': 15900.0,
                # 'lifetime': 4.0, # probably untrue
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'CHL': { # chlorophyll b
            'D' : 3.4,
            'S1': {
                'energy': 15300.0,
                'lifetime': 4.0,
                'reorg': 0.0, # actually will need to read this in
                },
            'S2': {
                'energy': 16700.0,
                # 'lifetime': 4.0, # probably untrue
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'KC1': {
            'S1': {
                # 'energy': 15970.0,
                'energy': 15600.0,
                'lifetime': 4.0,
                'reorg': 0.0,
                },
            'S2': {
                # 'energy': 17330.0,
                'energy': 16700.0,
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'KC2': {
            'S1': {
                'energy': 15600.0,
                'lifetime': 4.0,
                'reorg': 0.0,
                },
            'S2': {
                'energy': 16700.0,
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'A86': {
            'S1': {
                'energy': 15870.0,
                'lifetime': 0.01,
                'reorg': 0.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006, # 60fs
                'reorg': 0.0,
                }
            },
        'DD6': {
            'S1': {
                'energy': 14900.0,
                'lifetime': 0.01,
                'reorg': 0.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        'LUT': {
            'D' : 0.1, # kieran carotenoid pathway table 1
            'S1': {
                'energy': 15400.0,
                'lifetime': 0.01,
                'reorg': 450.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        'NEX': {
            'D' : 0.15,
            'S1': { # no real data for this or XAT
                'energy': 15200.0,
                'lifetime': 0.01,
                'reorg': 450.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        'XAT': {
            'D' : 0.18,
            'S1': { # no real data for this or XAT
                'energy': 15200.0,
                'lifetime': 0.01,
                'reorg': 450.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        }

# check where these are from!!!
# muh and renger
site_energies = {
        'CHL601' : 15405,
        'CLA602' : 14940,
        'CLA603' : 14850,
        'CLA604' : 14820,
        'CHL605' : 15465,
        'CHL606' : 15385,
        'CHL607' : 15225,
        'CHL608' : 15215,
        'CHL609' : 15475,
        'CLA610' : 14790,
        'CLA611' : 14950,
        'CLA612' : 14940,
        'CLA613' : 14840,
        'CLA614' : 14940,
        'LUT620' : 14940, # testing
        'LUT621' : 14940,
        # 'LUT620' : 14500, # guess from chris - both luteins 14500
        # 'LUT621' : 14500,
        # 'LUT620' : 19538, # 620/621 from kieran carotenoid paper
        # 'LUT621' : 19748,
        }

novod_energies = {
        'CHL601' : 15890,
        'CLA602' : 15160,
        'CLA603' : 15283,
        'CLA604' : 15460,
        'CHL605' : 15679,
        'CHL606' : 15851,
        'CHL607' : 15712,
        'CHL608' : 15763,
        'CHL609' : 15721,
        'CLA610' : 15073,
        'CLA611' : 15115,
        'CLA612' : 15097,
        'CLA613' : 15175,
        'CLA614' : 15264,
        'LUT620' : 16940, # testing
        'LUT621' : 16940,
        }

dict_from_vangelis = {
    'CHL601' : [262, 266, 270],
    'CLA602' : [232, 238, 244],
    'CLA603' : [233, 239, 245],
    'CLA604' : [250, 252, 254],
    'CHL605' : [256, 258, 260],
    'CHL606' : [257, 259, 261],
    'CHL607' : [263, 267, 271],
    'CHL608' : [264, 268, 272],
    'CHL609' : [265, 269, 273],
    'CLA610' : [234, 240, 246],
    'CLA611' : [235, 241, 247],
    'CLA612' : [236, 242, 248],
    'CLA613' : [237, 243, 249],
    'CLA614' : [251, 253, 255],
    'XAT622' : [283, 284, 285],
    'NEX623' : [280, 281, 282],
    'LUT620' : [274, 276, 278],
    'LUT621' : [275, 277, 279],
    }
