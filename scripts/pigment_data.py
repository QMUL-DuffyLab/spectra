#!/usr/bin/env python3

# note : usage is pigment_data[ligand_code][state][variable]
pigment_data = {
        'CLA': { # note that chlorophyll Qy := S1, Qx := S2
            'D' : 4.5,
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
            'D' : 3.6,
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
            }
        }

# check where these are from!!!
# muh and renger 
site_energies = {
        'CHL601' : 15405,
        'CLA602' : 14940,
        'CLA603' : 14850,
        'CLA604' : 14820,
        'CHL605' : 15465,
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
        }
