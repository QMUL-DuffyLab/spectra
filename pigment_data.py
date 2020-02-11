#!/usr/bin/env python3

# note : usage is pigment_data[ligand_code][state][variable]
pigment_data = {
        'CLA': { # note that chlorophyll Qy := S1, Qx := S2
            'S1': {
                'energy': 14900.0,
                'lifetime': 4.0,
                'reorg': 0.0, # actually will need to read this in
                },
            'S2': {
                'energy': 16300.0,
                'lifetime': 4.0, # probably untrue
                'reorg': 0.0,
                }
            },
        'CLC': {
            'S1': {
                'energy': 15970.0,
                'lifetime': 4.0,
                'reorg': 0.0,
                },
            'S2': {
                'energy': 17330.0,
                'lifetime': 4.0,
                'reorg': 0.0,
                }
            },
        'A86': {
            'S1': {
                'energy': 14900.0,
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
