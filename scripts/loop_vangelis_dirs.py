import os

input_dirs = [
        'structures/tLHCII_NpH/1/C',
        'structures/tLHCII_NpH/1/E',
        'structures/tLHCII_NpH/1/H',
        'structures/tLHCII_NpH/2/H',
        'structures/tLHCII_NpH/3/C',
        'structures/tLHCII_NpH/3/E',
        'structures/tLHCII_NpH/3/H',
        'structures/tLHCII_NpH/4/C',
        'structures/tLHCII_NpH/4/E',
        'structures/tLHCII_NpH/4/H',
        'structures/tLHCII_NpH/5/C',
        'structures/tLHCII_NpH/5/E',
        'structures/tLHCII_NpH/5/H',
        'structures/tLHCII_NpH/6/C',
        'structures/tLHCII_NpH/6/E',
        'structures/tLHCII_NpH/6/H',
        ]

output_dirs = [
        'out/tLHCII_NpH/1/C',
        'out/tLHCII_NpH/1/E',
        'out/tLHCII_NpH/1/H',
        'out/tLHCII_NpH/2/H',
        'out/tLHCII_NpH/3/C',
        'out/tLHCII_NpH/3/E',
        'out/tLHCII_NpH/3/H',
        'out/tLHCII_NpH/4/C',
        'out/tLHCII_NpH/4/E',
        'out/tLHCII_NpH/4/H',
        'out/tLHCII_NpH/5/C',
        'out/tLHCII_NpH/5/E',
        'out/tLHCII_NpH/5/H',
        'out/tLHCII_NpH/6/C',
        'out/tLHCII_NpH/6/E',
        'out/tLHCII_NpH/6/H',
        ]
average_dirs = [
        'out/tLHCII_NpH/1/C/C',
        'out/tLHCII_NpH/1/E/E',
        'out/tLHCII_NpH/1/H/H',
        'out/tLHCII_NpH/2/C/C',
        'out/tLHCII_NpH/2/E/E',
        'out/tLHCII_NpH/2/H/H',
        'out/tLHCII_NpH/3/C/C',
        'out/tLHCII_NpH/3/E/E',
        'out/tLHCII_NpH/3/H/H',
        'out/tLHCII_NpH/4/C/C',
        'out/tLHCII_NpH/4/E/E',
        'out/tLHCII_NpH/4/H/H',
        'out/tLHCII_NpH/5/C/C',
        'out/tLHCII_NpH/5/E/E',
        'out/tLHCII_NpH/5/H/H',
        'out/tLHCII_NpH/6/C/C',
        'out/tLHCII_NpH/6/E/E',
        'out/tLHCII_NpH/6/H/H',
        ]

for i in range(len(input_dirs)):
    # print("./scripts/run.py -i {} -o {} -f 0 -T 300.0 -c 2 -ps 1 -pc 0".format(input_dirs[i], output_dirs[i]))
    # print("python scripts/average_spectra.py -r 1 -i {}".format(output_dirs[i]))
    # os.system("./scripts/run.py -i {} -o {} -f 0 -T 300.0 -c 2 -ps 1 -pc 0".format(input_dirs[i], output_dirs[i]))
    os.system("python scripts/average_spectra.py -r 1 -i {}".format(average_dirs[i]))
