from astropy.table import Table
import clustr
import reglib
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
        'filename',
        type=str,
        help='Fits file used for scaling relation'
        )

def main():
    args = parser.parse_args()
    fits_file = args.filename
    catalog = Table.read(fits_file)

    xlabel = 'lambda'
    ylabel = 'lr500'

    # Make any necessary flag cuts, etc.
    # NOTE: You probably don't want this cut for the real run!!
    catalog = catalog[catalog['lambda'] > 30]

    x = catalog[xlabel]
    y = catalog[ylabel]

    # Scale it to log
    # ...

    #def run_lrgs(x, y, err_x, err_y, _xycov=None, nmc=500, dirichlet=True):
    intercept, slope, sigma = reglib.runlrgs(x, y, err_x, err_y)

    return

if __name__ == '__main__':
    main()
