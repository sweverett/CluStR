#-------------------------------------------------------
# We'll define useful classes here

class Config(object):
    '''
    Used for CluStR config processing
    Some options:
    - scale_luminosity: Divide luminosity columns by E(z)^-3/2
    '''
    _required_keys = []
    _default_run_name = 'clustr'
    def __init__(self, config_file, run_options):
        with open(config_file, 'r') as stream:

            self._config = yaml.safe_load(stream)
# TODO: was prefix (in arparse) changing to run_name
        # TODO: Any logical parsing here:
        #...
        self.run_options = run_options
        if run_options.run_name is None:
            self.run_name = _default_run_name
        else:
            self.run_name = run_options.run_name
        return

    # The following are so we can access the config
    # values similarly to a dict
    def __getitem__(self, key):
        return self._config.__dict__[key]

    def __setitem__(self, key, value):
        self._config.__dict__[key] = value

    def __delitem__(self, key):
        del self._config.__dict__[key]

    def __contains__(self, key):
        return key in self._config.__dict__

    def __len__(self):
        return len(self._config.__dict__)

    def __repr__(self):
        return repr(self._config.__dict__)
    def run_name()

class ArgumentParser:
    def __init__(self,)

    def parse_opts():
    ''' Parse command line arguments '''
    parser = argparse.ArgumentParser()
    # Required argument for catalog
    parser.add_argument('catalog', help='FITS catalog to open')
    # Required arguement for axes
    valid_axes = ['l500kpc', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500',
                  'tr500', 'tr500cc', 'lambda']
    parser.add_argument('y', help='what to plot on y axis', choices=valid_axes)
    parser.add_argument('x', help='what to plot on x axis', choices=valid_axes)
    # Optional argument for file prefix
    parser.add_argument('-p', '--prefix', help='prefix for output file')
    # Optional arguments for any flag cuts
    # FIX: in the future, make an allowed choices vector work!
    parser.add_argument(
        '-f',
        '--flags',
        nargs='+',
        type=str,
        help=(
            'Input any desired flag cuts as a list of flag names '
            '(with "" and no spaces!)'
        )
    )
    # Optional argument for which files to be saved
    # FIX: Implement!

    return parser.parse_args()


def Ez(z)
    Om = 0.3
    H_0 = 0.7
    h = H_0/100
    return np.sqrt(Om*(1.+z)**3 + h)


class Flag:
    def __init__(self,flag,boolean,cut_or_range):
        self.flag = flag
        self.boolean = boolean
        self.cut_or_range = cut_or_range
    def check_flag(self):
        if flag.lower() in boolean:
            try:
                #what would I get from this statement?
                bool_type = config.bools[flag+'_bool_type']
            except KeyError:
                raise TypeError()
            if bool_type is True or bool_type is False:
                return 'bool'
            else :
                raise TypeError()
        elif flag.lower() in cut_or_range:
            try:
                cut_type = config.bools[flag+'_bool_type']
            except KeyError:
                try:
                    range_type = config.ranges[flag+'_range_type']
                except KeyError:
                    raise TypeError()
                if range_type == 'inside' or range_type == 'outside':
                    return 'range'
                else:
                    raise TypeError()
        if cut_type == 'above' or cut_type == 'below':
                return 'cutoff'
        else:
            raise TypeError()
    raise TypeError()

#inheritance would be good here
#would inheritance be good here since the only common input is flags?
class Catalog(Flag):
    #no need for 2 labels anymore
    column_labels = {
        'lambda',
        'l500kpc',
        'lr2500',
        'lr500',
        'lr500cc',
        't500kpc',
        'tr2500',
        'tr500',
        'tr500cc'
    }
    def __init__(self,data,flags,flag,run_options):
        self.data = data
        self.flags = flags
        self.run_options = run_options
flag_class = Flag(INPUT)
    def create_cuts(self):
        mask = np.zeros(len(data), dtype=bool)
        for flag in flags:
            try:
                flag_type = flag_class.check_flag()
            except TypeError()
                continue
            if flag_type == 'bool':
                bool_type = config.bools[flag+'_bool_type']
                pass
                if isinstance(bool_type, bool):
                    #cut = data[flag] ==
                    pass
                else:
                    #print()
                    continue
            elif flag_type == 'cutoff': #why is this highlighted for if but not elif?
                cutoff = config.cutoffs[flag+'_cut']
                cut_type = config.cutoffs[flag+'_cut_type']

                if cut_type == 'above':
                    cut = data[flag] < cutoff
                elif cut_type == 'below':
                    cut = data[flag] > cutoff
                else:
                    print()
                    continue
            elif flag_type == 'range':
                fmin = config.ranges[flag+'_range_min']
                fmax = config.ranges[flag+'_range_max']
                range_type = config.ranges[flag+'_range_type']
                if range_type == 'inside':
                    cut = (data[flag] < fmax) & (data[flag > fmin])
                elif range_type == 'outside':
                    cut = (data[flag] > fmin) & (data[flag] < fmax)
                else:
                    print()
                    continue
            mask |= cut
            print()
    def get_data(self):
        #hdulist = fits.open(options.catalog)
        #data = hdulist[1].data
        #
        label_x = run_options.x
        label_y = run_options.y
        x = data[label_x]
        y = data[label_y]

        # Number of original data
        N = np.size(x)
        #where does the l come from? should it also be lambda?
        if label_x[0] == 'l' and label_y != 'lambda':
            x /= Ez(data['redshift'])
        if label_y[0] == 'l' and label_x != 'lambda':
            y /= Ez(data['redshift'])

        #how should error be more accurate here?

        flags = options.flags
        if flags is not None:
            # FIX: Should be more error handling than this!
            # FIX: Should write method to ensure all the counts are what we expect

            mask = create_cuts(self)
            x[mask] = -1
            y[mask] = -1

            print (
                'NOTE: `Removed` counts may be redundant, '
                'as some data fail multiple flags.'
            )

            # Take rows with good data, and all flagged data removed
            good_rows = np.all([x != -1, y != -1], axis=0)

            x = x[good_rows]
            y = y[good_rows]
            x_err = x_err[good_rows]
            y_err = y_err[good_rows]

            print 'Accepted {} data out of {}'.format(np.size(x), N)
        if np.size(x) == 0:
            print (
                '\nWARNING: No data survived flag removal. '
                'Suggest changing flag parameters in `param.config`.'
                '\n\nClosing program...\n'
            )
            raise SystemExit(2)

        print 'mean x error:', np.mean(x_err)
        print 'mean y error:', np.mean(y_err)

        hdulist.close()

        return (x, y, x_err, y_err)

class Data:
    def __init__(self,x, y, x_err, y_err, x_obs, y_obs, nmc):
        self.x = x
        self.y = y
        self.x_err = x_err
        self.y_err = y_err
        self.x_obs = x_obs
        self.y_obs = y_obs
        self.nmc = nmc

    def scale(self):
        log_x = np.log(self.x)
        x_piv = np.median(log_x)
        log_y = np.log(self.y)
        return (log_x-x_piv, log_y, x_err/self.x, y_err/self.y, x_piv)


#fit that accepts data/ should that go here or in the fitter class?

#replace check_prefix with run_name
#dont believe this is needed anymore?
def run_name:
    def __init__(self, run_options, parameters):
        self.options = options
        self.parameters = parameters
        pass
#class function of config
#inheritance here (from Data class)
class Fitter:
    pass
#inheritance here (from Data/Fitter class)
class SaveData:
    def __init__(self, options, parameters, , ,)
        pass
#class function in fitter class

#-------------------------------------------------------
# We'll write the main function here

def main():  # pylint: disable=missing-docstring

    # Parse all inputted options, regardless of param.config file
    options = parse_opts()

    # Set useful parameters from configure file
    config_file = 'param.config'
    set_parameters(config_file)

    # Set default prefix if none entered
    check_prefix(options)

    print '\nInputted options: {}'.format(options)
    print '\nGrabbing data...'

    # Grab and process data from catalog, including flag removal
    data_obs = get_data(options)

    # Scale for linear fitting
    scaled_data = scale(*data_obs)

    print '\nFitting data...'

    # Fit data using linmix, lrgs, or both
    kelly_scaled_fit, mantz_scaled_fit = fit(*scaled_data[:4])
    (x_min, x_max) = (np.min(scaled_data[0]), np.max(scaled_data[0]))

    print '\nMaking plots...'

    # Make all desired plots
    plotlib.make_plots(
        options, PARAMETERS, METHODS, data_obs, kelly_scaled_fit,
        mantz_scaled_fit, scaled_data[4], x_min, x_max
    )

    if PARAMETERS['save_data'] is True:
        print '\nSaving data...'
        save_data(options, PARAMETERS, METHODS, data_obs, kelly_scaled_fit,
                  mantz_scaled_fit, scaled_data[4], x_min, x_max)

    print '\nDone!'



if __name__ == '__main__':

    main()
