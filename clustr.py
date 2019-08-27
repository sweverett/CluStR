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


class Ez: #function
    def __init__(self,parameters,z):
        self.z = z
        self.Om = parameters['Om']
        self.H_0 = parameters['H_O']
        h = self.H_0/100
    def show(self):
        return np.sqrt(self.Om*(1.+self.z)**3 + h)


class Flag:
    def __init__(self,flag,data,options,boolean,cut_or_range):
        self.flag = flags
        self.data = data
        self.boolean = boolean
        self.cut_or_range = cut_or_range
        self.options = options
    def check_flag(self):
        pass
#inheritance would be good here
class Catalog:
    column_labels = {
        'lambda': 'lambda',
        'l500kpc': '500_kiloparsecs_band_lumin',
        'lr2500': 'r2500_band_lumin',
        'lr500': 'r500_band_lumin',
        'lr500cc': 'r500_core_cropped_band_lumin',
        't500kpc': '500_kiloparsecs_temperature',
        'tr2500': 'r2500_temperature',
        'tr500': 'r500_temperature',
        'tr500cc': 'r500_core_cropped_temperature'
    }
    def __init__():
        pass
    def create_cuts(self):
        pass

    def get_data(self):

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

    def fit(self):
        pass
#fit that accepts data
#replace check_prefix with run_name
def run_name:
    def __init__(self, options, parameters):
        self.options = options
        self.parameters = parameters
        pass
#class function of config
class Fitter:
    pass
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
