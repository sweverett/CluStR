#-------------------------------------------------------
# We'll define useful classes here

class Config(object):
    '''
    Used for CluStR config processing
    Some options:
    - scale_luminosity: Divide luminosity columns by E(z)^-3/2
    '''
    _required_keys = []

    def __init__(self, config_file):
        with open(config_file, 'r') as stream:

            self._config = yaml.safe_load(stream)

        # TODO: Any logical parsing here:
        #...

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

class CheckDependencies:
    '''
    check if all required packages are
    installed, and, if not, ask if these packages should be downloaded.
    Otherwise exists program.
    '''
    def __init__(self,):
        pass

class Ez: #function
    def __init__(self,parameters,z):
        self.z = z
        self.Om = parameters['Om']
        self.H_0 = parameters['H_O']
        h = self.H_0/100
    def show(self):
        return np.sqrt(self.Om*(1.+self.z)**3 + h)

class FitsLabel:
    def __init__(self,labels,axis_name):
        self.labels = labels
        self.axis_name = axis_name
    def show(self):
        return self.labels[self.axis_name]

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
#pass in a file name not all the variables ^
    def scale(self):
        log_x = np.log(self.x)
        x_piv = np.median(log_x)
        log_y = np.log(self.y)
        return (log_x-x_piv, log_y, x_err/self.x, y_err/self.y, x_piv)

    def fit(self):
        pass
#fit that accepts data
class CheckPrefix:
    def __init__(self, options, parameters):
        self.options = options
        self.parameters = parameters
        pass
#class function of config
class SaveData:
    def __init__(self, options, parameters, , ,)
        pass
#class function in fitter class

#-------------------------------------------------------
# We'll write the main function here

def main():
    pass

if __name__ == '__main__':

    main()
