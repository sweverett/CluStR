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

#-------------------------------------------------------
# We'll write the main function here

def main():
    pass

if __name__ == '__main__':
    main()
