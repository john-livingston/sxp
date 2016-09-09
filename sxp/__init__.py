import phot
import util
import vis


class SxpDataset(dict):

    """
    temporarily here until the output gets upgraded from current (bad) format.
    requires the user to supply a YAML config file with various config. values,
    see examples dir.
    """

    def __init__(self, setup, aor=None, geom='3x3', radius='setup'):

        self._config = setup['config']
        self['star'] = self._config['star']
        self['planet'] = self._config['planet']

        if aor:
            self['aor'] = aor
        else:
            self['aor'] = self._config['aor']
        if geom:
            self['geom'] = geom
        else:
            self['geom'] = self._config['geom']

        cube, t, f, r, s, c = util.load_data(self._config['data_dir'], self['aor'])

        self['pix'] = util.get_pix(cube, geom=self['geom'], normalize=True)
        self['t'] = t

        try:
            i = r.index(self._config['radius'])
        except:
            sys.exit('invalid radius selection')
        print("using radius: {}".format(r[i]))

        self['radius'] = r[i]
        self['f'] = f[i]
        self['fcor'] = f[i].copy() # initially same as f, but can be updated
        self['s'] = s[i]
