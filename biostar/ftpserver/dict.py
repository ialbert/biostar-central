
from collections.abc import MutableMapping




class UserTable(MutableMapping):


    def __getitem__(self, y):
        """ x.__getitem__(y) <==> x[y] """

        pass

    def __setitem__(self, *args, **kwargs):
        """ Set self[key] to value. """

        pass

    def __len__(self, *args, **kwargs):
        """ Return len(self). """

        return len(self.keys())


    def __iter__(self, *args, **kwargs):
        """ Implement iter(self). """

        return (k for k in self.keys())


    def __delitem__(self, *args, **kwargs):
        """ Delete self[key]. """
        # Take key, value pair out,
        pass

