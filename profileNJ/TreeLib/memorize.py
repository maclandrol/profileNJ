from collections import Hashable as hashable
from functools import partial


class memorize(object):
    """Cache function output when it's called and return it
    later when the same function is called with the same input,
    in this case, memorize use a hash to determine value to reevalute
    """

    def __init__(self, function):
        self.function = function
        self.cache = {}

    def __call__(self, hash, *args, **kwargs):
        """Call to memorize, (as decorator)"""

        if hash in self.cache:
            return self.cache[hash]

        elif not isinstance(hash, hashable) or hash is None:
            # hash is None or uncachable
            return self.function(*args, **kwargs)

        else:
            output = self.function(*args, **kwargs)
            self.cache[hash] = output
            return output

    def __repr__(self):
        """Return cached data"""
        for hash, data in self.cache:
            print hash, '=============>\n', data

    def __get__(self, obj, objtype):
        """Instance methods support"""
        return partial(self.__call__, obj)
