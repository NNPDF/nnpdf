"""
fileloader.py

Utilities to locate files in one or more paths.
These are useful to find templates in a way that suits us more than the jinja
loaders, where a lot of the functionality is either irrelevant or inadequate.

For example, we don't want to open the files that we find, and we want to
be able to recover the full path in every case.

The loader classes have two methods:

  - ``find`` that returns a tuple
``(parent, name)`` so that ``parent/name`` is an existing file. The split
depends on the particular finder class.

  - ``hint_files``  returns a list of entries in the same format as ``find``.
  The exact entries depend on the finder class.

"""
import itertools
import pathlib
import abc

class FinderError(Exception): pass
class FileNotInPaths(FileNotFoundError, FinderError): pass

class AbstractFinder(metaclass=abc.ABCMeta): # pragma: no cover
    @abc.abstractmethod
    def find(self, name): pass

    @abc.abstractmethod
    def hint_files(self): pass


class Finder(AbstractFinder):
    def __init__(self, path):
        """Locate files relative to ``path``.
        Note that all_files is not recursive."""
        self.path = pathlib.Path(path)
        if not self.path.is_dir():
            raise ValueError("Finder path must be a directory.")

    def hint_files(self):
        return ((self.path, p.name) for p in self.path.iterdir())

    def find(self, name):
        """Return a tuple ``(self.path, name)`` if ``self.path/name`` exists,
        otherwse raise `FileNotFoundError`. ``name`` must always be relative
        to the path."""
        #Dont' convert to path the name we are going to return, but keep it
        #as it was
        np = pathlib.Path(name)
        if np.is_absolute():
            raise FinderError(f'Invalid absolute path: {np}')

        lp = self.path / name
        if not lp.exists():
            raise FileNotInPaths(lp)
        return self.path, name

class ModuleFinder(Finder):
    def __init__(self, module):
        super().__init__(module.__path__[0])

class FallbackFinder(AbstractFinder):
    def __init__(self, paths):
        """Fallback finder from a list of Finders or paths"""
        self.finders = [p if isinstance(p, AbstractFinder) else Finder(p) for p in paths]

    def find(self, name):
        for f in self.finders:
            try:
                return f.find(name)
            except FileNotInPaths:
                pass
        raise FileNotInPaths(name)

    def hint_files(self):
        return itertools.chain.from_iterable(f.hint_files() for f in self.finders)
