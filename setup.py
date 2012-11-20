setup.py
Here is an example setup.py:

from distutils.core import setup
import sys

sys.path.append('googlemaps')
import googlemaps


setup(name='googlemaps',
      version='1.0',
      author='John Kleint',
      author_email='py-googlemaps-general@lists.sourceforge.net',
      url='http://sourceforge.net/projects/py-googlemaps/',
      download_url='https://sourceforge.net/projects/py-googlemaps/files/',
      description='Easy geocoding, reverse geocoding, driving directions, and local search in Python via Google.',
      long_description=googlemaps.GoogleMaps.__doc__,
      package_dir={'': 'googlemaps'},
      py_modules=['googlemaps'],
      provides=['googlemaps'],
      keywords='google maps local search ajax api geocode geocoding directions navigation json',
      license='Lesser Affero General Public License v3',
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Intended Audience :: Developers',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2',
                   'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
                   'License :: OSI Approved :: GNU Affero General Public License v3',
                   'Topic :: Internet',
                   'Topic :: Internet :: WWW/HTTP',
                   'Topic :: Scientific/Engineering :: GIS',
                  ],
     )
Don’t worry if you don’t have a download_url yet. Notice that we import our own module in order to pull in the __doc__ string for the long_description; you are free to use whatever text is appropriate. If you have a package instead of a module, you’d replace py_modules with packages and remove the package_dir. You can find the list of classifiers at PyPI’s list of trove classifiers. If your code depends on other third-party packages/modules, you can specify those with a required keyword argument.

This covers our code, but we need one more file to pull in the documentation: MANIFEST.in in the package root:

recursive-include doc/html *
prune doc/html/.doctrees/
exclude doc/html/.buildinfo
include LICENSE.txt
Now, to build your package, you just run setup.py sdist:

$ python setup.py sdist
If all goes well, this will create a tarball in dist/googlemaps-1.0.tar.gz. Let’s make sure there are no problems installing it and verify the presence of important files with cheesecake_index:

$ sudo easy_install cheesecake

$ cheesecake_index --path=dist/googlemaps-1.0.tar.gz

Ensure at the minimum that your package is able to be installed. Notice that the cheesecake index includes Pylint as one component, so you’re already ahead of the game. Personally I think the score is weighted a bit heavily toward installability and documentation, but a relative cheesecake index of at least 70% seems like a reasonable target.

Packaging Resources
distutils, the standard Python Distribution Utilities: http://docs.python.org/distutils/index.html
setuptools, an enhanced, extended distutils and home of easy_install: http://peak.telecommunity.com/DevCenter/setuptools
Cheesecake, package “kwalitee” checker: http://pycheesecake.org/
PyPI’s list of trove classifiers: http://pypi.python.org/pypi?%3Aaction=list_classifiers
