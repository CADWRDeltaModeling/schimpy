from setuptools import setup
import versioneer

VERSION = versioneer.get_version()

setup(
    name="schimpy", version=versioneer.get_version(), cmdclass=versioneer.get_cmdclass()
)
